#!/usr/bin/env python3
"""
Generate leaves.geo from an lsys-gen .out file: finds terminal branch tips
(degree-1 nodes in the branch graph, excluding the root) and attaches
floating leaf blades there, oriented using the same tangent + perpendicular
basis + golden-angle phyllotaxis approach as the MATLAB aspen reference.

No petiole, no boolean subtraction against the wood -- leaves float near
the tip. This is intentionally simpler than the MATLAB script since the
tree.geo wood is generated separately by lsys_to_geo.py.

Usage:
    python3 lsys_to_leaves_geo.py <input.out> <output.geo> \\
        [--scale 10] [--leaves-per-tip 3] [--leaf-radius 8] \\
        [--standoff 1.5] [--seed 42]
"""
import argparse
import json
import math
import os
import random

from geo_common import transform_lsys_to_hfss


def parse_lsys_output(path):
    """Parse groups of raw (unchained) line segments and their width."""
    groups = []
    current_group = None
    current_width = None

    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("Start Group"):
            current_group = []
            current_width = None
        elif line.startswith("Width:"):
            current_width = float(line.split(":")[1].strip())
        elif line.startswith("End Group"):
            if current_group:
                groups.append({"segments": current_group, "width": current_width})
            current_group = None
        elif line == "line":
            p1 = tuple(map(float, lines[i + 1].split()))
            p2 = tuple(map(float, lines[i + 2].split()))
            if current_group is not None:
                current_group.append((p1, p2))
            i += 2
        i += 1

    return groups


def round_point(p, ndigits=5):
    return tuple(round(c, ndigits) for c in p)


def build_graph(groups):
    """
    Build an undirected graph over all raw segment endpoints (rounded to
    merge floating-point duplicates). Returns:
        edges: list of (p1, p2, width)
        degree: dict point -> number of connected edges
        root: the very first point seen (start of the trunk)
    """
    edges = []
    degree = {}
    root = None

    for group in groups:
        width = group["width"] if group["width"] else 1.0
        for p1, p2 in group["segments"]:
            rp1, rp2 = round_point(p1), round_point(p2)
            if root is None:
                root = rp1
            edges.append((rp1, rp2, width))
            degree[rp1] = degree.get(rp1, 0) + 1
            degree[rp2] = degree.get(rp2, 0) + 1

    return edges, degree, root


def build_adjacency(edges):
    """point -> list of (neighbor_point, width) for the edge connecting them."""
    adj = {}
    for p1, p2, width in edges:
        adj.setdefault(p1, []).append((p2, width))
        adj.setdefault(p2, []).append((p1, width))
    return adj


def find_shoots(adj, root):
    """
    For each terminal tip (degree-1 node, not root), walk backward through
    the chain of degree-2 pass-through nodes until hitting a real branching
    point (degree != 2) or the root. This reconstructs the full terminal
    'shoot' polyline -- e.g. tip -> c2 -> c1 -> branch_junction -- instead
    of just the single last edge. Returns each shoot as ordered from base
    (junction) to tip, matching the MATLAB shoot.points convention, along
    with the per-segment widths.
    """
    shoots = []
    for node, neighbors in adj.items():
        if len(neighbors) == 1 and node != root:
            chain = [node]
            widths = []
            current_neighbor, w = neighbors[0]
            chain.append(current_neighbor)
            widths.append(w)
            prev, current = node, current_neighbor

            while current != root and len(adj[current]) == 2:
                n1, w1 = adj[current][0]
                n2, w2 = adj[current][1]
                nxt, w_next = (n2, w2) if n1 == prev else (n1, w1)
                chain.append(nxt)
                widths.append(w_next)
                prev, current = current, nxt

            # reverse so index 0 = base (junction/root), last index = tip
            points = list(reversed(chain))
            seg_widths = list(reversed(widths))
            shoots.append({"points": points, "widths": seg_widths})

    return shoots

def sample_polyline(points, widths, fraction):
    """Arc-length interpolation along a polyline: returns (point, tangent, width)."""
    seg_lengths = []
    for i in range(len(points) - 1):
        d = [points[i + 1][k] - points[i][k] for k in range(3)]
        seg_lengths.append(sum(c * c for c in d) ** 0.5)

    total = sum(seg_lengths)
    if total < 1e-9:
        return points[0], (0.0, 0.0, 1.0), widths[0]

    target = max(0.0, min(1.0, fraction)) * total
    cumulative = 0.0
    for i, seg_len in enumerate(seg_lengths):
        if cumulative + seg_len >= target or i == len(seg_lengths) - 1:
            local_t = 0.0 if seg_len < 1e-9 else (target - cumulative) / seg_len
            p0, p1 = points[i], points[i + 1]
            point = tuple(p0[k] + local_t * (p1[k] - p0[k]) for k in range(3))
            tangent = unit_vector([p1[k] - p0[k] for k in range(3)])
            width = widths[i]
            return point, tangent, width
        cumulative += seg_len

    return points[-1], unit_vector([points[-1][k] - points[-2][k] for k in range(3)]), widths[-1]


def unit_vector(v):
    mag = sum(c * c for c in v) ** 0.5
    if mag < 1e-9:
        return (0.0, 0.0, 1.0)
    return (v[0] / mag, v[1] / mag, v[2] / mag)

def point_to_line_distance(point, line_point, line_dir_unit):
    """Perpendicular distance from `point` to the infinite line through
    `line_point` with unit direction `line_dir_unit`."""
    w = [point[k] - line_point[k] for k in range(3)]
    proj_len = dot(w, line_dir_unit)
    proj = [line_point[k] + proj_len * line_dir_unit[k] for k in range(3)]
    return sum((point[k] - proj[k]) ** 2 for k in range(3)) ** 0.5

def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return sum(a[k] * b[k] for k in range(3))


def perpendicular_basis(direction):
    d = unit_vector(direction)
    reference = (0.0, 0.0, 1.0)
    if abs(dot(d, reference)) > 0.9:
        reference = (0.0, 1.0, 0.0)
    e1 = unit_vector(cross(d, reference))
    e2 = unit_vector(cross(d, e1))
    return e1, e2


def write_leaf_solid(f, tag, center, n, radius, thickness):
    """
    A 'circle with thickness' is just a short, wide cylinder -- use gmsh's
    analytic OpenCASCADE Cylinder() primitive directly instead of building
    a spline-based polygon + plane surface. This sidesteps floating-point
    planarity failures entirely (no user-supplied points/splines to go
    slightly out-of-plane) and meshes more predictably/cheaply in HFSS
    than a faceted NURBS surface would.

    Leaves are deliberately a fixed-radius circle: this keeps the written
    geometry consistent with the center-distance overlap guard used during
    placement.
    """
    base = [center[k] - 0.5 * thickness * n[k] for k in range(3)]
    vec = [thickness * n[k] for k in range(3)]
    f.write(
        f"Cylinder({tag}) = {{{base[0]:.6f}, {base[1]:.6f}, {base[2]:.6f}, "
        f"{vec[0]:.6f}, {vec[1]:.6f}, {vec[2]:.6f}, {radius:.6f}}};\n"
    )


class SpatialHashGrid:
    """
    Buckets points into cells of size `cell_size` so collision checks only
    look at nearby cells instead of every point ever inserted. Turns the
    leaf-placement check from O(n^2) into ~O(n) on average, since real
    trees spread leaves through 3D space rather than cramming them all
    into one tiny region.
    """
    def __init__(self, cell_size):
        self.cell_size = max(cell_size, 1e-6)
        self.cells = {}

    def _cell_index(self, point):
        return tuple(int(math.floor(point[k] / self.cell_size)) for k in range(3))

    def insert(self, point):
        self.cells.setdefault(self._cell_index(point), []).append(point)

    def any_within(self, point, min_distance):
        cx, cy, cz = self._cell_index(point)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in self.cells.get((cx + dx, cy + dy, cz + dz), []):
                        d = sum((point[k] - other[k]) ** 2 for k in range(3)) ** 0.5
                        if d <= min_distance:
                            return True
        return False


def write_leaves_geo(shoots, output_path, scale, leaves_per_tip, leaf_radius,
                      standoff, thickness, seed):
    random.seed(seed)
    golden_angle = math.pi * (3 - math.sqrt(5))
    manifest_records = []

    with open(output_path, "w") as f:
        f.write("// Auto-generated floating leaves for terminal branch tips\n")
        f.write("SetFactory(\"OpenCASCADE\");\n")
        f.write("Geometry.Tolerance = 1e-5;\n\n")

        cylinder_tag = 0
        leaf_count = 0
        unresolved_collisions = 0
        # Each leaf is a radius-``leaf_radius`` circular cylinder, so a
        # center-distance of 2 radii is a conservative non-overlap bound
        # regardless of its orientation.
        min_separation = leaf_radius * 2.0
        grid = SpatialHashGrid(cell_size=min_separation)

        shoot_records = []
        for shoot_index, shoot in enumerate(shoots):
            points = [transform_lsys_to_hfss(tuple(c * scale for c in p)) for p in shoot["points"]]
            widths = shoot["widths"]

            stem_start = points[0]
            stem_end = points[-1]
            axis_dir = unit_vector([stem_end[k] - stem_start[k] for k in range(3)])
            shoot_max_radius = 0.0

            # spread leaves along the shoot length, base to tip, with a
            # little jitter so they don't line up too mechanically
            fractions = [0.13 + (0.92 - 0.13) * i / max(1, leaves_per_tip - 1)
                         for i in range(leaves_per_tip)]
            fractions = [max(0.05, min(0.98, fr + random.uniform(-0.03, 0.03)))
                         for fr in fractions]

            for i, fraction in enumerate(fractions):
                point, tangent, width = sample_polyline(points, widths, fraction)
                e1, e2 = perpendicular_basis(tangent)
                radius_offset = max(0.5, width * scale / 2.0) + standoff

                # try several azimuth angles; accept the first one that
                # clears minimum separation from every already-placed leaf
                # (same strategy as the MATLAB reference's retry loop)
                center = None
                for attempt in range(8):
                    theta = (i * golden_angle + shoot_index * 0.47
                             + random.uniform(0, 2 * math.pi) + attempt * 0.55)
                    outward = unit_vector([
                        math.cos(theta) * e1[k] + math.sin(theta) * e2[k]
                        for k in range(3)
                    ])
                    u_try = unit_vector([outward[k] + 0.2 * tangent[k] for k in range(3)])
                    candidate = [
                        point[k] + radius_offset * outward[k] + leaf_radius * 0.5 * u_try[k]
                        for k in range(3)
                    ]

                    is_separated = not grid.any_within(candidate, min_separation)

                    if is_separated:
                        center = candidate
                        u, outward_final = u_try, outward
                        break

                if center is None:
                    # no separated placement found after 8 tries -- skip
                    # this leaf entirely rather than risk an intersection
                    unresolved_collisions += 1
                    continue

                grid.insert(center)
                v = unit_vector(cross(outward_final, u))
                n_axis = unit_vector(cross(u, v))

                # Keep every leaf at the requested radius.  Previously this
                # was randomized after collision placement, meaning the
                # fixed min_separation check did not match the geometry
                # actually written to the .geo file.
                size = leaf_radius

                cylinder_tag += 1
                write_leaf_solid(
                    f, cylinder_tag, center, n_axis, size, thickness
                )
                leaf_count += 1
                
                # NEW: track how far this leaf's edge sits from the shoot's own axis
                dist = point_to_line_distance(center, stem_start, axis_dir) + size
                if dist > shoot_max_radius:
                    shoot_max_radius = dist

                manifest_records.append({
                    "leaf_tag": cylinder_tag,
                    "shoot_index": shoot_index,
                    "center": [round(c, 6) for c in center],
                    "radius": round(size, 6),
                    "thickness": round(thickness, 6),
                })

            # NEW: only emit a shoot record if it actually got at least one leaf
            if shoot_max_radius > 0:
                shoot_records.append({
                    "shoot_index": shoot_index,
                    "stem_start": [round(c, 6) for c in stem_start],
                    "stem_end": [round(c, 6) for c in stem_end],
                    "max_radius": round(shoot_max_radius, 6),
                })

        print(f"Wrote {leaf_count} leaves (as cylinders) across {len(shoots)} shoots to {output_path}")
        if unresolved_collisions > 0:
            print(f"NOTE: {unresolved_collisions} leaf(ves) were skipped because no "
                  f"non-overlapping placement was found after 8 attempts. Consider "
                  f"increasing --standoff, reducing --leaves-per-tip, or a lower "
                  f"--leaf-radius if you need every requested leaf placed.")

    manifest_path = os.path.splitext(output_path)[0] + "_manifest.json"
    with open(manifest_path, "w") as mf:
        json.dump({"leaves": manifest_records, "shoots": shoot_records, "num_shoots": len(shoots)}, mf, indent=2)
    print(f"Wrote manifest ({len(manifest_records)} leaf records) to {manifest_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--scale", type=float, default=10.0,
                         help="must match the --scale used for lsys_to_geo.py")
    parser.add_argument("--leaves-per-tip", type=int, default=3)
    parser.add_argument("--leaf-radius", type=float, default=27.0,
                         help="approx leaf half-length in mm")
    parser.add_argument("--standoff", type=float, default=0.5,
                         help="extra gap (mm) beyond branch radius before leaf starts")
    parser.add_argument("--thickness", type=float, default=0.4)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    groups = parse_lsys_output(args.input)
    edges, degree, root = build_graph(groups)
    adj = build_adjacency(edges)
    shoots = find_shoots(adj, root)
    print(f"Found {len(shoots)} terminal shoots out of {len(degree)} nodes")

    write_leaves_geo(
        shoots, args.output, args.scale, args.leaves_per_tip,
        args.leaf_radius, args.standoff, args.thickness, args.seed
    )
