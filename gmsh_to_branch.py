#!/usr/bin/env python3
"""
Convert an lsys-gen 'generic' .out file into a Gmsh OpenCASCADE .geo file
containing tapered cone segments fused into one solid tree, following the
same pattern as the MATLAB aspen-branch generator (Cone() per segment,
overlap trick to guarantee clean boolean joins, single BooleanUnion at end).

Usage:
    python3 lsys_to_geo.py <input.out> <output.geo> [--scale 10] [--radius-scale 1.0] [--overlap 0.5]

    --scale         multiplies all XYZ coordinates (raw turtle units -> mm)
    --radius-scale  multiplies the group "Width" value to get segment radius
    --overlap       fraction of segment length each cone is extended past
                    its real endpoints, to guarantee overlap for BooleanUnion
"""
import sys
import argparse

from geo_common import transform_lsys_to_hfss


def parse_lsys_output(path):
    """Parse groups of connected line segments and their width."""
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
            p1 = list(map(float, lines[i + 1].split()))
            p2 = list(map(float, lines[i + 2].split()))
            if current_group is not None:
                current_group.append((p1, p2))
            i += 2

        i += 1

    return groups


def chain_segments_to_polylines(segments, tol=1e-6):
    """
    Merge a list of (p1, p2) line segments into ordered polylines, where
    consecutive segments share an endpoint. Falls back to independent
    2-point polylines if nothing chains.
    """
    def close(a, b):
        return all(abs(a[k] - b[k]) < tol for k in range(3))

    polylines = []
    used = [False] * len(segments)

    for i in range(len(segments)):
        if used[i]:
            continue
        used[i] = True
        chain = [segments[i][0], segments[i][1]]

        extended = True
        while extended:
            extended = False
            for j in range(len(segments)):
                if used[j]:
                    continue
                p1, p2 = segments[j]
                if close(chain[-1], p1):
                    chain.append(p2)
                    used[j] = True
                    extended = True
                elif close(chain[-1], p2):
                    chain.append(p1)
                    used[j] = True
                    extended = True

        polylines.append(chain)

    return polylines


def unit_vector(v):
    mag = sum(c * c for c in v) ** 0.5
    if mag < 1e-12:
        return (0.0, 0.0, 1.0)
    return (v[0] / mag, v[1] / mag, v[2] / mag)


def write_geo(groups, output_path, scale, radius_scale, overlap_fraction):
    cone_id = 0

    with open(output_path, "w") as f:
        f.write("// Auto-generated from lsys-gen .out file\n")
        f.write("// Tapered cone segments fused into one solid tree\n")
        f.write('SetFactory("OpenCASCADE");\n')
        f.write("Geometry.Tolerance = 1e-7;\n\n")

        for group in groups:
            width = group["width"] if group["width"] else 1.0
            radius = max(0.01, (width * radius_scale) / 2.0)

            polylines = chain_segments_to_polylines(group["segments"])

            for chain in polylines:
                n_segs = len(chain) - 1
                for k in range(n_segs):
                    # Apply transformation directly during scaling
                    p0 = transform_lsys_to_hfss([c * scale for c in chain[k]])
                    p1 = transform_lsys_to_hfss([c * scale for c in chain[k + 1]])

                    direction = unit_vector([p1[m] - p0[m] for m in range(3)])
                    seg_len = sum((p1[m] - p0[m]) ** 2 for m in range(3)) ** 0.5
                    ext = overlap_fraction * seg_len

                    # extend endpoints slightly so neighboring cones overlap
                    q0 = [p0[m] - ext * direction[m] if k == 0 else p0[m] for m in range(3)]
                    q1 = [p1[m] + ext * direction[m] if k == n_segs - 1 else p1[m] for m in range(3)]

                    vec = [q1[m] - q0[m] for m in range(3)]

                    cone_id += 1
                    # OpenCASCADE's Cone() rejects equal r0/r1 (degenerate --
                    # that's just a cylinder). Since our per-group width is
                    # constant, use Cylinder() here; Cone() is reserved for
                    # a future version that tapers radius across groups.
                    f.write(
                        f"Cylinder({cone_id}) = {{{q0[0]:.6f}, {q0[1]:.6f}, {q0[2]:.6f}, "
                        f"{vec[0]:.6f}, {vec[1]:.6f}, {vec[2]:.6f}, "
                        f"{radius:.6f}}};\n"
                    )

        if cone_id > 1:
            tags = ", ".join(str(t) for t in range(2, cone_id + 1))
            f.write("\n")
            f.write(
                f'tree() = BooleanUnion {{ Volume{{1}}; Delete; }} '
                f'{{ Volume{{{tags}}}; Delete; }};\n'
            )

    print(f"Wrote {cone_id} cone segments to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--scale", type=float, default=10.0,
                         help="multiplies XYZ coords (raw turtle units -> mm)")
    parser.add_argument("--radius-scale", type=float, default=1.0,
                         help="multiplies group Width to get segment radius")
    parser.add_argument("--overlap", type=float, default=0.1,
                         help="fraction of segment length to extend cone ends for clean unions")
    args = parser.parse_args()

    groups = parse_lsys_output(args.input)
    print(f"Parsed {len(groups)} groups")
    write_geo(groups, args.output, args.scale, args.radius_scale, args.overlap)