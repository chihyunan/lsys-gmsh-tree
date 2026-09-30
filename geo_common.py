# geo_common.py
import math


def unit_vector(v):
    mag = sum(c * c for c in v) ** 0.5
    if mag < 1e-9:
        return (0.0, 0.0, 1.0)
    return (v[0] / mag, v[1] / mag, v[2] / mag)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(a[k] * b[k] for k in range(3))


def sub(a, b):
    return tuple(a[k] - b[k] for k in range(3))


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def parse_frequency_hz(freq_val):
    if isinstance(freq_val, (int, float)):
        return float(freq_val)
    freq_str = str(freq_val).strip().upper()
    units = {"GHZ": 1e9, "MHZ": 1e6, "KHZ": 1e3, "HZ": 1.0}
    for unit, mult in units.items():
        if freq_str.endswith(unit):
            return float(freq_str.replace(unit, "").strip()) * mult
    return float(freq_str)


def boxes_overlap(b1, b2):
    return (b1[0] <= b2[3] and b1[3] >= b2[0] and
            b1[1] <= b2[4] and b1[4] >= b2[1] and
            b1[2] <= b2[5] and b1[5] >= b2[2])


def cylinder_aabb(p0, p1, radius, direction):
    box = [0.0] * 6
    for k in range(3):
        lateral = radius * math.sqrt(max(0.0, 1.0 - direction[k] ** 2))
        box[k] = min(p0[k], p1[k]) - lateral
        box[k + 3] = max(p0[k], p1[k]) + lateral
    return box


def closest_distance_segments(p1, q1, p2, q2):
    d1, d2, r = sub(q1, p1), sub(q2, p2), sub(p1, p2)
    a, e, f = dot(d1, d1), dot(d2, d2), dot(d2, r)
    eps = 1e-12
    if a <= eps and e <= eps:
        s, t = 0.0, 0.0
    elif a <= eps:
        s, t = 0.0, clamp(f / e, 0.0, 1.0)
    else:
        c = dot(d1, r)
        if e <= eps:
            t, s = 0.0, clamp(-c / a, 0.0, 1.0)
        else:
            b = dot(d1, d2)
            denom = a * e - b * b
            s = clamp((b * f - c * e) / denom, 0.0, 1.0) if denom > eps else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t, s = 0.0, clamp(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t, s = 1.0, clamp((b - c) / a, 0.0, 1.0)
    c1 = [p1[k] + s * d1[k] for k in range(3)]
    c2 = [p2[k] + t * d2[k] for k in range(3)]
    return sum((c1[k] - c2[k]) ** 2 for k in range(3)) ** 0.5


def cylinders_overlap(spec_i, spec_j):
    dist = closest_distance_segments(spec_i["base"], spec_i["tip"], spec_j["base"], spec_j["tip"])
    return dist <= (spec_i["radius"] + spec_j["radius"])


def group_overlapping_cylinders(specs):
    n = len(specs)
    parent = list(range(n))
    aabbs = [cylinder_aabb(s["base"], s["tip"], s["radius"], s["direction"]) for s in specs]

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            if boxes_overlap(aabbs[i], aabbs[j]) and cylinders_overlap(specs[i], specs[j]):
                union(i, j)

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())

# geo_common.py
def write_cylinder_solid(f, tag, base, vec, radius):
    """Write a single OpenCASCADE Cylinder() primitive: base point + axis
    vector (direction * height) + radius. Used identically by leaf disks
    and FE-BI shoot cylinders -- one code path, not two."""
    f.write(
        f"Cylinder({tag}) = {{{base[0]:.6f}, {base[1]:.6f}, {base[2]:.6f}, "
        f"{vec[0]:.6f}, {vec[1]:.6f}, {vec[2]:.6f}, {radius:.6f}}};\n"
    )
    
def write_boolean_fragments(f, tags):
    tags_str = ", ".join(str(t) for t in tags)
    # BooleanFragments computes conformal boundaries across all input volumes
    # and deletes the originals. Like BooleanUnion, no trailing semicolon is needed.
    f.write(f"BooleanFragments{{ Volume{{{tags_str}}}; Delete; }}{{}}\n")

def write_pairwise_boolean_union(f, group_members):
    """
    Captures the dynamically assigned output tag from Gmsh OpenCASCADE
    to prevent 'Unknown entity' errors on deleted input tags.
    """
    f.write(f"cur_tag = {group_members[0]};\n")
    for tool in group_members[1:]:
        f.write(
            f"tmp[] = BooleanUnion{{ Volume{{cur_tag}}; Delete; }}{{ Volume{{{tool}}}; Delete; }};\n"
        )
        f.write("cur_tag = tmp[0];\n")

def transform_lsys_to_hfss(p):
    """Maps Y-up L-system coordinates to HFSS Z-up space: (x, y, z) -> (x, -z, y)"""
    """Bakes in the +90deg X-axis rotation HFSS was applying post-import:
    (x, y, z) -> (x, -z, y). Doing it here means leaves (and shoot data)
    are written already in HFSS's frame -- no more group Rotate needed
    after import."""
    x, y, z = p
    return (x, -z, y)