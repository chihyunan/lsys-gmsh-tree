#!/usr/bin/env python3
"""
Generate a FE-BI cylinder .geo from a leaves manifest.
Individual cylinders are written to .geo without Gmsh Booleans.
Overlap groups are saved to a corresponding _groups.json file so HFSS 
can perform clean ACIS Unions per cluster.

Usage:
    python3 gmsh_to_febi.py   --frequency 1.26GHz
"""
import argparse
import json
import math
import os

from geo_common import (
    unit_vector,
    parse_frequency_hz,
    group_overlapping_cylinders,
    write_cylinder_solid,
)


def build_cylinder_specs(manifest, padding_mm):
    specs = []
    for shoot in manifest.get("shoots", []):
        p0, p1 = tuple(shoot["stem_start"]), tuple(shoot["stem_end"])
        d = unit_vector([p1[k] - p0[k] for k in range(3)])
        length = sum((p1[k] - p0[k]) ** 2 for k in range(3)) ** 0.5

        radius = shoot["max_radius"] + padding_mm
        height = length + 2.0 * padding_mm
        base = [p0[k] - padding_mm * d[k] for k in range(3)]
        tip = [p1[k] + padding_mm * d[k] for k in range(3)]

        specs.append({"base": base, "tip": tip, "radius": radius, "height": height, "direction": d})
    return specs


def write_febi_geo(specs, output_path):
    """Writes pure cylinder primitives to .geo with ZERO Gmsh boolean operations."""
    with open(output_path, "w") as f:
        f.write("// Auto-generated FE-BI domain cylinders (Unions handled in HFSS)\n")
        f.write('SetFactory("OpenCASCADE");\n\n')

        for idx, spec in enumerate(specs):
            tag = idx + 1
            vec = [spec["height"] * spec["direction"][k] for k in range(3)]
            write_cylinder_solid(f, tag, spec["base"], vec, spec["radius"])

    print(f"Wrote {len(specs)} cylinders to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    parser.add_argument("output")
    parser.add_argument(
        "--frequency",
        default="1.26GHz",
        help="must match the HFSS sim frequency -- drives lambda/10 padding",
    )
    args = parser.parse_args()

    with open(args.manifest, "r") as mf:
        manifest = json.load(mf)

    freq_hz = parse_frequency_hz(args.frequency)
    padding_mm = (299792458.0 / freq_hz) * 1000.0 / 10.0

    specs = build_cylinder_specs(manifest, padding_mm)
    groups = group_overlapping_cylinders(specs)

    # 1. Write clean .geo with no Gmsh OCC Booleans
    write_febi_geo(specs, args.output)

    # 2. Save group index lists using standard path convention
    groups_path = os.path.splitext(args.output)[0] + "_groups.json"
    with open(groups_path, "w") as gf:
        json.dump(groups, gf, indent=2)
    print(f"Wrote {len(groups)} union groups to {groups_path}")