#!/usr/bin/env python3
"""
run by python3 gmsh_main.py

run_pipeline.py -- single entry point for lsys .out -> .geo -> .step generation.

Runs, in order:
  1. run-example.sh              (lsys-gen -> .out)
  2. gmsh_to_branch.py + gmsh    (.out -> wood .geo -> wood .step)
  3. gmsh_to_leaves.py + gmsh    (.out -> leaves .geo -> leaves .step, + manifest.json)
  4. gmsh_to_febi.py + gmsh      (manifest.json -> febi .geo -> febi .step)

step_to_febi.py (HFSS import + FE-BI/plane-wave/analysis setup) is intentionally
NOT called here -- it's run separately once the .step files above exist.

Edit CONFIG below; nothing else should need to change for a normal run.
"""
import os
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CONFIG = {
    "RUN_EXAMPLE_SH": "./run-example.sh",
    "TREE_INPUT_DIR": "Tree_inputs",
    "TREE_NAME": "2.3ternary_A",

    "BRANCH_SCALE": 50, #default 50
    "RADIUS_SCALE": 30, #default 30
    "OVERLAP": 0.05, #small enough to avoid HFSS clash 

    "LEAF_SCALE": 50, #should be same as branch
    "LEAVES_PER_TIP": 20, #multiply by #of terminal branches to get leaf count 
    "LEAF_RADIUS": 27.1, #mm

    "FREQUENCY": "1.26GHz",  # drives lambda/10 FE-BI padding in gmsh_to_febi.py
}


def run(cmd):
    print("+ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, cwd=SCRIPT_DIR)


def main():
    base = os.path.join(CONFIG["TREE_INPUT_DIR"], CONFIG["TREE_NAME"])

    out_file = base + ".out"

    wood_geo = base + ".geo"
    wood_step = base + ".step"

    leaves_geo = base + "_leaves.geo"
    leaves_step = base + "_leaves.step"
    manifest_json = base + "_leaves_manifest.json"

    febi_geo = base + "_febi.geo"
    febi_step = base + "_febi.step"

    # 1. lsys-gen: produce the raw .out branch/segment data
    run([CONFIG["RUN_EXAMPLE_SH"], CONFIG["TREE_INPUT_DIR"], base])

    # 2. Wood geometry
    run([
        "python3", "gmsh_to_branch.py", out_file, wood_geo,
        "--scale", str(CONFIG["BRANCH_SCALE"]),
        "--radius-scale", str(CONFIG["RADIUS_SCALE"]),
        "--overlap", str(CONFIG["OVERLAP"]),
    ])
    run(["gmsh", wood_geo, "-0", "-o", wood_step])

    # 3. Leaf geometry (also writes the leaves manifest FE-BI depends on)
    run([
        "python3", "gmsh_to_leaves.py", out_file, leaves_geo,
        "--scale", str(CONFIG["LEAF_SCALE"]),
        "--leaves-per-tip", str(CONFIG["LEAVES_PER_TIP"]),
        "--leaf-radius", str(CONFIG["LEAF_RADIUS"]),
    ])
    run(["gmsh", leaves_geo, "-0", "-o", leaves_step])

    # 4. FE-BI domain cylinders, pre-unioned where shoots overlap
    run([
        "python3", "gmsh_to_febi.py", manifest_json, febi_geo,
        "--frequency", CONFIG["FREQUENCY"],
    ])
    run(["gmsh", febi_geo, "-0", "-o", febi_step])

    print("\nDone. STEP files ready:")
    print("  wood:  {}".format(wood_step))
    print("  leaves:{}".format(leaves_step))
    print("  febi:  {}".format(febi_step))
    print("\nRun step_to_febi.py separately with these paths to continue in HFSS.")


if __name__ == "__main__":
    main()