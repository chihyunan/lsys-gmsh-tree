# lsys-gmsh-tree

A procedural tree model generator that converts L-System outputs into 3D CAD geometries (`.step`) and FE-BI boundary domains via Gmsh for electromagnetic simulations (e.g., Ansys HFSS).

## Origin

This repository extends the [LSys](https://github.com/glk1001/LSys) engine, originally written by Jon Leech (UNC, 1991) based on *The Algorithmic Beauty of Plants* by P. Prusinkiewicz and A. Lindenmayer.

## Additions

Added an automated Python-to-Gmsh processing pipeline to convert raw L-System data (`.out`) into 3D STEP geometries:
- **`gmsh_to_branch.py`**: Generates 3D branch/wood geometry (`.geo` / `.step`).
- **`gmsh_to_leaves.py`**: Generates leaf geometries and outputs a leaf placement manifest (`manifest.json`).
- **`gmsh_to_febi.py`**: Generates pre-unioned FE-BI domain boundary cylinders around canopy clusters based on wavelength. (Currently not working)

## Usage

### Requirements
- Python 3
- [Gmsh](https://gmsh.info/) (must be added to system `PATH`)

### Running the Script

Execute the master entry-point script:

```bash
python3 gmsh_main.py 
```

- Default input is Tree_inputs/2.3ternary_A from tree example 
- modify or swap the input file for different tree 
- outputs also sent to Tree_inputs/ directory. 

### Output Files 

- input_leaves.step : leaves step 
- input.step : branches step 
