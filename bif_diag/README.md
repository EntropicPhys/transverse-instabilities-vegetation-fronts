# Bifurcation diagram

This folder contains the MATLAB/pde2path code used to compute the bifurcation diagram for the BWH vegetation model.

The main file is `cmds.m`. It should be executed in MATLAB after `pde2path` has been initialized and added to the MATLAB path.

## How to run

From MATLAB, move to this folder and run:

```matlab
cd bif_diag
cmds
```

Alternatively, if MATLAB is already open in this folder, simply run:

```matlab
cmds
```

## Workflow

The script `cmds.m` performs the following steps:

1. Defines the dimensional parameters of the BWH model.
2. Initializes the homogeneous vegetation branch.
3. Continues the homogeneous branch in the precipitation parameter.
4. Switches branch at the first bifurcation point.
5. Continues the patterned branch.
6. Computes the bare-soil branch.
7. Plots the resulting bifurcation diagram.

## Output

The continuation data are written to local pde2path folders, including

```text
fig2/hom
fig2/T1
fig2/hom1
```

## Modifying parameters

The main model parameters are defined near the beginning of `cmds.m`. To reproduce a different parameter regime, edit the parameter values in that file before running the script.

## Notes

This script is intended as a minimal reproduction script for the bifurcation diagram. It is not written as a general-purpose continuation package. Users who modify the continuation settings should be familiar with the basic pde2path workflow.
