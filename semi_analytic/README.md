# Semi-analytic transverse-instability calculation

This folder contains the MATLAB/pde2path code used to compute co-moving front solutions and the associated semi-analytic transverse-instability calculation for the BWH vegetation model.

The main file is `cmds.m`. This file is executed in MATLAB after `pde2path` has been initialized and added to the MATLAB path.

## How to run

From MATLAB, move to this folder and run:

```matlab
cd semi_analytic
cmds
```

Alternatively, if MATLAB is already open in this folder, simply run:

```matlab
cmds
```

## Workflow

The script `cmds.m` performs the following steps:

1. Defines the dimensional parameters of the BWH model.
2. Converts the dimensional parameters into the nondimensional parameters used in the continuation problem.
3. Initializes a one-dimensional co-moving front solution.
4. Uses time integration to obtain a frozen front state.
5. Continues the front branch using pde2path.
6. Loads the computed continuation points.
7. Computes the semi-analytic transverse-instability quantity.
8. Plots the resulting stability curve and diagnostic eigenvalues.

## Output

The continuation data are written to a local pde2path folder whose name is determined by the value of the parameter `ff`.

For example, if

```matlab
ff = 0.7;
```

then the output folder is

```text
0p7
```

The script also produces diagnostic plots for the spectrum slope, the approximate translational eigenvalue, and the approximate adjoint-kernel eigenvalue.

## Modifying parameters

The main model parameters are defined near the beginning of `cmds.m`. To study a different parameter regime, edit the parameter values in that file before running the script.

In particular, the parameter `ff` determines the output folder name and the front branch being computed.

## Notes

This script combines time integration, numerical continuation, and post-processing. It is intended as a minimal reproduction script for the semi-analytic transverse-instability calculation. Users who modify the continuation settings should be familiar with the basic pde2path workflow.
