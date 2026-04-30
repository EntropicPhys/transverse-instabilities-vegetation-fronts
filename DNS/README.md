# Direct numerical simulations

This folder contains the Python/Dedalus code used to run two-dimensional direct numerical simulations of perturbed vegetation fronts.

The main file is `cmds.py`. This file is executed from the terminal using Python, preferably with MPI.

## How to run

From a terminal, move to this folder and run:

```bash
cd DNS
mpiexec -n 4 python3 cmds.py
```

The number after `-n` is the number of MPI processes. For example, the command above runs the simulation using four processes.

For a small serial test, one may run:

```bash
python3 cmds.py
```

or equivalently:

```bash
mpiexec -n 1 python3 cmds.py
```

## Workflow

The script `cmds.py` performs the following steps:

1. Defines the dimensional parameters of the BWH model.
2. Converts the dimensional parameters into nondimensional form.
3. Computes the homogeneous vegetation state used in the initialization.
4. Sets up a two-dimensional Dedalus problem.
5. Uses a Fourier basis in the periodic direction and a Chebyshev basis in the bounded direction.
6. Evolves the system in a co-moving frame.
7. Monitors diagnostic quantities such as the residual and transverse modulation.
8. Writes simulation snapshots to a local output folder.

## Output

The simulation output is written to a Dedalus output folder whose name is generated from the parameter values `f` and `eta`.

The folder name has the form

```text
bwh_FI_f_<f>_eta_<eta>_var_2
```

with decimal points replaced by `p`.

## Notes

This script is intended as a minimal reproduction script for the direct numerical simulations. The default resolution may be computationally expensive, so it is recommended to test the script first with a reduced resolution or shorter final time before running full simulations on a workstation or computing cluster.
