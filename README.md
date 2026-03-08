# 🌱 Transverse instabilities of vegetation fronts: Biotic and abiotic drivers — codes 

This repository contains the minimal codes used to generate the main numerical results in the manuscript:

> **Transverse instabilities of vegetation fronts: Biotic and abiotic drivers**  
> (authors, journal status, year)

The repository contains the minimal code to reproduce the results of the manuscript: 
- The bifurcation diagram contains the calculation of the bifurcation diagram shown in Fig.~2 of the manuscript.
- Semi-analytic TI contains the calculation of data and a posteriori post-processing to determine the curvature of the dispersion relation
- DNS contains a minimal setup for a two-dimensional direct numerical simulation script based on the Dedalus library.  
---

## 🔧 Folder structure
```
TI/
├── Bifurcation diagram/
│   ├── sG.m               # RHS of the bwh-model in the co-moving frame
│   ├── sGjac.m            # Analytical Jacobian
│   ├── bwhinit.m          # Model initialization 
│   ├── oosetfemops.m      # FEM matrix setup
│   ├── sgbra.m            # continuation output
│   └── cmds.m             # bifurcation diagram
├── semi-analytic TI/
│   ├── sG.m               # RHS for single species
│   ├── sGjac.m            # Jacobian for single species
│   ├── oosetfemops.m      # FEM matrix setup 
│   ├── bwhinit.m          # Model initialization
│   ├── tintfreeze.m       # Direct numerical simulation based on semi-implicit integration,
│   ├── LAp.m              # numerical linear operator
│   ├── LOp.m              # numerical adjoint linear operator
│   ├── compute_lam_eta.m  # Compute the spectrum slope 
│   └── cmds.m             # co-frame continuation and stability curve
├── DNS/
│   └── cmds.m             # two-dimensional simulation with Dedalus

```
---

## Requirements
- MATLAB (R2020 or later recommended)
- `pde2path` (https://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/)
- Dedalus (https://github.com/DedalusProject/dedalus/tree/master)
---

## 📎 References

- Main paper: [], *submitted/in preparation*

---


