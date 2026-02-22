# 🌱 Transverse instabilities of vegetation fronts: Biotic and abiotic drivers — codes 

This repository contains code used to generate the main numerical results in the manuscript:

> **Transverse instabilities of vegetation fronts: Biotic and abiotic drivers**  
> (authors, journal status, year)

The repository contains the minimal codes to reproduce the results of the manuscript: 
- Bifurcation diagram contains the calculation of the bifurcation diagram shown in Fig.~2 of the manuscript.
- Semi-analytic TI contains the calculation of data and posteriori post-processing to determine the curvature of the dispersion relation
- DNS containst a minimal setpup for two-dimensional direct numerical simulations script based on Dedalus library.  
---

## 🔧 Folder structure
```
TI/
├── Bifurcation diagram/
│   ├── sG.m               # RHS of the bwh-model in the co-moving frame
│   ├── sGjac.m            # Analytical Jacobian
│   ├── bwhinit.m          # Model initialization 
│   ├── oosetfemops.m      # FEM matrix setup
│   └── cmds.m             # bifurcation diagram
├── semi-analytic TI/
│   ├── sG.m               # RHS for single species
│   ├── sGjac.m            # Jacobian for single species
│   ├── oosetfemops.m      # FEM matrix setup 
│   ├── bwhinit.m          # Model initialization
│   ├── L.m                # numerical linear operator
│   ├── Lap.m              # numerical adjoint linear operator
│   ├── tint.m             # Direct numerical simulation based on semi-implicit integration,
│   ├── cmds1.m            # co-frame continuation
│   └── cmds2.m            # stability curve
├── DNS/
│   └── cmds.m             # two-dimensional simulation

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

## Contact
- Name, affiliation
- Email (optional)


