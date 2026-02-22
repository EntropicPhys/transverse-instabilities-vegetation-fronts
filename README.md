# 🌱 Transverse instabilities of vegetation fronts: Biotic and abiotic drivers — codes 

This repository contains code used to generate the main numerical results in the manuscript:

> **Transverse instabilities of vegetation fronts: Biotic and abiotic drivers**  
> (authors, journal status, year)

The implementation is primarily based on **MATLAB + pde2path** of the one-dimensional planar-front and linearized eigenvalue problems.

---

## 🔧 Folder structure
```
bwhRS/
├── community_model/
│   ├── sG.m               # RHS of community model equations
│   ├── sGjac.m            # Analytical Jacobian
│   ├── bwhinit.m          # Model initialization
│   ├── bwhinit.m          # Set up starter pde2path struct. 
│   ├── sgbra.m            # Output quantifiers for cont.
│   ├── oosetfemops.m      # FEM matrix setup
│   ├── bfbb.m             # Stability range extraction from branches
│   ├── cswibra.m          # Modified version of cswibra for plotting purposes
│   ├── nplotbra.m         # Same as plotbra, with added normalization (1st argument) 
│   ├── sGdns.m            # Treat the chi-diffusion explicitly/ used for tintxs
│   ├── uplot1.m           # Special windows plot for continuation
│   ├── uplot2.m           # Only plot <B(.,chi)>  over chi 
│   ├── uplot3.m           # Also plot al(chi) for B_t=al(chi)*B+D_chi*B''
│   ├── tintxs.m           # Direct numerical simulation based on semi-implicit integration,
│   ├── bdbdt.m            # Transient Dynamics under Precipitation Changes
│   └── bbdns.m            # Time integration with continuously varying precipitation 
├── single_species_model/
│   ├── sG.m               # RHS for single species
│   ├── sGjac.m            # Jacobian for single species
│   ├── oosetfemops.m      # FEM matrix setup 
│   ├── bwhinit.m          # Model initialization
│   ├── sgbra.m            # Output quantifiers for cont.
│   ├── bpjac.m            # Matrix necessary for bif-cont.
│   ├── tint.m           # Direct numerical simulation based on semi-implicit integration,
│   ├── cmds1.m            # Continuation: homogeneous + Turing branches
│   └── cmds2.m            # Brute-force Busse-Balloon computation

```
---

## Requirements
- MATLAB (R2020 or later recommended)
- `pde2path` (https://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/)
- Dedalus (https://github.com/DedalusProject/dedalus/tree/master)
---

## Reproducibility map (fill in with your actual scripts)
- **Fig. 2 (planar front)**: `code/planar_front_1d/cont_planar_front.m`
- **Fig. 3 (kernel modes)**: `code/planar_front_1d/kernel_modes.m`
- **Fig. 4 (dispersion relation)**: `code/transverse_instability/dispersion_relation.m`
- **Fig. 5 (threshold curves)**: `code/transverse_instability/compute_thresholds.m`

---

## 📊 Parameter Table

| Index | Parameter | Description                | Symbol | Typical Value |
|-------|-----------|----------------------------|--------|----------------|
| 1     | `pp`      | Precipitation              | \(P\)  | 300            |
| 2     | `Lam0`    | Base growth rate           | \(\Lambda_0\) | 8        |
| 3     | `Ga`      | Mortality scale            | \(\Gamma\)     | 10       |
| 4     | `A`       | Max infiltration rate      | \(A\)  | 3000           |
| 5     | `R`       | Evaporation feedback       | \(R\)  | 0.7            |
| 6     | `L0`      | Base evaporation rate      | \(L_0\)| 200            |
| ...   | ...       | ...                        | ...    | ...            |
| 19    | `chimin`  | Min trait value            | \(\chi_{\min}\) | 0       |
| 20    | `chimax`  | Max trait value            | \(\chi_{\max}\) | 1       |

---

## 📎 References

- Main paper: [], *submitted/in preparation*

---

## Contact
- Name, affiliation
- Email (optional)


