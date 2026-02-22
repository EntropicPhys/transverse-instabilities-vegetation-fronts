# Transverse instabilities of vegetation fronts: Biotic and abiotic drivers — code & data

This repository contains code and (lightweight) data used to generate the main numerical results in the manuscript:

> **Transverse instabilities of vegetation fronts: Biotic and abiotic drivers**  
> (authors, journal status, year)

The implementation is primarily based on **MATLAB + pde2path** finite-element discretizations of the relevant 1D planar-front and 2D linearized eigenvalue problems.

---

## Folder structure
- `code/planar_front_1d/`: compute planar front profiles and the discrete translational/adjoint kernel modes used in the dispersion-relation formula.
- `code/transverse_instability/`: compute dispersion relations and instability thresholds as a function of transverse wavenumber and parameters.
- `code/simulations_2d/` (optional): direct numerical simulations used for verification/illustration.
- `data/processed/`: exported profiles, kernel modes, and spectra used in plotting.
- `results/`: heavy outputs (branches, logs, checkpoints). **Not tracked by git**.
- `figures/`: generated figures (can be tracked or regenerated).
- `manuscript/`: LaTeX source (or pointers to Overleaf) and bibliography.

---

## Quick start
1. Clone the repository.
2. Install requirements (below).
3. In MATLAB, run:
   - `code/utils/setup_pde2path.m` (adds pde2path paths, sets defaults)
   - `code/transverse_instability/compute_thresholds.m`

The scripts write outputs to `results/` and export lightweight artifacts to `data/processed/`.

---

## Requirements
- MATLAB (R2020 or later recommended)
- `pde2path` (finite-element continuation framework)
- OOPDE library (shipped with / linked from pde2path)

---

## Reproducibility map (fill in with your actual scripts)
- **Fig. 2 (planar front)**: `code/planar_front_1d/cont_planar_front.m`
- **Fig. 3 (kernel modes)**: `code/planar_front_1d/kernel_modes.m`
- **Fig. 4 (dispersion relation)**: `code/transverse_instability/dispersion_relation.m`
- **Fig. 5 (threshold curves)**: `code/transverse_instability/compute_thresholds.m`

---

## Parameter table (fill in)
| Symbol | Meaning | Typical value | Units |
|---|---|---:|---|
| `P` | precipitation |  |  |
| `D_W` | water diffusion |  |  |
| `D_H` | surface-water diffusion |  |  |
| ... | ... | ... | ... |

---

## Citation
If you use this repository, please cite the manuscript and/or the archived release (Zenodo DOI).

[![DOI](https://zenodo.org/badge/<YOUR_BADGE_ID>.svg)](https://doi.org/10.5281/zenodo.<YOUR_RECORD_ID>)

---

## License
Choose a license (e.g. MIT/BSD-3-Clause/GPL-3.0) and add `LICENSE`.

---

## Contact
- Name, affiliation
- Email (optional)


