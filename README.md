# Phase-locked flow decomposition and energy triad analysis

This repository contains MATLAB code for phase-locked flow decomposition of time-resolved velocity fields and for the computation of a plane-averaged energy budget separating mean, periodic (phase-locked), and residual contributions.

The method is designed for the analysis of CFD simulations of stirred tanks and single-use bioreactors, but is otherwise geometry-agnostic.

---

## Contents

- `src/`
  - `flow_decomposition.m` — main routine
  - `buildSMask.m` — construction of spatial validity masks
  - `plotModes.m` — visualization utility
- `examples/`
  - `run_demo.m` — example driver script
- `data/`
  - *(not included in this repository; see below)*

---

## Method overview

The decomposition is based on the identification of coherent harmonics of the impeller shaft frequency from the power spectral density of a plane-averaged kinetic-energy proxy. Selected harmonics are reconstructed via linear regression onto sine/cosine bases, yielding:

- mean flow component,
- phase-locked (periodic) component,
- residual (stochastic) component.

A plane-averaged energy triad is then computed from these contributions.

The implementation does **not** rely on phase averaging or blade-angle binning and operates directly in the time domain.

---

## Data layout and conventions

All spatial fields must follow the same canonical layout:

- rows correspond to the **z-direction** (`Nz`)
- columns correspond to the **x-direction** (`Nx`)

Required array sizes:

- `X, Z` : `[Nz × Nx]` coordinate grids
- `U, V, W` : `[Nz × Nx × T]` velocity components
- `S` : `[Nz × Nx]` or `[Nz × Nx × T]` validity mask (`1 = valid`, `0 = excluded`)
- `t` : `[T × 1]` uniformly sampled time vector

With this convention:

- the x-axis vector is `X(1,:)`
- the z-axis vector is `Z(:,1)`

The example dataset provided on Zenodo already satisfies these conventions.

---

## Demo dataset

The demonstration dataset used by `run_demo.m` is **not stored in this GitHub repository** due to size limitations.

It is archived separately on Zenodo:

**Demo dataset DOI:**  
*(insert Zenodo dataset DOI here once published)*

After downloading the dataset, place the file in the local `data/` folder so that the directory structure is:

```
data/
  dataset_demo.mat
```

---

## Running the demo

1. Download the demo dataset from Zenodo and place it in `data/`
2. Add the `src/` folder to your MATLAB path
3. Run the example script:

```matlab
run examples/run_demo.m
```

The script:

- loads the demo dataset,
- constructs the spatial mask,
- performs the phase-locked decomposition,
- generates diagnostic figures.

---

## Requirements

- MATLAB (tested with recent releases)
- No additional toolboxes beyond standard MATLAB functionality

---

## License

The MATLAB code in this repository is released under the **MIT License** (see `LICENSE`).

The demo dataset archived on Zenodo is released under **Creative Commons Attribution 4.0 International (CC-BY-4.0)**.

---

## Citation

If you use this code or the accompanying dataset in academic work, please cite:

- the associated journal article (once published),
- the Zenodo archive corresponding to this code and/or dataset.

Citation details will be added here once the Zenodo DOI for the code release is available.

---

## Contact

For questions or comments, please contact the corresponding author of the associated publication.

