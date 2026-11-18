# Median observers: simulations for the ECC paper

All figures and numbers of the paper come from one script, `median_ECC.m`.
Run it from this folder.

## Requirements

- MATLAB R2020a or later (`tiledlayout`, `exportgraphics`)
- Control System Toolbox (`place`)
- YALMIP and MOSEK

## Files

| File | Content |
| --- | --- |
| `median_ECC.m` | observer design, simulations, figures; all functions are local functions at the end of the file |
| `figures/fig_observers.pdf` | Fig. 1: (a) Luenberger, stubborn and median observers under bursts of q = 2 outliers; (b) median observer when Definition 1 holds, when it is violated, and with Gaussian noise |
| `figures/fig_limits.pdf` | Fig. 2: (a) q_max certified by the LMIs of Theorem 2 versus T_s; (b) largest mu T_mem versus q for a scalar plant, against the necessary condition of Proposition 1 |
| `data/fig2_lmi.mat` | LMI results behind Fig. 2 (`qS`, `qU`, `aS`, `aU`), written when `recompute = true` |

## Sections of the script

| Section | Produces |
| --- | --- |
| Plant, observer design | gains L_i (Ts = 25 ms, q = 2, lambda = 1), kappa and the certified floor |
| Disturbances, simulations | Fig. 1 and the numbers quoted in the text (marked "text only" in the code) |
| Fig. 2(a), Fig. 2(b) | q_max and a_max, computed or loaded from `data/fig2_lmi.mat` |
| Fig. 1, Fig. 2 | the two PDFs in `figures/` |

## Fig. 2 data

The LMI sweeps take about 15-30 min with MOSEK.

- `recompute = true`: recompute and save `data/fig2_lmi.mat`.
- `recompute = false`: load `data/fig2_lmi.mat`. If the file is missing, the script
  uses values stored in the code (last run, Oct 2026) and prints a warning.

Run once with `recompute = true` so that Fig. 2 comes from saved results.

## History

Older exploratory scripts (Lur'e and switched co-design, toy example) are in the
first commit (`7cb2d00`). The branches `median-observer-worst-case-*` (MIMO
worst case) and `python-codebase-refactor-*` (Python port) were created by an
automated agent and are not used for the paper.
