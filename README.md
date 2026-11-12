# annealed-bond-TAFM

Code and data for

> C. Won, *Annealed bonds double the cost of free vortices in the triangular Ising antiferromagnet* (in preparation).

## The model

Ising spins σ ∈ {0, 1} live on the sites of a triangular lattice. Each elementary triangle may pair with one neighbor into a lozenge, which switches off the bond the two triangles share. The pairings form a matching of the honeycomb lattice dual to the triangles, and they are annealed together with the spins. The Hamiltonian is

```
H = Σ_{active bonds} δ(σ_i, σ_j) + D n_1 + μ n_mon
```

- n_1 is the number of sites in state 1.
- n_mon is the number of unpaired ("frustrated") triangles.
- In the language of the TAFM, J = 1/2 and h/J = D.

The reference model is the fixed-bond triangular Ising antiferromagnet (TAFM) on the same lattice.

## Layout

| Path | Contents |
|---|---|
| `src/` | Simulation code: lattice and matching geometry, Metropolis sweeps (spin flips, dimer toggles and slides, hexagon flips), the fixed-bond TAFM, finite-size analysis, susceptibility, and the exact minimum-energy matching (`emin.py`). |
| `analysis/` | Analysis that turns raw runs into the summaries in `data/`. Subfolders: `ising_check`, `small_field_mc`, `susceptibility`, `defect_costs`, and `rg_line` (the small-field transition line predicted from the TAFM critical surface). |
| `figures/` | Scripts that draw every figure of the paper from `data/`. |
| `data/` | Processed data used by the figures and tables. |
| `results/susceptibility/raw/` | Raw outputs of the susceptibility runs. They are small enough to include. |
| `matplotlibrc` | Embeds TrueType (Type 42) fonts in the PDF figures. |

## Quick start

```bash
pip install -r requirements.txt
cd annealed-bond-tafm                 # run everything from the repository root
python3 figures/make_figs.py          # Figs. 1, 4, 6, 7, 8, 12, 13
python3 figures/make_chi_figs.py      # Figs. 9, 10
python3 figures/fig_rgline.py         # Fig. 11
python3 figures/phase_snapshots/figA.py   # Fig. 5
python3 figures/phase_snapshots/figB.py   # Fig. 14
python3 figures/make_snaps.py         # Figs. 2, 3 (runs short simulations, about 5 min)
```

Figures are written to `output/figs/`; set `OUT_DIR` to change this. Run from the root directory so that the bundled `matplotlibrc` is used. With the versions below, the scripts reproduce the published figures pixel for pixel.

## Figures and their sources

| Paper | File | Script | Data |
|---|---|---|---|
| Fig. 1 | `fig1_mechanism` | `figures/make_figs.py` | — (schematic) |
| Fig. 2 | `fig2_groundstates` | `figures/make_snaps.py` | generated |
| Fig. 3 | `fig3_snapshots` | `figures/make_snaps.py` | generated |
| Fig. 4 | `fig4_xi_defects` | `figures/make_figs.py` | `figdata/fig3a_xi.csv`, `fig3b_defects.csv`, `defect_classes_L32.csv` |
| Fig. 5 | `fig_vortex` | `figures/phase_snapshots/figA.py` | `snapshots/ann_T050.npz`, `tafm_T040.npz` |
| Fig. 6 | `fig5_binder` | `figures/make_figs.py` | `figdata/fig2_binder.csv` |
| Fig. 7 | `fig7_tc_mu` | `figures/make_figs.py` | `figdata/fig5_tc_mu.csv` |
| Fig. 8 | `fig6_heat` | `figures/make_figs.py` | `figdata/fig4_heat.csv` |
| Fig. 9 | `fig10_chi` | `figures/make_chi_figs.py` | `susceptibility/summary.txt` |
| Fig. 10 | `fig11_mH` | `figures/make_chi_figs.py` | `susceptibility/summary.txt` |
| Fig. 11 | `fig12_rgline` | `figures/fig_rgline.py` | `figdata/defect_classes_L32.csv` and TAFM critical points from the literature |
| Fig. 12 | `fig8_eta_cos3` | `figures/make_figs.py` | `figdata/fig6a_eta.csv`, `fig6b_cos3.csv` |
| Fig. 13 | `fig9_window` | `figures/make_figs.py` | `figdata/fig2_binder.csv` |
| Fig. 14 | `fig_window_snap` | `figures/phase_snapshots/figB.py` | `snapshots/win_*.npz` |

## Reproducing the data

Every script states its usage in its docstring. The main entry points are listed below.

- **Production Monte Carlo** (Potts code with a frozen third state, D3 = 10): `src/annealed_long.py`, `src/annealed_fss.py`, `src/annealed_tafm.py`, `src/annealed_tafm_xi.py`, `src/annealed_defects.py`.
  - `ISING=1` selects the strictly two-state sweep (`src/ising_sweep.py`); Appendix C of the paper compares the two.
  - Summaries of these runs are in `data/production_summaries/` and `data/figdata/`.
- **Ising validation and entropy:** `analysis/ising_check/`, with output in `data/ising_check/`.
- **Small-field Monte Carlo:** `analysis/small_field_mc/tc_small_d2.py`, with output in `data/small_field_mc/`.
- **Susceptibility and reduced-field scaling:** `src/chi.py` writes the runs, and `analysis/susceptibility/analyze.py` analyzes them. Output is in `results/susceptibility/` and `data/susceptibility/`.
- **Exact defect costs (Table III):**
  - `analysis/defect_costs/testA.py 6` and `testMu.py` compute the costs by maximum-weight matching.
  - `runC.py` measures the thermal densities of isolated and adjacent defects; `analyzeC.py` analyzes them.
- **Phase snapshots:** `figures/phase_snapshots/runs.py` produces the configurations, and `vort.py`, `adjcount.py`, `wn_check.py` count vortices and defects.

The raw Monte Carlo time series of the production runs (`.npz` checkpoints) are too large to include. They can be regenerated with the scripts above, and the development history is in [cywon-commits/RhombileLattice](https://github.com/cywon-commits/RhombileLattice), branch `partIX-ising-validation`.

## Requirements

Python ≥ 3.10 with numpy, numba, scipy, matplotlib, pandas and networkx. The versions used for the paper are pinned in `requirements.txt`.

## License

Code: MIT (see `LICENSE`). Data in `data/` and `results/`: CC BY 4.0.

## Citation

See `CITATION.cff`. A Zenodo DOI will be added with the first release.
