"""Seed-combined psi_s^2(L) and effective eta from weighted fits psi_s^2 ~ L^-eta (L >= LMIN),
plus seed-combined <cos 3theta>, <cos 6theta>; input data/window/raw_observables_d21.csv."""
import os, numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import DATA
LMIN = int(os.environ.get("LMIN", 32))
TC = {1: 0.696, 2: 0.800, 3: 0.950, 5: 1.050, 10: 1.141}
d = pd.read_csv(DATA / "window" / "raw_observables_d21.csv")
def combine(v, e):
    w = 1 / e**2; m = np.sum(w*v)/np.sum(w); err = np.sqrt(1/np.sum(w))
    if len(v) > 1:
        chi2 = np.sum(w*(v-m)**2)/(len(v)-1); err *= max(1.0, np.sqrt(chi2))
    return m, err
rows, fits = [], []
for mu, Tc in TC.items():
    for T, g in d[(d.mu == mu) & d.spin_psi2.notna()].groupby("T"):
        for L, h in g.groupby("L"):
            r = dict(mu=mu, T=T, TTc=T/Tc, L=L, nruns=len(h))
            for k in ("spin_psi2", "spin_U4", "c3", "c6"):
                r[k], r[k+"_e"] = combine(h[k].values, h[k+"_e"].values)
            rows.append(r)
P = pd.DataFrame(rows)
for (mu, T), g in P[P.L >= LMIN].groupby(["mu", "T"]):
    if len(g) < 3: continue
    x = np.log(g.L.values); y = np.log(g.spin_psi2.values); s = g.spin_psi2_e.values / g.spin_psi2.values
    W = np.diag(1/s**2); A = np.c_[np.ones_like(x), -x]
    cov = np.linalg.inv(A.T @ W @ A); p = cov @ A.T @ W @ y
    chi2 = float(((y - A @ p)**2 / s**2).sum()) / (len(x) - 2)
    fits.append(dict(mu=mu, T=T, TTc=T/TC[mu], eta=p[1], eta_e=np.sqrt(cov[1, 1]) * max(1, np.sqrt(chi2)),
                     chi2dof=chi2, nL=len(g), Lmax=int(g.L.max()), lnA=p[0]))
F = pd.DataFrame(fits)
out = DATA / "window"; P.to_csv(out / "psi2_by_L.csv", index=False, float_format="%.6g"); F.to_csv(out / "eta_fits.csv", index=False, float_format="%.5g")
pd.set_option("display.width", 200); print(F.round(3).to_string())
