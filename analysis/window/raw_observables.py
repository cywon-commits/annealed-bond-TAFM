"""Per-run observables from the long D=1 runs (results/long/long_L*_d21*.npz):
spin and hub order parameters (psi^2, Binder U4), <cos 3theta>, <cos 6theta> of the spin order parameter.
Discards the first DISCARD sweeps; errors from an NB-block jackknife.
Usage: python3 analysis/window/raw_observables.py <dir with npz> [out.csv]"""
import sys, os, re, glob
import numpy as np
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import DATA
DISCARD = int(os.environ.get("DISCARD", 100000)); NB = int(os.environ.get("NB", 8))
def jk(x, f):
    n = len(x) // NB; bl = [np.r_[x[:i*n], x[(i+1)*n:NB*n]] for i in range(NB)]
    v = np.array([f(b) for b in bl]); full = f(x[:NB*n])
    return full, np.sqrt((NB-1) * np.mean((v - v.mean())**2))
def U4(p2): return 1 - np.mean(p2**2) / (2 * np.mean(p2)**2)
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "long_L*_d21*.npz"))):
    z = np.load(f); thin = int(z["thin"]); k0 = DISCARD // thin
    if z["E"].shape[0] - k0 < 4 * NB: continue
    m = re.search(r"_s(\d+)(_(\w+))?\.npz$", f); kind = m.group(3) or "potts"
    rec = dict(file=os.path.basename(f), L=int(z["L"]), mu=float(z["mu"]), T=float(z["T"]), seed=int(m.group(1)), kind=kind,
               nsamp=z["E"].shape[0] - k0)
    for name, (a, b) in (("hub", ("re", "im")), ("spin", ("sre", "sim"))):
        if a not in z.files:
            continue
        x = z[a][k0:, 0].astype(float); y = z[b][k0:, 0].astype(float); p2 = x*x + y*y
        rec[f"{name}_psi2"], rec[f"{name}_psi2_e"] = jk(p2, np.mean)
        rec[f"{name}_U4"], rec[f"{name}_U4_e"] = jk(p2, U4)
        if name == "spin":
            th = np.arctan2(y, x)
            rec["c3"], rec["c3_e"] = jk(np.cos(3*th), np.mean)
            rec["c6"], rec["c6_e"] = jk(np.cos(6*th), np.mean)
    rows.append(rec)
import pandas as pd
df = pd.DataFrame(rows); out = sys.argv[2] if len(sys.argv) > 2 else DATA / "window" / "raw_observables_d21.csv"
_pl.Path(out).parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, index=False, float_format="%.6g")
print(len(df), "runs ->", out)
