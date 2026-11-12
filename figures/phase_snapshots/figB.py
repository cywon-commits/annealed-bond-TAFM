import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
exec(open(_pl.Path(__file__).with_name("plot.py")).read())
cases=[("win_mu0_T0.49_random.npz",r"(a) $\mu=0$, $T/T_c=0.93$"),("win_mu10_T0.8_random.npz",r"(b) $\mu=10$, $T/T_c=0.70$"),
       ("win_mu10_T1.06_random.npz",r"(c) $\mu=10$, $T/T_c=0.93$"),("win_mu10_T1.2_random.npz",r"(d) $\mu=10$, $T/T_c=1.05$")]
fig=plt.figure(figsize=(7.0,2.35)); W0=0.215
for k,(f,tt) in enumerate(cases):
    z=np.load(SNAP/f); ax=fig.add_axes([0.003+k*(W0+0.004),0.03,W0,0.74])
    d,fr=panel(ax,z["s"],z["mate"],"",RW=26,RH=26,amp=0.7,dsz=3.5,fsz=1.0)
    ax.set_title(tt+"\n"+rf"$\langle|\psi_s|^2\rangle={z['P'].mean():.3f}$, frustr. {fr:.3f}",fontsize=7.4)
    print(f,z["P"].mean(),d,fr)
axw=fig.add_axes([0.885,0.30,0.11,0.36]); wheel(axw,ticks=False)
for k in range(6):
    a=k*np.pi/3; ax_lw=1.4 if k%2==0 else 0.5
    axw.plot([0.7*np.cos(a),np.cos(a)],[0.7*np.sin(a),np.sin(a)],color="k",lw=ax_lw)
axw.set_title(r"arg $\psi_{\rm loc}$",fontsize=7.2,pad=2)
fig.text(0.94,0.27,"thick ticks:\nthe three\n$Z_3$ states",fontsize=6.3,ha="center",va="top")
for e in ("pdf","png"): fig.savefig(FIGS/f"fig_window_snap.{e}",dpi=220 if e=="png" else None)
