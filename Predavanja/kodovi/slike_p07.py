"""Slike za Predavanje 07 (Od sive atmosfere do teorijskih LTR spektara).

Svi modeli su izračunati kodom model_atmosfere_p07.py (LTR, ravnoteža zračenja, kontinuum bez linija):
  slajd_p07_milne_T_LD.png    - siva atmosfera: T/T_eff i potamnjenje ka rubu za Š-Š, Eddingtonovo i tačno rešenje
  slajd_p07_lambda_iteracija.png - Lambda iteracija Eddingtonovog rešenja -> Hopfovo rešenje (q(tau) i greška fluksa)
  slajd_p07_korekcija.png     - Unsöld-Lucy korekcija za T_eff = 10000 K, log g = 4: T(tau_R) po iteracijama i greška fluksa
  slajd_p07_sunce.png         - izlazni fluks modela Sunca, siva atmosfera i Plankova funkcija; dubina tau_nu = 2/3
  slajd_p07_neprozracnost.png - doprinosi neprozračnosti na tau_R = 1 za T_eff = 5777 K i 10000 K
  slajd_p07_teff_niz.png      - izlazni fluks za niz T_eff (log g = 4)
  slajd_p07_logg.png          - uticaj log g na Balmerov skok (T_eff = 8000 K)
  slajd_p07_balmer_skok.png   - Balmerov skok D u funkciji T_eff za log g = 3, 4, 5

Upotreba:  python slike_p07.py [izlazni_direktorijum] [--samo-milne]
Mreža modela se računa jednom (oko 25 min na dva jezgra) i čuva u modeli_p07.pkl pored slika.
"""
import os
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_atmosfere_p07 as M

plt.rcParams.update({"mathtext.fontset": "cm", "font.size": 13,
                     "axes.spines.top": False, "axes.spines.right": False})
C_PL, C_MOD, C_SIV = "#1b4f9c", "#d17a00", "0.45"
NU = M.mreza_frekvencija(300)
NUF = M.mreza_frekvencija(1000)
LAMF = M.c / NUF * 1e7          # nm

TEFF = (4500, 5000, 5500, 6000, 6500, 7000, 8000, 9000, 10000, 11000, 12000)
LOGG = (3.0, 4.0, 5.0)
PARAM = [(5777, 4.44), (8000, 2.0)] + [(T, g) for g in LOGG for T in TEFF]


def izracunaj(p):
    Teff, lg = p
    m = M.Model(Teff, lg, nd=80, nu=NU)
    m.unsold_lucy(10, verbose=False)
    m.spektar(NUF)
    d = dict(Teff=Teff, logg=lg, tauR=m.tauR, T=m.T, T_sivo=m.T_sivo, P=m.P, es=m.es,
             F=m.F_fino, tau_f=m.tau_fino, ist=m.istorija, maxerr=float(np.max(abs(m.istorija[-1][1]))))
    print(f"  T_eff = {Teff:5d} K, log g = {lg:.2f}: max |dH/H| = {d['maxerr']:.4f}", flush=True)
    return p, d


def mreza(kes):
    if os.path.exists(kes):
        return pickle.load(open(kes, "rb"))
    with Pool(2) as pool:
        res = dict(pool.map(izracunaj, PARAM))
    pickle.dump(res, open(kes, "wb"))
    return res


def nm_osa(ax, lo=100, hi=5000, tick=(100, 200, 300, 500, 1000, 2000, 3000, 5000)):
    ax.set_xscale("log"); ax.set_xlim(lo, hi)
    t = [t for t in tick if lo <= t <= hi]
    ax.set_xticks(t); ax.set_xticklabels([str(x) for x in t]); ax.minorticks_off()
    ax.set_xlabel(r"talasna dužina $\lambda$ [nm]")


def ivice(ax, y=None):
    for l, n in ((M.c / (M.NU_RH / 4) * 1e7, "Balmer"), (M.c / (M.NU_RH / 9) * 1e7, "Pašen"), (1642, r"H$^-$")):
        ax.axvline(l, color="0.75", lw=0.9, ls="--", zorder=0)


# ---------------------------------------------------------------- slike
def korekcija(R, out, kljuc=(10000, 4.0)):
    d = R[kljuc]
    lt = np.log10(d["tauR"])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    cm = plt.get_cmap("Blues")
    sel = [0, 1, 2, 4, len(d["ist"]) - 1]
    for j, it in enumerate(sel):
        T, e = d["ist"][it]
        lab = "siva (početak)" if it == 0 else ("konačno" if it == sel[-1] else f"iteracija {it}")
        col = C_SIV if it == 0 else (C_MOD if it == sel[-1] else cm(0.45 + 0.15 * j))
        a1.plot(lt, T, color=col, lw=2.0 if it in (0, sel[-1]) else 1.3, label=lab)
        a2.semilogy(lt, np.maximum(abs(e), 1e-6), color=col, lw=2.0 if it in (0, sel[-1]) else 1.3)
    a1.set_xlim(-5, 2); a1.set_ylim(None, 2.0 * d["Teff"])
    a1.set_xlabel(r"$\log\,\tau_R$"); a1.set_ylabel(r"$T$ [K]"); a1.legend(frameon=False, fontsize=11)
    a2.set_xlim(-5, 2); a2.set_ylim(1e-5, 1)
    a2.set_xlabel(r"$\log\,\tau_R$"); a2.set_ylabel(r"$|\,H - H_0\,|\,/\,H_0$")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def sunce(R, out):
    d = R[(5777, 4.44)]
    lam = LAMF
    # siva atmosfera (Hopf): F_nu(0) = 2 pi Int B_nu(T(t)) E2(t) dt
    t = np.logspace(-6, 2.5, 400)
    Tg = 5777 * (0.75 * (t + M.q_hopf(t)))**0.25
    Bt = M.planck(NUF[None, :], Tg[:, None])
    Fg = 2 * np.pi * np.trapezoid(Bt * expn(2, t)[:, None], t, axis=0)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.6, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.4, 1]})
    a1.plot(lam, np.pi * M.planck(NUF, 5777.) * 1e5, color=C_PL, lw=1.8, label=r"$\pi B_\nu(T_{\rm eff})$")
    a1.plot(lam, Fg * 1e5, color=C_SIV, lw=1.6, ls="--", label="siva atmosfera (Hopf)")
    a1.plot(lam, d["F"] * 1e5, color=C_MOD, lw=1.8, label="ne-siva LTR (naš model)")
    a1.set_ylabel(r"$F_\nu(0)$ [$10^{-5}$ erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]", fontsize=11)
    a1.set_ylim(0, None); a1.legend(frameon=False, fontsize=11, loc="upper left"); ivice(a1)
    # dubina na kojoj je tau_nu = 2/3, i temperatura tamo
    lt = np.log10(d["tauR"])
    l23 = np.array([np.interp(np.log10(2 / 3), np.log10(d["tau_f"][:, j]), lt) for j in range(NUF.size)])
    a2.plot(lam, l23, color=C_MOD, lw=1.8)
    a2.set_ylabel(r"$\log\,\tau_R\,(\tau_\nu = 2/3)$", fontsize=12)
    a2.set_ylim(-2.5, 0.5); ivice(a2)
    ax2 = a2.twinx(); ax2.spines["right"].set_visible(True)
    ylo, yhi = a2.get_ylim()
    tk = np.array([5000, 6000, 7000, 8000])
    mono = d["tauR"] > 1e-3                                 # T(tau) monotono raste ispod tau ~ 1e-3
    pos = np.interp(tk, d["T"][mono], lt[mono])
    m = (pos > ylo) & (pos < yhi)
    ax2.set_ylim(ylo, yhi); ax2.set_yticks(pos[m]); ax2.set_yticklabels([str(x) for x in tk[m]])
    ax2.set_ylabel(r"$T$ [K]", fontsize=12)
    nm_osa(a2, 100, 5000)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def neprozracnost(R, out):
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=False)
    boje = {"H I b-f": "#1b4f9c", "H I f-f": "#6fa8dc", "H$^-$ b-f": "#d17a00", "H$^-$ f-f": "#f1b05a",
            "Tomson (e$^-$)": "0.3", "Rejli (H I)": "0.65"}
    for ax, key in zip(axs, [(5777, 4.44), (10000, 4.0)]):
        d = R[key]
        i = np.argmin(abs(d["tauR"] - 1))
        es = {q: np.array([v[i]]) for q, v in d["es"].items()}
        comp = M.neprozracnost(NUF, np.array([d["T"][i]]), es, komponente=True)
        rho = es["rho"][0]
        tot = sum(v[0] for v in comp.values())
        for n, v in comp.items():
            ax.loglog(LAMF, v[0] / rho, color=boje[n], lw=1.5, label=n)
        ax.loglog(LAMF, tot / rho, color="k", lw=2.0, label="ukupno")
        ax.set_title(rf"$T_{{\rm eff}}$ = {key[0]} K:  $\tau_R = 1$,  $T$ = {d['T'][i]:.0f} K", fontsize=12)
        nm_osa(ax, 100, 5000, tick=(100, 200, 500, 1000, 2000, 5000))
        ivice(ax)
    axs[0].set_ylim(1e-3, 1e2); axs[1].set_ylim(1e-3, 1e4)
    axs[0].set_ylabel(r"$\chi_\nu/\rho$ [cm$^2$ g$^{-1}$]")
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, frameon=False, fontsize=11, loc="lower center", ncol=7)
    fig.tight_layout(rect=(0, 0.09, 1, 1)); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def teff_niz(R, out):
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    cm = plt.get_cmap("coolwarm_r")
    Ts = (4500, 5500, 6500, 8000, 10000, 12000)
    for j, T in enumerate(Ts):
        F = R[(T, 4.0)]["F"]
        ax.plot(LAMF, F / F.max(), color=cm(j / (len(Ts) - 1)), lw=1.8, label=f"{T} K")
    nm_osa(ax, 100, 5000); ivice(ax)
    ax.set_ylabel(r"$F_\nu(0)\,/\,\max F_\nu$"); ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=11, title=r"$T_{\rm eff}$ (log g = 4)", title_fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def logg(R, out):
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    cm = plt.get_cmap("Blues")
    j550 = np.argmin(abs(LAMF - 550))
    for j, g in enumerate((2.0, 3.0, 4.0, 5.0)):
        F = R[(8000, g)]["F"]
        ax.plot(LAMF, F / F[j550], color=cm(0.4 + 0.2 * j), lw=1.9, label=f"log g = {g:.0f}")
    ax.set_xlim(250, 1000); ax.set_ylim(0, None)
    ax.set_xlabel(r"talasna dužina $\lambda$ [nm]"); ax.set_ylabel(r"$F_\nu(0)\,/\,F_\nu(550\,{\rm nm})$")
    ivice(ax); ax.legend(frameon=False, fontsize=11, title=r"$T_{\rm eff}$ = 8000 K", title_fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def balmer_skok(R, out):
    lam2 = M.c / (M.NU_RH / 4) * 1e7                   # 364.7 nm
    jm = np.searchsorted(LAMF[::-1], lam2)            # LAMF opada
    lamr = LAMF[::-1]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    cm = plt.get_cmap("Blues")
    for j, g in enumerate(LOGG):
        D = []
        for T in TEFF:
            F = R[(T, g)]["F"][::-1]
            D.append(np.log10(F[jm] / F[jm - 1]))       # dugotalasna / kratkotalasna strana ivice
        ax.plot(TEFF, D, "o-", color=cm(0.45 + 0.25 * j), lw=1.9, ms=5, label=f"log g = {g:.0f}")
    ax.set_xlabel(r"$T_{\rm eff}$ [K]"); ax.set_ylabel(r"Balmerov skok $D$")
    ax.legend(frameon=False, fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)
    return lamr[jm], lamr[jm - 1]


# ---------------------------------------------------------------- siva atmosfera: rešenja Milneovog problema
C_SS, C_ED, C_HO = "#d17a00", "#1b4f9c", "black"          # iste boje kao u predavanju 06


def lambda_iteracija_eddington(n_max=3000, sacuvaj=(0, 1, 2, 10, 100, 1000, 3000)):
    """B^(n+1) = Lambda[B^(n)], polazeći od Eddingtonovog rešenja B^(0) = 3H (tau + 2/3), H = 1.
    Atmosfera je polubeskonačna: ispod najdublje tačke S se linearno ekstrapolira (M.momenti)."""
    t = np.concatenate([[0.0], np.logspace(-5, 3.3, 220)])[:, None]
    S = 3 * (t + 2 / 3)
    q, err = {}, []
    for n in range(n_max + 1):
        J, H = M.momenti(t, S)
        err.append(np.max(abs(H[t[:, 0] < 20, 0] - 1)))
        if n in sacuvaj:
            q[n] = (S / 3 - t)[:, 0].copy()
        S = J
    return t[:, 0], q, np.array(err)


def I_povrsina(t, S, mu):
    """I(0, mu) = Int S exp(-t/mu) dt/mu (trapez na finoj mreži)."""
    if mu == 0:
        return S[0]
    tt = np.linspace(0, 40 * mu, 4000)
    return np.trapezoid(np.interp(tt, t, S) * np.exp(-tt / mu) / mu, tt)


def milne_T_LD(t, out):
    """T/T_eff za Š–Š, Eddingtonovo i tačno (Hopfovo) rešenje + potamnjenje ka rubu."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    lt = np.linspace(-3, 1, 300); tau = 10**lt
    a1.plot(lt, (tau + 0.5)**0.25, color=C_SS, lw=2.2, label="Š–Š")
    a1.plot(lt, (0.75 * (tau + 2 / 3))**0.25, color=C_ED, lw=2.2, label="Eddington")
    a1.plot(lt, (0.75 * (tau + M.q_hopf(tau)))**0.25, color=C_HO, lw=1.6, ls="--", label="tačno (Hopf)")
    a1.axhline(1, color="0.6", lw=0.8, zorder=0)
    a1.set_xlabel(r"$\log\tau$"); a1.set_ylabel(r"$T\,/\,T_{\rm eff}$"); a1.set_xlim(-3, 1)
    a1.legend(frameon=False, loc="upper left")
    mu = np.linspace(0, 1, 60)
    a2.plot(mu, 2 / 3 * (mu + 0.5), color=C_SS, lw=2.2, label="Š–Š")
    a2.plot(mu, 3 / 5 * (mu + 2 / 3), color=C_ED, lw=2.2, label="Eddington")
    SH = 3 * (t + M.q_hopf(t))
    I1 = I_povrsina(t, SH, 1.0)
    a2.plot(mu, [I_povrsina(t, SH, m) / I1 for m in mu], color=C_HO, lw=1.6, ls="--", label="tačno (Hopf)")
    # posmatranja: vrednosti iz tabele sa slajda (predavanje 2025)
    a2.plot([1, 0.8, 0.6, 0.44, 0.2, 0], [1, 0.92, 0.81, 0.70, 0.49, 0.40], "o", ms=7,
            mfc="#b2182b", mec="white", label="Sunce (posmatranje)")
    a2.set_xlabel(r"$\mu$"); a2.set_ylabel(r"$I(0,\mu)\,/\,I(0,1)$")
    a2.set_xlim(1.02, -0.02); a2.set_ylim(0.25, 1.05); a2.legend(frameon=False, loc="lower left", fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def milne_lambda(t, q, err, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    cm = plt.get_cmap("Blues")
    ns = sorted(q)
    for j, n in enumerate(ns):
        lab = "Eddington ($n$ = 0)" if n == 0 else f"$n$ = {n}"
        col = C_ED if n == 0 else (C_HO if n == ns[-1] else cm(0.35 + 0.6 * j / len(ns)))
        a1.plot(np.log10(np.maximum(t, 1e-5)), q[n], color=col, lw=2.0 if n in (0, ns[-1]) else 1.4, label=lab)
    lt_ = np.linspace(-5, 2, 300)
    a1.plot(lt_, M.q_hopf(10**lt_), color="#b2182b", lw=1.6, ls="--", label="Hopf (tačno)")
    a1.axhline(1 / np.sqrt(3), color="0.6", lw=0.8, ls=":"); a1.axhline(0.7104, color="0.6", lw=0.8, ls=":")
    a1.text(-4.9, 1 / np.sqrt(3) - 0.006, r"$q(0) = 1/\sqrt{3}$", fontsize=10, color="0.4", va="top")
    a1.text(1.9, 0.7104 + 0.003, r"$q(\infty) = 0.7104$", fontsize=10, color="0.4", ha="right")
    a1.set_xlim(-5, 2); a1.set_ylim(0.555, 0.725)
    a1.set_xlabel(r"$\log\tau$"); a1.set_ylabel(r"$q^{(n)}(\tau) = S^{(n)}/3H - \tau$")
    h_, l_ = a1.get_legend_handles_labels()
    n = np.arange(err.size)
    a2.loglog(n[1:], err[1:], color=C_HO, lw=2.0)
    a2.set_xlabel(r"broj iteracija $n$"); a2.set_ylabel(r"max $|H^{(n)}/H - 1|$  ($\tau < 20$)")
    a2.legend(h_, l_, frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(d, exist_ok=True)
    t, q, err = lambda_iteracija_eddington()                 # ~1.5 min
    milne_T_LD(t, os.path.join(d, "slajd_p07_milne_T_LD.png"))
    milne_lambda(t, q, err, os.path.join(d, "slajd_p07_lambda_iteracija.png"))
    if "--samo-milne" in sys.argv:
        sys.exit(0)
    R = mreza(os.path.join(d, "modeli_p07.pkl"))
    korekcija(R, os.path.join(d, "slajd_p07_korekcija.png"))
    sunce(R, os.path.join(d, "slajd_p07_sunce.png"))
    neprozracnost(R, os.path.join(d, "slajd_p07_neprozracnost.png"))
    teff_niz(R, os.path.join(d, "slajd_p07_teff_niz.png"))
    logg(R, os.path.join(d, "slajd_p07_logg.png"))
    print("Balmerov skok meren na", balmer_skok(R, os.path.join(d, "slajd_p07_balmer_skok.png")), "nm")
