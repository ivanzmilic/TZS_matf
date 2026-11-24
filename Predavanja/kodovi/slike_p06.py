"""Slike za Predavanje 06 (Modelovanje atmosfera u ravnoteži zračenja).

  slajd_p06_izotermna.png      - izotermna atmosfera: S = 1 i J = Lambda[1] = 1 - E2(tau)/2
  slajd_p06_lambda.png         - transformacija S -> J = Lambda[S] za dve zadate funkcije izvora
  slajd_p06_provera.png        - provera rešenja: Lambda[S]/S i F(tau)/F za Š-Š, Eddingtonovo i Hopfovo rešenje
  slajd_p06_potamnjenje.png    - potamnjenje ka rubu: I(0,mu)/I(0,1) za tri rešenja + posmatranja (stara tabela)
  slajd_p06_hopf.png           - Hopfova funkcija q(tau) i T(tau)/T_eff za tri rešenja
  slajd_p06_spektar_sive.png   - izlazni fluks sive atmosfere i Plankova funkcija na T_eff (gore) i njihov odnos (dole)

Jedinice: H = F/(4 pi) = 1, pa je S_SS = 4(tau + 1/2), S_Edd = 3(tau + 2/3), S_Hopf = 3(tau + q(tau)).
Hopfova funkcija je data aproksimacijom q(tau) = 0.7104 - 0.1331 exp(-3.4488 tau)  (q(0) = 0.5773, q(inf) = 0.7104).

Upotreba: python slike_p06.py [izlazni_direktorijum]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import expn

plt.rcParams.update({"mathtext.fontset": "cm", "font.size": 13,
                     "axes.spines.top": False, "axes.spines.right": False})
C_SS, C_ED, C_HO, C_S, C_J = "#d17a00", "#1b4f9c", "black", "#1b4f9c", "#d17a00"
TMAX = 45.0   # gornja granica integracije po optičkoj dubini (umesto beskonačnosti)


def q_hopf(t):
    return 0.7104 - 0.1331 * np.exp(-3.4488 * t)


S_SS = lambda t: 4.0 * (t + 0.5)
S_ED = lambda t: 3.0 * (t + 2.0 / 3.0)
S_HO = lambda t: 3.0 * (t + q_hopf(t))


def Lambda(S, tau):
    """J(tau) = 1/2 Int_0^inf S(t) E1(|t - tau|) dt"""
    up = quad(lambda t: S(t) * expn(1, t - tau), tau, tau + TMAX, limit=200)[0]
    dn = quad(lambda t: S(t) * expn(1, tau - t), 0.0, tau, limit=200)[0] if tau > 0 else 0.0
    return 0.5 * (up + dn)


def Hflux(S, tau):
    """H(tau) = 1/2 [Int_tau^inf S E2(t - tau) dt - Int_0^tau S E2(tau - t) dt]"""
    up = quad(lambda t: S(t) * expn(2, t - tau), tau, tau + TMAX, limit=200)[0]
    dn = quad(lambda t: S(t) * expn(2, tau - t), 0.0, tau, limit=200)[0] if tau > 0 else 0.0
    return 0.5 * (up - dn)


def I_out(S, mu):
    """I(0, mu) = Int_0^inf S(t) exp(-t/mu) dt/mu"""
    if mu == 0:
        return S(0.0)
    return quad(lambda t: S(t) * np.exp(-t / mu) / mu, 0.0, TMAX * mu, limit=200)[0]


def izotermna(out):
    lt = np.linspace(-3, 3, 300)
    tau = 10**lt
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.plot(lt, np.ones_like(lt), color=C_S, lw=2.4, label=r"$S = B$")
    ax.plot(lt, 1 - 0.5 * expn(2, tau), color=C_J, lw=2.4, label=r"$J = \Lambda[S]$")
    ax.set_xlabel(r"$\log\tau$"); ax.set_ylim(0, 1.08); ax.set_xlim(-3, 3)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def lambda_transformacija(out):
    lt = np.linspace(-3, 2.5, 70)
    tau = 10**lt
    cl = lambda t: np.log10(np.maximum(t, 1e-3))        # S = const za tau < 10^-3
    funs = [(lambda t: 0.7 + 0.1 * cl(t), r"$S$ linearno po $\log\tau$"),
            (lambda t: 0.65 + 0.05 * (cl(t) + 1.0)**2, r"$S$ kvadratno po $\log\tau$")]
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.4))
    for ax, (S, ttl) in zip(axs, funs):
        ax.plot(lt, S(tau), color=C_S, lw=2.4, label=r"$S$")
        ax.plot(lt, [Lambda(S, t) for t in tau], color=C_J, lw=2.4, label=r"$J = \Lambda[S]$")
        ax.set_title(ttl, fontsize=13); ax.set_xlabel(r"$\log\tau$"); ax.set_ylim(0, 1.5); ax.set_xlim(-3, 2.5)
        ax.legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def provera(out):
    tau = np.linspace(0, 4, 60)
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.4))
    for S, c, lab in [(S_SS, C_SS, "Š–Š"), (S_ED, C_ED, "Eddington"), (S_HO, C_HO, "Hopf")]:
        ls = "--" if S is S_HO else "-"
        axs[0].plot(tau, [Lambda(S, t) / S(t) for t in tau], color=c, lw=2.2, ls=ls, label=lab)
        axs[1].plot(tau, [Hflux(S, t) for t in tau], color=c, lw=2.2, ls=ls, label=lab)
    axs[0].set_ylabel(r"$\Lambda_\tau[S]\,/\,S$"); axs[1].set_ylabel(r"$F(\tau)\,/\,F$")
    axs[0].set_ylim(0.85, 1.05); axs[1].set_ylim(0.9, 1.4)
    for ax in axs:
        ax.axhline(1, color="0.6", lw=0.8, zorder=0); ax.set_xlabel(r"$\tau$"); ax.set_xlim(0, 4)
        ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def potamnjenje(out):
    mu = np.linspace(0, 1, 60)
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    ax.plot(mu, 2 / 3 * (mu + 0.5), color=C_SS, lw=2.2, label="Š–Š")
    ax.plot(mu, 3 / 5 * (mu + 2 / 3), color=C_ED, lw=2.2, label="Eddington")
    I1 = I_out(S_HO, 1.0)
    ax.plot(mu, [I_out(S_HO, m) / I1 for m in mu], color=C_HO, lw=1.6, ls="--", label="Hopf")
    # posmatranja: vrednosti iz tabele sa slajda (predavanje 2025)
    ax.plot([1, 0.8, 0.6, 0.44, 0.2, 0], [1, 0.92, 0.81, 0.70, 0.49, 0.40], "o", ms=7,
            mfc="#b2182b", mec="white", label="posmatranje")
    ax.set_xlabel(r"$\mu$"); ax.set_ylabel(r"$I(0,\mu)\,/\,I(0,1)$")
    ax.set_xlim(1.02, -0.02); ax.set_ylim(0.25, 1.05); ax.legend(frameon=False, loc="lower left")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def hopf(out):
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.4))
    t = np.linspace(0, 3, 300)
    axs[0].plot(t, q_hopf(t), color=C_HO, lw=2.2, label=r"$q(\tau)$ (Hopf)")
    axs[0].axhline(2 / 3, color=C_ED, lw=2.0, label=r"Eddington: $q = 2/3$")
    axs[0].set_xlabel(r"$\tau$"); axs[0].set_ylabel(r"$q(\tau)$"); axs[0].set_ylim(0.55, 0.73); axs[0].set_xlim(0, 3)
    axs[0].legend(frameon=False, loc="lower right")
    lt = np.linspace(-3, 1, 300); tau = 10**lt
    axs[1].plot(lt, (tau + 0.5)**0.25, color=C_SS, lw=2.2, label="Š–Š")
    axs[1].plot(lt, (0.75 * (tau + 2 / 3))**0.25, color=C_ED, lw=2.2, label="Eddington")
    axs[1].plot(lt, (0.75 * (tau + q_hopf(tau)))**0.25, color=C_HO, lw=1.6, ls="--", label="Hopf")
    axs[1].axhline(1, color="0.6", lw=0.8, zorder=0)
    axs[1].set_xlabel(r"$\log\tau$"); axs[1].set_ylabel(r"$T\,/\,T_{\rm eff}$"); axs[1].set_xlim(-3, 1)
    axs[1].legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)


def spektar_sive(out):
    x = np.linspace(0.05, 12, 160)                       # x = h nu / k T_eff
    theta = lambda t: (0.75 * (t + q_hopf(t)))**0.25     # T / T_eff
    planck = lambda xx, th: xx**3 / np.expm1(xx / th)
    # F_nu(0)/pi = 2 Int_0^inf B_nu(T(t)) E2(t) dt
    Fg = np.array([2 * quad(lambda t: planck(xx, theta(t)) * expn(2, t), 0, TMAX, limit=200)[0] for xx in x])
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(5.6, 5.0), sharex=True, gridspec_kw={"height_ratios": [1.5, 1]})
    ax.plot(x, planck(x, 1.0), color=C_S, lw=2.2, label=r"$\pi B_\nu(T_{\rm eff})$")
    ax.plot(x, Fg, color=C_J, lw=2.2, label=r"$F_\nu(0)$, siva atmosfera")
    ax.set_ylabel(r"$F_\nu(0)$  [proizv. jed.]"); ax.set_ylim(0, None); ax.legend(frameon=False)
    # odnos: blizu maksimuma ~4 % manjak, u UV rastući višak
    ar.axhline(1.0, color="0.6", lw=1)
    ar.plot(x, Fg / planck(x, 1.0), color=C_J, lw=2.2)
    ar.set_xlabel(r"$h\nu\,/\,kT_{\rm eff}$"); ar.set_ylabel(r"$F_\nu(0)\,/\,\pi B_\nu(T_{\rm eff})$")
    ar.set_xlim(0, 12); ar.set_ylim(0.9, 2.0)
    fig.tight_layout(); fig.align_ylabels([ax, ar])
    fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)
    trap = lambda y: float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))   # radi i sa starim i sa novim numpy
    return trap(Fg) / trap(planck(x, 1.0))


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(d, exist_ok=True)
    izotermna(os.path.join(d, "slajd_p06_izotermna.png"))
    lambda_transformacija(os.path.join(d, "slajd_p06_lambda.png"))
    provera(os.path.join(d, "slajd_p06_provera.png"))
    potamnjenje(os.path.join(d, "slajd_p06_potamnjenje.png"))
    hopf(os.path.join(d, "slajd_p06_hopf.png"))
    r = spektar_sive(os.path.join(d, "slajd_p06_spektar_sive.png"))
    print("odnos integrala fluksa (siva / Plank na T_eff) na mreži x = 0.05-12: %.3f" % r)
