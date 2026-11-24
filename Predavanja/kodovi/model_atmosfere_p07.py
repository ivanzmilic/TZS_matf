"""Jednostavni LTR modeli zvezdanih atmosfera u ravnoteži zračenja (Predavanje 07).

Fizika (cgs jedinice):
  - neprozračnost: H I vezano-slobodni (n = 1..8) i slobodno-slobodni prelazi, H- vezano-slobodni (John 1988)
    i slobodno-slobodni prelazi (Bell & Berrington 1987, fit iz Gray 2005, jedn. 8.13), Tomsonovo rasejanje na
    elektronima, Rejlijevo rasejanje na H I; faktor stimulisane emisije je uključen za prave apsorpcije;
  - jednačina stanja: H, H-, He i nekoliko metala kao donora elektrona (Na, Mg, Al, Si, Ca, Fe), Saha + Bolcman,
    konstantne particione funkcije;
  - hidrostatička ravnoteža na Roselandovoj skali: dp/dtau_R = g rho / chi_R (pritisak zračenja zanemaren);
  - formalno rešenje JPZ preko eksponencijalnih integrala (J = Lambda[S], H = Phi[S]/4, S linearna po delovima);
    koherentno rasejanje: S_nu = (kappa_nu B_nu + sigma_nu J_nu) / chi_nu, rešeno ubrzanom Lambda iteracijom (ALI);
  - početna struktura: siva atmosfera, T^4 = 3/4 Teff^4 [tau_R + q(tau_R)];
  - temperaturska korekcija: Unsöld-Lucy postupak.

Ograničenja: nema spektralnih linija, metali nisu uključeni u neprozračnost (samo kao donori elektrona), nema
konvekcije, pritiska zračenja ni ne-LTR efekata. Modeli služe za ilustraciju kontinuuma (Balmerov skok, Pašenov
nagib, H- minimum na 1.64 um). Provera: za Hopfovu S daje J = S i H = const (1-2 %, tačnost fita za q).

Upotreba:  import model_atmosfere_p07 as M;  m = M.Model(5777, 4.44, nd=80, nu=M.mreza_frekvencija(300)).unsold_lucy(10)
           (slike za predavanje: slike_p07.py; samostalno pokretanje računa model Sunca)
"""
import os
import sys
import numpy as np

# ---------------------------------------------------------------- konstante (cgs)
h, c, k, me = 6.62607e-27, 2.99792e10, 1.380649e-16, 9.10938e-28
NU_RH = 3.28805e15                     # Ridbergova frekvencija za vodonik [Hz]: lambda_1 = 91.18 nm, lambda_2 = 364.7 nm
sigSB, eV, mH = 5.67037e-5, 1.602177e-12, 1.6735e-24
sigT = 6.6524e-25

# hemijski sastav (brojne zastupljenosti u odnosu na H), energije jonizacije [eV], U_II/U_I
METALI = {  # A, E_ion, U_I, U_II
    "Na": (1.7e-6, 5.139, 2.0, 1.0), "Mg": (4.0e-5, 7.646, 1.0, 2.0), "Al": (2.8e-6, 5.986, 6.0, 1.0),
    "Si": (3.2e-5, 8.152, 9.0, 6.0), "Ca": (2.2e-6, 6.113, 1.0, 2.0), "Fe": (3.2e-5, 7.902, 25.0, 30.0)}
A_HE = 0.085
MASA_PO_H = mH * (1 + 4 * A_HE + sum(v[0] * 30 for v in METALI.values()))   # masa po atomu H (gruba)


def saha_fak(T):
    """(2 pi m_e k T / h^2)^{3/2}"""
    return (2 * np.pi * me * k * T / h**2)**1.5


# ---------------------------------------------------------------- jednačina stanja
def jednacina_stanja(P, T):
    """Za dati gasni pritisak P i temperaturu T vraća n_H (ukupno), n_H0, n_H+, n_H-, n_e, rho."""
    s = saha_fak(T)
    phiH = 2 * 0.5 * s * np.exp(-13.598 * eV / (k * T))      # n_H+ n_e / n_H0  (U_H0 = 2, U_H+ = 1)
    phiHm = 4.0 * s * np.exp(-0.754 * eV / (k * T))            # n_H0 n_e / n_H-  (2 U_H0 / U_H- = 4)
    phiM = {el: 2 * v[3] / v[2] * s * np.exp(-v[1] * eV / (k * T)) for el, v in METALI.items()}

    def elektroni(ne):
        # broj elektrona po atomu H i broj čestica po atomu H, za dato n_e
        xH = phiH / (phiH + ne)
        e_po_H = xH + sum(v[0] * phiM[el] / (phiM[el] + ne) for el, v in METALI.items())
        cest_po_H = 1 + A_HE + sum(v[0] for v in METALI.values())
        return e_po_H, cest_po_H

    # rešavamo  P/kT = n_H * cest_po_H + n_e,  n_e = n_H * e_po_H   bisekcijom po log n_e
    Ntot = P / (k * T)
    lo, hi = np.full_like(T, -10.0), np.log10(Ntot)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        ne = 10**mid
        e_po_H, cest = elektroni(ne)
        nH = (Ntot - ne) / cest
        f = nH * e_po_H - ne          # > 0: treba više elektrona
        lo = np.where(f > 0, mid, lo)
        hi = np.where(f > 0, hi, mid)
    ne = 10**(0.5 * (lo + hi))
    e_po_H, cest = elektroni(ne)
    nH = (Ntot - ne) / cest
    nHp = nH * phiH / (phiH + ne)
    nH0 = nH - nHp
    nHm = nH0 * ne / phiHm
    rho = nH * MASA_PO_H
    return dict(nH=nH, nH0=nH0, nHp=nHp, nHm=nHm, ne=ne, rho=rho)


# ---------------------------------------------------------------- neprozračnost
def sigma_Hm_bf(lam_um):
    """H- vezano-slobodni presek [cm^2], John (1988), lam u mikronima."""
    lam0 = 1.6419
    C = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]
    out = np.zeros_like(lam_um)
    m = lam_um < lam0
    x = 1 / lam_um[m] - 1 / lam0
    f = sum(Cn * x**(n / 2) for n, Cn in enumerate(C))
    out[m] = 1e-18 * lam_um[m]**3 * x**1.5 * f
    return np.clip(out, 0, None)


def alpha_Hm_ff(lamA, T):
    """H- slobodno-slobodna neprozračnost po atomu H0 i po jedinici P_e [cm^2 / (dyn cm^-2)]; Gray (2005) jedn. 8.13."""
    lamc = np.clip(lamA, 1823.0, 151890.0)
    L = np.log10(lamc)
    f0 = -2.2763 - 1.6850 * L + 0.76661 * L**2 - 0.053346 * L**3
    f1 = 15.2827 - 9.2846 * L + 1.99381 * L**2 - 0.142631 * L**3
    f2 = -197.789 + 190.266 * L - 67.9775 * L**2 + 10.6913 * L**3 - 0.625151 * L**4
    th = np.log10(5040.0 / T)
    a = 1e-26 * 10**(f0 + f1 * th + f2 * th**2)
    return a * (lamA / lamc)**2          # van opsega fita: ~ lambda^2


def neprozracnost(nu, T, es, komponente=False):
    """Vraća (kappa_nu, sigma_nu) [cm^-1] za mrežu frekvencija nu (1D) i slojeve (T, es) -> oblik (ND, NF)."""
    T = T[:, None]
    lamA = c / nu[None, :] * 1e8
    stim = 1 - np.exp(-h * nu[None, :] / (k * T))
    nH0, nHm, ne, nHp = (es[q][:, None] for q in ("nH0", "nHm", "ne", "nHp"))
    # H I vezano-slobodni: sum po n, Kramers (g_bf = 1)
    kap = np.zeros((T.shape[0], nu.size))
    for n in range(1, 9):
        nu_n = NU_RH / n**2
        dE = 13.598 * (1 - 1 / n**2) * eV
        Nn = nH0 * 2 * n**2 / 2.0 * np.exp(-dE / (k * T))
        sig = np.where(nu[None, :] >= nu_n, 2.815e29 / (n**5 * nu[None, :]**3), 0.0)
        kap += Nn * sig
    # H I slobodno-slobodni
    kap += 3.692e8 * ne * nHp / (np.sqrt(T) * nu[None, :]**3)
    kap *= stim
    k_Hbf = kap - 3.692e8 * ne * nHp / (np.sqrt(T) * nu[None, :]**3) * stim
    # H- vezano-slobodni i slobodno-slobodni
    k_Hmbf = nHm * sigma_Hm_bf(lamA / 1e4) * stim
    k_Hmff = nH0 * alpha_Hm_ff(lamA, T) * ne * k * T          # P_e = n_e k T; stim. emisija već u fitu
    kap = kap + k_Hmbf + k_Hmff
    # rasejanje
    lamc = np.maximum(lamA, 1026.0)
    s_e = ne * sigT + 0 * lamA
    s_R = nH0 * (5.799e-13 / lamc**4 + 1.422e-6 / lamc**6 + 2.784 / lamc**8)
    sca = s_e + s_R
    if komponente:
        return {"H I b-f": k_Hbf, "H I f-f": kap - k_Hbf - k_Hmbf - k_Hmff, "H$^-$ b-f": k_Hmbf,
                "H$^-$ f-f": k_Hmff, "Tomson (e$^-$)": s_e, "Rejli (H I)": s_R}
    return kap, sca


def planck(nu, T):
    return 2 * h * nu**3 / c**2 / np.expm1(h * nu / (k * T))


def dBdT(nu, T):
    x = h * nu / (k * T)
    return planck(nu, T) * x / T * np.exp(x) / np.expm1(x)


# ---------------------------------------------------------------- frekvencijska mreža
def mreza_frekvencija(n=700):
    lam = np.geomspace(50e-7, 8e-4, n)              # 50 nm - 8 um
    ivice = [c / (NU_RH / n**2) * 1e7 for n in range(1, 7)] + [1642.0]   # nm (H I ivice, H- prag)
    extra = []
    for e in ivice:
        extra += [e * 1e-7 * (1 - 1e-4), e * 1e-7 * (1 + 1e-4)]
    lam = np.sort(np.concatenate([lam, extra]))
    return c / lam[::-1]                             # rastuće nu


# ---------------------------------------------------------------- struktura
def q_hopf(t):
    return 0.7104 - 0.1331 * np.exp(-3.4488 * t)


def hidrostatika(T, tauR, logg, nu):
    """Integracija dp/dtau_R = g rho / chi_R po ln tau (iterativno u svakom koraku). Vraća P i EOS."""
    g = 10**logg
    P = np.zeros_like(tauR)

    def chiR_rho(Pi, Ti):
        es = jednacina_stanja(np.array([Pi]), np.array([Ti]))
        kap, sca = neprozracnost(nu, np.array([Ti]), es)
        chi = (kap + sca)[0]
        w = dBdT(nu, Ti)
        dnu = np.gradient(nu)
        chiR = np.sum(w * dnu) / np.sum(w / chi * dnu)
        return chiR, es["rho"][0]

    # početni sloj: P = tau g / kappa_R (iteracija)
    Pi = 1.0
    for _ in range(30):
        chiR, rho = chiR_rho(Pi, T[0])
        Pi = tauR[0] * g * rho / chiR
    P[0] = Pi
    for i in range(1, tauR.size):
        dlt = np.log(tauR[i] / tauR[i - 1])
        chiR0, rho0 = chiR_rho(P[i - 1], T[i - 1])
        f0 = tauR[i - 1] * g * rho0 / chiR0
        Pg = P[i - 1] + f0 * dlt
        for _ in range(4):        # trapezno pravilo, iterirano
            chiR1, rho1 = chiR_rho(Pg, T[i])
            f1 = tauR[i] * g * rho1 / chiR1
            Pg = P[i - 1] + 0.5 * (f0 + f1) * dlt
        P[i] = Pg
    return P


def rosseland(nu, T, chi):
    dnu = np.gradient(nu)
    w = dBdT(nu[None, :], T[:, None])
    return np.sum(w * dnu, 1) / np.sum(w / chi * dnu, 1)


# ---------------------------------------------------------------- formalno rešenje
from scipy.special import expn, exp1


def _En(n, x):
    out = np.zeros_like(x)
    m = x < 60.0
    out[m] = expn(n, x[m])
    return out


def _E234(x):
    """E2, E3, E4 iz E1 rekurzijom E_{n+1} = (exp(-x) - x E_n) / n (za x < 60; iznad su zanemarljivi)."""
    x = np.minimum(x, 60.0)
    ex = np.exp(-x)
    E1 = exp1(np.maximum(x, 1e-300))
    E2 = np.where(x > 0, ex - x * E1, 1.0)
    E3 = (ex - x * E2) / 2
    E4 = (ex - x * E3) / 3
    return E2, E3, E4


def momenti(tau, S):
    """J = 1/2 Int S E1(|t - tau|) dt  i  H = 1/2 [Int_tau^inf S E2(t - tau) dt - Int_0^tau S E2(tau - t) dt],
    sa S linearnom po delovima između tačaka mreže (tačni integrali preko E_n). tau, S: (ND, NF)."""
    ND, NF = tau.shape
    J = np.zeros((ND, NF))
    H = np.zeros((ND, NF))
    dS_dno = (S[-1] - S[-2]) / (tau[-1] - tau[-2])
    for i in range(ND):
        ti = tau[i]
        a = np.abs(tau[:-1] - ti)          # krajevi segmenata
        b = np.abs(tau[1:] - ti)
        gore = tau[:-1] >= ti              # segment dublje od tačke i (t > tau_i)
        lo, hi = np.where(gore, a, b), np.where(gore, b, a)     # lo < hi
        S_lo = np.where(gore, S[:-1], S[1:])
        S_hi = np.where(gore, S[1:], S[:-1])
        D = np.maximum(hi - lo, 1e-300)
        E2l, E3l, E4l = _E234(lo)
        E2h, E3h, E4h = _E234(hi)
        lo_c, hi_c = np.minimum(lo, 60.0), np.minimum(hi, 60.0)
        # J (jezgro E1): A0 = -E2, A1 = -x E2 - E3 ;  H (jezgro E2): A0 = -E3, A1 = -x E3 - E4
        for out, (El1, Eh1, El2, Eh2), sgn in ((J, (E2l, E2h, E3l, E3h), 1.0),
                                              (H, (E3l, E3h, E4l, E4h), np.where(gore, 1.0, -1.0))):
            dA0 = El1 - Eh1
            dA1 = (lo_c * El1 + El2) - (hi_c * Eh1 + Eh2)
            integ = (S_lo * (hi * dA0 - dA1) + S_hi * (dA1 - lo * dA0)) / D
            out[i] = 0.5 * np.sum(sgn * integ, 0)
        # ispod najdublje tačke: linearna ekstrapolacija S
        e2, e3, e4 = _E234(tau[-1] - ti)
        J[i] += 0.5 * (S[-1] * e2 + dS_dno * e3)
        H[i] += 0.5 * (S[-1] * e3 + dS_dno * e4)
        # iznad najviše tačke (0 < t < tau_0): S = S_0
        a2, a3, _ = _E234(ti - tau[0])
        b2, b3, _ = _E234(ti)
        J[i] += 0.5 * S[0] * (a2 - b2)
        H[i] -= 0.5 * S[0] * (a3 - b3)
    return J, H


def lambda_dijag(tau):
    """Dijagonala Lambda operatora za S linearnu po delovima: Lambda*_ii = 1/2 [w(D+) + w(D-)],
    w(D) = 1 - (1/2 - E3(D)) / D, gde su D+- optički razmaci do susednih tačaka."""
    L = np.zeros_like(tau)
    D = np.diff(tau, axis=0)
    w = 1 - (0.5 - _En(3, D)) / np.maximum(D, 1e-12)
    w = np.where(D < 1e-6, 0.5 * D, w)             # mali razmaci: w ~ D/2
    L[:-1] += 0.5 * w
    L[1:] += 0.5 * w
    return L


def rasejanje(tau, B, eps, S0, tol=1e-3, nmax=40):
    """Ubrzana Lambda iteracija (ALI) za S = eps B + (1 - eps) J, sa dijagonalnim Lambda*."""
    Ld = lambda_dijag(tau)
    S = S0.copy()
    for _ in range(nmax):
        J, H = momenti(tau, S)
        S_novo = (eps * B + (1 - eps) * (J - Ld * S)) / (1 - (1 - eps) * Ld)
        # konvergenciju proveravamo tamo gde prava apsorpcija nije zanemarljiva (eps > 0.01)
        gotovo = np.max(np.where(eps > 0.01, abs(S_novo / S - 1), 0.0)) < tol
        S = S_novo
        if gotovo:
            break
    J, H = momenti(tau, S)
    return S, J, H


# ---------------------------------------------------------------- model
class Model:
    def __init__(self, Teff, logg, nd=70, nu=None):
        self.Teff, self.logg = Teff, logg
        self.nu = mreza_frekvencija() if nu is None else nu
        self.tauR = np.logspace(-6, 2.5, nd)
        self.T = Teff * (0.75 * (self.tauR + q_hopf(self.tauR)))**0.25
        self.T_sivo = self.T.copy()
        self.istorija = []
        self.struktura()

    def struktura(self):
        self.P = hidrostatika(self.T, self.tauR, self.logg, self.nu)
        self.es = jednacina_stanja(self.P, self.T)
        self.kap, self.sca = neprozracnost(self.nu, self.T, self.es)
        self.chi = self.kap + self.sca
        self.chiR = rosseland(self.nu, self.T, self.chi)
        # tau_nu = Int chi_nu / chi_R dtau_R
        r = self.chi / self.chiR[:, None]
        dtR = np.diff(self.tauR)[:, None]
        self.tau = np.vstack([r[0] * self.tauR[0], r[0] * self.tauR[0] + np.cumsum(0.5 * (r[1:] + r[:-1]) * dtR, 0)])
        # visina: dz = -dtau_R / chi_R
        dz = -np.diff(self.tauR) / (0.5 * (self.chiR[1:] + self.chiR[:-1]))
        self.z = np.concatenate([[0.0], np.cumsum(dz)])
        self.z -= self.z[np.argmin(abs(self.tauR - 1))]      # z = 0 na tau_R = 1

    def zracenje(self):
        B = planck(self.nu[None, :], self.T[:, None])
        eps = self.kap / self.chi                          # udeo prave apsorpcije
        # polazno S: odnos S/B iz prethodne iteracije temperature (ako postoji)
        S0 = B * self.S_B if getattr(self, "S_B", None) is not None else B.copy()
        S, self.J, self.H = rasejanje(self.tau, B, eps, S0)  # koherentno rasejanje (ALI)
        self.S, self.S_B = S, S / B
        self.B = B
        dnu = np.gradient(self.nu)
        self.Htot = np.sum(self.H * dnu, 1)
        self.H0 = sigSB * self.Teff**4 / (4 * np.pi)
        self.Fnu = 4 * np.pi * self.H[0]
        return self.Htot / self.H0 - 1

    def spektar(self, nu_fino):
        """Izlazni fluks F_nu(0) na finoj mreži frekvencija, za konvergiranu strukturu (sa rasejanjem)."""
        kap, sca = neprozracnost(nu_fino, self.T, self.es)
        chi = kap + sca
        r = chi / self.chiR[:, None]                       # Roselandova skala ostaje ista
        dtR = np.diff(self.tauR)[:, None]
        tau = np.vstack([r[0] * self.tauR[0], r[0] * self.tauR[0] + np.cumsum(0.5 * (r[1:] + r[:-1]) * dtR, 0)])
        B = planck(nu_fino[None, :], self.T[:, None])
        eps = kap / chi
        # polazno S/B: interpolacija sa grube mreže, pa nekoliko Lambda iteracija
        SB = np.array([np.interp(np.log(nu_fino), np.log(self.nu), self.S_B[d]) for d in range(self.T.size)])
        S, J, H = rasejanje(tau, B, eps, B * SB)
        self.nu_fino, self.tau_fino, self.chi_fino = nu_fino, tau, chi
        self.F_fino = 4 * np.pi * H[0]
        return self.F_fino

    def unsold_lucy(self, n_iter=8, verbose=True):
        """Unsöld-Lucy korekcija temperature; RZ je uslov za pravu apsorpciju (kappa), tau_F je skala fluks-srednje neprozračnosti."""
        for it in range(n_iter):
            err = self.zracenje()
            self.istorija.append((self.T.copy(), err.copy()))
            if verbose:
                print(f"  iter {it}: max |dH/H| = {np.max(abs(err)):.4f}")
            dnu = np.gradient(self.nu)
            Btot = np.sum(self.B * dnu, 1)
            Jtot = np.sum(self.J * dnu, 1)
            kP = np.sum(self.kap * self.B * dnu, 1) / Btot      # RZ: samo prava apsorpcija
            kJ = np.sum(self.kap * self.J * dnu, 1) / Jtot
            kH = np.sum(self.chi * self.H * dnu, 1) / self.Htot
            dH = self.H0 - self.Htot
            # tau po fluks-srednjoj neprozračnosti
            dtF = np.diff(self.tauR) * 0.5 * ((kH / self.chiR)[1:] + (kH / self.chiR)[:-1])
            integ = np.concatenate([[0.0], np.cumsum(0.5 * (dH[1:] + dH[:-1]) * dtF)])
            # Eddingtonova aproksimacija: J ~ 3 Int_0^tau H dtau_H + 2 H(0)  ->  potrebna promena J:
            dJ = 3 * integ + 2 * dH[0]
            # nova B tako da kappa_P B' = kappa_J (J + dJ)
            dB = (kJ / kP) * (Jtot + dJ) - Btot
            dT = dB / (4 * sigSB * self.T**3 / np.pi)
            dT = np.clip(dT, -0.15 * self.T, 0.15 * self.T)
            self.T = self.T + dT
            self.struktura()
        err = self.zracenje()
        self.istorija.append((self.T.copy(), err.copy()))
        if verbose:
            print(f"  kraj : max |dH/H| = {np.max(abs(err)):.4f}")
        return self


if __name__ == "__main__":
    import time
    t0 = time.time()
    m = Model(5777, 4.44, nd=80, nu=mreza_frekvencija(300))
    m.unsold_lucy(10)
    i1 = np.argmin(abs(m.tauR - 1))
    print("tau_R=1: T=%.0f P=%.3e ne=%.3e rho=%.3e" % (m.T[i1], m.P[i1], m.es['ne'][i1], m.es['rho'][i1]))
    print("time", time.time() - t0)
