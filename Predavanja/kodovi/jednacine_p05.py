"""Jednačine za Predavanje 05 (Atmosfere u RZ), renderovane LaTeX-om u PNG.

Konvencija: frekvencija nu (ne lambda); fizički fluks je F,
astrofizički \\mathcal{F} = F/pi, Eddingtonov H = F/(4 pi); emisivnost eta_nu.

Upotreba:  python jednacine_p05.py [izlazni_direktorijum]
Potrebno:  pdflatex (sa standalone, amsmath, mathtools, tikz) i pdftocairo (poppler).
Slike se u .pptx ubacuju uglavnom sa 1.8-2x prirodnom veličinom (10pt -> 18-20pt),
a LaTeX izvor stoji i u alt-tekstu slike u prezentaciji.
"""
import json
import os
import subprocess
import sys
import tempfile

DPI = 600

EQ = {
    # slajd 2: formalno rešenje
    "s02_formalno_resenje": r"""\begin{aligned}
I_\nu(\tau_\nu,\mu>0) &= \int_{\tau_\nu}^{\infty} S_\nu(t)\,e^{-(t-\tau_\nu)/\mu}\,dt/\mu\\
I_\nu(\tau_\nu,\mu<0) &= \int_{0}^{\tau_\nu} S_\nu(t)\,e^{-(\tau_\nu-t)/|\mu|}\,dt/|\mu|
\end{aligned}""",
    # slajd 5: održanje energije
    "s05_luminoznost": r"""\begin{aligned}
L(r) &= \frac{dE^{\rm out}-dE^{\rm in}}{dt} = \mathrm{const}\\
4\pi r^2\,F(r) &= \mathrm{const}
\end{aligned}""",
    "s05_fluks_const": r"""\begin{aligned}
F(r) &\approx \mathrm{const}\\
\int_0^\infty F_\nu(z)\,d\nu &= \mathrm{const}
\end{aligned}""",
    # slajd 6: definicije fluksa
    "s06_fluks_def": r"""\begin{aligned}
F_\nu &= \oint I_\nu\cos\theta\,\sin\theta\,d\theta\,d\phi\\
\mu &= \cos\theta\\
F_\nu &= 2\pi\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu
\end{aligned}""",
    "s06_fluks": r"\mathmakebox[1.3em][l]{F_\nu} = 2\pi\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu",
    "s06_astrofizicki_fluks": r"\mathmakebox[1.3em][l]{\mathcal{F}_\nu} = 2\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu",
    "s06_eddingtonov_fluks": r"\mathmakebox[1.3em][l]{H_\nu} = \frac{1}{2}\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu",
    # prelazni slajd posle definicije fluksa: šta posmatramo
    "s06b_F_povrsina": r"F_\nu(0) = 2\pi\int_0^1 I_\nu(0,\mu)\,\mu\,d\mu",
    "s06b_I_povrsina": r"I_\nu(0,\mu) = \int_0^\infty S_\nu(t)\,e^{-t/\mu}\,dt/\mu",
    # slajd 7: linearna zavisnost intenziteta od mu
    "s07_I_linearno": r"I_\nu^+(0,\mu) = a_\nu + b_\nu\mu",
    # slajd 8
    "s08_mapiranje": r"I_\nu^+(0,\mu)\;\rightarrow\;S_\nu(\tau_\nu)\;\rightarrow\;T(\tau_\nu)",
    # slajd 12: ME aproksimacija za fluks
    "s12_ME_fluks": r"""\begin{aligned}
S_\nu &= a_\nu + b_\nu\tau_\nu\\
F_\nu(0) &= 2\pi\int_{-1}^{1} I_\nu(0,\mu)\,\mu\,d\mu = 2\pi\int_{0}^{1}(a_\nu+b_\nu\mu)\,\mu\,d\mu\\
F_\nu(0) &= \pi\left(a_\nu+\tfrac{2}{3}b_\nu\right) = \pi S_\nu(\tau_\nu=2/3)
\end{aligned}""",
    # slajd 13: siva atmosfera
    "s13_siva_ME": r"""\begin{aligned}
F &= \pi S(\tau=2/3)\\
T_{\rm eff} &= T(\tau=2/3)
\end{aligned}""",
    # slajd 16: ravnoteža zračenja
    "s16_RZ": r"""\begin{gathered}
dE_{\rm abs} = dE_{\rm em}\\
dA\,dt\,dl\oint\int_0^\infty \chi_\nu I_\nu\,d\nu\,d\Omega
 = dA\,dt\,dl\oint\int_0^\infty \eta_\nu\,d\nu\,d\Omega\\
\int_0^\infty \chi_\nu J_\nu\,d\nu = \int_0^\infty \chi_\nu S_\nu\,d\nu
\end{gathered}""",
    # slajd 17: srednji intenzitet i Eddingtonov fluks
    "s17_J_ugao": r"\mathmakebox[1.3em][l]{J_\nu} = \frac{1}{4\pi}\oint I_\nu(\theta,\phi)\,\sin\theta\,d\theta\,d\phi",
    "s17_J_mu": r"\mathmakebox[1.3em][l]{J_\nu} = \frac{1}{2}\int_{-1}^{1} I_\nu(\mu)\,d\mu",
    "s17_H_mu": r"\mathmakebox[1.3em][l]{H_\nu} = \frac{1}{2}\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu",
    # slajd 18 (ranije 19): ravnoteža zračenja u unutrašnjosti zvezde
    # z = visina (raste ka spolja), pa je d tau = -chi dz i F = -(16 sigma T^3 / 3 chi) dT/dz
    "s19_RTE": r"""\begin{aligned}
&\mu\frac{dI_\nu(z,\mu)}{dz} = -\chi_\nu(z)\,I_\nu(z,\mu) + \eta_\nu(z)\\[4pt]
&\frac{d}{dz}\int_0^\infty d\nu\int I_\nu(z,\mu)\,\mu\,d\Omega = \frac{d}{dz}\int_0^\infty F_\nu(z)\,d\nu = 0
\end{aligned}""",
    "s19_F_const": r"F(z) = \int_0^\infty F_\nu(z)\,d\nu = \pi\mathcal{F}(z) = \mathrm{const}",
    "s19_F_zvezda": r"F_* = \pi\mathcal{F}_* = \sigma T_{\rm eff}^4",
    "s19_RTE_tau": r"-\mu\frac{1}{\chi_\nu(z)}\frac{dI_\nu(z,\mu)}{dz} = I_\nu(z,\mu) - S_\nu(z)\qquad\qquad \int d\nu\int\mu\,d\Omega",
    "s19_F_K": r"F = -4\pi\int_0^\infty\frac{1}{\chi_\nu}\frac{dK_\nu}{dz}\,d\nu = -c\int_0^\infty\frac{1}{\chi_\nu}\frac{dp_\nu}{dz}\,d\nu = \mathrm{const}",
    "s19_F_p": r"F = -c\,\frac{1}{\chi}\frac{dp^{\rm zr}}{dz}\qquad\quad p^{\rm zr} = \frac{4\pi}{c}K\qquad\quad K\approx\frac{1}{3}J\approx\frac{1}{3}B = \frac{1}{3}\frac{\sigma T^4}{\pi}",
    "s19_F_gradT": r"F = -\frac{c}{\chi}\,\frac{4\pi}{c}\,\frac{4}{3\pi}\,\sigma T^3\frac{dT}{dz} = -\frac{16}{3}\frac{\sigma}{\chi}T^3\,\tikz[baseline=(X.base)]\node[draw=red,thick,ellipse,inner sep=0.5pt](X){$\dfrac{dT}{dz}$};",
    # slajd 20: K integral
    "s20_K_def": r"K = \int_0^\infty K_\nu\,d\nu = \frac{1}{4\pi}\int_0^\infty\!\oint I_\nu\cos^2\theta\,d\Omega\,d\nu",
    # slajd 21: momenti
    "s21_momenti": r"""\begin{aligned}
J &= \frac{1}{2}\int_{-1}^{1} I(\mu)\,d\mu\\
H &= \frac{1}{2}\int_{-1}^{1} I(\mu)\,\mu\,d\mu\\
K &= \frac{1}{2}\int_{-1}^{1} I(\mu)\,\mu^2\,d\mu
\end{aligned}""",
    "s21_izotropno": r"""\begin{aligned}
J &= I = 3K\\
H &= 0
\end{aligned}""",
    # slajd 22
    "s22_S_J": r"S = J",
    "s22_S_B_J": r"S = B = J",
    # slajd 23
    "s23_Milne": r"\mu\frac{dI}{d\tau} = I - S = I - J",
    # slajd 23: Lambda operator (dodata integracija po mu)
    "s23_lambda": r"""\begin{aligned}
J &= \frac{1}{2}\int_{-1}^{1} I(\mu)\,d\mu\\
&= \frac{1}{2}\int_0^1 d\mu\int_\tau^\infty S(t)\,e^{(\tau-t)/\mu}\,dt/\mu
 + \frac{1}{2}\int_{-1}^0 d\mu\int_0^\tau S(t)\,e^{(t-\tau)/|\mu|}\,dt/|\mu|\\
&= \frac{1}{2}\int_0^\infty S(t)\,@E1@(|t-\tau|)\,dt\\
J(\tau) &= \Lambda_\tau[S(t)]
\end{aligned}""",
    # slajd 26: srednji intenzitet kroz Lambda operator (znak eksponenta, faktor 1/2)
    "s26_J_lambda": r"""\begin{aligned}
J(\tau) &= \frac{1}{2}\int_{-1}^{1} I(\tau,\mu)\,d\mu\\
&= \frac{1}{2}\int_0^1\!\!\int_\tau^\infty S(t)\,e^{-(t-\tau)/\mu}\,\frac{dt}{|\mu|}\,d\mu
 + \frac{1}{2}\int_{-1}^0\!\!\int_0^\tau S(t)\,e^{(\tau-t)/\mu}\,\frac{dt}{|\mu|}\,d\mu\\
&= \frac{1}{2}\int_0^\infty S(t)\,E_1(|t-\tau|)\,dt
\end{aligned}""",
    # slajd 27: Milneov problem preko Lambda operatora (Phi daje astrofizicki fluks)
    "s27_milne_lambda": r"""\begin{aligned}
J &= \Lambda[S]\\
J &= S = B\\
B &= \Lambda[B]\\
F &= \pi\,\Phi[B] = \sigma T_{\rm eff}^4
\end{aligned}""",
    # slajd 26: Mihalas: Phi = 2 int S E2 (astrofizicki fluks), X = 2 int S E3
    "s26_operatori": r"""\begin{aligned}
\mathcal{F}_\tau &= \Phi_\tau[S(t)]\\
K_\tau &= \frac{1}{4}X_\tau[S(t)]
\end{aligned}""",
}

TEMPLATE = r"""\documentclass[10pt,border=3pt]{standalone}
\usepackage{amsmath,amssymb,mathtools,tikz}
\usetikzlibrary{shapes.geometric}
\begin{document}
$\displaystyle %s$
\end{document}
"""


def render(name, tex, outdir):
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "eq.tex"), "w") as f:
            f.write(TEMPLATE % tex)
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "eq.tex"],
                       cwd=tmp, check=True, stdout=subprocess.DEVNULL)
        info = subprocess.run(["pdfinfo", "eq.pdf"], cwd=tmp, check=True,
                              capture_output=True, text=True).stdout
        w_pt, h_pt = [float(v) for v in
                      next(l for l in info.splitlines() if l.startswith("Page size")).split()[2:5:2]]
        subprocess.run(["pdftocairo", "-png", "-r", str(DPI), "-transp", "-singlefile",
                        "eq.pdf", os.path.join(os.path.abspath(outdir), name)], cwd=tmp, check=True)
    return {"w_in": w_pt / 72.0, "h_in": h_pt / 72.0, "tex": tex}


if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "jednacine_p05")
    os.makedirs(outdir, exist_ok=True)
    meta = {name: render(name, tex, outdir) for name, tex in EQ.items()}
    with open(os.path.join(outdir, "jednacine_p05.json"), "w") as f:
        json.dump(meta, f, indent=1)
    for k, v in meta.items():
        print(f"{k:24s} {v['w_in']:.2f} x {v['h_in']:.2f} in (prirodno)")
