"""Jednačine za Predavanje 07 (Od sive atmosfere do teorijskih LTR spektara) - samo nove jednačine;
ostale su preuzete iz predavanja 06 (jednacine_p06.py), renderovane LaTeX-om u PNG.

Konvencija (ista kao u jednacine_p05.py): frekvencija nu (ne lambda); fizički fluks je F,
astrofizički \\mathcal{F} = F/pi, Eddingtonov H = F/(4 pi); emisivnost eta_nu; visina z (raste ka spolja),
pa je d tau_nu = -chi_nu dz.

Upotreba:  python jednacine_p06.py [izlazni_direktorijum]
Potrebno:  pdflatex (sa standalone, amsmath, mathtools, tikz) i pdftocairo (poppler).
LaTeX izvor svake jednačine stoji i u alt-tekstu slike u prezentaciji.
"""
import json
import os
import subprocess
import sys
import tempfile

DPI = 600

EQ = {
    # slajd 2: podsetnik (Milneov problem)
    "p07_podsetnik": r"""\begin{aligned}
\mu\frac{dI}{d\tau} &= I - S, \qquad S = J = B = \frac{\sigma T^4}{\pi}\\
J(\tau) &= \Lambda_\tau[S(t)] \quad\Rightarrow\quad B(\tau) = \Lambda_\tau[B(t)]\\
F &= \pi\mathcal{F} = 4\pi H = \sigma T_{\rm eff}^4 = \mathrm{const}
\end{aligned}""",
    # slajd 2: podsetnik na približna rešenja
    "p07_momenti": r"\frac{dH}{d\tau} = J - S = 0 \;\Rightarrow\; H = \frac{F}{4\pi} = \mathrm{const}, \qquad \frac{dK}{d\tau} = H \;\Rightarrow\; K = H\tau + \mathrm{const}",
    "p07_resenja": r"""\begin{array}{lcc}
 & \text{\v{S}varc\v{s}ild--\v{S}uster} & \text{Eddington}\\[2pt] \hline
\text{pretpostavka} & I^{\pm}\ \text{izotropni u polusferama} & K = J/3,\;\; J(0) = 2H(0)\\[4pt]
S(\tau) = J(\tau) & \dfrac{F}{\pi}\left(\tau + \dfrac{1}{2}\right) & \dfrac{3F}{4\pi}\left(\tau + \dfrac{2}{3}\right)\\[10pt]
T^4(\tau) & T_{\rm eff}^4\left(\tau + \dfrac{1}{2}\right) & \dfrac{3}{4}\,T_{\rm eff}^4\left(\tau + \dfrac{2}{3}\right)\\[10pt]
T = T_{\rm eff} & \tau = 1/2 & \tau = 2/3
\end{array}""",
    # slajd 4: Lambda iteracija
    "p07_lambda_iter": r"B^{(n+1)}(\tau) = \Lambda_\tau\!\left[B^{(n)}(t)\right], \qquad B^{(0)}(\tau) = \frac{3F}{4\pi}\left(\tau + \frac{2}{3}\right), \qquad q^{(n)}(\tau) = \frac{4\pi}{3F}\,B^{(n)}(\tau) - \tau",
    # slajdovi po Grayu (2005), pogl. 9
    "p07g_skaliranje": r"T(\tau_0) = S_0\,T_\odot(\tau_0), \qquad S_0 = \frac{T_{\rm eff}}{T_{\rm eff}^{\odot}}",
    "p07g_doprinos": r"F_\nu(0) = \int_{-\infty}^{\infty} \underbrace{2\pi\ln 10\; S_\nu(\tau_\nu)\,E_2(\tau_\nu)\,\frac{\chi_\nu}{\chi_0}\,\tau_0}_{\text{funkcija doprinosa}}\; d\log\tau_0",
    # Unsöld-Lucy korekcija
    "p07_UL_1": r"J(\tau) \approx 3\int_0^{\tau} H\,d\tau'_F + 2H(0) \quad\Rightarrow\quad \Delta J(\tau) \approx 3\int_0^{\tau}\Delta H\,d\tau'_F + 2\,\Delta H(0)",
    "p07_UL_dH": r"\Delta H(\tau) = \frac{\sigma T_{\rm eff}^4}{4\pi} - H(\tau)",
    "p07_UL_2": r"""\begin{aligned}
\kappa_P\,(B + \Delta B) &= \kappa_J\,(J + \Delta J)\\
\Delta B &= \frac{\kappa_J}{\kappa_P}\,(J + \Delta J) - B, \qquad \Delta T = \frac{\pi\,\Delta B}{4\sigma T^3}
\end{aligned}""",
    "p07_UL_def": r"\kappa_J = \frac{\displaystyle\int_0^\infty \kappa_\nu J_\nu\,d\nu}{J}, \qquad \kappa_P = \frac{\displaystyle\int_0^\infty \kappa_\nu B_\nu\,d\nu}{B}, \qquad d\tau_F = -\chi_F\,dz",
    # naš model
    "p07_model": r"""\begin{aligned}
\chi_\nu &= \underbrace{\kappa_\nu^{\rm bf,ff}({\rm H\,I}) + \kappa_\nu^{\rm bf,ff}({\rm H}^-)}_{\kappa_\nu} + \underbrace{n_e\,\sigma_T + n_{\rm H\,I}\,\sigma^{\rm R}_{\nu}}_{\sigma_\nu}\\
S_\nu &= \frac{\kappa_\nu B_\nu(T) + \sigma_\nu J_\nu}{\chi_\nu}\\
\frac{dp}{d\tau_R} &= \frac{g\,\rho}{\chi_R}\\
T^4(\tau_R) &= \tfrac{3}{4}\,T_{\rm eff}^4\left[\tau_R + q(\tau_R)\right] \quad \text{(po\v{c}etna struktura)}
\end{aligned}""",
    "p07_fluks": r"F_\nu(0) = 2\pi\int_0^\infty B_\nu\!\left(T(t_\nu)\right)E_2(t_\nu)\,dt_\nu, \qquad \tau_\nu = \int_0^{\tau_R}\frac{\chi_\nu}{\chi_R}\,d\tau_R'",
    "p07_D": r"D = \log_{10}\frac{F_\nu(\nu_2^-)}{F_\nu(\nu_2^+)}, \qquad h\nu_2 = \frac{13.6\,{\rm eV}}{4}\;\;(\lambda_2 = 364.7\,{\rm nm})",
}


TEMPLATE = r"""\documentclass[10pt,border=3pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,mathtools,tikz}
\usetikzlibrary{shapes.geometric,positioning}
\begin{document}
%s
\end{document}
"""


def render(name, body, tex, outdir):
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "eq.tex"), "w") as f:
            f.write(TEMPLATE % body)
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
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "jednacine_p07")
    os.makedirs(outdir, exist_ok=True)
    meta = {name: render(name, r"$\displaystyle %s$" % tex, tex, outdir) for name, tex in EQ.items()}
    with open(os.path.join(outdir, "jednacine_p07.json"), "w") as f:
        json.dump(meta, f, indent=1)
    for k, v in meta.items():
        print(f"{k:28s} {v['w_in']:.2f} x {v['h_in']:.2f} in (prirodno)")
