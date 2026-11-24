"""Jednačine za Predavanje 06 (Modelovanje atmosfera u ravnoteži zračenja), renderovane LaTeX-om u PNG.

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
    # slajd 2: podsetnik
    "p06_podsetnik": r"""\begin{aligned}
\mu\frac{dI}{d\tau} &= I - S, \qquad S = J = B = \frac{\sigma T^4}{\pi}\\
J(\tau) &= \Lambda_\tau[S(t)] = \frac{1}{2}\int_0^\infty S(t)\,E_1(|t-\tau|)\,dt\\
F &= \pi\mathcal{F} = 4\pi H = \sigma T_{\rm eff}^4 = \mathrm{const}
\end{aligned}""",
    # slajd 3: pretpostavke modela
    "p06_pretpostavke": r"""\begin{aligned}
\mu\frac{dI_\nu}{dz} &= -\chi_\nu I_\nu + \eta_\nu\\[3pt]
\int_0^\infty \chi_\nu S_\nu\,d\nu &= \int_0^\infty \chi_\nu J_\nu\,d\nu\\[3pt]
S_\nu &= B_\nu(T)\\[3pt]
\chi_\nu &= f(T,p;\nu)
\end{aligned}""",
    # slajd 4: RZ u sivoj i ne-sivoj atmosferi
    "p06_RZ_siva": r"S = J",
    "p06_RZ_nesiva": r"\int_0^\infty \chi_\nu J_\nu\,d\nu = \int_0^\infty \chi_\nu S_\nu\,d\nu",
    # slajd 5: izotermna atmosfera
    "p06_izotermna": r"""\begin{aligned}
S &= B = 1\\
J(\tau) &= \Lambda_\tau[1] = 1 - \tfrac{1}{2}E_2(\tau)\\
J(0) &= \tfrac{1}{2}
\end{aligned}""",
    # slajd 7: Milneov problem preko Lambda operatora
    "p06_milne_lambda": r"J(\tau) = \Lambda_\tau[J(t)] = \frac{1}{2}\int_0^\infty J(t)\,E_1(|t-\tau|)\,dt",
    # slajd 8: momenti JPZ (siva atmosfera, RZ)
    "p06_momenti_JPZ": r"""\begin{aligned}
\mu\frac{dI}{d\tau} &= I - S\\
\frac{dH}{d\tau} &= J - S = 0 \;\;\rightarrow\;\; H = \frac{F}{4\pi} = \mathrm{const}\\
\frac{dK}{d\tau} &= H \;\;\rightarrow\;\; K = H\tau + C
\end{aligned}""",
    # slajd 9: Švarcšild-Šusterovo rešenje
    "p06_SS_def": r"""\begin{aligned}
J^+ &= \int_0^1 I(\mu)\,d\mu, \qquad J^- = \int_{-1}^0 I(\mu)\,d\mu\\
J &= \tfrac{1}{2}\left(J^+ + J^-\right)\\
\int_0^1 I(\mu)\,\mu\,d\mu &\approx \tfrac{1}{2}J^+, \qquad \int_{-1}^0 I(\mu)\,\mu\,d\mu \approx -\tfrac{1}{2}J^-\\
F &= 2\pi\int_{-1}^{1} I(\mu)\,\mu\,d\mu = \pi\left(J^+ - J^-\right)
\end{aligned}""",
    "p06_SS_izvodjenje": r"""\begin{aligned}
\mu\frac{dI}{d\tau} &= I - J \qquad \Big|\int_0^1 d\mu,\;\; \int_{-1}^0 d\mu\\
\frac{1}{2}\frac{dJ^+}{d\tau} &= J^+ - J\\
-\frac{1}{2}\frac{dJ^-}{d\tau} &= J^- - J\\
\frac{1}{2}\frac{d}{d\tau}\left(J^+ + J^-\right) &= J^+ - J^- = \frac{F}{\pi}\\
\frac{dJ}{d\tau} &= \frac{F}{\pi} \;\;\rightarrow\;\; J = \frac{F}{\pi}\,\tau + C
\end{aligned}""",
    # slajd 10: veza sa astrofizičkim fluksom
    "p06_SS_astro": r"J(\tau) = \mathcal{F}\left(\tau + \tfrac{1}{2}\right), \qquad \mathcal{F} = \frac{F}{\pi}\;\;\text{(astrofizi\v{c}ki fluks)}",
    "p06_SS_resenje": r"""\begin{aligned}
C &= J(0) = \tfrac{1}{2}J^+(0) = \frac{F}{2\pi} \qquad \left(J^-(0) = 0\right)\\
J(\tau) &= S(\tau) = \frac{F}{\pi}\left(\tau + \frac{1}{2}\right)
\end{aligned}""",
    # slajd 10: Eddingtonovo rešenje
    "p06_Edd_1": r"\frac{dK}{d\tau} = H \;\;\rightarrow\;\; K = H\tau + C",
    "p06_Edd_2": r"K(\tau) \approx \frac{J(\tau)}{3} = \frac{S(\tau)}{3}",
    "p06_Edd_3": r"S(\tau) = 3H\tau + C'",
    "p06_Edd_4": r"J(0) = 2H \;\;\rightarrow\;\; C' = 2H",
    "p06_Edd_5": r"S(\tau) = 3H\left(\tau + \frac{2}{3}\right) = \frac{3}{4}\frac{F}{\pi}\left(\tau + \frac{2}{3}\right)",
    # slajd 11: uslovi koje proveravamo
    "p06_provera": r"\Lambda_\tau[S] = S, \qquad F(\tau) = \pi\,\Phi_\tau[S] = \mathrm{const}",
    # slajd 12: poređenje
    "p06_por_SS": r"""\begin{aligned}
S(\tau) &= \frac{F}{\pi}\left(\tau + \frac{1}{2}\right)\\[2pt]
T^4(\tau) &= T_{\rm eff}^4\left(\tau + \frac{1}{2}\right)\\[2pt]
\frac{I(0,\mu)}{I(0,1)} &= \frac{2}{3}\left(\mu + \frac{1}{2}\right)\\[2pt]
I(0,0) &= 0.33\,I(0,1)
\end{aligned}""",
    "p06_por_Edd": r"""\begin{aligned}
S(\tau) &= \frac{3}{4}\frac{F}{\pi}\left(\tau + \frac{2}{3}\right)\\[2pt]
T^4(\tau) &= \frac{3}{4}\,T_{\rm eff}^4\left(\tau + \frac{2}{3}\right)\\[2pt]
\frac{I(0,\mu)}{I(0,1)} &= \frac{3}{5}\left(\mu + \frac{2}{3}\right)\\[2pt]
I(0,0) &= 0.40\,I(0,1)
\end{aligned}""",
    # slajd 13: posmatrački test (vrednosti iz stare tabele, 2025)
    "p06_tabela_LD": r"""\begin{array}{cccc}
\mu & \text{posmatranje} & \text{R.Z. (\v{S}--\v{S})} & \text{R.Z. (Edd)}\\ \hline
1 & 1 & 1 & 1\\
0.8 & 0.92 & 0.87 & 0.88\\
0.6 & 0.81 & 0.73 & 0.76\\
0.44 & 0.70 & 0.63 & 0.66\\
0.2 & 0.49 & 0.47 & 0.52\\
0 & 0.40 & 0.33 & 0.40
\end{array}""",
    # slajd 14: Hopfovo rešenje
    "p06_Hopf": r"""\begin{aligned}
T^4(\tau) &= \frac{3}{4}\,T_{\rm eff}^4\left[\tau + q(\tau)\right]\\[2pt]
q(0) &= \frac{1}{\sqrt{3}} \approx 0.577 \;\le\; q(\tau) \;\le\; q(\infty) \approx 0.710
\end{aligned}""",
    # slajd 16: izlazni fluks sive atmosfere
    "p06_fluks_siva": r"""\begin{aligned}
F_\nu(0) &= 2\pi\int_0^\infty B_\nu\!\left(T(t)\right)E_2(t)\,dt\\
&\approx \pi B_\nu\!\left(T(\tau = 2/3)\right) = \pi B_\nu(T_{\rm eff})
\end{aligned}""",
    # slajd 18: ne-siva atmosfera, kvalitativno
    "p06_EB_nu": r"F_\nu(0) \approx \pi B_\nu\!\left(T(\tau_\nu = 2/3)\right)",
    "p06_Balmer_plus": r"""\begin{aligned}
\chi_\nu &< \bar\chi \;\;\Rightarrow\;\; z(\tau_\nu = 2/3) < z(\bar\tau = 2/3)\\
&\Rightarrow\;\; T(\tau_\nu = 2/3) > T_{\rm eff}
\end{aligned}""",
    "p06_Balmer_minus": r"""\begin{aligned}
\chi_\nu &> \bar\chi \;\;\Rightarrow\;\; z(\tau_\nu = 2/3) > z(\bar\tau = 2/3)\\
&\Rightarrow\;\; T(\tau_\nu = 2/3) < T_{\rm eff}
\end{aligned}""",
    # slajd 19: hidrostatička ravnoteža
    "p06_HR": r"""\begin{aligned}
\frac{dp}{dz} &= -\rho g, \qquad d\tau = -\chi\,dz\\[2pt]
\frac{dp}{d\tau} &= \frac{\rho}{\chi}\,g\\[2pt]
\chi &= \chi(p,T), \qquad \rho = \rho(p,T)\\[2pt]
T(z) &= \;?\,, \qquad p(z) = \;?
\end{aligned}""",
    # slajd 20: računanje strukture
    "p06_struktura_ulaz": r"T(\tau),\qquad g,\qquad p_0 = p(\tau_0)",
    "p06_struktura": r"""\begin{aligned}
\frac{dp}{d\tau} &= \frac{\rho}{\chi}\,g\\
\rho_0 &= f(T_0,p_0)\\
\chi_0 &= f(T_0,p_0)\\
p_1 &= p_0 + \left(\frac{dp}{d\tau}\right)_{\!0}(\tau_1-\tau_0) \;\;\rightarrow\;\; \rho_1,\;\chi_1,\;\dots
\end{aligned}""",
    # slajd 21: mapiranje tau -> z
    "p06_tau_z": r"\tau_i - \tau_{i-1} = \left(z_{i-1} - z_i\right)\frac{\chi_i + \chi_{i-1}}{2}",
    # slajd 22: srednje neprozračnosti
    "p06_sr_moment": r"\frac{dK_\nu}{dz} = -\chi_\nu H_\nu \qquad\rightarrow\qquad \frac{dK}{dz} = -\int_0^\infty \chi_\nu H_\nu\,d\nu",
    "p06_sr_fluks": r"\chi_F = \frac{\displaystyle\int_0^\infty \chi_\nu H_\nu\,d\nu}{H} \qquad\Rightarrow\qquad \frac{dp^{\rm zr}}{d\tau_F} = \frac{4\pi}{c}H = \frac{F}{c} = \frac{\sigma T_{\rm eff}^4}{c}",
    "p06_sr_Plank": r"\chi_P = \frac{\displaystyle\int_0^\infty \chi_\nu B_\nu(T)\,d\nu}{B(T)}, \qquad B(T) = \frac{\sigma T^4}{\pi}",
    "p06_sr_Roseland": r"H = -\int_0^\infty \frac{1}{\chi_\nu}\frac{dK_\nu}{dz}\,d\nu \equiv -\frac{1}{\chi_R}\frac{dK}{dz} \qquad\Rightarrow\qquad \frac{1}{\chi_R} = \frac{\displaystyle\int_0^\infty \frac{1}{\chi_\nu}\frac{\partial B_\nu}{\partial T}\,d\nu}{\displaystyle\int_0^\infty \frac{\partial B_\nu}{\partial T}\,d\nu}",
    # slajd 23: ne-siva atmosfera
    "p06_nesiva": r"""\begin{aligned}
\mu\frac{dI_\nu}{dz} &= -\chi_\nu I_\nu + \eta_\nu\\[3pt]
2\pi\int_0^\infty\!\!\int_{-1}^{1} I_\nu(\mu)\,\mu\,d\mu\,d\nu &= \sigma T_{\rm eff}^4 = \mathrm{const}\\[3pt]
\int_0^\infty \chi_\nu B_\nu\,d\nu &= \int_0^\infty \chi_\nu J_\nu\,d\nu\\[3pt]
\frac{dp}{dz} &= -\rho g
\end{aligned}""",
    # slajd 24: temperaturska korekcija
    "p06_korekcija": r"""\begin{aligned}
\int_0^\infty \chi_\nu B_\nu(T_0)\,d\nu &\neq \int_0^\infty \chi_\nu J_\nu\,d\nu\\
\int_0^\infty \chi_\nu B_\nu(T_0+\Delta T)\,d\nu &= \int_0^\infty \chi_\nu J_\nu\,d\nu\\
B_\nu(T_0+\Delta T) &\approx B_\nu(T_0) + \left(\frac{\partial B_\nu}{\partial T}\right)_{T_0}\Delta T\\
\int_0^\infty \chi_\nu B_\nu(T_0)\,d\nu + \Delta T\int_0^\infty \chi_\nu\frac{\partial B_\nu}{\partial T}\,d\nu &= \int_0^\infty \chi_\nu J_\nu\,d\nu
\end{aligned}""",
    "p06_dT": r"\Delta T(\tau) = \frac{\displaystyle\int_0^\infty \chi_\nu\left[J_\nu - B_\nu(T_0)\right]d\nu}{\displaystyle\int_0^\infty \chi_\nu\,\frac{\partial B_\nu}{\partial T}\,d\nu}",
    # slajd 30: numeričko rešavanje JPZ
    "p06_kratke_karakteristike": r"""\begin{aligned}
\Delta &= \left(\tau_\nu^{\,i+1} - \tau_\nu^{\,i}\right)/\mu, \qquad \bar S_\nu = \tfrac{1}{2}\left(S_\nu^{\,i+1} + S_\nu^{\,i}\right)\\[2pt]
I_\nu^{\,i} &= I_\nu^{\,i+1}\,e^{-\Delta} + \bar S_\nu\left(1 - e^{-\Delta}\right)
\end{aligned}""",
    "p06_lambda_iteracija": r"S^{(n+1)}(\tau) = \Lambda_\tau\!\left[S^{(n)}\right], \qquad \text{uz}\;\; F = \mathrm{const}",
}

# Šema modelovanja (slajd 25) - TikZ, renderuje se kao slika
SEMA = r"""\begin{tikzpicture}[node distance=5.5mm, font=\small,
  kutija/.style={draw, fill=cyan!12, rounded corners=1pt, minimum width=62mm, minimum height=7.5mm, align=center},
  odluka/.style={draw, fill=cyan!12, diamond, aspect=3.6, inner sep=1pt, align=center},
  >=stealth]
\node[kutija, trapezium, trapezium left angle=110, trapezium right angle=110, minimum width=40mm] (ul)
  {$T_{\rm eff},\;\log g,\;$ hem.\ sastav;\quad $T^{(0)}(\tau)$};
\node[kutija, below=of ul] (mod) {model:\quad $T,\;p,\;p_e,\;\chi_\nu$};
\node[kutija, below=of mod] (tau) {$\displaystyle \tau_\nu(\tau) = \int_0^{\tau}\frac{\chi_\nu(t)}{\bar\chi(t)}\,dt,\qquad B_\nu(\tau_\nu)$};
\node[kutija, below=of tau] (jpz) {jedna\v{c}ina prenosa zra\v{c}enja:\quad $I_\nu(\tau_\nu,\mu),\;\; J_\nu(\tau_\nu)$};
\node[odluka, below=of jpz] (rz) {$\displaystyle\int_0^\infty\!\chi_\nu B_\nu\,d\nu = \int_0^\infty\!\chi_\nu J_\nu\,d\nu\;?$};
\node[kutija, below=of rz] (dT) {$\displaystyle \Delta T(\tau) = \frac{\int_0^\infty \chi_\nu\,[J_\nu - B_\nu(T^{(0)})]\,d\nu}{\int_0^\infty \chi_\nu\,(\partial B_\nu/\partial T)\,d\nu}$};
\node[kutija, below=of dT, minimum width=30mm] (T) {$T(\tau) = T^{(0)}(\tau) + \Delta T(\tau)$};
\node[draw, fill=cyan!12, rounded corners=6pt, right=14mm of rz, minimum height=7mm, minimum width=14mm] (kraj) {kraj};
\draw[->] (ul) -- (mod); \draw[->] (mod) -- (tau); \draw[->] (tau) -- (jpz); \draw[->] (jpz) -- (rz);
\draw[->] (rz) -- node[right]{ne} (dT); \draw[->] (dT) -- (T);
\draw[->] (rz) -- node[above]{da} (kraj);
\draw[->] (T.west) -- ++(-26mm,0) |- (mod.west);
\end{tikzpicture}"""

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
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "jednacine_p06")
    os.makedirs(outdir, exist_ok=True)
    meta = {name: render(name, r"$\displaystyle %s$" % tex, tex, outdir) for name, tex in EQ.items()}
    meta["p06_sema"] = render("p06_sema", SEMA, SEMA, outdir)
    with open(os.path.join(outdir, "jednacine_p06.json"), "w") as f:
        json.dump(meta, f, indent=1)
    for k, v in meta.items():
        print(f"{k:28s} {v['w_in']:.2f} x {v['h_in']:.2f} in (prirodno)")
