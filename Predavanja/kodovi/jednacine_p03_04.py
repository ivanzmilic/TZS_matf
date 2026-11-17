"""Jednačine za Predavanje 03_04 (Neprozračnost), renderovane LaTeX-om u PNG.

Konvencija: frekvencija nu (ne lambda); fizički fluks je F,
astrofizički \\mathcal{F} = F/pi, Eddingtonov H = F/(4 pi); emisivnost eta_nu.

Upotreba:  python jednacine_p03_04.py [izlazni_direktorijum]
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
    # Konvencija (03_04): E_i = energija nivoa (E_1 = -13.6 eV), Delta E_i = E_i - E_1 = energija ekscitacije,
    # E_ion = energija jonizacije; sigma = presek po čestici (bf, ff); i = donji, k = gornji nivo;
    # j = stepen jonizacije (Saha). chi je SAMO neprozračnost.
    # slajd 2
    "s02_formalno": r"I_\nu^+ = I_\nu^0\,e^{-\tau_\nu} + \int_0^{\tau_\nu} S_\nu(t)\,e^{-t}\,dt",
    # slajd 10
    "s10_tau_ref": r"\tau(z) = \int_z^{z_{\max}} \chi_0(z')\,dz'",
    # slajd 14: Borov model
    "s14_E_bohr": r"E_i = -\frac{Ze^2}{2r_i} < 0",
    "s14_bohr": r"""\begin{aligned}
m v_i r_i &= i\,\frac{h}{2\pi}\\
\frac{m v_i^2}{r_i} &= \frac{Ze^2}{r_i^2}
\end{aligned}""",
    "s14_E_C": r"E_i = -\frac{C}{i^2}",
    "s14_E1": r"""\begin{aligned}
E_1 &= -13.6\,\mathrm{eV}\\
E_\infty &= 0
\end{aligned}""",
    "s14_dE": r"\Delta E_i = E_i - E_1 = C\left(1-\frac{1}{i^2}\right) > 0",
    "s14_dE_box": r"\Delta E_i = E_{\rm ion}\left(1-\frac{1}{i^2}\right)",
    "s14_Eion": r"""\begin{aligned}
\Delta E_1 &= 0\\
E_{\rm ion} &= C = -E_1 = 13.6\,\mathrm{eV}
\end{aligned}""",
    "s14_tabela": r"""\begin{array}{l}
\text{Granice serija}\\[10pt]
\begin{array}{l@{\hspace{3.2em}}c@{\hspace{3.2em}}r}
\text{Serija} & i & \lambda_g\,[\mathrm{nm}]\\[5pt]
\text{Lyman} & 1 & 91.2\\[4pt]
\text{Balmer} & 2 & 364.6\\[4pt]
\text{Paschen} & 3 & 820.4\\[4pt]
\text{Brackett} & 4 & 1459.0\\[4pt]
\text{Pfund} & 5 & 2279.0
\end{array}
\end{array}""",
    "s14_nu_g": r"h\nu_{g,i} = E_{\rm ion} - \Delta E_i = E_{\rm ion} - E_{\rm ion}\left(1-\frac{1}{i^2}\right) = \frac{E_{\rm ion}}{i^2}",
    "s14_nu_g_box": r"\nu_{g,i} = \frac{E_{\rm ion}}{h\,i^2}",
    "oznaka_dE": r"\Delta E",
    # slajd 15
    "s15_dE_osa": r"\Delta E\ [\mathrm{eV}]",
    "s15_prelaz": r"h\nu = E_k - E_i",
    "s15_rydberg": r"h\nu_{ik} = E_{\rm ion}\left(\frac{1}{i^2} - \frac{1}{k^2}\right)",
    "s15_bf": r"h\nu = (E_{\rm ion} - \Delta E_i) + \frac{1}{2}mv^2",
    # slajd 16
    "s16_kramers": r"\sigma_{i,\nu}^{\rm bf} = C\,g_{\rm bf}\,\frac{1}{\nu^3 i^5},\quad (\nu \ge \nu_{g,i})",
    "s16_prag": r"\sigma_{i,\nu_{g,i}}^{\rm bf} \propto i",
    "s16_ni_sigma": r"n_i\,\sigma_{i,\nu}^{\rm bf} = C\,\frac{g_{\rm bf}}{\nu^3 i^5}\,\frac{n_1}{g_1}\,g_i\,e^{-\Delta E_i/kT}",
    "s16_chi_bf": r"\chi_\nu^{\rm bf} = \sum_{i=i_g}^{\infty} n_i\,\sigma_{i,\nu}^{\rm bf} = C\,\frac{n_1}{g_1}\,\frac{1}{\nu^3}\sum_{i=i_g}^{\infty}\frac{g_i}{i^5}\,g_{\rm bf}\,e^{-\Delta E_i/kT}",
    "s16_chi_f": r"\chi_\nu^{\rm bf} = f(n_e, n^+, T, \nu)",
    # slajd 17
    "s17_chi_bf": r"\chi_\nu^{\rm bf} = \sum_{i=i_g}^{\infty} n_i\,\sigma_{i,\nu}^{\rm bf}",
    # slajd 18
    "s18_dsigma": r"d\sigma_\nu^{\rm ff} = C''\,Z_{\rm ef}^2\,g_{\rm ff}\,\frac{1}{v}\,\frac{1}{\nu^3}\,dv \qquad\qquad (v,\ v+dv)",
    "s18_sigma_ff": r"\sigma_\nu^{\rm ff} = C''\,g_{\rm ff}\,\frac{Z_{\rm ef}^2}{\nu^3}\int\frac{1}{v}\,f(v)\,dv = C^*\,\frac{1}{\nu^3}\,\frac{1}{\sqrt{T}}\,g_{\rm ff}",
    "s18_chi_ff": r"\chi_\nu^{\rm ff} = C^*\,\frac{n_e n^+}{\sqrt{T}}\,g_{\rm ff}\,\frac{1}{\nu^3} = f(n_e, n^+, T, \nu)",
    # slajd 19
    "s19_Hminus": r"\chi_\nu(H^-) = n_{H^-}\,\sigma_{\nu,H^-}^{\rm bf} + n_{H0}\,P_e\,\sigma_{\nu,H^-}^{\rm ff} = f(T, n_e, n_{H0}, \nu)",
    "oznaka_sigma": r"\sigma",
    # slajdovi 21, 22
    "s21_HI": r"\chi_\nu^{\rm bf}(H\,I) = \sum_i n_i\,\sigma_{i,\nu}^{\rm bf}",
    "s21_Hminus": r"\chi_\nu(H^-) = n_{H^-}\,\sigma_{\nu,H^-}^{\rm bf} + n_{H0}\,P_e\,\sigma_{\nu,H^-}^{\rm ff}",
    # slajd 23
    "s23_boltzmann": r"\frac{n_k}{n_i} = \frac{g_k}{g_i}\,e^{-(E_k - E_i)/kT}",
    "s23_saha": r"\frac{n_{j+1}}{n_j} = \frac{1}{n_e}\,\Phi(T) = \frac{1}{n_e}\left(\frac{2\pi m_e kT}{h^2}\right)^{3/2}\frac{2U_{j+1}}{U_j}\,e^{-E_{{\rm ion},j}/kT}",
    # slajd 29
    "s29_vodonik": r"""\begin{aligned}
p &= (n_e + n_p + n_{H0})\,kT\\
n_e &= n_p\\
\frac{n_p n_e}{n_{H0}} &= \left(\frac{2\pi m_e kT}{h^2}\right)^{3/2}\frac{2g_p}{g_{H0}}\,e^{-E_{\rm ion}/kT}
\end{aligned}""",
    # slajd 31
    "s31_boltzmann": r"n_i = n_{H0}\,\frac{g_i\,e^{-\Delta E_i/kT}}{U_{H0}}",
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
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "jednacine_p03_04")
    os.makedirs(outdir, exist_ok=True)
    meta = {name: render(name, tex, outdir) for name, tex in EQ.items()}
    with open(os.path.join(outdir, "jednacine_p03_04.json"), "w") as f:
        json.dump(meta, f, indent=1)
    for k, v in meta.items():
        print(f"{k:24s} {v['w_in']:.2f} x {v['h_in']:.2f} in (prirodno)")
