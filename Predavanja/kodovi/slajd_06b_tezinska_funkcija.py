"""Slika za prelazni slajd posle definicije fluksa (Predavanje 05):
gore  - Sunčev disk sa tri položaja (mu = 1, 0.5, 0.2),
dole  - težinska funkcija formalnog rešenja na površini, e^{-t/mu}/mu,
        čija je srednja vrednost <tau> = mu  (=> Eddington-Barbier).

Upotreba: python slajd_06b_tezinska_funkcija.py [izlazni.png]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"mathtext.fontset": "cm", "font.size": 13,
                     "axes.spines.top": False, "axes.spines.right": False})

MUS = [1.0, 0.5, 0.2]
COLS = ["#1b4f9c", "#d17a00", "#b2182b"]

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(4.6, 5.6), gridspec_kw={"height_ratios": [1.05, 1]})

# --- disk sa potamnjenjem ka rubu (prosti linearni zakon, samo ilustracija)
n = 600
x = np.linspace(-1, 1, n)
X, Y = np.meshgrid(x, x)
R2 = X**2 + Y**2
mu = np.sqrt(np.clip(1 - R2, 0, 1))
img = np.where(R2 <= 1, 1 - 0.6 * (1 - mu), np.nan)
ax0.imshow(img, extent=(-1, 1, -1, 1), cmap="afmhot", vmin=-0.15, vmax=1.25, origin="lower")
ANG = {1.0: 0.0, 0.5: np.deg2rad(140), 0.2: np.deg2rad(25)}   # smer na disku, samo radi preglednosti
LAB = {1.0: (0.0, 0.42), 0.5: (-0.95, 0.95), 0.2: (1.05, 0.78)}
for m, c in zip(MUS, COLS):
    r = np.sqrt(1 - m**2)                       # r/R = sin(theta)
    px, py = r * np.cos(ANG[m]), r * np.sin(ANG[m])
    ax0.plot(px, py, "o", ms=9, mfc=c, mec="white", mew=1.5, zorder=3)
    ax0.annotate(rf"$\mu={m:g}$", (px, py), xytext=LAB[m], color=c, ha="center", fontsize=13,
                 arrowprops=dict(arrowstyle="-", color=c, lw=1))
ax0.text(0, -1.17, r"centar: $\mu=1$  $\longrightarrow$  rub: $\mu\to 0$", ha="center", fontsize=12)
ax0.set_xlim(-1.25, 1.25); ax0.set_ylim(-1.3, 1.12); ax0.set_aspect("equal"); ax0.axis("off")

# --- težinska funkcija
t = np.linspace(0, 3, 600)
for m, c in zip(MUS, COLS):
    ax1.plot(t, np.exp(-t / m) / m, color=c, lw=2.2, label=rf"$\mu={m:g}$")
    ax1.axvline(m, color=c, ls="--", lw=1)
ax1.set_xlim(0, 3); ax1.set_ylim(0, 5.3)
ax1.set_xlabel(r"$t = \tau_\nu$")
ax1.set_ylabel(r"$e^{-t/\mu}/\mu$")
ax1.text(1.25, 2.1, r"$\langle\tau_\nu\rangle=\mu$", fontsize=14)
ax1.text(1.25, 1.55, "(isprekidane linije)", fontsize=10, color="0.35")
ax1.legend(frameon=False, loc="upper right", fontsize=11)

fig.tight_layout(h_pad=0.6)
out = sys.argv[1] if len(sys.argv) > 1 else "slajd_06b_tezinska_funkcija.png"
fig.savefig(out, dpi=300, transparent=False, facecolor="white")
print("saved", out)
