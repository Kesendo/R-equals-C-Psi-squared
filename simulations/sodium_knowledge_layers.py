"""Sodium in agarose, layer by layer: what we know about a sample stacked as transparent filters.

The plane is the motion of one sodium species: horizontally w0*tau of a slow motion, vertically a / J2, the share of
J2 carried by a fast part that is flat in frequency (extreme narrowing; a <= J2 since b >= 0). In that two-motion model
(the two-step picture) J_m = a + b / (1 + m^2 (w0 tau)^2) with b = J0 - a, so every point of the plane predicts J1 and J2, and every piece of knowledge becomes a layer that lets
light through where the point agrees with it:
  J1, J2   the measured spectral densities (a band of +-10 % relative, our choice: the source gives no error bars);
  tau      an independent rotational correlation time, 12 +- 2 ns in 2.4 % agarose (the width is our choice);
  free ion a >= the flat density of free Na+ in the same apparatus, 1/(2 * 56 ms) (an assumption: the gel slows water).
The bottom edge a = 0 is the single-motion model, the one Lorentzian of experiments/THE_NUCLEUS_AS_A_KEY.md.
Light is the product of the layers' transmissions, so the brightness-only panel carries everything; the hue only says
which layer let light through, and the palette is the viewer's choice. The script also prints the w0*tau interval that
J1 and J2 alone allow (product of their bands above 1/2): within +-10 % they bound tau from below only.

Data: agarose, Nimerovsky, Sieme, Rezaei-Ghaleh, Methods 228 (2024) 55, Supplementary Table S4 (J in s^-1, 23Na at
105.84 MHz); tissue, the review 23Na relaxometry (PMC12406536) Table 1, T1 read as the weighted mean 0.4 J1 + 1.6 J2.
Writes visualizations/sodium_knowledge_layers_{neon,okabe,brightness}.png.
"""
import math, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "visualizations")

NU_NA = 105.84e6                      # 23Na Larmor frequency of the agarose study, Hz
W0 = 2 * math.pi * NU_NA
FREE = 1 / (2 * 0.056)                # free Na+ in the same apparatus: T1 = T2 = 56 ms, flat J = R/2

AGAROSE = {"3 % (1)": (183, 13.3, 11.5), "3 % (2)": (153, 11.8, 10.4), "6 %": (270, 15.1, 12.5),
           "9 %": (359, 19.1, 15.1), "12 %": (447, 21.6, 16.1), "15 %": (515, 25.8, 20.1)}

# the same samples' other layers: the T1 the authors fit in common to both mobile populations (Table S6, ms; mono-
# exponential, so J1 = J2 = 1/(2 T1), the density at w0 and 2w0, not a flat density: the larger mobile population has
# T2 << T1) and the water's 17O T1 (Fig. 1E, ms, read off the figure by eye, one 3 % point; the reference is 6.86 ms)
T1_MOBILE = {"3 % (1)": 38.6, "3 % (2)": 42.2, "6 %": 35.2, "9 %": 30.2, "12 %": 26.2, "15 %": 22.8}
T1_O17 = {"3 % (1)": 5.3, "3 % (2)": 5.3, "6 %": 4.6, "9 %": 4.05, "12 %": 3.55, "15 %": 3.1}

def tissue(t1, t2s, t2f):             # ms -> (J0, J1, J2): R2f = J0 + J1, R2s = J1 + J2, R1 = 0.4 J1 + 1.6 J2
    r1, r2s, r2f = 1e3 / t1, 1e3 / t2s, 1e3 / t2f
    j1 = (1.6 * r2s - r1) / 1.2
    return (r2f - j1, j1, r2s - j1)
TISSUE = {"cartilage, 3 T": tissue(14.45, 12.62, 0.348), "calf muscle, 1.5 T": tissue(32.3, 12.27, 0.46)}

def two_motion(J0, J1, J2):
    """The closed reading: the (a, b, x) with J_m = a + b/(1 + m^2 x) for m = 0, 1, 2, or None."""
    f = lambda x: (J0 - J1) * (1 - 1 / (1 + 4 * x)) - (J0 - J2) * (1 - 1 / (1 + x))
    lo, hi = 1e-9, 1e9
    if f(lo) * f(hi) > 0: return None
    for _ in range(400):
        mid = math.sqrt(lo * hi)
        if f(lo) * f(mid) <= 0: hi = mid
        else: lo = mid
    x = math.sqrt(lo * hi); b = (J0 - J1) / (1 - 1 / (1 + x)); a = J0 - b
    return (a, b, x) if a >= 0 and b > 0 else None

LX = np.linspace(math.log10(0.05), math.log10(200), 520)       # log10(w0 tau)
FY = np.linspace(0, 1, 300)                                    # a / J2
X, F = np.meshgrid(10 ** LX, FY)

def band(pred, meas, rel=0.10):
    return np.exp(-0.5 * ((pred - meas) / (rel * abs(meas))) ** 2)

def layers(J, agarose):
    J0, J1, J2 = J
    a = F * max(J2, 0.0); b = J0 - a; x = X ** 2
    out = [("J₁", band(a + b / (1 + x), J1)), ("J₂", band(a + b / (1 + 4 * x), J2))]
    if J2 <= 0:                                  # no positive density gives a negative J2: nothing passes
        out = [(n, np.zeros_like(X)) for n, _ in out]
    if agarose:
        tau = X / W0
        out.append(("τ = 12 ns", np.exp(-0.5 * ((tau - 12e-9) / 2e-9) ** 2)))
        out.append(("free ion", 1 / (1 + np.exp(-(a - FREE) / 0.4))))
    return out

PALETTES = {
    "neon": dict(bg="#080b12", text="#7ef9ff", light=(1, 1, 1),
                 tint={"J₁": "#08f7fe", "J₂": "#fe53bb", "τ = 12 ns": "#ffe700", "free ion": "#b14aed"}),
    "okabe": dict(bg="#14181d", text="#e8edf2", light=(1, 1, 1),
                  tint={"J₁": "#0072B2", "J₂": "#E69F00", "τ = 12 ns": "#009E73", "free ion": "#CC79A7"}),
    "brightness": dict(bg="#101010", text="#e0e0e0", light=(1, 1, 1), tint=None),
}

def rgb(h): return np.array([int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])

def render(ls, P):
    every = np.prod([t for _, t in ls], axis=0)[..., None]    # the light that passes every layer
    if P["tint"] is None:                                     # brightness only: nothing but that product
        return np.repeat(every, 3, axis=2)
    img = np.zeros(X.shape + (3,))
    for name, t in ls:                                        # each layer's own window, faint, in its colour (screen)
        img = 1 - (1 - img) * (1 - 0.45 * rgb(P["tint"][name]) * t[..., None])
    return np.clip(img * (1 - every) + every, 0, 1)           # where all pass: white

def check():
    for n, J in AGAROSE.items():
        sol = two_motion(*J)
        assert sol is not None, n
        a, b, x = sol
        tau = math.sqrt(x) / W0
        assert 11e-9 < tau < 14e-9 and 0 < a < J[2], (n, tau, a)
        ls = layers(J, True)
        light = np.prod([t for _, t in ls], axis=0)
        iy, ix = np.unravel_index(np.argmax(light), light.shape)
        # the brightest pixel sits on the closed reading, within the grid spacing
        assert abs(LX[ix] - math.log10(math.sqrt(x))) < 0.05 and abs(FY[iy] - a / J[2]) < 0.05, n
        # the single-motion edge (a = 0) never lets J1 and J2 through together
        assert (ls[0][1][0] * ls[1][1][0]).max() < 1e-6, n
        p = ls[0][1] * ls[1][1]; ok = np.where(p.max(axis=0) > 0.5)[0]
        lo, hi = 10 ** LX[ok[0]], 10 ** LX[ok[-1]]
        assert 8e-9 < lo / W0 < 10e-9, n           # J1 and J2 alone: tau above about 9 ns, the upper end loose
        print(f"   J1 and J2 alone: w0tau {lo:5.2f} to {hi:6.1f}, tau {lo / W0 * 1e9:4.1f} to {hi / W0 * 1e9:6.1f} ns")
        print(f"{n:9s} J0={J[0]:4d} J1={J[1]:5.1f} J2={J[2]:5.1f}  two motions: a={a:5.2f} b={b:6.1f} "
              f"w0tau={math.sqrt(x):5.2f} tau={tau * 1e9:5.2f} ns; edge a=0 dark")
    ratios = []
    for n, J in AGAROSE.items():
        a = two_motion(*J)[0]; jm = 1e3 / (2 * T1_MOBILE[n]); ratios.append(a / jm)
        print(f"{n:9s} a/free={a / FREE:4.2f}  water slowing={6.86 / T1_O17[n]:4.2f}  a/mobile flat={a / jm:5.3f}")
    assert sum(0.815 < r < 0.845 for r in ratios) == 5, ratios
    for n, J in TISSUE.items():
        print(f"{n:18s} J0={J[0]:7.1f} J1={J[1]:6.1f} J2={J[2]:6.1f}  two motions: {two_motion(*J)}")

def draw(pal):
    P = PALETTES[pal]
    rows = list(AGAROSE.items()) + list(TISSUE.items())
    fig, axes = plt.subplots(2, 4, figsize=(22, 9.6), facecolor=P["bg"])
    for ax, (name, J) in zip(axes.flat, rows):
        ag = name in AGAROSE
        ax.imshow(render(layers(J, ag), P), origin="lower", aspect="auto",
                  extent=[LX[0], LX[-1], FY[0], FY[-1]], interpolation="bilinear")
        ax.set_facecolor(P["bg"])
        for s in ax.spines.values(): s.set_color(P["text"])
        ax.tick_params(colors=P["text"], labelsize=9)
        ax.set_xticks([-1, 0, 1, 2]); ax.set_xticklabels(["0.1", "1", "10", "100"])
        ax.set_title(f"{'agarose ' if ag else ''}{name}\nJ₀, J₁, J₂ = {J[0]:.0f}, {J[1]:.1f}, {J[2]:.1f} s⁻¹",
                     color=P["text"], fontsize=11)
        ax.set_xlabel("ω₀τ of the slow motion", color=P["text"], fontsize=10)
        ax.set_ylabel("a / J₂ (flat part)", color=P["text"], fontsize=10)
        ax.text(LX[0] + 0.05, 0.02, "a = 0: one motion", color=P["text"], fontsize=8, va="bottom")
    fig.suptitle("Sodium, layer by layer: each thing we know lets light through where it agrees; "
                 "white = every layer passes", color=P["text"], fontsize=15)
    fig.text(0.30, 0.015, "layers: J₁, J₂ (every panel), τ = 12 ns and free ion (agarose only)",
             ha="center", color=P["text"], fontsize=11)
    if P["tint"] is not None:
        for i, (k, c) in enumerate(P["tint"].items()):
            fig.text(0.58 + 0.08 * i, 0.015, f"■ {k}", color=c, fontsize=12)
    else:
        fig.text(0.58, 0.015, "brightness only: the product of the layers, no hue", color=P["text"], fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    path = os.path.join(OUT, f"sodium_knowledge_layers_{pal}.png")
    fig.savefig(path, dpi=110, facecolor=P["bg"]); plt.close(fig)
    print("wrote", os.path.relpath(path, os.path.join(HERE, "..")))

if __name__ == "__main__":
    check()
    for p in PALETTES: draw(p)
