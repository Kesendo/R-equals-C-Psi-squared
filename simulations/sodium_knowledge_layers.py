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

# three views of the same gels (Fig. 1A, 1C and 2D, read off by eye; the 0 % reference first, one 3 % point):
# water diffusion (1e-9 m^2/s), water 1H T1 (s), sodium diffusion (1e-9 m^2/s). Diffusion sees motion with no coupling
# (and the gel's obstruction), 1H T1 motion times a dipolar coupling fixed inside the molecule, 17O T1 and a motion
# times a quadrupolar coupling that can change.
D_WATER = {"0": 2.5, "3": 2.08, "6": 1.96, "9": 1.83, "12": 1.72, "15": 1.62}
T1_H1 = {"0": 3.0, "3": 2.6, "6": 2.22, "9": 2.2, "12": 1.9, "15": 1.75}
D_NA = {"0": 1.26, "3": 1.15, "6": 1.07, "9": 1.02, "12": 0.90, "15": 0.90}
PCT = {"3 % (1)": "3", "3 % (2)": "3", "6 %": "6", "9 %": "9", "12 %": "12", "15 %": "15"}

# the larger mobile population (Table S5, m = 1): mono-exponential T2 (ms), beside the common T1 of Table S6
POPS = {"3 % (1)": (29, 14), "3 % (2)": (33, 12), "6 %": (33, 11), "9 %": (45, 8), "12 %": (32, 12), "15 %": (39, 9)}  # Table S5: p3b, p2 (%)
T2_MOBILE1 = {"3 % (1)": 5.45, "3 % (2)": 6.9, "6 %": 3.55, "9 %": 2.9, "12 %": 2.55, "15 %": 1.85}
# the rest of the fit, for the line width: Table S5 m = 2 (T2, ms), Table S4 T2s and T2l (ms) of the slow population
T2_MOBILE2 = {"3 % (1)": 33, "3 % (2)": 35, "6 %": 26, "9 %": 21, "12 %": 19, "15 %": 18}
T2_SLOW = {"3 % (1)": (5.09, 40.28), "3 % (2)": (6.08, 45.05), "6 %": (3.51, 36.34), "9 %": (2.64, 29.26),
           "12 %": (2.14, 26.59), "15 %": (1.85, 21.82)}
# the 1D 23Na line width without any echo (Fig. 2B, FWHM in Hz, read off by eye; the 0 % reference first)
FWHM_NA = {"0": 7.3, "3": 12.5, "6": 17.0, "9": 17.9, "12": 21.7, "15": 21.8}

def central_rate_is_mean_longitudinal():
    """Exact, from the repo's sodium generator (simulations/sodium_two_mirrors.py): the central-transition coherence
    decays at 3(J1 + J2), the mean of the two longitudinal rates 6 J1 and 6 J2, at every J0, splitting wQ and offset w."""
    import sympy as sp
    from sodium_two_mirrors import generator
    J0, J1, J2, wQ, w, lam = sp.symbols("J0 J1 J2 wQ w lam")
    _, _, L = generator(sp.Rational(3, 2), J0, J1, J2, wQ, w)
    blk = [0 * 4 + 1, 1 * 4 + 2, 2 * 4 + 3]                     # coherence order +1: (3/2,1/2), (1/2,-1/2), (-1/2,-3/2)
    cp = (lam * sp.eye(3) - L.extract(blk, blk)).det()
    central = -3 * (J1 + J2) - sp.I * w
    pops = [4 * k + k for k in range(4)]
    pp = sp.factor((lam * sp.eye(4) - L.extract(pops, pops)).det())
    return sp.expand(cp.subs(lam, central)) == 0 and sp.expand(pp - lam * (lam + 6 * J1) * (lam + 6 * J2) * (lam + 6 * J1 + 6 * J2)) == 0

def central_blind_to_splitting():
    """Exact, on the same generator: the central coherence is coupled to no other element, its entry is
    -3(J1 + J2) - i w with no wQ at all (the splitting T20 takes the same value on m = 1/2 and m = -1/2), and the hard
    pi pulse about x (the antidiagonal flip m -> -m, fixed up to phase by P Ix P^-1 = Ix, P Iy P^-1 = -Iy,
    P Iz P^-1 = -Iz) keeps T20 and turns Iz: an echo refocuses every static term odd in m and none even in m."""
    import sympy as sp
    from sodium_two_mirrors import generator, spin_ops
    J0, J1, J2, wQ, w = sp.symbols("J0 J1 J2 wQ w")
    I = sp.Rational(3, 2)
    _, _, L = generator(I, J0, J1, J2, wQ, w)
    c = 1 * 4 + 2                                                # |1/2><-1/2|
    alone = all(L[c, j] == 0 and L[j, c] == 0 for j in range(16) if j != c)
    entry = sp.expand(L[c, c] - (-3 * (J1 + J2) - sp.I * w)) == 0
    ms, Iz, T = spin_ops(I)
    Ip = sp.zeros(4)
    for k in range(1, 4): Ip[k - 1, k] = sp.sqrt(I * (I + 1) - ms[k] * (ms[k] + 1))
    Ix, Iy = (Ip + Ip.T) / 2, (Ip - Ip.T) / (2 * sp.I)
    P = sp.Matrix(4, 4, lambda a, b: sp.I if a + b == 3 else 0)
    Pi = P.inv()
    flip = all(sp.simplify(P * A * Pi - s * A) == sp.zeros(4) for A, s in [(Ix, 1), (Iy, -1), (Iz, -1), (T[0], 1)])
    return alone and entry and flip

def echo_on_the_letters():
    """The same echo read by the repo's letter instruments (framework.letter_cube, the Theodolite's letter layer
    framework.diagnostics.cube_letters), in the sodium page's labelling: bit a = outer/inner, bit b = sign of m, so
    Iz = Z_b + Z_a Z_b / 2 and the first-order splitting is Z_a. The hard pi pulse m -> -m flips bit b alone; it
    reverses the offset and keeps the splitting, no plain site-local monomial sandwich reverses both, and Z_a moves no
    string of the central transition (a kept, b flipped) and every string of the satellites (a flipped)."""
    import sympy as sp
    from framework import cube_tables as ct, letter_cube as lc
    from framework.diagnostics import cube_letters as cl
    h = sp.Rational(1, 2)
    level = {(0, 0): 3 * h, (1, 0): h, (1, 1): -h, (0, 1): -3 * h}       # (a, b) -> m, Z = +1 on bit value 0
    z = lambda bit: 1 - 2 * bit
    labels = all(z(b) + h * z(a) * z(b) == m and 3 * m * m - 15 * h * h == 3 * z(a) for (a, b), m in level.items())
    offset = ct.Setting(2, [('Z', (1,), 1), ('ZZ', (0, 1), h)], [])
    split = ct.Setting(2, [('Z', (0,), 1)], [])
    both = ct.Setting(2, offset.terms + split.terms, [])
    pulse = [cl.site_map('d', 1, 'd', 1, False), cl.site_map('a', 1, 'a', 1, False)]
    echo = cl.verify(offset, 'pi', pulse) and cl.verify(split, 'copy', pulse) and not cl.verify(split, 'pi', pulse)
    # plain sandwiches A·rho·B hold the pulses: none reverses both; the transposed branch finds rho -> rho^T itself,
    # which reverses every real diagonal H and is no pulse
    t = cl.read(both, 'pi', transposed=True)
    none_for_both = (cl.read(both, 'pi', transposed=False) is None and t is not None
                     and t.kinds == [('d', 'd'), ('d', 'd')] and t.x == [1, 1] and t.y == [1, 1])
    central, satellite = ['IX', 'IY', 'ZX', 'ZY'], ['XI', 'YI', 'XZ', 'YZ']
    blind = all(lc.commutator_coefficient('ZI', s) is None for s in central)
    moved = all(lc.commutator_coefficient('ZI', s) is not None for s in satellite)
    return labels and echo and none_for_both and blind and moved

def fwhm_of(comps, lor, sigma=0.0):
    """FWHM (Hz) of a sum of lines centred at zero, comps = [(weight, T2 in ms)]: each a Lorentzian of width
    1/(pi T2) + lor, convolved with a Gaussian of standard deviation sigma (a Voigt profile) when sigma > 0."""
    from scipy.special import voigt_profile
    f = np.linspace(-1000, 1000, 400001)
    s = sum(p * voigt_profile(f, sigma, (1e3 / (math.pi * t2) + lor) / 2) if sigma > 0 else
            p * (g / 2) / (f ** 2 + (g / 2) ** 2) / math.pi
            for p, t2, g in ((p, t2, 1e3 / (math.pi * t2) + lor) for p, t2 in comps))
    above = f[s >= s.max() / 2]
    return above[-1] - above[0]

def tissue(t1, t2s, t2f):            # ms -> (J0, J1, J2): R2f = J0 + J1, R2s = J1 + J2, R1 = 0.4 J1 + 1.6 J2
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
    print("three views (ratio against 0 %): water D, 1H, 17O, 17O/D, 1H/D, 17O/1H | sodium D, a (vs free ion), a/D")
    for n, J in AGAROSE.items():
        c = PCT[n]; dw = D_WATER["0"] / D_WATER[c]; h = T1_H1["0"] / T1_H1[c]; o = 6.86 / T1_O17[n]
        dn = D_NA["0"] / D_NA[c]; aa = two_motion(*J)[0] / FREE
        print(f"{n:9s} {dw:4.2f} {h:4.2f} {o:4.2f}  {o / dw:4.2f} {h / dw:4.2f} {o / h:4.2f} | {dn:4.2f} {aa:4.2f} {aa / dn:4.2f}")
    assert central_rate_is_mean_longitudinal()
    print("central transition decays at the mean of the two longitudinal rates (exact, every J0, wQ, w)")
    page_ratio = [7.1, 6.1, 9.9, 10.4, 10.3, 12.3]; page_share = [25.6, 25.2, 24.2, 26.0, 24.8, 24.6]
    for k, n in enumerate(AGAROSE):
        r = T1_MOBILE[n] / T2_MOBILE1[n]; p3b, p2 = POPS[n]; share = 0.4 * p3b + p2
        assert round(r, 1) == page_ratio[k] and round(share, 1) == page_share[k], n
        print(f"{n:9s} large mobile population: T1/T2 = {r:4.1f}; share decaying near T1 in the fit 0.4*p3b + p2 = {share:4.1f} %")
    # a single fitted longitudinal rate ~ 0.4 J1 + 1.6 J2 against the central rate J1 + J2: at most 2.5, as J2 -> 0
    assert max((1.0 + t) / (0.4 + 1.6 * t) for t in np.linspace(0, 10, 10001)) == 2.5
    # the fit's readout: a spectral point at the line's centre weighs each component by amplitude times its T2; the
    # authors' inverse Laplace transforms at the shortest filter (Figs. S16 simulated, S5 measured, peaks below 7 ms)
    S16 = [0.31, 0.34, 0.26, 0.22, 0.23, 0.20]; S5 = [0.31, 0.33, 0.24, 0.17 + 0.08, 0.21, 0.19]
    for k, n in enumerate(AGAROSE):
        p3b, p2 = POPS[n]; p1 = 100 - p3b - p2; t2s, t2l = T2_SLOW[n]
        fast, slow = [(0.6 * p3b, t2s), (p1, T2_MOBILE1[n])], [(0.4 * p3b, t2l), (p2, T2_MOBILE2[n])]
        amp = sum(p for p, _ in fast) / 100
        w = sum(p * t for p, t in fast) / sum(p * t for p, t in fast + slow)
        print(f"{n:9s} fast share: amplitude {amp:4.2f}, weighted by T2 {w:4.2f}; Fig. S16 {S16[k]:4.2f}, Fig. S5 {S5[k]:4.2f}")
        # the printed integrals carry two decimals and an inverse Laplace transform's smearing; what is asserted is
        # the separation: the T2-weighted share sits ten times closer to Fig. S16 than the amplitude share does
        assert abs(w - S16[k]) * 10 < abs(amp - S16[k]), n
    assert central_blind_to_splitting()
    print("central coherence alone, entry -3(J1+J2) - i w with no wQ; the pi pulse keeps T20 and turns Iz (exact)")
    assert echo_on_the_letters()
    print("on the letters: the pulse flips bit b, the splitting Z_a sits on bit a; no site-local pulse reverses both")
    # the 1D line without an echo against the line the echo fit's populations predict, a reading and no gate: the
    # difference moves by several Hz with the shape assumed for the instrument's part (Lorentzian, Gaussian, none),
    # each matched to the 0 % reference (T2 = 56.3 ms, Table S7), so the line width decides nothing about a static spread
    from scipy.optimize import brentq
    ref = 56.3
    lor = FWHM_NA["0"] - 1e3 / (math.pi * ref)
    sig = brentq(lambda g: fwhm_of([(1, ref)], 0, g) - FWHM_NA["0"], 1e-3, 5)
    print(f"1D line width (Hz), measured minus predicted; instrument as Lorentzian {lor:.2f}, "
          f"as Gaussian {2 * math.sqrt(2 * math.log(2)) * sig:.2f}, or none")
    for n in AGAROSE:
        p3b, p2 = POPS[n]; p1 = 100 - p3b - p2; t2s, t2l = T2_SLOW[n]
        comps = [(0.6 * p3b, t2s), (0.4 * p3b, t2l), (p1, T2_MOBILE1[n]), (p2, T2_MOBILE2[n])]
        meas = FWHM_NA[PCT[n]]
        for label, cs in (("fit's amplitudes", comps), ("amplitudes ~ p*T2", [(p * t, t) for p, t in comps])):
            d = [meas - fwhm_of(cs, lor), meas - fwhm_of(cs, 0, sig), meas - fwhm_of(cs, 0)]
            print(f"{n:9s} Fig. 2B {meas:5.1f}  {label:18s} {d[0]:+5.1f} {d[1]:+5.1f} {d[2]:+5.1f}")
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
