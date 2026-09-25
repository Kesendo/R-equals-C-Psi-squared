# -*- coding: utf-8 -*-
"""envelope_n2_rises.py

Do the successive local maxima of CΨ fall for every two-qubit state under local Z-dephasing?
Read literally, no, in every class. This producer computes why, and the structural question that
survives.

Book: CΨ = Tr(ρ²)·L₁(ρ)/(d − 1) with L₁ = Σ_{i≠j}|ρ_ij| in the computational basis (d = 4 at N = 2),
Lindblad dρ/dt = −i[H, ρ] + Σ_l γ(Z_l ρ Z_l − ρ) (jump √γ Z_l, so a one-site coherence decays at 2γ).
A local maximum is a grid sample above its left neighbour and not below its right one, refined by
golden-section search inside the two neighbouring grid steps. No prominence threshold is used anywhere:
a threshold would be an error code here.

(A) Local fields on a pure product state: H = X⊗I + 0.37·I⊗Y from |01⟩. The state stays a product, a
    product state has P = P₁P₂ and L₁ = (1 + ℓ₁)(1 + ℓ₂) − 1, and at γ = 0
    CΨ(t) = [(1 + |sin 2t|)(1 + |sin 0.74t|) − 1]/3. The two clocks are incommensurate; the MAIN peaks
    (smooth, nondegenerate maxima) rise, and continuity in γ carries the rise into damping.
(B) Micro-maxima under number-conserving H. Damping moves the zeros of ρ₀₁,₁₀, which follows the damped
    one-excitation clock, against the zeros of the (0,1)- and (1,2)-block coherences, which keep the
    undamped clock exactly (their dissipator is the uniform rate 2γ). Where two such corners of L₁ bracket
    a slope that changes sign, a micro-maximum can grow just above the trough, below the next main peak:
    a literal rise. In a pure state the corners coincide at γ = 0 and split; this needs an entry that
    passes through zero, which is not generic. For ψ₀ = a|00⟩ + b|01⟩ + c|11⟩ (a, b, c real) under XXX
    the whole trajectory is in closed form (xxx_pure_cpsi), with the corners exactly at nπ/4 and nπ/ω,
    ω = 2√(4 − γ²). In the mixed example (|i0⟩⟨i0| + |i+⟩⟨i+|)/2 under XXX at γ = 0.05 the two corners sit
    0.035 apart at γ = 0, and damping turns the fixed one, a kink on a rising slope, into a trough.
(C) The structural question. At γ = 0 a number-conserving H gives U(T) = a phase on each excitation
    sector, T = 2π/Ω the one-excitation block period, so CΨ(t + T) = CΨ(t) for EVERY state; for a pure
    state CΨ = [(c + √p + √(s − p))² − 1]/3 with p(t) = |ψ₀₁(t)|² one sinusoid, so every γ = 0 maximum
    is equal. The main peak of a period is the largest value of CΨ in it (windows anchored at the
    lowest point of the γ = 0 period); for pure states also the largest value between consecutive γ = 0
    minima (a half-period when p crosses s/2). The scans read whether damping ever makes a main peak
    rise.
(D) Beside the finite N ≥ 3 atlas (Symphony conventions: XY = Σ (J/2)(XX + YY), Heisenberg =
    Σ (J/4)(XX + YY + ZZ), Bell+ on the two leading sites): the carrier's occupied energies, its
    γ = 0 maxima, and the N = 3, Q = 2000 row on both chains with the atlas reader.
(E) The local-unitary example: |0⟩ ⊗ (I + sin 60°·X + cos 60°·Z)/2 under Z-dephasing alone crosses 1/4
    downward; an R_y on qubit 2 that turns its Bloch vector onto the equator lifts CΨ back above 1/4.

Run: python simulations/envelope_n2_rises.py   (about two minutes; writes
simulations/results/envelope_n2_rises.txt)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

RESULTS = Path(__file__).resolve().parent / "results" / "envelope_n2_rises.txt"

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
GOLD = (np.sqrt(5.0) - 1.0) / 2.0
H_XXX = np.kron(X, X) + np.kron(Y, Y) + np.kron(Z, Z)


# ---------------------------------------------------------------- primitives
def cpsi(rho):
    """CΨ = Tr(ρ²)·L₁/(d − 1)."""
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return float(np.real(np.trace(rho @ rho)) * l1 / (rho.shape[0] - 1))


def site_op(n, site, m):
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, m if k == site else I2)
    return out


def liouvillian(H, gamma):
    """Row-major vec (vec[a·d + b] = ρ[a, b]); jumps √γ Z_l on every site."""
    d = H.shape[0]
    n = int(round(np.log2(d)))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(n):
        Zk = site_op(n, k, Z)
        L = L + gamma * (np.kron(Zk, Zk.T) - np.eye(d * d))
    return L


def golden(f, a, b, sign=1.0, iterations=90):
    """Golden-section search for the maximum (sign = +1) or minimum (sign = −1) of a one-humped f."""
    g = lambda x: sign * f(x)
    c, d = b - GOLD * (b - a), a + GOLD * (b - a)
    gc, gd = g(c), g(d)
    for _ in range(iterations):
        if gc > gd:
            b, d, gd = d, c, gc
            c = b - GOLD * (b - a)
            gc = g(c)
        else:
            a, c, gc = c, d, gd
            d = a + GOLD * (b - a)
            gd = g(d)
    t = 0.5 * (a + b)
    return t, f(t)


class Trajectory:
    """ρ(t) = R e^{Λt} R⁻¹ ρ₀ from one eigendecomposition of L; CΨ(t) by calling it."""

    def __init__(self, H, gamma, rho0):
        lam, R = np.linalg.eig(liouvillian(H, gamma))
        self.lam, self.R, self.d = lam, R, H.shape[0]
        self.c0 = np.linalg.solve(R, rho0.reshape(-1).astype(complex))

    def rho(self, t):
        return (self.R @ (np.exp(self.lam * t) * self.c0)).reshape(self.d, self.d)

    def __call__(self, t):
        return cpsi(self.rho(t))


def local_maxima(f, tmax, dt, t0=0.0):
    """All local maxima of f on [t0, tmax], refined; no prominence threshold."""
    ts = np.arange(t0, tmax + 0.5 * dt, dt)
    vals = np.array([f(t) for t in ts])
    return [golden(f, ts[i - 1], ts[i + 1]) for i in range(1, len(ts) - 1)
            if vals[i] > vals[i - 1] and vals[i] >= vals[i + 1]], ts, vals


def local_minima(ts, vals, f):
    return [golden(f, ts[i - 1], ts[i + 1], sign=-1.0) for i in range(1, len(ts) - 1)
            if vals[i] < vals[i - 1] and vals[i] <= vals[i + 1]]


def rising_pairs(mx):
    return [(mx[k - 1], mx[k]) for k in range(1, len(mx)) if mx[k][1] > mx[k - 1][1]]


def ket(label):
    single = {"0": np.array([1, 0]), "1": np.array([0, 1]),
              "+": np.array([1, 1]) / np.sqrt(2), "-": np.array([1, -1]) / np.sqrt(2),
              "i": np.array([1, 1j]) / np.sqrt(2)}
    return np.kron(single[label[0]], single[label[1]]).astype(complex)


def fmt(mx, k):
    return ", ".join(f"{c:.6f} (t {t:.4f})" for t, c in mx[:k])


# ---------------------------------------------------------------- (A) local fields
H_A = np.kron(X, I2) + 0.37 * np.kron(I2, Y)


def closed_form_a(t):
    return ((1.0 + abs(np.sin(2.0 * t))) * (1.0 + abs(np.sin(0.74 * t))) - 1.0) / 3.0


def part_a():
    lines = ["(A) local fields, pure product state: H = X⊗I + 0.37·I⊗Y, ρ₀ = |01⟩⟨01|"]
    rho0 = np.outer(ket("01"), ket("01").conj())
    L0 = liouvillian(H_A, 0.0)
    v0 = rho0.reshape(-1)
    ts = np.linspace(0.0, 40.0, 4001)
    dev = max(abs(cpsi((expm(L0 * t) @ v0).reshape(4, 4)) - closed_form_a(t)) for t in ts)
    lines.append(f"  γ = 0 closed form [(1 + |sin 2t|)(1 + |sin 0.74t|) − 1]/3 against exact "
                 f"propagation, 4001 times in [0, 40]: max |Δ| = {dev:.2e}")
    t1, c1 = golden(closed_form_a, 0.0, np.pi / 2)
    t2, c2 = golden(closed_form_a, np.pi / 2, np.pi)
    lines.append(f"  closed form, first two maxima (golden-section on [0, π/2] and [π/2, π]): "
                 f"{c1:.6f} (t {t1:.4f}) → {c2:.6f} (t {t2:.4f}), rise {c2 - c1:+.6f}")
    for gamma in (0.0, 0.002, 0.05):
        L = liouvillian(H_A, gamma)
        f = lambda t, L=L: cpsi((expm(L * t) @ v0).reshape(4, 4))
        step = expm(L * 0.002)
        v = v0.astype(complex).copy()
        n = 20000
        vals = np.empty(n + 1)
        for i in range(n + 1):
            vals[i] = cpsi(v.reshape(4, 4))
            v = step @ v
        tg = 0.002 * np.arange(n + 1)
        mx = [golden(f, tg[i - 1], tg[i + 1]) for i in range(1, n)
              if vals[i] > vals[i - 1] and vals[i] >= vals[i + 1]]
        lines.append(f"  γ = {gamma:g}: {len(mx)} maxima in t < 40, first two {fmt(mx, 2)}; "
                     f"rising successive pairs {len(rising_pairs(mx))}")
    return lines


# ---------------------------------------------------------------- (B) micro-maxima
def xxx_pure_cpsi(a, b, c, gamma, t):
    """Exact CΨ(t) for ψ₀ = a|00⟩ + b|01⟩ + c|11⟩ (a, b, c ≥ 0) under H = XX + YY + ZZ, jumps √γ Z_l.

    |00⟩ and |11⟩ are eigenstates; the block (|01⟩, |10⟩) precesses at 4 about x with its coherence
    damped at 4γ, so r_z = s·e^{−2γt}[cos ωt + (2γ/ω) sin ωt] and r_y = −(4s/ω)·e^{−2γt} sin ωt,
    s = b², ω = 2√(4 − γ²), r_x ≡ 0. The (0,1)- and (1,2)-block coherences are e^{−2γt} times their
    undamped values, of modulus b·a|cos 2t|, b·a|sin 2t| and b·c|cos 2t|, b·c|sin 2t|; ρ₀₀,₁₁ has modulus
    e^{−4γt}·a·c."""
    s = b * b
    w = 2.0 * np.sqrt(4.0 - gamma * gamma)
    e2 = np.exp(-2.0 * gamma * t)
    rz = s * e2 * (np.cos(w * t) + (2.0 * gamma / w) * np.sin(w * t))
    ry = -(4.0 * s / w) * e2 * np.sin(w * t)
    l1 = 2.0 * (e2 * b * (a + c) * (abs(np.cos(2.0 * t)) + abs(np.sin(2.0 * t)))
                + e2 * e2 * a * c + abs(ry) / 2.0)
    purity = (a ** 4 + c ** 4 + (s * s + ry * ry + rz * rz) / 2.0
              + 2.0 * e2 * e2 * s * (a * a + c * c) + 2.0 * e2 ** 4 * a * a * c * c)
    return purity * l1 / 3.0


def pure_micro_example(label, a, b, c, gamma, n):
    """The n-th pair of corners (nπ/4, nπ/ω), the micro-maximum between them and the next main peak."""
    w = 2.0 * np.sqrt(4.0 - gamma * gamma)
    f = lambda t: xxx_pure_cpsi(a, b, c, gamma, t)
    c1, c2 = n * np.pi / 4.0, n * np.pi / w
    grid = np.linspace(c1, c2, 2001)[1:-1]
    vals = np.array([f(t) for t in grid])
    k = int(np.argmax(vals))
    tm, vm = golden(f, grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)])
    tn, vn = golden(f, c2, (n + 1) * np.pi / 4.0)
    psi = np.array([a, b, 0.0, c], dtype=complex)
    L = liouvillian(H_XXX, gamma)
    v0 = np.outer(psi, psi.conj()).reshape(-1)
    check = max(abs(cpsi((expm(L * t) @ v0).reshape(4, 4)) - f(t)) for t in (c1, c2, tm, tn))
    return [f"  {label}, γ = {gamma:g}: corners at {n}π/4 = {c1:.6f} (CΨ {f(c1):.9f}) and {n}π/ω = "
            f"{c2:.6f} (CΨ {f(c2):.9f}); micro-maximum {vm:.9f} at t {tm:.6f}, {vm - max(f(c1), f(c2)):.2e} "
            f"above the higher corner; next main peak {vn:.7f} at t {tn:.5f}: a rising successive pair; "
            f"the matrix exponential at these four times agrees with the closed form to {check:.1e}"]


def rho_mixed_b():
    a, b = ket("i0"), ket("i+")
    return 0.5 * (np.outer(a, a.conj()) + np.outer(b, b.conj()))


def part_b():
    lines = ["(B) micro-maxima under number-conserving H (literal counting)"]
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(5):
        v = np.abs(rng.normal(size=3))
        a, b, c = v / np.linalg.norm(v)
        gamma = rng.uniform(0.0, 0.6)
        psi = np.array([a, b, 0.0, c], dtype=complex)
        L = liouvillian(H_XXX, gamma)
        v0 = np.outer(psi, psi.conj()).reshape(-1)
        for t in np.linspace(0.0, 12.0, 97):
            worst = max(worst, abs(cpsi((expm(L * t) @ v0).reshape(4, 4)) - xxx_pure_cpsi(a, b, c, gamma, t)))
    lines.append(f"  pure a|00⟩ + b|01⟩ + c|11⟩ under XXX: closed form against exact propagation, 5 random "
                 f"(a, b, c, γ), 97 times each: max |Δ| = {worst:.2e}")
    lines += pure_micro_example("|0⟩ ⊗ R_y(π/4)|0⟩ = cos(π/8)|00⟩ + sin(π/8)|01⟩",
                                np.cos(np.pi / 8), np.sin(np.pi / 8), 0.0, 0.5, 1)
    f8 = lambda t: xxx_pure_cpsi(np.cos(np.pi / 8), np.sin(np.pi / 8), 0.0, 0.5, t)
    w8 = 2.0 * np.sqrt(4.0 - 0.25)
    count = 0
    for n in range(1, 13):
        c1, c2 = n * np.pi / 4.0, n * np.pi / w8
        grid = np.linspace(c1, c2, 2001)[1:-1]
        count += max(f8(t) for t in grid) > max(f8(c1), f8(c2))
    lines.append(f"    the same state: a micro-maximum between the corners of {count} of the first 12 corner pairs")
    v = np.array([0.61997, 0.75128, 0.22631])
    a, b, c = v / np.linalg.norm(v)
    lines += pure_micro_example("ψ₀ ∝ 0.61997|00⟩ + 0.75128|01⟩ + 0.22631|11⟩ (normalized)", a, b, c, 0.05, 5)

    rho0 = rho_mixed_b()
    tr = Trajectory(H_XXX, 0.05, rho0)
    lines.append(f"  mixed ρ₀ = (|i0⟩⟨i0| + |i+⟩⟨i+|)/2 under XXX, γ = 0.05: purity "
                 f"{np.real(np.trace(rho0 @ rho0)):.6f}, CΨ(0) = {cpsi(rho0):.6f}")
    mx, ts, vals = local_maxima(tr, 12.0, 5e-5)
    mn = local_minima(ts, vals, tr)
    rp = rising_pairs(mx)
    lines.append(f"    literal maxima on a dt = 5·10⁻⁵ grid in t < 12: {len(mx)}; first six {fmt(mx, 6)}")
    lines.append(f"    rising successive pairs {len(rp)}; the first: {rp[0][0][1]:.6f} (t {rp[0][0][0]:.4f}) → "
                 f"{rp[0][1][1]:.6f} (t {rp[0][1][0]:.4f})")
    tr0 = Trajectory(H_XXX, 0.0, rho0)
    z0 = brentq(lambda t: tr0.rho(t)[1, 2].imag, 2.150, 2.165, xtol=1e-15)
    z5 = brentq(lambda t: tr.rho(t)[1, 2].imag, 2.150, 2.159, xtol=1e-15)
    zf = (np.pi + np.arctan(3.0)) / 2.0
    h = 1e-6
    slopes = lambda f: ((f(zf) - f(zf - h)) / h, (f(zf + h) - f(zf)) / h)
    s0, s5 = slopes(tr0), slopes(tr)
    lines.append(f"    the first micro-maximum between two corners (ρ₀₁,₁₀ is imaginary here, so |ρ₀₁,₁₀| passes "
                 f"through zero): Im ρ₀₁,₁₀ = 0 at t {z0:.6f} = 3π/16 + π/2 at γ = 0 (the only minimum nearby, CΨ "
                 f"{tr0(z0):.6f}) and at t {z5:.6f} at γ = 0.05; ρ₀₀,₁₀ = 0 at (π + arctan 3)/2 = {zf:.6f} at both "
                 f"(|ρ₀₀,₁₀| there {abs(tr0.rho(zf)[0, 2]):.1e} and {abs(tr.rho(zf)[0, 2]):.1e})")
    lines.append(f"    at γ = 0 the fixed zero is a kink on a rising slope (slopes {s0[0]:+.3f} → {s0[1]:+.3f}); at "
                 f"γ = 0.05 it is a trough (slopes {s5[0]:+.3f} → {s5[1]:+.3f}, CΨ {tr(zf):.6f}); the micro-maximum "
                 f"sits {rp[0][0][1] - tr(z5):.2e} above the moved corner")
    heights, pairs = [], []
    for lower, _ in rp:
        left = max((m for m in mn if m[0] < lower[0]), key=lambda m: m[0])
        right = min((m for m in mn if m[0] > lower[0]), key=lambda m: m[0])
        heights.append(lower[1] - max(left[1], right[1]))
        pairs.append((left[0], right[0]))
    lines.append(f"    every rising pair starts at a micro-maximum between two minima: its height above the higher "
                 f"one is {min(heights):.1e} to {max(heights):.1e} (a main peak rises ~0.1 above its troughs); the "
                 f"minimum pairs sit {min(b - a for a, b in pairs):.3f} to {max(b - a for a, b in pairs):.3f} apart: "
                 + ", ".join(f"{a:.3f}/{b:.3f}" for a, b in pairs[:4]) + ", …")
    return lines


# ---------------------------------------------------------------- (C) main peaks
def hopping(delta):
    """(XX + YY) + δ(Z₁ − Z₂): hopping 2 and a staggered field; ZZ and a uniform field are invisible to CΨ."""
    return np.kron(X, X) + np.kron(Y, Y) + delta * (np.kron(Z, I2) - np.kron(I2, Z))


def block_clock(H, rho0):
    """γ = 0 one-excitation population p(t) = A + 2|z|cos(Ωt + φ); returns (A, |z|, φ, Ω, s)."""
    E, V = np.linalg.eigh(H[1:3, 1:3])
    Rb = rho0[1:3, 1:3]
    Rp = V.conj().T @ Rb @ V
    z = V[0, 0] * Rp[0, 1] * np.conj(V[0, 1])
    A = abs(V[0, 0]) ** 2 * Rp[0, 0].real + abs(V[0, 1]) ** 2 * Rp[1, 1].real
    return A, abs(z), np.angle(z), abs(E[1] - E[0]), np.trace(Rb).real


def window_maxima(f, bounds, points=400):
    """Largest value of f in each window; also whether it sits on the window's left or right edge."""
    out = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        ts = np.linspace(a, b, points)
        vals = np.array([f(t) for t in ts])
        k = int(np.argmax(vals))
        t, v = golden(f, ts[max(k - 1, 0)], ts[min(k + 1, points - 1)])
        out.append((t, v, "left" if k == 0 else ("right" if k == points - 1 else "inside")))
    return out


def period_windows(H, rho0, periods):
    """Windows of one block-clock period T = 2π/Ω, anchored at the lowest point of the γ = 0 period."""
    T = 2.0 * np.pi / block_clock(H, rho0)[3]
    f0 = Trajectory(H, 0.0, rho0)
    ts = np.linspace(0.0, T, 2001)
    k = int(np.argmin([f0(t) for t in ts]))
    t0, _ = golden(f0, ts[max(k - 1, 0)], ts[min(k + 1, 2000)], sign=-1.0)
    return [t0 + j * T for j in range(periods + 1)]


def hump_windows(H, psi, tmax):
    """For a pure state: windows between consecutive γ = 0 minima, the turning points of p(t) (all of
    them when p crosses s/2, every other one otherwise); each window holds one γ = 0 maximum."""
    rho0 = np.outer(psi, psi.conj())
    A, zabs, phi, Om, s = block_clock(H, rho0)
    turning = [(n * np.pi - phi) / Om for n in range(-2, int(tmax * Om / np.pi) + 3)]
    turning = [t for t in turning if 0.0 <= t <= tmax]
    if A - 2 * zabs < s / 2 < A + 2 * zabs:
        return turning
    keep_low = s / 2 >= A + 2 * zabs          # f increasing on the range of p: minima where p is lowest
    return [t for t in turning if (np.cos(Om * t + phi) < 0) == keep_low]


def named_single_block_kets():
    """Kets whose one-excitation part is one basis ket: p(0) ∈ {0, s}, so p crosses s/2 iff δ < 1
    (hopping 2, staggered splitting 4δ); δ ∈ {0, 0.5, 2, 3} keeps every case off that boundary."""
    labels = ["01", "10"] + [a + b for a in "01" for b in "+-i"] + [b + a for a in "01" for b in "+-i"]
    return {k: ket(k) for k in labels}


def part_c():
    lines = ["(C) main peaks under number-conserving H"]
    rng = np.random.default_rng(20260925)
    worst = 0.0
    for _ in range(5):
        A_ = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
        rho = A_ @ A_.conj().T
        rho /= np.trace(rho)
        H = hopping(rng.uniform(0.0, 2.5)) + rng.uniform(-1.0, 1.0) * np.kron(Z, Z)
        T = 2.0 * np.pi / block_clock(H, rho)[3]
        f = Trajectory(H, 0.0, rho)
        worst = max(worst, max(abs(f(t + T) - f(t)) for t in np.linspace(0.0, 3.0, 31)))
    lines.append(f"  γ = 0 periodicity: |CΨ(t + T) − CΨ(t)| for 5 random rank-2 mixed states under random "
                 f"hopping, staggered field and ZZ: max {worst:.1e}")
    rng_half = np.random.default_rng(5)
    for delta in (0.0, 0.5, 2.0):
        worst = 0.0
        for _ in range(4):
            A_ = rng_half.normal(size=(4, 2)) + 1j * rng_half.normal(size=(4, 2))
            rho = A_ @ A_.conj().T
            rho /= np.trace(rho)
            H = hopping(delta)
            T = 2.0 * np.pi / block_clock(H, rho)[3]
            f = Trajectory(H, 0.0, rho)
            worst = max(worst, max(abs(f(t + T / 2) - f(t)) for t in np.linspace(0.0, 3.0, 61)))
        lines.append(f"  γ = 0, staggered field δ = {delta:g}: |CΨ(t + T/2) − CΨ(t)| for 4 random rank-2 states: "
                     f"max {worst:.1e}" + (" (at δ = 0, U(T/2) is a diagonal phase times the qubit swap)"
                                           if delta == 0.0 else ""))
    worst_spread = worst_pred = 0.0
    for _ in range(6):
        H = hopping(rng.uniform(0.0, 2.5))
        v = rng.normal(size=4) + 1j * rng.normal(size=4)
        psi = v / np.linalg.norm(v)
        mx, _, _ = local_maxima(Trajectory(H, 0.0, np.outer(psi, psi.conj())), 12.0, 0.002)
        vals = np.array([c for _, c in mx])
        A, zabs, _, _, s = block_clock(H, np.outer(psi, psi.conj()))
        p = s / 2 if A - 2 * zabs <= s / 2 <= A + 2 * zabs else (A + 2 * zabs if A + 2 * zabs < s / 2 else A - 2 * zabs)
        cc = abs(psi[0]) + abs(psi[3])
        pred = ((cc + np.sqrt(p) + np.sqrt(s - p)) ** 2 - 1.0) / 3.0
        worst_spread = max(worst_spread, float(vals.max() - vals.min()))
        worst_pred = max(worst_pred, float(np.max(np.abs(vals - pred))))
    lines.append(f"  γ = 0 lemma, 6 random pure states: spread of the maxima {worst_spread:.1e}; largest "
                 f"|maximum − [(c + √p + √(s − p))² − 1]/3 at the p nearest s/2| {worst_pred:.1e}")
    psi = rng.normal(size=4) + 1j * rng.normal(size=4)
    psi /= np.linalg.norm(psi)
    collapse = 0.0
    for rho0 in (np.outer(psi, psi.conj()), rho_mixed_b()):
        ref = Trajectory(H_XXX, 0.05, rho0)
        for H in (hopping(0.0), hopping(0.0) + 0.5 * np.kron(Z, Z),
                  H_XXX + 0.4 * (np.kron(Z, I2) + np.kron(I2, Z))):
            cur = Trajectory(H, 0.05, rho0)
            collapse = max(collapse, max(abs(cur(t) - ref(t)) for t in np.linspace(0.0, 10.0, 201)))
    lines.append(f"  collapse: XY, XXZ Δ = 0.5 and XXX with a uniform field against XXX, one pure and one "
                 f"mixed start, γ = 0.05, 201 times in [0, 10]: max |Δ| = {collapse:.1e}")

    fb = rho_mixed_b()
    A, zabs, phi, Om, s = block_clock(H_XXX, fb)
    half = [t for t in ((n * np.pi - phi) / Om for n in range(-2, 12)) if 0.0 <= t <= 6.0]
    for gamma in (0.0, 0.05):
        wm = window_maxima(Trajectory(H_XXX, gamma, fb), period_windows(H_XXX, fb, 8))
        lines.append(f"  (B)'s main peaks per period, γ = {gamma:g}: " + ", ".join(f"{v:.6f}" for _, v, _ in wm))
        wh = window_maxima(Trajectory(H_XXX, gamma, fb), half)
        lines.append(f"  (B)'s main peaks per half-period (turning points of its block population; at δ = 0 already "
                     f"a period of CΨ), γ = {gamma:g}: " + ", ".join(f"{v:.6f}" for _, v, _ in wh))

    named = named_single_block_kets()
    for gamma in (0.05, 0.2, 0.5):
        for scheme in ("period", "hump"):
            rng = np.random.default_rng(20260925)
            runs = rising = 0
            closest = -np.inf
            edges = {"left": 0, "right": 0}
            windows = interior = 0
            for delta in (0.0, 0.5, 2.0, 3.0):
                H = hopping(delta)
                cases = [("named", psi) for psi in named.values()]
                for _ in range(8):
                    v = rng.normal(size=4) + 1j * rng.normal(size=4)
                    cases.append(("pure", v / np.linalg.norm(v)))
                if scheme == "period":
                    for _ in range(8):
                        A_ = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
                        r = A_ @ A_.conj().T
                        cases.append(("mixed", r / np.trace(r)))
                for kind, st in cases:
                    rho0 = st if st.ndim == 2 else np.outer(st, st.conj())
                    bounds = (period_windows(H, rho0, 8) if scheme == "period"
                              else hump_windows(H, st, 12.0))
                    wm = window_maxima(Trajectory(H, gamma, rho0), bounds)
                    vals = [v for _, v, _ in wm]
                    windows += len(wm)
                    for _, _, e in wm:
                        if e == "inside":
                            interior += 1
                        else:
                            edges[e] += 1
                    diff = max(vals[k + 1] - vals[k] for k in range(len(vals) - 1))
                    runs += 1
                    rising += diff > 0
                    closest = max(closest, diff)
            what = ("per period (8 periods; 14 named, 8 random pure, 8 random rank-2 mixed states per H)"
                    if scheme == "period" else
                    "per hump between γ = 0 minima (t ≤ 12; 14 named and 8 random pure states per H)")
            lines.append(f"  γ = {gamma:g}, main peaks {what}: {runs} runs, {windows} windows ({interior} with an "
                         f"interior maximum), rising trajectories {rising}, largest successive difference "
                         f"{closest:+.2e}; window maxima on the left edge {edges['left']} (such a window cannot "
                         f"rise: its left edge closes the previous window), on the right edge {edges['right']}")
    return lines


# ---------------------------------------------------------------- (D) beside the atlas
def chain_h(n, kind, J=1.0):
    letters, coef = ((X, Y), J / 2.0) if kind == "XY" else ((X, Y, Z), J / 4.0)
    return sum(coef * site_op(n, b, P) @ site_op(n, b + 1, P) for b in range(n - 1) for P in letters)


def bell_carrier(n):
    psi = np.zeros(2 ** n, dtype=complex)
    psi[0] = psi[(1 << (n - 1)) | (1 << (n - 2))] = 1.0 / np.sqrt(2.0)
    return psi


def atlas_reader_rises(curve, bar=1e-3):
    """The atlas reader: local maxima (left endpoint if descending; interior rise-or-flat in, strict
    fall out), three-point apex heights, predecessor rises above the reporting bar."""
    idx = ([0] if curve[0] > curve[1] else []) + [i for i in range(1, len(curve) - 1)
                                                  if curve[i] >= curve[i - 1] and curve[i] > curve[i + 1]]
    apex = []
    for i in idx:
        if i == 0:
            apex.append(curve[0])
            continue
        y0, y1, y2 = curve[i - 1], curve[i], curve[i + 1]
        den = y0 - 2 * y1 + y2
        v = y1 - (y0 - y2) ** 2 / (8 * den) if den < 0 else y1
        apex.append(v if v >= y1 else y1)
    diffs = [apex[k] - apex[k - 1] for k in range(1, len(apex))]
    return len(apex), sum(d > bar for d in diffs), max(diffs)


def part_d():
    lines = ["(D) beside the finite N ≥ 3 atlas (Symphony conventions; energies in units of J)"]
    for n, kind in ((3, "XY"), (4, "XY"), (5, "XY"), (3, "Heisenberg")):
        H = chain_h(n, kind)
        w, V = np.linalg.eigh(H)
        c = V.conj().T @ bell_carrier(n)
        # the carrier's weight on each eigenvector is either O(1) or zero up to rounding (~1e-16)
        occupied = sorted({round(float(e), 6) + 0.0 for e, a in zip(w, c) if abs(a) > 1e-9})
        f = lambda t, V=V, c=c, w=w, n=n: (np.abs(V @ (c * np.exp(-1j * w * t))).sum() ** 2 - 1.0) / (2 ** n - 1)
        t = np.linspace(0.0, 125.0, 250001)
        amps = np.abs((V * c) @ np.exp(-1j * np.outer(w, t)))
        cp = (amps.sum(axis=0) ** 2 - 1.0) / (2 ** n - 1)        # pure: L₁ = (Σ|ψ_i|)² − 1
        m = [golden(f, t[i - 1], t[i + 1])[1] for i in range(1, len(t) - 1)
             if cp[i] >= cp[i - 1] and cp[i] > cp[i + 1]]
        lines.append(f"  {kind} N = {n}: occupied energies {occupied}; γ = 0, J·t ≤ 125: {len(m)} maxima, "
                     f"spread {max(m) - min(m):.2e}, first eight " + ", ".join(f"{x:.6f}" for x in m[:8]))
    for n, kind, J in ((3, "XY", 20.0), (3, "Heisenberg", 20.0), (4, "XY", 5.0), (5, "XY", 5.0)):
        tr = Trajectory(chain_h(n, kind, J), 0.01, np.outer(bell_carrier(n), bell_carrier(n).conj()))
        curve = np.array([tr(25.0 * i / 1599) for i in range(1600)])
        count, rises, largest = atlas_reader_rises(curve)
        raw = atlas_reader_rises(curve, bar=0.0)[1]
        lines.append(f"  {kind} N = {n}, J = {J:g}, γ = 0.01 (Q = {J / 0.01:g}), t ≤ 25, 1600 points, atlas reader: "
                     f"{count} maxima, rises above the 10⁻³ bar {rises} ({raw} raw), largest successive "
                     f"difference {largest:+.2e}")
    return lines


# ---------------------------------------------------------------- (E) the local-unitary lift
def part_e():
    lines = ["(E) a local unitary lifts CΨ back above 1/4"]
    th = np.deg2rad(60.0)
    rz = np.cos(th)
    cp = lambda rx, rzz: ((1.0 + rx * rx + rzz * rzz) / 2.0) * rx / 3.0   # |0⟩ ⊗ Bloch (rx, 0, rz)
    rx0 = np.sin(th)
    rxc = brentq(lambda r: cp(r, rz) - 0.25, 0.1, rx0, xtol=1e-15)
    rho = lambda rx, rzz: np.kron(np.diag([1.0, 0.0]),
                                  0.5 * (I2 + rx * X + rzz * Z)).astype(complex)
    lines.append(f"  |0⟩ ⊗ (I + sin 60°·X + cos 60°·Z)/2: CΨ(0) = {cpsi(rho(rx0, rz)):.6f}; Z-dephasing "
                 f"shrinks r_x only, and CΨ = 1/4 at r_x = {rxc:.6f} (CΨ there {cpsi(rho(rxc, rz)):.6f})")
    angle = np.degrees(np.arctan2(rz, rxc))
    r = np.hypot(rxc, rz)
    lines.append(f"  R_y on qubit 2 by {angle:.1f}° turns its Bloch vector onto the equator, (|r|, 0, 0) with "
                 f"|r| = {r:.6f}: CΨ = {cpsi(rho(r, 0.0)):.6f}, the largest any unitary on qubit 2 reaches "
                 f"at fixed purity")
    return lines


def main():
    out = ["Successive local maxima of CΨ at N = 2 under local Z-dephasing "
           "(simulations/envelope_n2_rises.py)", ""]
    for part in (part_a, part_b, part_c, part_d, part_e):
        out.extend(part())
        out.append("")
    text = "\n".join(out)
    print(text)
    RESULTS.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
