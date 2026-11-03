"""Second route for the F86 envelope Q_peak, independent of the C# ResonanceScan on the way to max_t.

The observable is K_b(Q, t) = 2·Re⟨ρ(t)| S |∂ρ/∂J_b⟩ on the (n, n+1) coherence block with the Dicke
probe and the spatial-sum kernel S, at γ₀ = 0.05, J = Q·γ₀; the F86 curve is the envelope
max_t |K_b(Q, t)|, and Q_peak, K_max and HWHM_left/Q_peak are read off it per bond class (Endpoint:
the two end bonds; Interior: the rest) by a transcription of the C# ParabolicPeakFinder on the C#
default Q grid (153 points on [0.2, 4]).

What is independent of C#:
  * the block. c = 2 goes through the F89 eigendecomposition route of
    f89_to_f86_kbond_via_eigendecomp.py (F89 J = F86 J / 2); c >= 3 through a numpy build of the
    (n, n+1) block from the bit strings;
  * the kernel, used through its rank-N factor S = 2 Σ_l w_l w_lᵀ;
  * Duhamel in the closed form f(t) = A e(t) − e(t)∘rowsum(A) + t e(t)∘B, vectorised over t
    (A_rc = X_rc c0_c / (λ_c − λ_r) off the degenerate pairs, B the degenerate diagonal sum);
  * the maximiser in t: a dense scan on [0.3, 40] and a bounded Brent step around its winner
    (C# seeds a 45-point grid on [0.2, 2.4]·t_peak and refines by golden section).

Usage:
  python simulations/f86_envelope_second_route.py N c          # per bond and per class
  python simulations/f86_envelope_second_route.py --fixture    # all fixture blocks
  python simulations/f86_envelope_second_route.py --fixture N c  # one block, merged into it
  python simulations/f86_envelope_second_route.py N c Q_HI     # per bond on [0.2, Q_HI], dQ = 0.025,
                                                               # global and canonical (first local max)
  python simulations/f86_envelope_second_route.py --fourmode N # four-mode model, class peaks
  python simulations/f86_envelope_second_route.py --bare       # bare doubled-PTF 2x2 envelope
  python simulations/f86_envelope_second_route.py --continuous N  # c=2: canonical peak per bond,
                                                               # grid reading vs continuous search in Q
  python simulations/f86_envelope_second_route.py --lineshape N c  # class-mean y(x) = K(Q*(1+x))/K_max at
                                                               # the ShapeFunctionWitnesses x points
The fixture is compute/RCPsiSquared.Core.Tests/Fixtures/python_envelope_route_summary.json, read by
ValidateAgainstPythonEnvelopeRouteTests. --fourmode is the second route of FourModeResonanceScanTests
(the R-even vector in a degenerate top singular pair, as InterChannelSvd picks it); --bare that of
C2BareDoubledPtfClosedFormTests (brute Duhamel by a block exponential, no closed form).
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f89_to_f86_kbond_via_eigendecomp as f89  # noqa: E402

G0 = 0.05
Q_GRID = np.linspace(0.2, 4.0, 153)
FIXTURE_BLOCKS = [(2, 5), (2, 6), (2, 7), (2, 8), (3, 5), (3, 6), (3, 7), (3, 8), (4, 7)]
# (4, 8) is in the fixture too, but a default run does not recompute it: its block (dim 3920) takes hours
# here; re-read it explicitly with `--fixture 8 4`, which merges it.
FIXTURE = (Path(__file__).resolve().parent.parent / "compute" / "RCPsiSquared.Core.Tests"
           / "Fixtures" / "python_envelope_route_summary.json")


def _states(N, k):
    return sorted(sum(1 << (N - 1 - s) for s in c) for c in combinations(range(N), k))


def _block(N, n):
    """(n, n+1) block: dephasing diagonal D, per-bond M_H, Dicke probe, kernel rows W (S = 2 WᵀW)."""
    P, Qs = _states(N, n), _states(N, n + 1)
    pi = {p: i for i, p in enumerate(P)}
    qi = {q: i for i, q in enumerate(Qs)}
    Mq = len(Qs)
    M = len(P) * Mq
    nb = N - 1
    D = np.zeros(M)
    Mh = np.zeros((nb, M, M), complex)

    def flips(s):
        for b in range(nb):
            hi, lo = (s >> (N - 1 - b)) & 1, (s >> (N - 2 - b)) & 1
            if hi != lo:
                yield b, s ^ ((1 << (N - 1 - b)) | (1 << (N - 2 - b)))

    for p in P:
        for q in Qs:
            i = pi[p] * Mq + qi[q]
            D[i] = -2 * G0 * bin(p ^ q).count("1")
            for b, pf in flips(p):
                Mh[b, pi[pf] * Mq + qi[q], i] += -1j
            for b, qf in flips(q):
                Mh[b, pi[p] * Mq + qi[qf], i] += 1j
    probe = np.full(M, 1 / (2 * np.sqrt(len(P) * Mq)), complex)
    W = np.zeros((N, M))
    for site in range(N):
        m = 1 << (N - 1 - site)
        for p in P:
            if not p & m and (p | m) in qi:
                W[site, pi[p] * Mq + qi[p | m]] = 1.0
    return D, Mh, probe, W


def _k_functions(ev, R, Ri, c0, W, Xs, kernel_scale):
    """|K_b(t)| for each bond matrix X_b = R⁻¹ M_b R, vectorised over t."""
    WR = W @ R
    Dm = ev[None, :] - ev[:, None]
    deg = np.abs(Dm) <= 1e-10
    out = []
    for X in Xs:
        Xc = X * c0[None, :]
        A = np.where(deg, 0, Xc / np.where(deg, 1, Dm))
        rs, B = A.sum(1), np.where(deg, Xc, 0).sum(1)

        def K(t, A=A, rs=rs, B=B):
            t = np.atleast_1d(np.asarray(t, float))
            E = np.exp(np.outer(ev, t))
            rho = WR @ (E * c0[:, None])
            dr = WR @ (A @ E - E * rs[:, None] + t[None, :] * E * B[:, None])
            return np.abs(kernel_scale * 2 * np.sum(np.conj(rho) * dr, axis=0).real)
        out.append(K)
    return out


def _k_at_q_c2(q, N):
    L, basis, de = f89.build_se_de_full((q / 2) * G0, G0, N)
    n = len(basis)
    ev, R = np.linalg.eig(L)
    Ri = np.linalg.inv(R)
    c0 = Ri @ np.full(n, np.sqrt(2 / (N * N * (N - 1))) / 2, dtype=complex)
    W = np.zeros((N, n))
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l != i and (jk[1] if jk[0] == l else jk[0]) == i:
                W[l, idx] = 1.0
    Xs = [Ri @ f89.build_per_bond_M(b, N, basis, de) @ R for b in range(N - 1)]
    return _k_functions(ev, R, Ri, c0, W, Xs, 1.0)


def _k_at_q_block(q, model):
    D, Mh, probe, W = model
    L = np.diag(D).astype(complex) + q * G0 * Mh.sum(0)
    ev, R = np.linalg.eig(L)
    Ri = np.linalg.inv(R)
    Xs = [Ri @ Mh[b] @ R for b in range(Mh.shape[0])]
    return _k_functions(ev, R, Ri, Ri @ probe, W, Xs, 2.0)


def _envelope(K, ts):
    ks = K(ts)
    i = int(np.argmax(ks))
    r = minimize_scalar(lambda t: -K(t)[0], bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                        method="bounded", options={"xatol": 1e-12})
    return max(-r.fun, ks[i])


def envelope_curves(N, c, q_grid=Q_GRID):
    """bonds × Q array of max_t |K_b(Q, t)| on q_grid (default: the C# Q grid)."""
    if c == 2:
        ts, kq = np.linspace(0.3, 40.0, 1600), (lambda q: _k_at_q_c2(q, N))
    else:
        model = _block(N, c - 1)
        ts, kq = np.linspace(0.3, 40.0, 800), (lambda q: _k_at_q_block(q, model))
    k = np.zeros((N - 1, len(q_grid)))
    for j, q in enumerate(q_grid):
        for b, K in enumerate(kq(q)):
            k[b, j] = _envelope(K, ts)
    return k


def peak(grid, k, canonical=False):
    """Transcription of the C# ParabolicPeakFinder: (Q_peak, K_max, HWHM_left). Default: the
    global maximum on the grid (Find); canonical=True: the first local maximum from below
    (FindFirstLocalMaximum), falling back to the global one when there is none."""
    i = int(np.argmax(k))
    if canonical:
        i = next((j for j in range(1, len(k) - 1) if k[j] > k[j - 1] and k[j] >= k[j + 1]), i)
    kmax, qp = k[i], grid[i]
    if 0 < i < len(grid) - 1:
        x0, x1, x2 = grid[i - 1:i + 2]
        y0, y1, y2 = k[i - 1:i + 2]
        det = (x0 - x1) * (x0 - x2) * (x1 - x2)
        a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / det
        b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / det
        cc = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / det
        if a < 0:
            qs = -b / (2 * a)
            if abs(qs - grid[i]) <= grid[i + 1] - grid[i - 1]:
                qp, kmax = qs, cc - b * b / (4 * a)
    half, hl = kmax / 2, None
    for j in range(i, 0, -1):
        if k[j] < half:
            x0, x1, y0, y1 = grid[j], grid[j + 1], k[j], k[j + 1]
            hl = qp - (x0 + (half - y0) * (x1 - x0) / (y1 - y0))
            break
    return qp, kmax, hl


def class_peaks(k):
    nb = k.shape[0]
    out = {}
    for name, cls in (("interior", list(range(1, nb - 1))), ("endpoint", [0, nb - 1])):
        qp, km, hl = peak(Q_GRID, k[cls].mean(0))
        out[name] = (float(qp), float(km), float(hl / qp))
    return out


def fixture_entry(N, c):
    e = {"label": f"c{c}_N{N}", "c": c, "N": N, "n": c - 1, "gamma_0": G0}
    k = envelope_curves(N, c)
    for b in range(N - 1):
        qp, km, hl = peak(Q_GRID, k[b])
        print(f"c={c} N={N} bond {b}: Q_peak={qp!r} K_max={km!r} HWHM_left/Q_peak={hl / qp if hl else None!r}", flush=True)
    for name, (qp, km, r) in class_peaks(k).items():
        e[f"Q_peak_{name}"], e[f"K_max_{name}"], e[f"hwhm_left_over_q_peak_{name}"] = qp, km, r
    return e


def write_fixture(blocks):
    """Computes the given (c, N) blocks and merges them into the fixture, ordered by (c, N)."""
    entries = {(e["c"], e["N"]): e for e in json.loads(FIXTURE.read_text())} if FIXTURE.exists() else {}
    for c, N in blocks:
        entries[(c, N)] = fixture_entry(N, c)
        print(entries[(c, N)], flush=True)
        ordered = [entries[key] for key in sorted(entries)]
        FIXTURE.write_bytes((json.dumps(ordered, indent=2) + "\n").encode())


def _rev(s, N):
    return int(format(s, f"0{N}b")[::-1], 2)


def fourmode_curves(N, q_grid=Q_GRID, n=1):
    """bonds x Q envelope of the four-mode model: the HD-1 and HD-3 channel vectors and the top
    singular pair of the HD-1 -> HD-3 block of the total hopping; in a degenerate top pair the
    R-even direction (R the chain reversal). Max over t by a dense scan on [0.3, 40] + Brent."""
    D, Mh, probe, W = _block(N, n)
    P, Qs = _states(N, n), _states(N, n + 1)
    pi = {p: i for i, p in enumerate(P)}
    qi = {q: i for i, q in enumerate(Qs)}
    M = len(P) * len(Qs)
    R = np.zeros((M, M))
    for p_ in P:
        for q_ in Qs:
            R[pi[_rev(p_, N)] * len(Qs) + qi[_rev(q_, N)], pi[p_] * len(Qs) + qi[q_]] = 1
    hd = np.array([bin(p_ ^ q_).count("1") for p_ in P for q_ in Qs])

    def proj(h):
        idx = np.where(hd == h)[0]
        Pm = np.zeros((M, len(idx)), complex)
        Pm[idx, np.arange(len(idx))] = 1
        return Pm

    P1, P3 = proj(1), proj(3)
    V = P1.conj().T @ Mh.sum(0) @ P3
    U, sv, _ = np.linalg.svd(V)
    deg = int(np.sum(np.abs(sv - sv[0]) / sv[0] < 1e-8))
    us = [P1 @ U[:, k] for k in range(deg)]
    u0 = np.linalg.svd(np.column_stack([u + R @ u for u in us]), full_matrices=False)[0][:, 0]
    u0 /= np.linalg.norm(u0)
    v0 = P3 @ (V.conj().T @ (P1.conj().T @ u0))
    v0 /= np.linalg.norm(v0)
    chan = lambda h: (hd == h).astype(complex) / np.sqrt((hd == h).sum())
    B = np.column_stack([chan(1), chan(3), u0, v0])
    pr = lambda m: B.conj().T @ m @ B
    De, Mb, Mt, pe, We = pr(np.diag(D)), [pr(m) for m in Mh], pr(Mh.sum(0)), B.conj().T @ probe, W @ B
    ts = np.linspace(0.3, 40.0, 1600)
    k = np.zeros((N - 1, len(q_grid)))
    for j, q in enumerate(q_grid):
        L = De + q * G0 * Mt
        ev, Rv = np.linalg.eig(L)
        Ri = np.linalg.inv(Rv)
        Ks = _k_functions(ev, Rv, Ri, Ri @ pe, We, [Ri @ m @ Rv for m in Mb], 2.0)
        for b, K in enumerate(Ks):
            k[b, j] = _envelope(K, ts)
    return k


def bare_envelope():
    """Bare doubled-PTF 2x2 model, L = [[-2, 2ix], [2ix, -6]], V = [[0, 2i], [2i, 0]], rho0 = (1, 0),
    gamma0 = 1, tau = 4t: K by a block exponential (Duhamel), maximised over t, then over x.
    Returns (x_peak, tau*, x_half, ratio)."""
    V = np.array([[0, 2j], [2j, 0]])

    def K(x, t):
        L = np.array([[-2, 2j * x], [2j * x, -6]])
        Bm = np.zeros((4, 4), complex)
        Bm[:2, :2] = L
        Bm[2:, 2:] = L
        Bm[:2, 2:] = V
        E = expm(Bm * t)
        r0 = np.array([1, 0], complex)
        return 2 * np.real(np.vdot(E[2:, 2:] @ r0, E[:2, 2:] @ r0))

    ts = np.linspace(0.01, 3, 600)

    def env(x):
        ks = np.abs([K(x, t) for t in ts])
        i = int(np.argmax(ks))
        r = minimize_scalar(lambda t: -abs(K(x, t)), bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                            method="bounded", options={"xatol": 1e-13})
        return -r.fun, r.x

    xs = np.linspace(0.05, 6, 400)
    ev = np.array([env(x)[0] for x in xs])
    i = int(np.argmax(ev))
    r = minimize_scalar(lambda x: -env(x)[0], bounds=(xs[i - 1], xs[i + 1]), method="bounded",
                        options={"xatol": 1e-12})
    xp, kp = r.x, -r.fun
    xh = brentq(lambda x: env(x)[0] - kp / 2, xs[1], xp, xtol=1e-15)
    return xp, 4 * env(xp)[1], xh, (xp - xh) / xp


def continuous_canonical(N):
    """c=2: per bond, the canonical peak read on the default Q grid against a continuous search
    (bounded Brent in Q around the grid reading, the half-maximum by brentq), with t* at the peak.
    Prints Q_peak and HWHM_left/Q_peak both ways and their differences."""
    ts = np.linspace(0.3, 40.0, 1600)
    k = envelope_curves(N, 2)

    def env_at(q, b):
        K = _k_at_q_c2(q, N)[b]
        ks = K(ts)
        i = int(np.argmax(ks))
        r = minimize_scalar(lambda t: -K(t)[0], bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                            method="bounded", options={"xatol": 1e-12})
        return max(-r.fun, ks[i]), r.x

    for b in range((N - 1 + 1) // 2):
        qg, kg, hg = peak(Q_GRID, k[b], canonical=True)
        h = Q_GRID[1] - Q_GRID[0]
        r = minimize_scalar(lambda q: -env_at(q, b)[0], bounds=(qg - h, qg + h), method="bounded",
                            options={"xatol": 1e-10})
        qc, kc = r.x, -r.fun
        xh = brentq(lambda q: env_at(q, b)[0] - kc / 2, Q_GRID[0], qc * 0.9, xtol=1e-12)
        print(f"N={N} bond {b}: grid Q_peak={qg!r} ratio={hg / qg!r} | continuous Q_peak={qc!r} "
              f"ratio={(qc - xh) / qc!r} t*={env_at(qc, b)[1]!r} | grid-continuous dQ={qg - qc:.3e} "
              f"dratio={hg / qg - (qc - xh) / qc:.3e}", flush=True)


LINESHAPE_X = (-0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 1.00)


def lineshape(N, c):
    """Class-mean envelope y(x) = K(Q*(1+x))/K_max, x = (Q - Q*)/Q*, the relative-Q lineshape of
    F86b Statement 2. The curve is read on [0.2, 6] (step 0.025) so that x = +1 stays on the grid
    for the Endpoint class (Q* near 2.6); Q* and K_max are the class peak on the default window
    [0.2, 4], the reading of the HWHM band; y between grid points is linear, as in the C#
    ShapeFunctionWitnesses."""
    grid = np.linspace(0.2, 6.0, int(round((6.0 - 0.2) / 0.025)) + 1)
    k = envelope_curves(N, c, grid)
    nb, inside = k.shape[0], grid <= 4.0 + 1e-12
    out = {}
    for name, cls in (("interior", list(range(1, nb - 1))), ("endpoint", [0, nb - 1])):
        curve = k[cls].mean(0)
        qp, km, _ = peak(grid[inside], curve[inside])
        out[name] = (qp, [float(np.interp(qp * (1 + x), grid, curve) / km) for x in LINESHAPE_X])
    return out


if __name__ == "__main__":
    if sys.argv[1:2] == ["--lineshape"]:
        N, c = int(sys.argv[2]), int(sys.argv[3])
        for name, (qp, ys) in lineshape(N, c).items():
            print(f"lineshape c={c} N={N} {name}: Q*={qp!r} y=" + " ".join(f"{x:+.2f}:{y!r}" for x, y in zip(LINESHAPE_X, ys)), flush=True)
        sys.exit()
    if sys.argv[1:2] == ["--bare"]:
        xp, taus, xh, ratio = bare_envelope()
        print(f"bare envelope: x_peak={xp!r} tau*={taus!r} x_half={xh!r} HWHM_left/x_peak={ratio!r}")
        sys.exit()
    if sys.argv[1:2] == ["--continuous"]:
        continuous_canonical(int(sys.argv[2]))
        sys.exit()
    if sys.argv[1:2] == ["--fourmode"]:
        N = int(sys.argv[2])
        cp = class_peaks(fourmode_curves(N))
        for name, (qp, km, r) in cp.items():
            print(f"four-mode N={N} {name}: Q_peak={qp!r} K_max={km!r} HWHM_left/Q_peak={r!r}")
        sys.exit()
    if sys.argv[1:2] == ["--fixture"]:
        # --fixture: every block of FIXTURE_BLOCKS; --fixture N c: that one block, merged in.
        write_fixture(FIXTURE_BLOCKS if len(sys.argv) == 2 else [(int(sys.argv[3]), int(sys.argv[2]))])
    else:
        N, c = int(sys.argv[1]), int(sys.argv[2])
        grid = Q_GRID if len(sys.argv) < 4 else np.linspace(0.2, float(sys.argv[3]),
                                                             int(round((float(sys.argv[3]) - 0.2) / 0.025)) + 1)
        k = envelope_curves(N, c, grid)
        for b in range(N - 1):
            for label, canon in (("global", False), ("canonical", True)):
                qp, km, hl = peak(grid, k[b], canonical=canon)
                print(f"c={c} N={N} bond {b} {label}: Q_peak={qp!r} K_max={km!r} "
                      f"HWHM_left/Q_peak={hl / qp if hl else None!r}")
        for name, (qp, km, r) in (class_peaks(k).items() if len(sys.argv) < 4 else []):
            print(f"c={c} N={N} {name}: Q_peak={qp!r} K_max={km!r} HWHM_left/Q_peak={r!r}")
