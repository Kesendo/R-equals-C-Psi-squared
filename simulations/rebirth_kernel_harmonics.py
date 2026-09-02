"""Which spatial harmonic q carries the rebirth kernel C_a(tau) on the cracked ring, and which one the
first rebirth carries. A READING, not a gate: nothing here passes or fails.

Verifier for one paragraph of experiments/THE_CRACKED_BELL.md (2026-09-02, the section "The (1,1) block").
Reuses ring_renewal_current.kernels (imported silently: that script runs its own gates at import time,
and their PASS lines are not this script's). Two reads:

  1. RMS over the run of |FFT_a C_a(tau)|_q, the kernel's own spatial weight by harmonic q; and the
     seed-blindness of the kernel, max|C(m=1) - C(m=2)| compared as zero (C_a contains no psi_0).
  2. The Parseval split of the first rebirth, V(T_0) = int_0^T0 sum_a C_a(T0 - s) n_free,a(s) ds = sum_q V_q,
     V_q = (1/N) int Ctilde_q(T0 - s) conj(ntilde_q(s)) ds, with T_0 the FIRST ZERO CROSSING OF THE FREE
     CURRENT OF THAT GEOMETRY (found on the grid and refined linearly; N = 8, m = 1, delta = 0.15 gives the
     page's 20.2953). A first version of this script integrated the N = 12 geometry to N = 8's T_0.

Run: python simulations/rebirth_kernel_harmonics.py
"""
import contextlib
import io
import sys

import numpy as np

sys.path.insert(0, "simulations")
with contextlib.redirect_stdout(io.StringIO()):
    import ring_renewal_current as rr


def kernels_for(n, m, delta, ts):
    rr.M_SEED = m
    return rr.kernels(n, delta, ts)


def weights(n, m, delta, tmax=58.0, steps=5800):
    ts = np.linspace(0.0, tmax, steps + 1)
    K, C, n_free, I_free = kernels_for(n, m, delta, ts)
    F = np.abs(np.fft.fft(C, axis=1))            # [t, q]
    return np.sqrt(np.mean(F ** 2, axis=0)), C


def first_zero(ts, f):
    """First zero crossing of f after t = 0, linear refinement between grid points."""
    for i in range(1, len(ts) - 1):
        if f[i] == 0.0:
            return ts[i]
        if f[i] * f[i + 1] < 0:
            return ts[i] + (ts[i + 1] - ts[i]) * f[i] / (f[i] - f[i + 1])
    raise RuntimeError("no zero crossing on the grid")


def v_split(n, m, delta, tmax=80.0, steps=8000):
    ts_scan = np.linspace(0.0, tmax, steps + 1)
    _, _, _, I_free = kernels_for(n, m, delta, ts_scan)
    T0 = first_zero(ts_scan, I_free)
    ts = np.linspace(0.0, T0, 4001)
    K, C, n_free, I_free = kernels_for(n, m, delta, ts)
    Cq = np.fft.fft(C, axis=1)                   # [t, q]
    nq = np.fft.fft(n_free, axis=1)              # [t, q]
    integrand = (Cq[::-1] * np.conj(nq)).real / n   # sum_a C_a(T0 - s) n_a(s) = (1/N) sum_q Cq(T0 - s) conj(nq(s))
    Vq = np.trapezoid(integrand, ts, axis=0)
    V_direct = np.trapezoid(np.sum(C[::-1] * n_free, axis=1), ts)
    return T0, Vq, V_direct


if __name__ == "__main__":
    print("Read 1: the kernel's own spatial weight, RMS over the run of |FFT_a C_a(tau)|_q")
    C_by_m = {}
    for n, m, delta in [(8, 1, 0.15), (8, 1, 0.05), (8, 2, 0.15), (12, 1, 0.1), (11, 1, 0.15), (16, 1, 0.1)]:
        rms, C = weights(n, m, delta)
        C_by_m[(n, m, delta)] = C
        order = np.argsort(-rms)
        print(f"  N={n} m={m} delta={delta}: RMS |C_q| by q = " + ", ".join(f"q{q}:{rms[q]:.4f}" for q in range(n)))
        nz = rms[1:]
        print(f"     largest q = {order[0]} (N/2 = {n/2}, 2m = {2*m}), second q = {order[1]}, ratio second/first = {rms[order[1]]/rms[order[0]]:.3f}; "
              f"q = 0 exactly {rms[0]:.1e}; over q >= 1: min {nz.min():.4f}, max {nz.max():.4f}, spread (max-min)/min = {(nz.max()-nz.min())/nz.min()*100:.1f}%")
    diff = np.max(np.abs(C_by_m[(8, 1, 0.15)] - C_by_m[(8, 2, 0.15)]))
    print(f"  seed-blindness: max|C_a(tau; m=1) - C_a(tau; m=2)| at N=8, delta=0.15 = {diff:.1e} (compared as zero: {diff == 0.0})")
    print("Read 2: the Parseval split of the first rebirth V(T_0) = sum_q V_q, T_0 the free current's first zero of THAT geometry")
    for n, m, delta in [(8, 1, 0.15), (12, 1, 0.1)]:
        T0, Vq, Vd = v_split(n, m, delta)
        tot = np.abs(Vq).sum()
        print(f"  N={n} m={m} delta={delta}: T_0 = {T0:.4f}, V(T0) direct = {Vd:.4f}, sum_q V_q = {Vq.sum():.4f}")
        print("     V_q by q: " + ", ".join(f"q{q}:{Vq[q]:+.4f}" for q in range(n)))
        print("     share of |V_q|: " + ", ".join(f"q{q}:{abs(Vq[q])/tot:.4f}" for q in range(n)))
        print(f"     share at q = +-2m: {(abs(Vq[2*m % n]) + abs(Vq[(-2*m) % n]))/tot*100:.1f}%; share at q = N/2: "
              f"{(abs(Vq[n//2])/tot*100 if n % 2 == 0 else float('nan')):.3f}%")
