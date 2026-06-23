#!/usr/bin/env python3
"""SUBSYSTEM_CROSSING Step-2 Case-C gate, my own (deep proof review 2026-06-22).

Claim under test (proof Step 2, Case C + headline): for ANY primitive CPTP map on 2 qubits,
the unique fixed point has CΨ(ρ*) < 1/4, and every trajectory from CΨ>1/4 eventually < 1/4.

CΨ definition (proof line 157): CΨ(ρ) = Tr(ρ²) · L1(ρ)/(d-1), L1 = off-diagonal l1-coherence.

Counterexample to test: depolarize-toward-σ, ε(ρ) = (1-p)ρ + p·Tr(ρ)·σ,
σ = 0.95·|Φ+⟩⟨Φ+| + 0.05·I/4. Is ε PRIMITIVE, and is CΨ(σ) > 1/4? Does the Bell+
trajectory ever cross below 1/4?
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

d = 4
I4 = np.eye(d, dtype=complex)
phi = np.zeros(d, dtype=complex); phi[0] = 1/np.sqrt(2); phi[3] = 1/np.sqrt(2)  # |Φ+⟩=(|00⟩+|11⟩)/√2
Phi = np.outer(phi, phi.conj())
bell = Phi.copy()  # Bell+ as the initial state (CΨ=1/3)


def L1_offdiag(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def CPsi(rho):
    return np.trace(rho @ rho).real * L1_offdiag(rho) / (d - 1)


def channel(rho, p, sigma):
    return (1 - p) * rho + p * np.trace(rho) * sigma


def superop(p, sigma):
    """Build the 16x16 superoperator of ε in column-stacking, to test primitivity."""
    S = np.zeros((d * d, d * d), dtype=complex)
    for k in range(d * d):
        E = np.zeros(d * d, dtype=complex); E[k] = 1
        rho = E.reshape(d, d, order='F')
        S[:, k] = channel(rho, p, sigma).reshape(-1, order='F')
    return S


def main():
    sigma = 0.95 * Phi + 0.05 * I4 / 4
    print(f"σ Hermitian? {np.allclose(sigma, sigma.conj().T)}, Tr σ = {np.trace(sigma).real:.6f}, "
          f"PSD min-eig = {np.linalg.eigvalsh(sigma).min():.6f}")
    print(f"CΨ(σ) = {CPsi(sigma):.6f}   (threshold 1/4 = 0.25)")
    print(f"CΨ(Bell+) = {CPsi(bell):.6f}   (expect 1/3 ≈ 0.3333)")
    print("-" * 64)
    for p in (1.0, 0.5, 0.1, 0.01):
        # fixed point check
        fp_resid = np.max(np.abs(channel(sigma, p, sigma) - sigma))
        # primitivity: superoperator spectrum
        ev = np.linalg.eigvals(superop(p, sigma))
        mods = np.sort(np.abs(ev))[::-1]
        lam1, lam2 = mods[0], mods[1]
        primitive = (abs(lam1 - 1) < 1e-9) and (lam2 < 1 - 1e-9)
        print(f"p={p}: ε(σ)=σ resid={fp_resid:.2e}; |λ|_1={lam1:.4f}, |λ|_2={lam2:.4f} "
              f"→ primitive(unique FP, gap)? {primitive}")
    print("-" * 64)
    # trajectory from Bell+ at p=0.1: does CΨ ever cross below 1/4?
    p = 0.1
    rho = bell.copy()
    below = False; mn = 1.0
    traj = []
    for n in range(1, 201):
        rho = channel(rho, p, sigma)
        c = CPsi(rho); mn = min(mn, c); traj.append(c)
        if c < 0.25:
            below = True
    print(f"Bell+ trajectory under ε (p={p}), n=1..200: CΨ_final={traj[-1]:.6f}, "
          f"CΨ_min={mn:.6f}, ever < 1/4? {below}")
    print("-" * 64)
    cs = CPsi(sigma)
    if cs > 0.25:
        print(f"VERDICT: a PRIMITIVE CPTP map has CΨ(ρ*)={cs:.4f} > 1/4 → Case C 'CΨ(ρ*)<1/4 for")
        print("all primitive maps' is FALSE; the headline 'eventual crossing' fails for this channel.")
    else:
        print("VERDICT: counterexample did NOT exceed 1/4 — Case C survives this probe.")


if __name__ == '__main__':
    main()
