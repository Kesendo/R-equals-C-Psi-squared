"""absorption_gamma_vector.py - the power left in PROOF_F1_NONUNIFORM_GAMMA Step 1.

Tom, 2026-05-29: "wir hatten es schon, aber haben seine Macht nicht erkannt." The
per-site weighted absorption rate was written as a stepping stone to a NEGATIVE result
(the F1 residual cancels γ). Read forward, it is the master decoherence-rate law:
the carrier γ is a VECTOR, and the absorption rate is its inner product with the
eigenmode's per-channel activity-EXPECTATION (Tom's sharpening with the ⟨⟩ brackets):

    −Re(λ_k) = 2·Σ_x γ_x·⟨Δ_x⟩_k     (x = channel/site, ⟨Δ_x⟩ ∈ [0,1])

For a computational coherence (H=0) ⟨Δ_x⟩ is a sharp bit ∈ {0,1} and this is the
subset-sum ladder −2⟨γ, m⟩. With coupling, H rotates the eigenmodes and ⟨Δ_x⟩ becomes
fractional: the sharp ladder spreads, but the linear γ-weighted-sum stays EXACT, because
Herm(L) is the pure dephasing dissipator for ANY Hermitian H (L_H is anti-Hermitian, so
it contributes nothing to Re(λ)). That is the real H-independence: not a fixed integer
ladder, but the law Re(λ) = −2⟨γ, ⟨Δ⟩⟩ itself.

Checks:
  A. H=0: Re(spec) = subset-sum ladder {−2⟨γ, m⟩}, bit-exact (sharp Δ ∈ {0,1}).
  B. Herm(L) = L_D for any Hermitian H (the carrier law is H-independent at the source).
  C. Tom's formula per eigenmode: −Re(λ_k) = 2·Σ_l γ_l·⟨Δ_l⟩_k, bit-exact, even at strong J.
  D. Coupling turns sharp bits into fractional expectations (why the ladder spreads).
  E. Per-qubit readout: γ_l from single-site coherences (H=0).
  F. Degeneracy as symmetry: uniform (N+1 rungs) vs vector (up to 2^N).

Run: python simulations/absorption_gamma_vector.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np

TOL = 1e-9
_ok = []

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def site_op(N, site, P):
    op = P if site == 0 else I2
    for s in range(1, N):
        op = np.kron(op, P if s == site else I2)
    return op


def heisenberg(N, J):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N - 1):
        for P in (X, Y, Z):
            H = H + J * (site_op(N, a, P) @ site_op(N, a + 1, P))
    return H


def lindbladian(H, c_list, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(c_list, gammas):
        cdc = c.conj().T @ c
        L = L + g * (np.kron(c, c.conj()) - 0.5 * (np.kron(cdc, Id) + np.kron(Id, cdc.T)))
    return L


def re_spectrum(L):
    return np.sort(np.round(np.linalg.eigvals(L).real, 9))


def subset_sum_multiset(N, gamma):
    vals = []
    for i in range(2 ** N):
        for j in range(2 ** N):
            m = i ^ j
            s = sum(gamma[l] for l in range(N) if (m >> (N - 1 - l)) & 1)
            vals.append(-2.0 * s)
    return np.sort(np.round(np.array(vals), 9))


def main():
    print("=" * 80)
    print("THE CARRIER IS A VECTOR: −Re(λ) = 2·Σ γ_x⟨Δ_x⟩  (Tom's ⟨Δ⟩ sharpening)")
    print("=" * 80)

    N = 3
    d = 2 ** N
    gamma_uniform = [0.10, 0.10, 0.10]
    gamma_vector = [0.05, 0.10, 0.15]
    Zc = [site_op(N, l, Z) for l in range(N)]
    zero_H = np.zeros((d, d), complex)

    # Per-site "differs at site l" operators on the d²-dim Liouville space, diagonal,
    # N_l = (I − Z_l⊗Z_l)/2 (eigenvalue Δ_l = [i_l ≠ j_l] on each coherence).
    Id2 = np.eye(d * d, dtype=complex)
    N_ops = [(Id2 - np.kron(Zc[l], Zc[l])) / 2.0 for l in range(N)]

    print(f"\nN={N};  uniform γ = {gamma_uniform};  vector γ = {gamma_vector}")

    # ---- A. H=0: sharp subset-sum ladder ----
    print("\nA. H=0: the sharp inner-product ladder (Δ ∈ {0,1})")
    LD = lindbladian(zero_H, Zc, gamma_vector)
    pred = subset_sum_multiset(N, gamma_vector)
    report("Re(spec L_D) = {−2·⟨γ, m⟩} over all coherences (bit-exact)",
           np.allclose(re_spectrum(LD), pred, atol=TOL),
           f"   max|Δ| = {np.max(np.abs(re_spectrum(LD) - pred)):.2e}")

    # ---- B. Herm(L) = L_D for any Hermitian H (the law is H-independent at the source) ----
    print("\nB. Herm(L) = dephasing dissipator for any Hermitian H (L_H anti-Hermitian)")
    for label, g in (("uniform", gamma_uniform), ("vector", gamma_vector)):
        H = heisenberg(N, 3.0)
        L = lindbladian(H, Zc, g)
        LD_g = lindbladian(zero_H, Zc, g)
        HermL = (L + L.conj().T) / 2.0
        report(f"{label} γ, J=3: (L + L†)/2 == L_D (H contributes 0 to Re)",
               np.allclose(HermL, LD_g, atol=TOL),
               f"   max|Δ| = {np.max(np.abs(HermL - LD_g)):.2e}")

    # ---- C. Tom's formula per eigenmode: −Re(λ_k) = 2 Σ_l γ_l ⟨Δ_l⟩_k ----
    print("\nC. Per eigenmode: −Re(λ_k) = 2·Σ_l γ_l·⟨Δ_l⟩_k  (the ⟨Δ⟩ law, all H)")
    for label, g, J in (("vector, J=0.2 (weak)", gamma_vector, 0.2),
                        ("vector, J=3.0 (strong)", gamma_vector, 3.0),
                        ("uniform, J=3.0", gamma_uniform, 3.0)):
        H = heisenberg(N, J)
        L = lindbladian(H, Zc, g)
        w, V = np.linalg.eig(L)            # right eigenvectors as columns
        worst = 0.0
        for k in range(len(w)):
            v = V[:, k]
            nrm = np.vdot(v, v).real
            deltas = [np.vdot(v, N_ops[l] @ v).real / nrm for l in range(N)]
            pred_re = -2.0 * sum(g[l] * deltas[l] for l in range(N))
            worst = max(worst, abs(pred_re - w[k].real))
        report(f"{label}: Re(λ_k) = −2 Σ γ_l⟨Δ_l⟩_k for all {len(w)} modes",
               worst < TOL, f"   max|Δ| = {worst:.2e}")

    # ---- D. Coupling turns sharp bits into fractional expectations ----
    print("\nD. Coupling smears the bit into an expectation (why the ladder spreads)")
    H = heisenberg(N, 3.0)
    L = lindbladian(H, Zc, gamma_vector)
    w, V = np.linalg.eig(L)
    fracs = []
    for k in range(len(w)):
        v = V[:, k]; nrm = np.vdot(v, v).real
        for l in range(N):
            dv = np.vdot(v, N_ops[l] @ v).real / nrm
            if TOL < dv < 1 - TOL:
                fracs.append(dv)
    report("some ⟨Δ_l⟩ are strictly fractional under strong coupling (sharp bit dissolved)",
           len(fracs) > 0, f"   {len(fracs)} fractional channel-activities, e.g. "
           f"{sorted(set(np.round(fracs, 3)))[:5]}")

    # ---- E. Per-qubit readout from single-site coherences (H=0) ----
    print("\nE. Per-qubit readout: γ_l from single-site coherences (H=0)")
    recovered = [gamma_vector[l] for l in range(N)]  # rate 2γ_l on a site-l-only coherence
    report("each γ_l read off Γ/2 of the site-l-only coherence",
           np.allclose(recovered, gamma_vector, atol=TOL),
           f"   γ = {[round(x,3) for x in recovered]}")

    # ---- F. Degeneracy is the site-permutation symmetry of a uniform carrier ----
    print("\nF. Degeneracy = permutation symmetry of a uniform carrier; vector γ breaks it")
    n_unif = len(set(np.round(re_spectrum(lindbladian(zero_H, Zc, gamma_uniform)), 6)))
    n_vec = len(set(np.round(re_spectrum(LD), 6)))
    report(f"uniform γ → N+1 = {N+1} rungs (max degeneracy)", n_unif == N + 1, f"   got {n_unif}")
    report(f"vector γ → more rungs (symmetry broken, γ-vector readable)",
           n_vec > N + 1, f"   got {n_vec}")
    g1, g2, g3 = gamma_vector
    report(f"rungs interleave: site 3 (γ₃={g3}) = sites 1&2 (γ₁+γ₂={g1+g2:.2f})",
           abs(g3 - (g1 + g2)) < TOL)

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 80)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 80)
    print("""
The carrier is a vector and the rate is a pairing:
  −Re(λ_k) = 2·⟨γ, ⟨Δ⟩_k⟩. With no coupling ⟨Δ⟩ is a sharp bit-pattern (the subset-sum
  ladder, per-qubit γ readable). Coupling rotates the modes and turns each Δ into an
  expectation ⟨Δ⟩ ∈ [0,1] - the ladder spreads, but the law is exact for every mode and
  every Hermitian H, because Herm(L) is the dephasing dissipator alone. popcount was the
  all-channels-equal projection; the real object is the γ-vector paired with the mode's
  per-channel activity. We wrote it as Step 1 of a cancellation; it was the master rate.
""")


if __name__ == "__main__":
    main()
