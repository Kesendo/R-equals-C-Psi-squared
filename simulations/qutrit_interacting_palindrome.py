#!/usr/bin/env python3
"""
_qutrit_interacting_palindrome.py - WIP scout: the F121 interacting follow-up.

The dissipator's partial palindrome at d>2 has a closed-form value about the
physical palindrome center -N*gamma (F121: paired = Σ d^N C(N,k)(d-1)^min(k,N-k),
54/81 at d=3,N=2, full iff d=2). This scout settles the FULL Liouvillian
L = L_H + L_D, correcting an earlier center-mismatch error and finding the
mechanism.

FINDINGS THIS SCOUT PINS (self-validating):

  [A] Re(lambda) = -2*gamma * <Q>_mode exactly for every Hermitian H, where
      <Q> is the right Hilbert-Schmidt Rayleigh quotient v^dag Qv/v^dag v.
      SU(3) symmetry quantizes these values; a random Hermitian H spreads them.

  [B] For SU(3) Heisenberg the symmetry QUANTIZES <Q> into {0,1,1.5,2} (rungs
      {0:6, 1:36, 1.5:12, 2:27}), clean real parts {0,-2,-3,-4}*g; the -3g rung is
      exactly <Q>=1.5 (a 50/50 Hamming-1/Hamming-2 mix). A random H spreads <Q>
      continuously -> no rungs.

  [C] The sampled SU(3)-Heisenberg H reduces pairing at both checked centers.
      About the physical center -N*gamma: dissipator 54 -> full L 48.
      About -3g (the two big rungs): dissipator 72 -> full L 60. The earlier
      "full L exceeds the ceiling (60>54)" was a center mismatch (full-L-best -3g
      vs dissipator-physical-center -2g). This is not universal over H: H=cI
      has L_H=0 and leaves every dissipator count unchanged.

  [D] The interacting paired count is H-DEPENDENT (SU(3) Heisenberg 60 robustly
      across J=0.05..10; off-diagonal-only 48; single bilinear 52; random ~0).
      So there is NO H-independent closed form for the interacting palindrome;
      the dissipator's 54 (about -N*gamma) belongs to the invariant dissipator
      skeleton; it is not an interacting paired-count invariant.

Reuses the committed verifier's infrastructure.
"""

import numpy as np
import sys, os
from itertools import product as iprod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qutrit_partial_palindrome import (  # noqa: E402
    gm_raw, GAMMA, site_op, H_su3_heisenberg, L_dephasing, L_hamiltonian,
    palindrome_pairs,
)


def hamming_diag(N, d=3):
    """Q on the d^{2N} coherence space: diagonal, Q[(i,j)] = Hamming(i,j).
    L_D = -2*gamma*Q for the full-Cartan {lam3,lam8} dephasing."""
    states = list(iprod(range(d), repeat=N))
    diag = []
    for i in states:
        for j in states:
            diag.append(sum(1 for a, b in zip(i, j) if a != b))
    return np.array(diag, dtype=float)


def q_expectations(Lf, Q):
    """Right Hilbert-Schmidt Rayleigh quotient v^dag Q v / v^dag v per mode."""
    ev, V = np.linalg.eig(Lf)
    qexp = np.sum(V.conj() * (Q[:, None] * V), axis=0) / np.sum(V.conj() * V, axis=0)
    return ev, qexp


def best_pairing(evals, gamma):
    best, bc = 0, None
    for cc in np.linspace(0, 1.2, 1201):
        p = palindrome_pairs(evals, cc / 2, tol=1e-4)
        if p > best:
            best, bc = p, cc
    return best, (-bc / 2 / gamma if bc is not None else None)


def paired_at_center(evals, center_in_g, gamma):
    """paired count for reflection about center = -center_in_g * gamma.
    palindrome_pairs(ev, Sg) reflects lambda -> -2*Sg - lambda (center -Sg),
    so to centre at -center_in_g*gamma we pass Sg = center_in_g*gamma."""
    return palindrome_pairs(evals, center_in_g * gamma, tol=1e-4)


def rung_hist(evals, gamma):
    h = {}
    for r in np.round(evals.real / gamma, 2):
        h[r] = h.get(r, 0) + 1
    return {k: h[k] for k in sorted(h)}


def main():
    N, d, g = 2, 3, GAMMA
    Q = hamming_diag(N, d)
    LD = L_dephasing(N, g, [gm_raw[2], gm_raw[7]])

    print("=" * 70)
    print("F121 interacting follow-up: the full L = L_H + L_D")
    print("=" * 70)

    # ---- [A] Re(lambda) = -2g<Q> exact, for SU(3) Heisenberg AND random H ----
    print("\n[A] Re(lambda) = -2*gamma*<Q>_mode  (universal right-HS Rayleigh reading):")
    np.random.seed(11)
    H_su3 = H_su3_heisenberg(N, [(0, 1)], J=1.0)
    A = np.random.randn(9, 9) + 1j * np.random.randn(9, 9)
    H_rand = (A + A.conj().T) / 2
    devs = {}
    for label, H in [("SU(3) Heisenberg", H_su3), ("random Hermitian", H_rand)]:
        ev, qexp = q_expectations(L_hamiltonian(H) + LD, Q)
        dev = np.max(np.abs(ev.real - (-2 * g * qexp.real)))
        devs[label] = dev
        print(f"    {label:<20} max|Re(lambda) - (-2g*Re<Q>)| = {dev:.2e}")
    assert max(devs.values()) < 1e-10, f"Absorption-Theorem right-HS identity should be exact: {devs}"
    print("    -> exact for both Hamiltonians; symmetry quantizes the SU(3) values, not the identity. OK")

    # ---- [B] symmetry quantizes <Q>; -3g rung = <Q>=1.5 ----
    print("\n[B] Symmetry quantizes <Q>. SU(3) Heisenberg <Q> distribution:")
    ev, qexp = q_expectations(L_hamiltonian(H_su3) + LD, Q)
    qh = {}
    for q in np.round(qexp.real, 3):
        qh[q] = qh.get(q, 0) + 1
    print(f"    <Q> = {qh}   (rungs {rung_hist(ev, g)})")
    assert qh.get(0.0) == 6 and qh.get(1.0) == 36 and qh.get(1.5) == 12 and qh.get(2.0) == 27, \
        f"SU(3) <Q> distribution off: {qh}"
    # random H spreads <Q> continuously
    ev_r, qexp_r = q_expectations(L_hamiltonian(H_rand) + LD, Q)
    n_distinct_rand = len(set(np.round(qexp_r.real, 3)))
    print(f"    random H: {n_distinct_rand} distinct <Q> values (continuous spread, no rungs)")
    assert n_distinct_rand > 30, "random H should spread <Q>, not quantize it"
    print("    -> the -3g rung is exactly <Q>=1.5 (Hamming-1/Hamming-2 mix); "
          "symmetry quantizes, it does not create. OK")

    # ---- [C] sampled SU(3)-Heisenberg H changes pairing at two fixed centers ----
    print("\n[C] Sampled SU(3)-Heisenberg H reduces pairing at two checked centers:")
    evD = np.linalg.eigvals(LD)
    evF = np.linalg.eigvals(L_hamiltonian(H_su3) + LD)
    print(f"    {'center':>10}{'dissipator':>13}{'full L':>9}{'H effect':>12}")
    rows = []
    for cg in (2.0, 3.0):
        pD = paired_at_center(evD, cg, g)
        pF = paired_at_center(evF, cg, g)
        rows.append((cg, pD, pF))
        print(f"    {-cg:>9.1f}g{pD:>13}{pF:>9}{'degrades' if pF < pD else 'NO':>12}")
    # physical center -N*g = -2g: dissipator 54, full L 48
    assert rows[0] == (2.0, 54, 48), f"physical-center row wrong: {rows[0]}"
    # -3g: dissipator 72, full L 60
    assert rows[1] == (3.0, 72, 60), f"-3g row wrong: {rows[1]}"
    bD = best_pairing(evD, g)
    bF = best_pairing(evF, g)
    print(f"    best-over-centers: dissipator {bD[0]} @ {bD[1]:.1f}g ; full L {bF[0]} @ {bF[1]:.1f}g")
    assert bD[0] == 72 and bF[0] == 60, f"best-over-centers off: {bD}, {bF}"
    print("    -> at both checked centers this H reduces pairing (54->48, 72->60). The old")
    print("       '60>54 exceeds' compared full-L@-3g vs dissipator@-2g: a center mismatch. OK")

    # ---- [D] the interacting count is H-dependent (no closed form) ----
    print("\n[D] The interacting paired count is H-DEPENDENT (no closed form):")
    np.random.seed(7)
    offdiag = [m for i, m in enumerate(gm_raw) if i not in (2, 7)]
    cases = []
    for J in (0.0, 0.05, 1.0, 10.0):
        cases.append((f"SU(3) Heisenberg J={J}",
                      sum(J * site_op(m, 0, N, 3) @ site_op(m, 1, N, 3) for m in gm_raw)))
    cases.append(("SU(3) off-diagonal gens", sum(site_op(m, 0, N, 3) @ site_op(m, 1, N, 3) for m in offdiag)))
    cases.append(("single bilinear lam1xlam1", site_op(gm_raw[0], 0, N, 3) @ site_op(gm_raw[0], 1, N, 3)))
    for s in range(3):
        B = np.random.randn(9, 9) + 1j * np.random.randn(9, 9)
        cases.append((f"random Hermitian seed {s}", (B + B.conj().T) / 2))
    counts = {}
    print(f"    {'Hamiltonian':<28}{'best paired':>12}{'center(g)':>11}")
    for label, H in cases:
        b, c = best_pairing(np.linalg.eigvals(L_hamiltonian(H) + LD), g)
        counts[label] = b
        print(f"    {label:<28}{b:>12}{(c if c else 0):>11.2f}")
    assert counts["SU(3) Heisenberg J=0.0"] == 72
    su3 = [counts[k] for k in counts if k.startswith("SU(3) Heisenberg") and "J=0.0" not in k]
    assert all(x == 60 for x in su3), f"tested nonzero SU(3) couplings not at 60: {su3}"
    rand = [counts[k] for k in counts if k.startswith("random")]
    assert max(rand) < 20, f"random H should barely pair: {rand}"
    print("    -> J=0 returns the dissipator best count 72; tested nonzero J gives 60; random H ~0.")
    print("       The count floats")
    print("       with H's symmetry; no H-independent closed form. 54 (about -N*g) belongs")
    print("       to the dissipator skeleton; it is not an interacting paired-count invariant.")

    print("\n" + "=" * 70)
    print("ALL CHECKS PASSED. For the sampled SU(3)-Heisenberg H, pairing decreases")
    print("54->48 about the physical center -N*gamma and 72->60 about -3g. This is")
    print("not universal over H: H=cI leaves the dissipator unchanged. The old '60>54'")
    print("compared full-L@-3g against dissipator@-2g, hence mismatched centers. For")
    print("the symmetric SU(3) Heisenberg the real parts sit at -2g<Q> exactly (<Q> in")
    print("{0,1,1.5,2}, the -3g rung = <Q>=1.5); this QUANTIZATION is a symmetry effect")
    print("(a random H spreads the values, while retaining the Rayleigh identity). The interacting count is H-dependent (60 for the tested")
    print("nonzero SU(3)-Heisenberg couplings, ~0 generic): no H-independent closed form.")
    print("The dissipator's 54 about -N*gamma is a skeleton count, not an interacting invariant.")
    print("=" * 70)


if __name__ == "__main__":
    main()
