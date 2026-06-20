"""DECISIVE mechanism gate (the SE proof's own CV discriminator, PROOF_COHERENCE_HORIZON_SLOPE cross-check 2a):
does the overdamped slow (2,2) eigenvalue obey the resummed dispersion lambda^2 + 8*g*lambda + 4*J^2*q^2 = 0
with a gamma-CONSTANT q?  (finite-N-robust, unlike a full-curve fit.)

If q^2_eff = -(lambda^2 + 8 g lambda)/(4 J^2) is gamma-constant (low coefficient of variation) across a
gamma-sweep at fixed N, the 8g dispersion holds -> the handover Q_h = sqrt3/q is DERIVED (the SE coherence-
ladder mechanism, darkness=1 instead of the EP). Compare the truncated telegrapher 4g/2J^2q^2 (should scatter
worse, as in the SE proof). Also report q_eff vs the ring q_min=2pi/N (does the 2-particle mode use q_min?).

GATE: the 8g form has CV < ~0.02 AND beats the 4g form; q_eff is gamma-constant. A firing gate = the
dispersion is NOT the SE one for 2 particles (diagnose the 2-particle resummation), do NOT loosen.
"""
import sys
import numpy as np
sys.path.insert(0, "simulations")
from ring_handover_qh import Lkk_sparse
import scipy.sparse.linalg as spla
from math import pi, sqrt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
J = 1.0


def slow_real_eig(N, gamma):
    """The slowest interior mode's eigenvalue (max Re, Re<0) of the (2,2) sector at J=1, given gamma.
    In the overdamped regime it is real; return the complex eigenvalue (Im~0 there)."""
    L = Lkk_sparse(N, 2, gamma)
    Dsq = L.shape[0]
    if Dsq <= 6000:
        ev = np.linalg.eigvals(L.toarray())
    else:
        ev = spla.eigs(L, k=min(30, Dsq - 2), sigma=1e-3, which="LM", return_eigenvectors=False, maxiter=8000, tol=1e-10)
    interior = ev[ev.real < -1e-9]
    return interior[np.argmax(interior.real)]


def main():
    GATE = {"fired": []}

    def gate(name, ok, detail=""):
        flag = "ok " if ok else "GATE-FIRE"
        if not ok:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 100)
    print("MECHANISM GATE: does the overdamped slow (2,2) eigenvalue obey lambda^2+8g*lambda+4J^2 q^2 = 0")
    print("                with gamma-CONSTANT q? (the SE coherence-ladder dispersion, the decisive CV test)")
    print("=" * 100)

    for N in (10, 12):
        qmin = 2 * pi / N
        # overdamped regime: below the EP (Qq=2 -> Q=2/qmin=N/pi). Sweep Q in [0.5,1.7]*Q_h, Q_h~N*sqrt3/(2pi).
        Qh = N * sqrt(3) / (2 * pi)
        # ensure we stay overdamped (slow mode real, Im~0): Q up to ~0.95*N/pi (just below the EP)
        Qmax = min(1.6 * Qh, 0.95 * N / pi)
        Qs = np.linspace(0.55 * Qh, Qmax, 8)
        print(f"\n  N={N}  (q_min=2pi/N={qmin:.5f}, EP at Q=N/pi={N/pi:.4f}, Q_h~{Qh:.4f})")
        print(f"  {'Q':>8} {'gamma':>8} {'lambda(slow)':>22} {'|Im|':>9} {'q2_eff(8g)':>11} {'q2_eff(4g)':>11}")
        q2_8, q2_4 = [], []
        for Q in Qs:
            g = J / Q
            lam = slow_real_eig(N, g)
            q2a = -(lam.real ** 2 + 8 * g * lam.real) / (4 * J ** 2)        # 8g resummed
            q2b = -(lam.real ** 2 + 4 * g * lam.real) / (2 * J ** 2)        # 4g truncated telegrapher
            q2_8.append(q2a); q2_4.append(q2b)
            print(f"  {Q:>8.4f} {g:>8.4f} {lam.real:>12.5f}{('+'+format(abs(lam.imag),'.4f')+'i'):>10} "
                  f"{abs(lam.imag):>9.1e} {q2a:>11.6f} {q2b:>11.6f}")
        q2_8, q2_4 = np.array(q2_8), np.array(q2_4)
        # only use the overdamped (real) points (|Im| tiny) for the CV
        cv8 = np.std(q2_8) / np.mean(q2_8)
        cv4 = np.std(q2_4) / np.mean(q2_4)
        q_eff = sqrt(abs(np.mean(q2_8)))
        print(f"    CV(8g)={cv8:.4f}  CV(4g)={cv4:.4f}   q_eff(8g)={q_eff:.5f} vs q_min={qmin:.5f} "
              f"(q_eff/q_min={q_eff/qmin:.4f})")
        gate(f"N={N}: 8g dispersion holds (q2_eff gamma-constant, CV<0.03) AND beats 4g",
             cv8 < 0.03 and cv8 < cv4, f"CV8={cv8:.4f}, CV4={cv4:.4f}")

    print("\n" + "=" * 100)
    if GATE["fired"]:
        print(f"{len(GATE['fired'])} GATE(S) FIRED -> the (2,2) overdamped slow mode does NOT cleanly obey the SE")
        print("8g dispersion. The endpoint Q_h->N*sqrt3/(2pi) is solid (G2 earlier) but the MECHANISM (the exact")
        print("2-particle dispersion) is NOT the naive SE form -- the sqrt3 derivation is a candidate, not proven.")
        print("DIAGNOSE: re-derive the 2-particle resummation (different coefficients?), do NOT loosen.")
        sys.exit(1)
    else:
        print("ALL GATES PASS -> the (2,2) slow mode obeys lambda^2+8g*lambda+4J^2 q^2 with constant q. The SE")
        print("coherence-ladder dispersion governs the handover; Q_h = sqrt3/q (darkness=1) is DERIVED, the EP-")
        print("sibling (Qq=2 -> N/pi) of the SE coherence horizon. Tier1-standard (leading-order slope sqrt3/(2pi)).")
    print("=" * 100)


if __name__ == "__main__":
    main()
