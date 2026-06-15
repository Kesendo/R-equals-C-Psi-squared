"""The off-gap argument for the two clocks (clock_hand_ladder, the open piece of ClockHandLadderClaim).

THE QUESTION. ClockHandLadderClaim (Tier1Candidate) says the coherence hand ω_mem = the F2b band edge
E1 = 2J·cos(π/(N+1)) for N≥3. The mechanism (the |vac⟩⟨ψ_k| modes are simultaneous eigenoperators of L_D at
−2γ and L_H at iE_k) is DERIVED; what is only VERIFIED N=3-5 is the GAP-DOMINANCE: that no mode at the exact
gap rate 2γ oscillates faster than the band edge. That is the open remainder.

WHAT IS CONFIRMED HERE (gate-first, full 4^N Liouvillian, N=3,4,5):
  (ii) max|Im| among the modes at EXACTLY Re=−2γ equals the band edge E1 (√2 / φ / √3 at N=3/4/5);
  the vacuum-anchored |vac⟩⟨ψ_k| ladder are exact eigenoperators at −2γ+iE_k (residual ~1e-16);
  (i)  the Q-FLOOR is the COHERENCE HORIZON Q*(N): the band edge is the gap (the coherence hand) only
       above Q*(N) (1.879 / 2.374 at N=4/5, → 2N/π); the numerical floor brackets Q*(N) exactly. Below it
       the overdamped {0,2}-coherence diffusion mode is slower than the −2γ band edge AND non-oscillating,
       so it takes the gap and the clock freezes. So regime (i) is RESOLVED by thread (a)
       (PROOF_COHERENCE_HORIZON_SLOPE); only (ii) remains the open proof.

THE REFRAMING (the advance over the first sketch). The first guess — that the only exact −2γ eigenoperators
are the vacuum/full-anchored ladders — is INCOMPLETE: there are 22 / 32 / 50 exact −2γ modes at N=3/4/5, more
than the 4N (12/16/20) the ladders give. The correct frame is the Absorption Theorem (PROOF_ABSORPTION_THEOREM,
Tier-1 derived): Re(λ) = −2γ·⟨n_XY⟩, so the exact −2γ eigenoperators are precisely the ⟨n_XY⟩=1 subspace
(operators with exactly one X/Y Pauli factor). The vacuum/full ladders are a SUBSET of it; the full exact-2γ
set is the L_H-invariant part of the n_XY=1 Pauli subspace (dim N·2^N).

THE OPEN TARGET (precise). Prove that the maximum oscillation frequency inside the protected n_XY=1 subspace
is exactly E1 = 2J·cos(π/(N+1)). The subtlety: L_H = [H, ·] moves a lone X/Y along the chain like a free
particle (n_XY=1 preserved, frequency a single-particle band energy E_k ≤ E1), BUT it can also turn a
background Z into an X/Y pair (n_XY → 3, e.g. [X_iX_{i+1}, Z_i] = −2i·Y_i·X_{i+1}), so the n_XY=1 subspace is
NOT closed under L_H and the dephased chain is genuinely INTERACTING (no clean free-fermion shortcut). The
vacuum/full anchoring is one way to be protected (H annihilates |vac⟩/|full⟩, so no background to activate);
characterizing ALL the protected n_XY=1 modes — and bounding their frequency by E1 — is the remaining proof.

Self-validating: the gates below assert. Reuses coherence_horizon_se_block._L_full / sine_U."""
import numpy as np
from coherence_horizon_se_block import _L_full, sine_U, qstar_se


def band_edge(N, J=1.0):
    return 2.0 * J * np.cos(np.pi / (N + 1))


def single_exc_state(N, amps):
    """sum_i amps[i]·|site i excited>, in _L_full's basis (site 0 most significant, excite = |1> = Z eig −1)."""
    v = np.zeros(2 ** N, complex)
    for i in range(N):
        v[1 << (N - 1 - i)] = amps[i]
    return v


def exact_2g_max_im(N, J=1.0, g=0.1):
    """(#exact−2γ modes, max|Im| over them). The gate: max|Im| must equal the band edge E1."""
    ev = np.linalg.eigvals(_L_full(N, J, g))
    at = ev[np.abs(ev.real + 2.0 * g) < 1e-6]
    return len(at), (float(np.abs(at.imag).max()) if len(at) else 0.0)


def vac_ladder_residual(N, J=1.0, g=0.1):
    """Worst residual of |vac⟩⟨ψ_k| as an eigenoperator at −2γ+iE_k (should be ~machine zero)."""
    L = _L_full(N, J, g)
    vac = np.zeros(2 ** N, complex); vac[0] = 1.0
    U = sine_U(N)
    worst = 0.0
    for k in range(N):
        rho = np.outer(vac, single_exc_state(N, U[:, k]).conj())   # |vac⟩⟨ψ_k|
        v = rho.flatten()                                          # row-stacking vec (the _L_full convention)
        lam = -2.0 * g + 1j * 2.0 * J * np.cos((k + 1) * np.pi / (N + 1))
        worst = max(worst, float(np.linalg.norm(L @ v - lam * v) / np.linalg.norm(v)))
    return worst


def q_floor_bracket(N, J=1.0):
    """(Q_fail, Q_hold): descending in Q, the band edge is the gap mode at Q_hold but not at the next-lower
    grid point Q_fail. The floor sits between them (and equals the Coherence Horizon Q*(N))."""
    be = band_edge(N, J)
    last_ok = None
    for Q in [20, 10, 6, 4, 3, 2.5, 2.25, 2, 1.85, 1.7, 1.5, 1.3, 1.1]:
        g = J / Q
        ev = np.linalg.eigvals(_L_full(N, J, g))
        nz = ev[ev.real < -1e-9]
        gap = nz.real.max()
        clock = float(np.abs(nz[np.abs(nz.real - gap) < 1e-6].imag).max())
        if abs(gap + 2 * g) < 1e-6 and abs(clock - be) < 1e-5:
            last_ok = Q
        elif last_ok is not None:
            return Q, last_ok         # (fail, hold): the floor is between
    return None, last_ok


def main():
    print("(ii) exact-(−2γ) max|Im| == band edge, and the vacuum ladder is an exact eigenoperator:")
    for N in (3, 4, 5):
        cnt, mx = exact_2g_max_im(N)
        be = band_edge(N)
        res = vac_ladder_residual(N)
        print(f"  N={N}: {cnt:3d} exact-(−2γ) modes; max|Im|={mx:.6f} = band edge {be:.6f}; "
              f"vac-ladder residual {res:.1e}")
        assert abs(mx - be) < 1e-5, f"N={N}: max|Im| at −2γ is not the band edge"
        assert res < 1e-10, f"N={N}: the vacuum ladder is not an exact −2γ eigenoperator"

    print("(i) the Q-floor = the Coherence Horizon Q*(N) (band edge is the gap only above it):")
    for N in (4, 5):
        qfail, qhold = q_floor_bracket(N)
        qstar = qstar_se(N)
        print(f"  N={N}: floor in ({qfail}, {qhold}); Coherence Horizon Q*(N) = {qstar:.4f} sits inside "
              f"-> regime (i) is thread (a)")
        assert qfail < qstar < qhold, f"N={N}: the gap-dominance floor must equal the Coherence Horizon Q*(N)"

    print("\nOK: (i) the gap-dominance Q-floor = the Coherence Horizon Q*(N) (resolved, thread a). "
          "(ii) max|Im| at the gap = band edge confirmed N=3,4,5; the general proof (max frequency in the "
          "protected n_XY=1 subspace = E1) is the open remainder -- see the module docstring.")


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")   # Windows console: allow the Unicode in the report
    main()
