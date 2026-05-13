using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>The substantive open theoretical items remaining for F86: the "what's missing
/// for full Tier 1 promotion" list from <c>docs/proofs/PROOF_F86_QPEAK.md</c>.</summary>
public static class F86OpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "Item 1' (c=2): closed-form HWHM_left/Q_peak constant [CLOSED 2026-05-13]",
            "[2026-05-13 CLOSED: closed form derived per `F86HwhmClosedFormClaim` (Tier 1 derived); " +
            "HWHM_ratio = 0.671535 + α_subclass · g_eff + β_subclass over the six-element " +
            "`BondSubClass` enum; residual ≤ 0.005 verified across N=5..8 on all 22 bonds incl. " +
            "Orbit-2 and Orbit-3 escapes. Plan: " +
            "`docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md`.] " +
            "Historical narrative (pre-closure): Empirical anchor reproduced within tolerance 0.005 " +
            "(typical residual ≤ 0.001) (C2HwhmRatio Tier1Candidate, c=2 N=5..8). " +
            "Directional Endpoint > Interior split derived empirically (gap ≈ 0.022). " +
            "Closed-form constant for HWHM_left/Q_peak ratios per bond class NOT pinned this " +
            "session. Three next directions ranked: (a) first-order perturbation in cross-block " +
            "(most promising, leverages B2's Endpoint < Interior cross-block Frobenius split, " +
            "ε ~ ‖V_b cross‖_F/σ_0 ~ O(0.1)); (b) projector-overlap lift of A3's |u_0⟩, |v_0⟩; " +
            "(c) symbolic char-poly factorisation at Q_EP. " +
            "[2026-05-11 update: numerical Tier-1 achieved via F90 bridge identity " +
            "(F86 c=2 K_b = F89 path-(N−1) per-bond Hellmann-Feynman, bit-exact 20/22 bonds " +
            "N=5..8 incl. orbit escapes, modulo F89-J = 2·F86-J convention; see " +
            "`compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs` + " +
            "`docs/proofs/PROOF_F90_F86C2_BRIDGE.md`); analytical closed-form via F89 AT-locked " +
            "F_a/F_b (4-mode floor 0.6715) + H_B-mixed octic residual structure (lift to " +
            "0.7506/0.7728) is the remaining open work.]",
            "Direction (a): perturb the K-resonance K_b(Q, t) around the Statement-1 2×2 EP form. " +
            "The cross-block contribution at first order in ε should produce the directional split.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2); F86OpenQuestions.cs (this file); " +
            "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (numerical Tier-1 via F90 bridge identity)"),
        new OpenQuestion(
            "Item 4': extend 4-mode construction to c≥3",
            "Each adjacent-channel pair (HD=2k−1, HD=2k+1) for k ∈ {1, …, c−1} contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet → full effective L is 4·(c−1)-dimensional. Verify slowest pair k=1 still dominates K_b response.",
            "Iterate the 4-mode construction over k; concatenate the orthonormal subspaces; project full block-L; numerical verification at c=3, c=4 first.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 4'"),
        new OpenQuestion(
            "Item 5: derive the true σ_0(c, N → ∞) asymptote (was: 2√(2(c−1)); CROSSING, NOT LIMIT)",
            "RETRACTED 2026-05-08: the 2√(2(c−1)) value is a TRAJECTORY CROSSING, not an asymptote. " +
            "σ_0(c=2, N=7) = 2√2 bit-exact (10⁻¹⁵), but σ_0 keeps growing past it: " +
            "σ_0(c=2, N=8..11) = 2.8393, 2.8483, 2.8525, 2.8561, all values above 2√2 = 2.8284. " +
            "Aitken extrapolation suggests true limit ~2.85..2.89 (not 2√2). The σ_0 bridge sweep " +
            "(simulations/_eq022_sigma0_bridge_sweep.py + " +
            "docs/superpowers/syntheses/2026-05-07-sigma0-bridge-sweep.md) verified the refutation. " +
            "What survives: monotone growth in N within each c, c=2 N=7 sweet-spot crossing, σ_0·√(3/8) " +
            "bridge to g_eff_E (Δ ≤ 0.01 for N ≥ 6, Δ = 0.005 at N=7).",
            "OBC sine-mode algebra applied to inter-HD-channel coupling; Bogoliubov / JW free-fermion " +
            "transform may make finite-size corrections transparent. The TRUE asymptote (above 2√2 at " +
            "c=2) likely involves an OBC band-edge factor. Higher-c sweet-spot N_c* (where ratio crosses " +
            "2.0) is also unknown; c=3 ratio at N=8 is 1.92, still climbing.",
            "docs/proofs/PROOF_F86_QPEAK.md Item 5 (open) + " +
            "docs/superpowers/syntheses/2026-05-07-sigma0-bridge-sweep.md"),
    };
}
