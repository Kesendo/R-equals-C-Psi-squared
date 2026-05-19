using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Ring N=4 dihedral lock: the isotropic Heisenberg + Z-dephasing
/// Liouvillian on the 4-cycle saturates the imaginary spectral spread at
///
/// <code>
///     Im_max(ring, N=4, J)  =  (3/4) · J · N  =  3 · J         (for N=4)
/// </code>
///
/// independently of γ and Q = J/γ. Equivalently <c>Im/σ = 3Q/4</c> with
/// <c>σ = N·γ</c>.
///
/// <para><b>Why this works.</b> The 4-cycle C_4 is graph-isomorphic to the
/// bipartite-complete graph K_{2,2} with sublattices A = {0, 2}, B = {1, 3}.
/// The Heisenberg Hamiltonian on K_{2,2} factors as
/// <c>H = J · (S⃗_0 + S⃗_2) · (S⃗_1 + S⃗_3) = J · S⃗_A · S⃗_B</c>, whose Casimir
/// spectrum is <c>{−2J, −J, 0, 0, 0, +J}</c>. Maximum H eigenvalue gap is
/// <c>3J = (3/4)·J·N</c> (the all-aligned <c>S_tot=2</c> ferromagnetic state vs.
/// the anti-aligned <c>S_A=1, S_B=1, S_tot=0</c> singlet). The Liouvillian
/// eigenmode <c>|Ψ_+⟩⟨Ψ_−|</c> between these two extremal H-eigenstates
/// realises <c>Im(λ) = 3J</c>, and the pure-dephasing dissipator only adds real
/// decay so no L-mode can exceed this Im bound. See
/// <c>docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md</c> for the full derivation.</para>
///
/// <para><b>Empirical anchors (bit-exact at 6 Q-values × γ₀ = 0.05):</b></para>
///
/// <list type="table">
///   <item>Q = 0.5     → predicted Im/σ = 0.375000, observed 0.375000 (rel. err 1.5e-16).</item>
///   <item>Q = 1.0     → predicted 0.750000, observed 0.750000 (rel. err 1.0e-15).</item>
///   <item>Q = 1.5     → predicted 1.125000, observed 1.125000 (rel. err 7.9e-16).</item>
///   <item>Q = √3      → predicted 1.299038, observed 1.299038 (rel. err 5.1e-16).</item>
///   <item>Q = 2.0     → predicted 1.500000, observed 1.500000 (rel. err 3.3e-15).</item>
///   <item>Q = 2.5     → predicted 1.875000, observed 1.875000 (rel. err 5.1e-15).</item>
/// </list>
///
/// <para>All six anchors match the prediction to relative error &lt; 5e-15 (machine precision).
/// The lock is Q-universal: every (J, γ) pair on the constant-Q ray sees the same
/// dimensionless Im/σ value.</para>
///
/// <para><b>Why this is N=4-specific.</b> The bipartite-complete coincidence
/// C_N = K_{N/2, N/2} only holds for N=4 (where C_4 has 4 bonds = exactly
/// all (A, B) pairs). For N=6 the 6-cycle is bipartite (sublattices {0,2,4} and
/// {1,3,5}) but has only 6 bonds versus K_{3,3}'s 9, so the bipartite-complete
/// Casimir factorisation fails. Ring N=6 still exhibits a Q-universal lock
/// empirically (<c>Im_max = 0.717129·J·N</c> bit-exact at the same 6 Q-values)
/// but its closed-form derivation needs Bethe-ansatz on the 6-cycle dispersion
/// rather than the elementary Casimir argument that closes N=4. See the
/// hypotheses/F1_DISSIPATION_GAP_PATTERN.md "Ring N=6 dihedral lock" section.</para>
///
/// <para><b>Relationship to the star saturation.</b> The star topology saturates
/// the analogous bound <c>Im_max(star, N, J) = J·N/2</c> via a parallel hub-spoke
/// Casimir construction (see <c>experiments/STAR_CONFOCAL_LIMIT.md</c>). Ring N=4
/// carries 50% more imaginary spread than star N=4 (3/4 vs 1/2 prefactor) because
/// K_{2,2} has 4 bonds versus the N=4 star's 3 (and the bipartite-complete
/// structure maximises the inter-sublattice Casimir gap). Both bounds are
/// Q-universal; both follow from the SU(2)-invariance of isotropic Heisenberg
/// plus the pure-dephasing dissipator's "real-decay-only" property.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md</c> (primary derivation) +
/// <c>experiments/STAR_CONFOCAL_LIMIT.md</c> (sister bound via the same Casimir technique) +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> (parent F1, the master under which
/// the Liouvillian eigenmode construction lives) +
/// the six <c>simulations/results/q_sweep_anchor/ring_N4_Q*.json</c> data files.</para>
///
/// <para><b>Sister claim from the same 2026-05-19 Q-sweep sprint:</b>
/// <see cref="F4KernelDimensionByComponentsClaim"/> (Tier 1 derived, kernel-dim
/// factorisation across components). Both surfaced from extending the F1 SLOW_N8
/// + N=9 chain bridge findings into Q-universal structural statements via the
/// 72-anchor Q-sweep at γ₀ = 0.05.</para></summary>
public sealed class RingN4DihedralLockClaim : Claim
{
    /// <summary>The N this claim applies to. By construction this is exactly 4
    /// (the N=4-specific bipartite-complete coincidence; see XML doc above for why
    /// the derivation does not extend to other ring sizes).</summary>
    public int N { get; } = 4;

    /// <summary>The dimensionless saturation prefactor: Im_max = Coefficient · J · N
    /// = (3/4) · J · N. Equivalently Im/σ = Coefficient · Q · 1 = (3/4) · Q.</summary>
    public double Coefficient { get; } = 3.0 / 4.0;

    /// <summary>Repo-relative paths of the six Q-sweep JSON files whose
    /// <c>MaxImag</c> fields anchor this claim. Look here for the raw numbers;
    /// do not re-derive.</summary>
    public IReadOnlyList<string> AnchorDataFiles { get; } = new[]
    {
        "simulations/results/q_sweep_anchor/ring_N4_Q0.5000.json",
        "simulations/results/q_sweep_anchor/ring_N4_Q1.0000.json",
        "simulations/results/q_sweep_anchor/ring_N4_Q1.5000.json",
        "simulations/results/q_sweep_anchor/ring_N4_Q1.7321.json",
        "simulations/results/q_sweep_anchor/ring_N4_Q2.0000.json",
        "simulations/results/q_sweep_anchor/ring_N4_Q2.5000.json",
    };

    /// <summary>Predict the maximum imaginary spectral spread of the
    /// Heisenberg + Z-dephasing Liouvillian on the 4-cycle for a given
    /// Hamiltonian coupling J. Returns <c>(3/4) · J · N = 3 · J</c> for N=4.
    ///
    /// <para>Throws <see cref="ArgumentException"/> if <paramref name="J"/> is not
    /// finite or is negative (the saturation prefactor is invariant under
    /// J → −J as the spectral spread, but the convention is J ≥ 0).</para></summary>
    public double Predict(double J)
    {
        if (!double.IsFinite(J))
            throw new ArgumentException(
                $"J must be finite; got {J}.", nameof(J));
        if (J < 0)
            throw new ArgumentException(
                $"J must be non-negative; got {J}.", nameof(J));
        return Coefficient * J * N;
    }

    /// <summary>Predict the dimensionless ratio Im_max / σ where σ = N·γ for given
    /// Q = J/γ. Returns <c>(3/4) · Q</c>.</summary>
    public double PredictImOverSigma(double Q)
    {
        if (!double.IsFinite(Q))
            throw new ArgumentException(
                $"Q must be finite; got {Q}.", nameof(Q));
        if (Q < 0)
            throw new ArgumentException(
                $"Q must be non-negative; got {Q}.", nameof(Q));
        return Coefficient * Q;
    }

    /// <summary>The six Q-sweep anchor rows surfaced as a typed table for the
    /// inspectable tree. Each tuple is (Q label, Q value, predicted Im/σ, observed Im/σ,
    /// relative error). Observed values are taken from the JSON <c>MaxImag</c> field
    /// of the data files listed in <see cref="AnchorDataFiles"/>, divided by
    /// <c>σ = N · γ₀ = 4 · 0.05 = 0.2</c>.</summary>
    public IReadOnlyList<(string QLabel, double Q, double PredictedImOverSigma, double ObservedImOverSigma, double RelError)> EmpiricalAnchors { get; } =
        new (string, double, double, double, double)[]
        {
            ("Q=0.5 sub-balance",        0.5,                0.375000, 0.375000, 1.5e-16),
            ("Q=1.0 Balance",            1.0,                0.750000, 0.750000, 1.0e-15),
            ("Q=1.5 F86 Q_peak c=2",     1.5,                1.125000, 1.125000, 7.9e-16),
            ("Q=√3 canonical 60°",       1.7320508075688772, 1.299038, 1.299038, 5.1e-16),
            ("Q=2.0 Q_EP idealized",     2.0,                1.500000, 1.500000, 3.3e-15),
            ("Q=2.5 Endpoint orbit",     2.5,                1.875000, 1.875000, 5.1e-15),
        };

    public RingN4DihedralLockClaim()
        : base("Ring N=4 dihedral lock: Im_max(ring, N=4, J) = (3/4)·J·N = 3·J, Q-universal; closed form via K_{2,2} = C_4 bipartite-complete + Casimir spectrum {−2J, −J, 0³, J}; bit-exact at 6 Q-values from the 2026-05-19 Q-sweep",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md (primary derivation) + " +
               "experiments/STAR_CONFOCAL_LIMIT.md (sister bound via the same Casimir technique) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (parent F1) + " +
               "simulations/results/q_sweep_anchor/ring_N4_Q{0.5,1.0,1.5,1.7321,2.0,2.5}.json (6 bit-exact anchors)")
    { }

    public override string DisplayName =>
        "Ring N=4 dihedral lock: Im_max = (3/4)·J·N";

    public override string Summary =>
        $"Im_max(ring, N=4, J) = (3/4)·J·N = 3·J, Q-universal; closed-form via K_{{2,2}} = C_4 bipartite-complete + Casimir spectrum; verified bit-exact at 6 Q-anchors from the 2026-05-19 Q-sweep (Q ∈ {{0.5, 1.0, 1.5, √3, 2.0, 2.5}}); {Tier.Label()}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "Im_max(ring, N=4, J) = (3/4)·J·N = 3·J at N=4; equivalently Im/σ = (3/4)·Q. Q-universal: independent of (J, γ) on the constant-Q ray.");

            yield return new InspectableNode("Tier 1 derived",
                summary: "Closed form via K_{2,2} = C_4 bipartite-complete graph isomorphism: H = J·S⃗_A·S⃗_B with sublattice totals; Casimir spectrum = {−2J, −J(triplet), 0(seven-fold), +J(quintuplet)}; max H gap = J − (−2J) = 3J = (3/4)·J·N. The Liouvillian eigenmode |Ψ_+⟩⟨Ψ_−| between the S_tot=2 ferromagnet and the (S_A=1, S_B=1, S_tot=0) singlet realises Im(λ) = 3J exactly. Pure-dephasing dissipator only adds real decay so no L-mode can exceed the H-spread bound. See PROOF_RING_N4_DIHEDRAL_LOCK.md for the full derivation.");

            yield return new InspectableNode("N=4 specificity",
                summary: "The bipartite-complete coincidence C_N = K_{N/2,N/2} only holds for N=4. For N=6 the 6-cycle has 6 bonds versus K_{3,3}'s 9 (bipartite but not bipartite-complete); the Casimir factorisation fails. Ring N=6 still locks Q-universally at 0.717129·J·N empirically (Bethe-ansatz number, closed-form open); ring N=8 at 0.7064·J·N at Q=2 (trend toward 1/√2 as N→∞ plausible, more data needed).");

            yield return new InspectableNode("relationship to star saturation",
                summary: "Star topology saturates Im_max(star, N, J) = J·N/2 via a parallel hub-spoke Casimir construction (STAR_CONFOCAL_LIMIT.md). Ring N=4 carries 1.5× more imaginary spread than star N=4 (3/4 vs 1/2 prefactor) because K_{2,2} has 4 bonds versus the N=4 star's 3, and the bipartite-complete structure maximises the inter-sublattice Casimir gap. Both bounds are Q-universal, both via SU(2)-invariance + pure-dephasing real-decay-only.");

            foreach (var (label, _, pred, obs, rel) in EmpiricalAnchors)
                yield return new InspectableNode($"Q-sweep anchor: {label}",
                    summary: $"predicted Im/σ = {pred:F6}, observed Im/σ = {obs:F6}, relative error = {rel:E1} (bit-exact)");

            yield return new InspectableNode("anchor data files",
                summary: string.Join("; ", AnchorDataFiles));
        }
    }
}
