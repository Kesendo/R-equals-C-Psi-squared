using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Ring N=4 dihedral lock: the isotropic Heisenberg + Z-dephasing
/// Liouvillian on the 4-cycle saturates the imaginary spectral spread at
///
/// <code>
///     Im_max(ring, N=4, J)  =  3 · J  =  (3/4) · J · N         (at N=4)
/// </code>
///
/// independently of γ and Q = J/γ. Equivalently <c>Im/σ = 3Q/4</c> with
/// <c>σ = N·γ</c>.
///
/// <para><b>Why this works.</b> The 4-cycle C_4 is graph-isomorphic to the
/// bipartite-complete graph K_{2,2} with sublattices A = {0, 2}, B = {1, 3}.
/// The Heisenberg Hamiltonian on K_{2,2} factors as
/// <c>H = J · (S_0 + S_2) · (S_1 + S_3) = J · S_A · S_B</c>, whose Casimir
/// spectrum is <c>{−2J, −J³, 0⁷, J⁵}</c> (multiplicities sum to 16 = 2⁴).
/// Maximum H eigenvalue gap is
/// <c>3J = (3/4)·J·N</c> (the all-aligned <c>S_tot=2</c> ferromagnetic state vs.
/// the anti-aligned <c>S_A=1, S_B=1, S_tot=0</c> singlet). The realising
/// Liouvillian eigenmode is <c>|0000⟩⟨Ψ_−|</c> with <c>λ = −4γ − 3i·J</c>, and
/// its conjugate: the ket must be one of the two Z-product members of the
/// quintuplet, since only then does the coherence carry a constant
/// <c>popcount(i⊕j)</c> and uniform dephasing act on it as a scalar. The upper
/// bound holds because D is Hilbert-Schmidt self-adjoint, the jump operators
/// Z_l being Hermitian; "adds only real decay" is not sufficient. See
/// <c>docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md</c> for the full derivation.</para>
///
/// <para><b>Empirical anchors (6 Q-values × γ₀ = 0.05, machine precision):</b></para>
///
/// <list type="table">
///   <item>Q = 0.5     → predicted Im/σ = 0.375000, observed 0.375000 (rel. err 1.2e-15).</item>
///   <item>Q = 1.0     → predicted 0.750000, observed 0.750000 (rel. err 1.5e-16).</item>
///   <item>Q = 1.5     → predicted 1.125000, observed 1.125000 (rel. err 7.9e-16).</item>
///   <item>Q = √3      → predicted 1.299038, observed 1.299038 (rel. err 1.0e-15).</item>
///   <item>Q = 2.0     → predicted 1.500000, observed 1.500000 (rel. err 3.3e-15).</item>
///   <item>Q = 2.5     → predicted 1.875000, observed 1.875000 (rel. err 5.1e-15).</item>
/// </list>
///
/// <para>All six anchors match the prediction to relative error &lt; 1e-14 (machine
/// precision; the largest is 5.1e-15, the smallest 1.5e-16, none exactly zero). The
/// errors above are computed from the two-sided <c>max|Im λ|</c> the law is stated
/// for, read off the stored sweep as <c>max(MaxImag, |MinImag|)</c>; at Q = 0.5, 1.0
/// and √3 it is <c>|MinImag|</c> that is marginally the larger, so those three differ
/// from the sweep's one-sided <c>MaxImag</c> reading. Both readings agree with (3/4)·Q
/// to machine precision, so the distinction lives in the last bit and moves nothing
/// else. The lock is Q-universal: every (J, γ) pair on the constant-Q ray sees the
/// same dimensionless Im/σ value.</para>
///
/// <para><b>Why this is N=4-specific.</b> The bipartite-complete coincidence
/// C_N = K_{N/2, N/2} only holds for N=4 (where C_4 has 4 bonds = exactly
/// all (A, B) pairs). For N=6 the 6-cycle is bipartite (sublattices {0,2,4} and
/// {1,3,5}) but has only 6 bonds versus K_{3,3}'s 9, so the bipartite-complete
/// Casimir factorisation fails. Ring N=6 still exhibits a Q-universal lock
/// (<c>Im_max = 0.717129·J·N</c> at the same 6 Q-values, machine precision)
/// and its constant is still elementary, <c>c_6 = (5+√13)/12</c> from factoring the
/// S_z=0 characteristic polynomial over ℚ, though not by the Casimir argument that
/// closes N=4. Bethe ansatz is what the N → ∞ limit needs (<c>c_∞ = ln 2</c>),
/// not the individual N. See the
/// hypotheses/F1_DISSIPATION_GAP_PATTERN.md "Ring N=6 dihedral lock" section.</para>
///
/// <para><b>Relationship to the star saturation.</b> The star topology saturates
/// the analogous bound <c>Im_max(star, N, J) = J·N/2</c> via a parallel hub-spoke
/// Casimir construction (see <c>experiments/STAR_CONFOCAL_LIMIT.md</c>). Ring N=4
/// carries 50% more imaginary spread than star N=4 (3/4 vs 1/2 prefactor), and the
/// two sit at the two ends of one range: J·N/2 is the MINIMUM of ΔE_max over
/// connected graphs with the star its unique minimiser (F147, searched exhaustively
/// at N ≤ 6, not proved past it), while at N=4 the 4-cycle sits at the MAXIMUM. Neither end is caused by bipartite-completeness or
/// by bond count. Ten of the 38 connected labelled graphs on four vertices reach 3J,
/// among them K_4, which is not bipartite and whose Casimir spectrum
/// <c>{−1.5J², −0.5J⁹, +1.5J⁵}</c> is a different multiplet structure altogether;
/// and the maximality does not survive N=5 (ring 3.1180339887 against K_5's 4).
/// Both closed forms are ΔE_max values; what turns either into an Im_max is F148,
/// which holds on every graph. "Adds only real decay" is not the reason and is not
/// sufficient: the bound rests on D being Hilbert-Schmidt self-adjoint, which is
/// what the Hermiticity of the Z_l buys.</para>
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
    /// = (3/4) · J · N. Equivalently Im/σ = Coefficient · Q = (3/4) · Q.</summary>
    public double Coefficient { get; } = 3.0 / 4.0;

    /// <summary>Repo-relative paths of the six Q-sweep JSON files whose
    /// <c>MaxImag</c> and <c>MinImag</c> fields anchor this claim (the law is the
    /// two-sided maximum of the two). Look here for the raw numbers;
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
        CasimirBoundClaimHelpers.RequireFiniteNonNegative(J, nameof(J));
        return Coefficient * J * N;
    }

    /// <summary>Predict the dimensionless ratio Im_max / σ where σ = N·γ for given
    /// Q = J/γ. Returns <c>(3/4) · Q</c>.</summary>
    public double PredictImOverSigma(double Q)
    {
        CasimirBoundClaimHelpers.RequireFiniteNonNegative(Q, nameof(Q));
        return Coefficient * Q;
    }

    /// <summary>The six Q-sweep anchor rows surfaced as a typed table for the
    /// inspectable tree. Each tuple is (Q label, Q value, predicted Im/σ, observed Im/σ,
    /// relative error). Q labels and Q values come from
    /// <see cref="CasimirBoundClaimHelpers.QSweepAnchorLabels"/> (shared with
    /// <see cref="StarImMaxBoundClaim.EmpiricalAnchors"/>). Observed values are taken
    /// from the JSON <c>MaxImag</c> field of the data files listed in
    /// <see cref="AnchorDataFiles"/>, as the two-sided <c>max(MaxImag, |MinImag|)</c>
    /// the claim is actually about, divided by <c>σ = N · γ₀ = 4 · 0.05 = 0.2</c>.
    /// They are the measured numbers, not the predicted ones, and the relative error
    /// is computed from them rather than stored.</summary>
    public IReadOnlyList<(string QLabel, double Q, double PredictedImOverSigma, double ObservedImOverSigma, double RelError)> EmpiricalAnchors { get; }
        = BuildEmpiricalAnchors();

    // Observed Im/σ per anchor, read off the stored Q-sweep JSON at γ₀ = 0.05 as
    // max(MaxImag, |MinImag|) / 0.2 (one per Q in CasimirBoundClaimHelpers.QSweepAnchorLabels,
    // in order). All at machine precision against (3/4)·Q, none exactly equal to it. Reproduced by
    // simulations/ring_n4_lock_gate.py §8.
    private static readonly double[] _observedImOverSigma =
        new[] { 0.37500000000000044, 0.7499999999999999, 1.1250000000000009,
                1.2990381056766593, 1.5000000000000049, 1.8750000000000095 };

    private static (string, double, double, double, double)[] BuildEmpiricalAnchors()
    {
        var labels = CasimirBoundClaimHelpers.QSweepAnchorLabels;
        var rows = new (string, double, double, double, double)[labels.Count];
        for (int i = 0; i < labels.Count; i++)
        {
            var (label, q) = labels[i];
            double pred = (3.0 / 4.0) * q;
            double obs = _observedImOverSigma[i];
            rows[i] = (label, q, pred, obs, Math.Abs(obs - pred) / pred);
        }
        return rows;
    }

    public RingN4DihedralLockClaim()
        : base("Ring N=4 dihedral lock: Im_max(ring, N=4, J) = 3·J = (3/4)·J·N at N=4, Q-universal; closed form via K_{2,2} = C_4 bipartite-complete + Casimir spectrum {−2J, −J³, 0⁷, J⁵} (multiplicities sum to 16 = 2⁴); 6 Q-anchors at machine precision from the 2026-05-19 Q-sweep. The prefactor is fixed at N=4, not a family law: the ring c_N falls to ln 2 along even N and rises to it along odd N.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md (primary derivation) + " +
               "simulations/ring_n4_lock_gate.py (gate) + " +
               "experiments/STAR_CONFOCAL_LIMIT.md (sister bound via the same Casimir technique) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (parent F1) + " +
               "simulations/results/q_sweep_anchor/ring_N4_Q{0.5,1.0,1.5,1.7321,2.0,2.5}.json (6 anchors)")
    { }

    public override string DisplayName =>
        "Ring N=4 dihedral lock: Im_max = 3·J";

    public override string Summary =>
        $"Im_max(ring, N=4, J) = 3·J = (3/4)·J·N at N=4, Q-universal; closed-form via K_{{2,2}} = C_4 bipartite-complete + Casimir spectrum; verified at 6 Q-anchors to machine precision from the 2026-05-19 Q-sweep (Q ∈ {{0.5, 1.0, 1.5, √3, 2.0, 2.5}}); {Tier.Label()}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "Im_max(ring, N=4, J) = 3·J, which at N=4 is (3/4)·J·N; equivalently Im/σ = (3/4)·Q. Q-universal: independent of (J, γ) on the constant-Q ray. The ·N form is the star's normalisation borrowed for comparison, not a ring family law: ΔE_max(ring)/(J·N) is 0.75 at N=4, 0.623607 at N=5, 0.717129 at N=6.");

            yield return new InspectableNode("Tier 1 derived",
                summary: "Closed form via K_{2,2} = C_4 bipartite-complete graph isomorphism: H = J·S_A·S_B with sublattice totals; Casimir spectrum = {−2J, −J(triplet), 0(seven-fold), +J(quintuplet)}, multiplicities 1+3+7+5 = 16; max H gap = J − (−2J) = 3J. The realising Liouvillian eigenmode is |0000⟩⟨Ψ_−| with λ = −4γ − 3i·J, and its conjugate: the ket must be one of the two Z-product members of the quintuplet, since only a constant popcount(i⊕j) across the operator's support lets uniform dephasing act as a scalar (the M = 0, ±1 members give eigenoperator residuals 1.155 and 1.000). The upper bound holds because D is Hilbert-Schmidt self-adjoint, the Z_l being Hermitian; 'adds only real decay' is not sufficient. See PROOF_RING_N4_DIHEDRAL_LOCK.md.");

            yield return new InspectableNode("N=4 specificity",
                summary: "The bipartite-complete coincidence C_N = K_{N/2,N/2} only holds for N=4. For N=6 the 6-cycle has 6 bonds versus K_{3,3}'s 9 (bipartite but not bipartite-complete); the Casimir factorisation fails. Ring N=6 still locks Q-universally at 0.717129·J·N, and that constant is elementary, c_6 = (5+√13)/12 from factoring the S_z=0 characteristic polynomial over ℚ; ring N=8 at 0.7064·J·N at Q=2, the largest root of 512c³ − 640c² + 232c − 25. Along EVEN N the sequence c_4 = 0.75, c_6 = 0.7171, c_8 = 0.7064 decreases past 1/√2 ≈ 0.7071 between N=6 and N=8 and continues below it, to c_∞ = ln 2 = 0.693147 from the Hulthén per-bond ground energy, NOT 1/√2. The ODD rings are frustrated and approach the same limit from BELOW, rising: c_5 = 0.623607, c_7 = 0.657883, c_9 = 0.671922. Two branches, not one descent. Bethe ansatz is what the LIMIT needs, not the individual N.");

            yield return new InspectableNode("relationship to star saturation",
                summary: "The star has Im_max(star, N, J) = J·N/2 via a parallel hub-spoke Casimir construction (STAR_CONFOCAL_LIMIT.md). Ring N=4 carries 1.5× more imaginary spread than star N=4. The two sit at opposite ends of one range: J·N/2 is the MINIMUM of ΔE_max over connected graphs with the star its unique minimiser (F147, exhaustive at N ≤ 6 and open past it), and at N=4 the 4-cycle sits at the MAXIMUM. Neither end is caused by bipartite-completeness: ten of the 38 connected labelled graphs on four vertices reach 3J, K_4 among them, which is not bipartite and carries the different multiplet structure {−1.5J², −0.5J⁹, +1.5J⁵}; and the maximality does not survive N=5 (ring 3.1180339887 vs K_5's 4). What turns either ΔE_max into an Im_max is F148, which holds on every graph. Its attainment half needs uniform γ over the FULL single-site {Z_l} set plus a fully polarised H extreme: Hermitian jumps alone buy only the bound: at rate 0.5, Z_0 alone gives 2.8961214544 against 3J, and γ = (1, 0.5, 0.5, 0.5) gives 2.7247448714. Neither is γ-universal, so saturation and γ-independence fail together.");

            foreach (var (label, _, pred, obs, rel) in EmpiricalAnchors)
                yield return new InspectableNode($"Q-sweep anchor: {label}",
                    summary: $"predicted Im/σ = {pred:F6}, observed Im/σ = {obs:F6}, relative error = {rel:E1} (machine precision, non-zero)");

            yield return new InspectableNode("anchor data files",
                summary: string.Join("; ", AnchorDataFiles));
        }
    }
}
