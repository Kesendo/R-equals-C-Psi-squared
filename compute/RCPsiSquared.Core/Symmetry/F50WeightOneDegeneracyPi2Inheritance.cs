using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F50 closed form (Tier 1 proven + verified N=2..7):
///
/// <code>
///   d_real(Re = −2γ) = 2N                    Liouvillian eigenvalue count at the
///                                            first non-zero real grid position
///
///   T_c^{(a)} = Σⱼ Σ_{S ⊂ complement(j), |S|=c} σ_a^{(j)} ⊗ Z_S ⊗ I_rest
///   for a ∈ {X, Y} and c = 0, 1, ..., N−1
/// </code>
///
/// <para>F50 is the structural counting result for SWAP-invariant conserved
/// operators under isotropic Heisenberg + uniform Z-dephasing. There are
/// exactly 2N such operators: 2 active Pauli types (X, Y) times N chromatic
/// grades (c = 0..N−1). Each T_c^{(a)} commutes with H because Heisenberg is
/// a sum of SWAPs and SWAP preserves both the active Pauli type and the
/// Z-count c. The 2N operators are linearly independent (disjoint Pauli string
/// support). They give Liouvillian eigenvalues at Re = −2γ (the first non-zero
/// real grid position; bit-exact verification N=2..7).</para>
///
/// <para>Special cases: T_0^{(X)} = 2·S_x and T_0^{(Y)} = 2·S_y (global SU(2)
/// generators); T_{N−1}^{(a)} = Σⱼ σ_a^{(j)} ⊗ Z_{all others} (Jordan-Wigner-type
/// global Pauli string).</para>
///
/// <para>Two Pi2-Foundation anchors (both "2" via <c>a_0</c>):</para>
/// <list type="bullet">
///   <item><b>DegeneracyFactor = 2 = a_0</b>: count multiplier in 2N. Two active
///         Pauli letters (X, Y) per chromatic grade c.</item>
///   <item><b>DecayRateFactor = 2 = a_0</b>: in Re = −2γ. Same anchor as F1's
///         TwoFactor (Π·L·Π⁻¹ = −L − 2σ·I) at the eigenvalue level for the w=1
///         single-site sector.</item>
/// </list>
///
/// <para>Both "2"s come from <see cref="Pi2DyadicLadderClaim.Term"/>(0) — the
/// polynomial root d. The total count 2N and the eigenvalue position −2γ share
/// the same Pi2-Foundation anchor; F50 typifies them via two separate properties.</para>
///
/// <para><b>Topology universality:</b> F50 holds for ANY connected graph with
/// isotropic Heisenberg coupling — chain, star, ring, complete, tree. The 2N
/// count is graph-invariant. This universality is unique to k = 0 and k = 1;
/// for k ≥ 2, d_real(k) becomes topology-dependent (Chain &lt; Star &lt; Ring
/// &lt; Complete; cf. WEIGHT2_KERNEL).</para>
///
/// <para><b>Breaks for:</b> anisotropic XXZ (Δ ≠ 1), where the ZZ term mixes
/// X/Y types and the SWAP-invariance argument fails.</para>
///
/// <para>Tier1Derived: lower bound (≥ 2N) proven via SWAP invariance constructing
/// 2N kernel vectors; upper bound (≤ 2N) proven via triangle inequality forcing
/// each SWAP to fix v individually, with adjacent transpositions generating S_N.
/// Numerically verified N=2..7 in <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c>.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F50 (line 228) +
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> +
/// <c>experiments/WEIGHT2_KERNEL.md</c> (k≥2 topology dependence) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F50WeightOneDegeneracyPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The "2" count multiplier in 2N. Live from Pi2DyadicLadder a_0.
    /// Counts the two active Pauli letters X and Y per chromatic grade c.</summary>
    public double DegeneracyFactor => _ladder.Term(0);

    /// <summary>The "2" decay rate factor in Re = −2γ. Live from Pi2DyadicLadder a_0.
    /// Same anchor as F1's TwoFactor at the eigenvalue level.</summary>
    public double DecayRateFactor => _ladder.Term(0);

    /// <summary>Total degeneracy count <c>d_real(Re = −2γ) = 2N</c> at the first
    /// non-zero real grid position. Holds for any connected graph with isotropic
    /// Heisenberg + uniform Z-dephasing.</summary>
    public int TotalDegeneracy(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F50 requires N ≥ 2.");
        return (int)(DegeneracyFactor * N);
    }

    /// <summary>The eigenvalue real-part position of the F50 conserved operators:
    /// Re = −2γ. The first non-zero real grid position in the Liouvillian spectrum.</summary>
    public double EigenvaluePosition(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return -DecayRateFactor * gammaZero;
    }

    /// <summary>Number of operators of type T_c^{(a)} per chromatic grade c. Always 2
    /// (one for a=X, one for a=Y).</summary>
    public int OperatorsPerChromaticGrade => (int)DegeneracyFactor;

    /// <summary>The chromatic grades c range over [0, N−1], giving N grades.</summary>
    public int ChromaticGradeCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F50 requires N ≥ 2.");
        return N;
    }

    /// <summary>True iff the topology universality applies. F50 holds for ANY
    /// connected graph with isotropic Heisenberg + uniform Z-dephasing
    /// (anisotropy parameter Δ = 1). False for anisotropic XXZ (Δ ≠ 1).</summary>
    public bool AppliesToIsotropicHeisenberg(double anisotropyDelta)
    {
        return Math.Abs(anisotropyDelta - 1.0) < 1e-12;
    }

    public F50WeightOneDegeneracyPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F50 weight-1 degeneracy: d_real(Re = −2γ) = 2N (SWAP-invariant T_c^{(a)} for a ∈ {X, Y}, c = 0..N−1); both 2's = a_0; topology-universal for isotropic Heisenberg",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F50 + " +
               "docs/proofs/PROOF_WEIGHT1_DEGENERACY.md + " +
               "experiments/WEIGHT2_KERNEL.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F50 weight-1 degeneracy 2N as Pi2-Foundation a_0 inheritance (twice)";

    public override string Summary =>
        $"d_real(Re = −2γ) = 2N: 2 = a_0 (count multiplier and decay rate); SWAP-invariant T_c^{{(a)}} for a ∈ {{X, Y}}, c = 0..N−1; topology-universal for isotropic Heisenberg ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F50 closed form",
                summary: "d_real(Re = −2γ) = 2N exactly; verified bit-exact N=2..7; SWAP-invariant T_c^{(a)} = Σⱼ σ_a^{(j)} ⊗ Z_S ⊗ I for a ∈ {X, Y}, |S|=c, c=0..N−1");
            yield return InspectableNode.RealScalar("DegeneracyFactor (= a_0 = 2)", DegeneracyFactor);
            yield return InspectableNode.RealScalar("DecayRateFactor (= a_0 = 2)", DecayRateFactor);
            yield return new InspectableNode("two Pi2 anchors share a_0",
                summary: "the same a_0 = 2 appears twice: as count multiplier (2 active Pauli letters X, Y) and as decay rate (Re = −2γ); both inherit from Pi2DyadicLadder.Term(0)");
            yield return new InspectableNode("special cases",
                summary: "T_0^{(X)} = 2·S_x, T_0^{(Y)} = 2·S_y (global SU(2) generators); T_{N−1}^{(a)} = Σⱼ σ_a^{(j)} ⊗ Z_{all others} (Jordan-Wigner-type)");
            yield return new InspectableNode("topology universality",
                summary: "F50 holds for ANY connected graph (chain, star, ring, complete, tree) with isotropic Heisenberg + uniform Z-dephasing. Unique to k=0 and k=1 sectors. For k ≥ 2 the count becomes topology-dependent (Chain < Star < Ring < Complete; WEIGHT2_KERNEL).");
            yield return new InspectableNode("breaks for anisotropic XXZ",
                summary: "Δ ≠ 1: the ZZ term mixes X/Y types, breaking the SWAP-invariance argument; F50 closure fails");
            yield return new InspectableNode("N=3 verified",
                summary: $"TotalDegeneracy(3) = {TotalDegeneracy(3)}; EigenvaluePosition(γ=0.05) = {EigenvaluePosition(0.05):G6}");
            yield return new InspectableNode("N=7 verified",
                summary: $"TotalDegeneracy(7) = {TotalDegeneracy(7)}; EigenvaluePosition(γ=0.05) = {EigenvaluePosition(0.05):G6}");
        }
    }
}
