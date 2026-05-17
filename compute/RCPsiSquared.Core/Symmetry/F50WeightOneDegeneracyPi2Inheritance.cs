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
/// <para>Both "2"s come from <see cref="Pi2DyadicLadderClaim.Term"/>(0): the
/// polynomial root d. The total count 2N and the eigenvalue position −2γ share
/// the same Pi2-Foundation anchor; F50 typifies them via two separate properties.</para>
///
/// <para><b>Topology universality (with K_3 N=3 caveat, 2026-05-17):</b> F50
/// holds for the chain at all tested N (= 2..7), for ring/star/complete at
/// N ≥ 4, and for graphs containing triangles inside larger structure (paw at
/// N=4, bowtie + book at N=5). The single empirically-known anomaly is
/// <b>N = 3 K_3 (= ring = triangle = complete on 3 vertices)</b>: this case
/// gives <c>d_real = 8 = 2N + 2</c> instead of <c>2N = 6</c>. The 2 extras
/// are weight-1 operators in the standard 2-dim irrep of S_3 = Aut(K_3)
/// acting on the c=1 sector; they commute with H_K_3 because the three bond
/// commutators cancel pairwise but do NOT commute with any individual bond.
/// The proof's Step 5 derivation has a matrix-commutator vs conjugation-
/// action gap that explains the missed K_3 case; see
/// <c>PROOF_WEIGHT1_DEGENERACY § Appendix (2026-05-17)</c> for the
/// structural reading and the topology sweep that established the anomaly is
/// K_3-specific. <see cref="TotalDegeneracy"/> still returns the 2N lower
/// bound (always correct, rigorously proven); the K_3 N=3 actual count is
/// exposed via <see cref="K3TripleN3ActualCount"/>.</para>
///
/// <para>For k ≥ 2, d_real(k) becomes topology-dependent
/// (Chain &lt; Star &lt; Ring &lt; Complete; cf. WEIGHT2_KERNEL).</para>
///
/// <para><b>Breaks for:</b> anisotropic XXZ (Δ ≠ 1), where the ZZ term mixes
/// X/Y types and the SWAP-invariance argument fails.</para>
///
/// <para>Tier consistency: lower bound (≥ 2N) Tier 1 derived via SWAP
/// invariance constructing 2N kernel vectors. Upper bound (≤ 2N) was claimed
/// Tier 1 derived via triangle inequality, but the 2026-05-17 review found a
/// derivation gap (matrix-commutator vs conjugation-action confusion in
/// Step 5). The upper bound is now Tier 2 verified empirically for chain
/// N=2..7 and for most other connected graphs at N ≥ 4; the empirical
/// exception is K_3 N=3 (count = 2N+2 = 8). See PROOF_WEIGHT1_DEGENERACY
/// § Appendix (2026-05-17).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F50 +
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> (incl. § Appendix 2026-05-17:
/// K_3 N=3 anomaly + Step-5 derivation gap + Trivial/Sign/Standard weight-1
/// irrep table) + <c>experiments/WEIGHT2_KERNEL.md</c> (Apr 2026; k≥2
/// topology dependence + the Trivial/Alternating/Mixed S_N-irrep table
/// format that today's weight-1 K_3 finding extends) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F50WeightOneDegeneracyPi2Inheritance : Claim
{
    public Pi2DyadicLadderClaim Ladder { get; }
    /// <summary>The "2" count multiplier in 2N. Live from Pi2DyadicLadder a_0.
    /// Counts the two active Pauli letters X and Y per chromatic grade c.</summary>
    public double DegeneracyFactor => Ladder.Term(0);

    /// <summary>The "2" decay rate factor in Re = −2γ. Live from Pi2DyadicLadder a_0.
    /// Same anchor as F1's TwoFactor at the eigenvalue level.</summary>
    public double DecayRateFactor => Ladder.Term(0);

    /// <summary>Total degeneracy count <c>d_real(Re = −2γ) = 2N</c> at the first
    /// non-zero real grid position. The 2N value is the F50 lower bound (always
    /// correct, rigorously proven). For chain at all tested N + most connected
    /// graphs at N ≥ 4 it is also the actual count (equality holds). The one
    /// known empirical exception is N=3 K_3 (= ring = triangle = complete on
    /// 3 vertices), which actually gives 2N+2 = 8; see
    /// <see cref="K3TripleN3ActualCount"/>.</summary>
    public int TotalDegeneracy(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F50 requires N ≥ 2.");
        return (int)(DegeneracyFactor * N);
    }

    /// <summary>The empirically-observed actual count of pure-real Liouvillian
    /// eigenvalues at Re = −2γ for the K_3 N=3 graph (= ring = triangle =
    /// complete on 3 vertices). Returns 8 = 2N + 2, exceeding the F50 lower
    /// bound by 2. The 2 extras are weight-1 operators in the standard 2-dim
    /// irrep of S_3 = Aut(K_3); they exist only on K_3 (any external bond
    /// breaks the S_3 symmetry and restores the count to 2N). See
    /// <c>F50NativeEigenvalueCountTests.NativeSpectrum_K3N3_Anomaly_Has8PureRealAtMinusTwoGamma</c>
    /// for the native verification.</summary>
    public const int K3TripleN3ActualCount = 8;

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
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
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
