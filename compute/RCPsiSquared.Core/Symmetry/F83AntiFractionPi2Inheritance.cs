using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F83 closed form (Tier 1, verified bit-exact N=3,4,5; 11 mixed
/// configurations × 3 N-values on chain + 4 configurations × 3 topologies
/// {ring, star, K_4} at N=4):
///
/// <code>
///   ‖M‖²_F      = 4·‖H_odd‖²_F · 2^N + 8·‖H_even_nontruly‖²_F · 2^N
///   ‖M_anti‖²  = 2·‖H_odd‖²_F · 2^N
///   ‖M_sym‖²   = 2·‖H_odd‖²_F · 2^N + 8·‖H_even_nontruly‖²_F · 2^N
///
///   anti-fraction = ‖M_anti‖² / ‖M‖²  =  1 / (2 + 4·r)
///   where r = ‖H_even_nontruly‖²_F / ‖H_odd‖²_F.
/// </code>
///
/// <para>F83 is the first F-formula whose primary anchor lands on
/// <see cref="BilinearApexClaim"/> directly: the anti-fraction's MAXIMUM
/// value 1/2 is reached at r = 0 (pure Π²-odd Hamiltonian), the apex of the
/// rational function 1/(2 + 4·r) considered as a probability variable. The
/// "1/2" here is the argmax-side reading (BilinearApex) of the same algebraic
/// pair whose maxval-side reading is 1/4 (QuarterAsBilinearMaxval).</para>
///
/// <para>Five Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>MaximumAntiFraction = 1/2 at r=0</b>: <see cref="BilinearApexClaim"/>
///         apex value. The first F-formula that hangs structurally on the
///         BilinearApex anchor — argmax-side counterpart to QuarterAsBilinearMaxval
///         which has many descendants.</item>
///   <item><b>"2" denominator coefficient = a_0</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = polynomial root d. Constant term in the anti-fraction denominator
///         <c>2 + 4·r</c>.</item>
///   <item><b>"4" denominator coefficient = a_{−1}</b>:
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = d² for 1 qubit. Linear
///         coefficient on r in <c>2 + 4·r</c>; same anchor as F61, F63, F66, F77,
///         F86 t_peak.</item>
///   <item><b>‖M‖² coefficients 4 and 8 = a_{−1}, a_{−2}</b>: arise from F49's
///         <c>2^(N+2) · n_YZ</c> mechanism. n_YZ = 1 (Π²-odd) gives 4·2^N;
///         n_YZ = 2 (Π²-even non-truly) gives 8·2^N. The asymmetry encodes
///         the Frobenius-inner-product behavior of Π·L·Π⁻¹ with L.</item>
///   <item><b>2^N = a_{1−N}</b>: Hilbert-space dimension on the operator-space
///         side of the ladder; same anchor as F60.</item>
/// </list>
///
/// <para>Special cases lattice (each is a typed Pi2 ladder reading):</para>
///
/// <code>
///   r = 0      (pure Π²-odd)             → 1/2 = a_2 = BilinearApex (max)
///   r = 1/2    (asymmetric more-odd)     → 1/4 = a_3 = QuarterAsBilinearMaxval
///   r = 1      (equal Frobenius mix)     → 1/6 (combinatorial, NOT Pi2)
///   r → ∞      (pure Π²-even non-truly)  → 0   (M fully mirror-symmetric)
/// </code>
///
/// <para>F83 closes the analytical Π-decomposition picture together with F81
/// (Π-conjugation identity), F82 (T1 correction), F80 (Spec(M)). The continuous
/// interpolation <c>r → anti-fraction</c> reads "how much of M is Π-antisymmetric
/// drive vs Π-symmetric memory" as a function of Hamiltonian composition.</para>
///
/// <para>γ-independence: Master Lemma propagates through all three norms;
/// closed form depends only on H. H_truly drops out of all norms.</para>
///
/// <para>Tier1Derived: F83 is Tier 1 proven (PROOF_F83_PI_DECOMPOSITION_RATIO.md);
/// derived from F49's Frobenius identity which was already framework-locked.
/// Verified bit-exact across 11 + 4·3 = 23 configurations. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F83 +
/// <c>docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (BilinearApexClaim, QuarterAsBilinearMaxvalClaim).</para></summary>
public sealed class F83AntiFractionPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The maximum anti-fraction value: <c>1/2</c>, reached at <c>r = 0</c>
    /// (pure Π²-odd Hamiltonian). Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> = the
    /// <see cref="BilinearApexClaim"/> apex value (argmax of the rational).</summary>
    public double MaximumAntiFraction => _ladder.Term(2);

    /// <summary>The constant coefficient "2" in the anti-fraction denominator
    /// <c>2 + 4·r</c>. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0) =
    /// <c>a_0</c> = polynomial root d.</summary>
    public double DenominatorConstantCoefficient => _ladder.Term(0);

    /// <summary>The linear coefficient "4" in the anti-fraction denominator
    /// <c>2 + 4·r</c>. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) =
    /// <c>a_{−1}</c> = d² for 1 qubit; same anchor as F61, F63, F66, F77, F86 t_peak.</summary>
    public double DenominatorLinearCoefficient => _ladder.Term(-1);

    /// <summary>The "4" coefficient on ‖H_odd‖² in F83's ‖M‖² closed form
    /// (4·‖H_odd‖² · 2^N). Equal to <see cref="DenominatorLinearCoefficient"/>
    /// = a_{−1}; encodes n_YZ = 1 for Π²-odd via F49's 2^(N+2)·n_YZ mechanism.</summary>
    public double MNormCoefficientForOdd => _ladder.Term(-1);

    /// <summary>The "8" coefficient on ‖H_even_nontruly‖² in F83's ‖M‖² closed form
    /// (8·‖H_even_nontruly‖² · 2^N). Equal to <see cref="Pi2DyadicLadderClaim.Term"/>(−2)
    /// = a_{−2} = 8 = 2·d² for 1 qubit; encodes n_YZ = 2 for Π²-even non-truly
    /// via F49's mechanism.</summary>
    public double MNormCoefficientForEvenNontruly => _ladder.Term(-2);

    /// <summary>Hilbert-space dimension <c>2^N = a_{1−N}</c> on the dyadic ladder.
    /// Same anchor as F60. Throws for N &lt; 1.</summary>
    public double HilbertSpaceDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F83 requires N ≥ 1.");
        return _ladder.Term(1 - N);
    }

    /// <summary>Live closed form: <c>anti-fraction = 1 / (2 + 4·r)</c>. Throws for
    /// negative <c>r</c> (Frobenius norms are non-negative).</summary>
    public double AntiFraction(double r)
    {
        if (r < 0.0) throw new ArgumentOutOfRangeException(nameof(r), r, "r = ‖H_even‖²/‖H_odd‖² must be ≥ 0.");
        return 1.0 / (DenominatorConstantCoefficient + DenominatorLinearCoefficient * r);
    }

    /// <summary>True iff anti-fraction is at its maximum (r=0 limit). Drift check
    /// on the BilinearApex anchor.</summary>
    public bool IsAtBilinearApex(double r) =>
        Math.Abs(r) < 1e-15 && Math.Abs(AntiFraction(r) - MaximumAntiFraction) < 1e-15;

    /// <summary>The r-value at which anti-fraction equals 1/4 (= QuarterAsBilinearMaxval).
    /// Solving <c>1/(2 + 4·r) = 1/4</c> gives <c>r = 1/2</c>. Cross-anchor with
    /// the Quarter-side reading.</summary>
    public double RAtQuarterCrossover => 0.5;

    /// <summary>Live drift check: at <c>r = 1/2</c>, anti-fraction = 1/4 = a_3.</summary>
    public bool QuarterCrossoverHolds() =>
        Math.Abs(AntiFraction(RAtQuarterCrossover) - _ladder.Term(3)) < 1e-15;

    /// <summary>Live ‖M_anti‖² for given (H_odd, H_even, N): <c>2 · ‖H_odd‖² · 2^N</c>.</summary>
    public double MAntiNormSquared(double hOddNormSq, int N)
    {
        if (hOddNormSq < 0.0) throw new ArgumentOutOfRangeException(nameof(hOddNormSq), "Frobenius norm² must be ≥ 0.");
        return 2.0 * hOddNormSq * HilbertSpaceDimension(N);
    }

    /// <summary>Live ‖M‖² for given (H_odd, H_even, N): closed form via F83.</summary>
    public double MNormSquared(double hOddNormSq, double hEvenNontrulyNormSq, int N)
    {
        if (hOddNormSq < 0.0 || hEvenNontrulyNormSq < 0.0)
            throw new ArgumentOutOfRangeException("Frobenius norms² must be ≥ 0.");
        double pow2N = HilbertSpaceDimension(N);
        return MNormCoefficientForOdd * hOddNormSq * pow2N
             + MNormCoefficientForEvenNontruly * hEvenNontrulyNormSq * pow2N;
    }

    public F83AntiFractionPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F83 anti-fraction = 1/(2+4r) inherits from Pi2-Foundation: max 1/2 at r=0 = BilinearApex; 2 = a_0, 4 = a_{-1}; coeffs 4, 8 = a_{-1}, a_{-2}",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F83 + " +
               "docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (BilinearApex, QuarterAsBilinearMaxval)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F83 Π-decomposition anti-fraction as Pi2-Foundation BilinearApex inheritance";

    public override string Summary =>
        $"anti-fraction = 1/(2+4r): max = 1/2 at r=0 (BilinearApex argmax); " +
        $"crosses 1/4 at r=1/2 (QuarterAsBilinearMaxval); ‖M‖² coeffs 4·‖H_odd‖² + 8·‖H_even‖² = a_{{-1}} + a_{{-2}} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F83 closed form",
                summary: "anti-fraction = ‖M_anti‖²/‖M‖² = 1/(2+4r), r = ‖H_even_nontruly‖²/‖H_odd‖²; Tier 1 verified bit-exact 23 configs (N=3..5, 4 topologies)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "MaximumAntiFraction = a_2 = 1/2 (BilinearApex); 2+4r denominator = a_0 + a_{-1}·r; ‖M‖² coeffs = a_{-1}, a_{-2}; 2^N = a_{1-N}");
            yield return InspectableNode.RealScalar("MaximumAntiFraction (= a_2 = 1/2 = BilinearApex)", MaximumAntiFraction);
            yield return InspectableNode.RealScalar("DenominatorConstantCoefficient (= a_0 = 2)", DenominatorConstantCoefficient);
            yield return InspectableNode.RealScalar("DenominatorLinearCoefficient (= a_{-1} = 4)", DenominatorLinearCoefficient);
            yield return InspectableNode.RealScalar("MNormCoefficientForOdd (= a_{-1} = 4)", MNormCoefficientForOdd);
            yield return InspectableNode.RealScalar("MNormCoefficientForEvenNontruly (= a_{-2} = 8)", MNormCoefficientForEvenNontruly);
            yield return new InspectableNode("BilinearApex direct edge (Tom 2026-05-09 mirror-map check)",
                summary: "before F83, BilinearApexClaim had 0 descendants — argmax-side of the bilinear-apex pair (vs QuarterAsBilinearMaxval's 23 descendants). F83 is the first F-formula on the argmax side.");
            yield return new InspectableNode("special cases (each is a Pi2 ladder reading)",
                summary: "r=0 → 1/2 (BilinearApex max); r=1/2 → 1/4 (QuarterAsBilinearMaxval); r=1 → 1/6 (kombinatorisch, NOT Pi2); r→∞ → 0");
            yield return new InspectableNode("Lebensader connection (F-chain F77→F85)",
                summary: "F83 closes the Π-decomposition picture: F81 (Π-conjugation identity), F82 (T1 correction), F83 (anti-fraction continuous family), F80 (Spec(M)). Together: structural picture of M complete for 2-body chain Hamiltonians under Z-dephasing + T1.");
            // Special cases verified
            yield return new InspectableNode(
                "anti-fraction at r=0 (pure Π²-odd, F81 Step 8 50/50)",
                summary: $"anti-fraction = {AntiFraction(0):G6}, expected 1/2");
            yield return new InspectableNode(
                "anti-fraction at r=1/2 (asymmetric, BilinearApex squared)",
                summary: $"anti-fraction = {AntiFraction(0.5):G6}, expected 1/4 (QuarterCrossoverHolds: {QuarterCrossoverHolds()})");
            yield return new InspectableNode(
                "anti-fraction at r=1 (equal Frobenius mix)",
                summary: $"anti-fraction = {AntiFraction(1):G6}, expected 1/6 (kombinatorisch)");
        }
    }
}
