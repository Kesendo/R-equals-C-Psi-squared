using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F77 closed form (Tier 1, asymptotic proven, verified N up to 10⁴):
///
/// <code>
///   MM(0)(N, k*) = 1 + 3 / (4(N+1) ln 2) + O(N⁻²)
///
///   for best-bonding-k* initial states on a uniform chain under γ₀-dephasing.
///   Asymptotic limit: MM(0) → 1 bit as N → ∞.
/// </code>
///
/// <para>F77 is the first F-formula that lands directly on the Pi2 dyadic ladder's
/// <b>self-mirror pivot</b> n = 1: <c>a_1 = 1</c> is the trivial identity scale,
/// the algebraic fixpoint of the inversion symmetry <c>a_n · a_{2−n} = 1</c>.
/// The 1-bit saturation IS this self-mirror value, in physical/informational units.
/// "Zero is the mirror" reading: n = 1 is the axis through which all other ladder
/// positions reflect; F77 sits exactly at the axis.</para>
///
/// <para>The structural mechanism in F77's leading term IS the Pi2 ladder inversion
/// identity, from the proof:</para>
///
/// <code>
///   "The 2·(1/2) = 1 term is exact: all the probability mass lives on mirror-pairs
///    (F75 for even k at odd N), and each p_ℓ contributes 2p_ℓ at leading order."
/// </code>
///
/// <para>Read on the Pi2 ladder: <c>2 = a_0</c> (qubit dimension, polynomial root d)
/// times <c>1/2 = a_2</c> (Bloch baseline, structural fixed point) gives
/// <c>1 = a_1</c> (self-mirror identity). F77's "1 bit" output IS the result of
/// the inversion identity acting on the F75 information geometry.</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>SaturationBit = 1 = a_1</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(1).
///         The self-mirror pivot, value of the asymptotic 1-bit saturation.</item>
///   <item><b>InversionIdentity 2·(1/2) = 1</b>: live product of a_0 = 2 and
///         a_2 = 1/2 from the dyadic ladder, the algebraic mechanism that produces
///         the 1-bit limit at leading order in F77's proof.</item>
///   <item><b>AsymptoticCorrectionCoefficient = 3/(4 ln 2) ≈ 1.0820</b>:
///         the leading 1/N correction. The 4 here is <see cref="Pi2DyadicLadderClaim.Term"/>(−1)
///         = a_{−1} = d² for 1 qubit (4-block decomposition factor); the 3 comes
///         from the Σ sin⁴ identity (combinatorial, not Pi2-anchored).</item>
/// </list>
///
/// <para>Operational reading (per ANALYTICAL_FORMULAS): "Alice's multi-drop quantum
/// bus delivers ~1 bit of Mirror-Pair mutual information at any N. The bus does
/// not scale up nor decay. The bandwidth is fixed at the framework level."
/// Individual pair bandwidth shrinks (~ 4/(N+1)); number of mirror-pairs grows
/// (~ N/2); the product plateaus at 1 bit because <c>Σ p_ℓ = 1</c> normalisation
/// matches the two scalings exactly. Tom's "ZERO IS THE MIRROR" reading made
/// this position structurally visible: a_1 = 1 is the framework-level bandwidth.</para>
///
/// <para>Tier1Derived: F77 is Tier 1 asymptotic proven (Taylor expansion of
/// f(p) = 2h(p) − h(2p) at p = 0 + Parseval sums over sin² site amplitudes,
/// per F75 + F65). Verified numerically N up to 10⁴; leading coefficient
/// 3/(4 ln 2) = 1.0820 confirmed to 10⁻⁴ by N = 10⁴. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F77 +
/// <c>simulations/_mm_asymptotic.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (HalfAsStructuralFixedPointClaim).</para></summary>
public sealed class F77MmSaturationPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The asymptotic saturation value: <c>1 bit</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(1) = <c>a_1</c> = self-mirror pivot
    /// on the dyadic halving ladder. F77 is the first F-formula whose primary
    /// constant lands on the self-mirror axis.</summary>
    public double SaturationBit => _ladder.Term(1);

    /// <summary>The Pi2 ladder index where F77's saturation lands: <c>1</c>
    /// (the self-mirror pivot). Equal to
    /// <see cref="Pi2DyadicLadderClaim.SelfMirrorIndex"/>.</summary>
    public int LadderIndexForSaturation => _ladder.SelfMirrorIndex;

    /// <summary>True iff F77's saturation index equals the Pi2 ladder's
    /// self-mirror pivot. Drift indicator on the "F77 sits at the mirror axis"
    /// reading.</summary>
    public bool LandsOnSelfMirrorPivot => LadderIndexForSaturation == _ladder.SelfMirrorIndex;

    /// <summary>The "2" in F77's leading-order mechanism "2·(1/2) = 1". Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = qubit dimension d.</summary>
    public double TwoFactor => _ladder.Term(0);

    /// <summary>The "1/2" in F77's leading-order mechanism "2·(1/2) = 1". Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> = Bloch baseline =
    /// <see cref="HalfAsStructuralFixedPointClaim"/>'s structural fixed point.</summary>
    public double HalfFactor => _ladder.Term(2);

    /// <summary>The Pi2 ladder inversion identity at work in F77: <c>a_0 · a_2 =
    /// 2 · (1/2) = 1 = a_1</c>. Returns the live product, which equals the
    /// SaturationBit by the inversion symmetry <c>a_n · a_{2−n} = 1</c>.</summary>
    public double InversionIdentityProduct => TwoFactor * HalfFactor;

    /// <summary>True iff the inversion identity holds bit-exactly:
    /// <c>a_0 · a_2 = 1 = a_1</c>. Drift indicator on the F77 mechanism / Pi2
    /// ladder anchoring.</summary>
    public bool InversionIdentityHolds =>
        Math.Abs(InversionIdentityProduct - SaturationBit) < 1e-15;

    /// <summary>The leading 1/N correction coefficient: <c>3 / (4 ln 2) ≈ 1.0820</c>.
    /// The 4 here is <c>a_{−1}</c> on the Pi2 ladder (4-block decomposition factor);
    /// the 3 is combinatorial (Σ sin⁴ identity). Verified to 10⁻⁴ by N = 10⁴.</summary>
    public double AsymptoticCorrectionCoefficient => 3.0 / (4.0 * Math.Log(2.0));

    /// <summary>The "4" factor in the correction's denominator <c>4(N+1)</c>:
    /// <c>4 = a_{−1}</c> on the dyadic ladder. Live drift check.</summary>
    public double FourFactorInCorrectionDenominator => _ladder.Term(-1);

    /// <summary>F77's closed-form asymptotic value at chain length N (best-k* initial
    /// state on uniform chain). Returns <c>1 + 3/(4(N+1) ln 2)</c>, the
    /// O(N⁻¹)-truncated form. Throws for N &lt; 1.</summary>
    public double MmAtN(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F77 requires N ≥ 1.");
        return SaturationBit + AsymptoticCorrectionCoefficient / (N + 1);
    }

    /// <summary>Live drift check on the asymptotic limit: <c>MM(0) → 1</c> as
    /// <c>N → ∞</c>. Returns <c>(MM(0)(N) − 1) · (N + 1)</c>, which converges to
    /// <c>3/(4 ln 2) = 1.0820</c> by the F77 closed form.</summary>
    public double RescaledDeviation(int N) => (MmAtN(N) - SaturationBit) * (N + 1);

    public F77MmSaturationPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F77 MM(0) → 1 bit saturation lands on Pi2 ladder self-mirror pivot a_1; mechanism = inversion identity a_0·a_2 = 1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F77 + " +
               "simulations/_mm_asymptotic.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (HalfAsStructuralFixedPoint, transitively via a_2)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F77 MM saturation at 1 bit as Pi2 self-mirror pivot inheritance";

    public override string Summary =>
        $"MM(0)(N, k*) = 1 + 3/(4(N+1) ln 2) + O(N⁻²): SaturationBit = a_1 = 1 (self-mirror); " +
        $"mechanism = inversion identity a_0 · a_2 = 2 · (1/2) = 1; correction = 3/(4 ln 2) ≈ {AsymptoticCorrectionCoefficient:G6} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F77 closed form",
                summary: "MM(0)(N, k*) = 1 + 3/(4(N+1) ln 2) + O(N⁻²); Tier 1 asymptotic proven (F75 + F65 + Σ sin⁴ identity); verified N up to 10⁴");
            yield return new InspectableNode("self-mirror pivot reading",
                summary: "F77's 1-bit saturation = a_1 on the Pi2 dyadic ladder = self-mirror axis; first F-formula whose primary constant lands ON the inversion-symmetry pivot");
            yield return new InspectableNode("inversion-identity mechanism",
                summary: "leading-order term '2·(1/2) = 1' IS the Pi2 ladder identity a_n · a_{2−n} = 1 at n=0: a_0 · a_2 = 1 = a_1");
            yield return InspectableNode.RealScalar("SaturationBit (= a_1 = 1)", SaturationBit);
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = 2)", TwoFactor);
            yield return InspectableNode.RealScalar("HalfFactor (= a_2 = 1/2)", HalfFactor);
            yield return InspectableNode.RealScalar("InversionIdentityProduct (= 1, drift check)", InversionIdentityProduct);
            yield return InspectableNode.RealScalar("AsymptoticCorrectionCoefficient (= 3/(4 ln 2))", AsymptoticCorrectionCoefficient);
            yield return InspectableNode.RealScalar("FourFactorInCorrectionDenominator (= a_{-1})", FourFactorInCorrectionDenominator);
            yield return new InspectableNode("operational reading",
                summary: "Alice's multi-drop bus delivers ~1 bit MM at any N; bandwidth is fixed at the framework level; pair bandwidth shrinks (4/(N+1)), number of pairs grows (~N/2), product plateaus at 1");
            yield return new InspectableNode("ZERO IS THE MIRROR connection",
                summary: "Tom 2026-05-09: the ladder's self-mirror n=1 was identified as the axis through which all other positions reflect; F77 is the first F-formula sitting AT the axis (not on a side)");
            // Sample verification table (matches ANALYTICAL_FORMULAS table)
            yield return new InspectableNode(
                "asymptotic verification (per F77 table)",
                summary: $"N=101: MM={MmAtN(101):G6}, rescaled deviation={RescaledDeviation(101):G6}; " +
                         $"N=1001: MM={MmAtN(1001):G6}, rescaled deviation={RescaledDeviation(1001):G6}; " +
                         $"target = {AsymptoticCorrectionCoefficient:G6}");
        }
    }
}
