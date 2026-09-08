using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F44 algebraic pair-rate log identity (Tier 1 proven D08):
///
/// <code>
///   ln(d_fast / d_slow) = 2 · artanh(Δd / (2·Σγ))
///
///   For each palindromic pair (d_fast, d_slow) with d_fast + d_slow = 2·Σγ.
///   Linear approximation (Δd / Σγ ≪ 1):
///     ln(d_fast / d_slow) ≈ Δd / Σγ
/// </code>
///
/// <para>F44 is an algebraic identity for positive decay-rate pairs. It defines
/// neither process probabilities nor work, energy, free-energy, or temperature
/// variables. No probability/work ensemble or dimensionless thermodynamic
/// exponent is part of this typed claim.</para>
///
/// <para><b>Relationship to F1 + F68.</b> F44 is the general algebraic
/// identity for ANY palindromic pair; F68 is a specific application to the
/// bonding-mode pair (α_b, α_p) under endpoint Z-dephasing where Σγ = γ₀.
/// The pair sum d_fast + d_slow = 2Σγ is F1's palindrome identity at the
/// eigenvalue level (same as F68's "α_b + α_p = 2γ₀"). What F44 adds is the
/// artanh-closed-form ratio reading: the RATIO (not just the sum) of the pair
/// rates has a closed expression in Δd / (2Σγ).</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>SumCoefficient = 2 = a_0</b>: in d_fast + d_slow = 2·Σγ
///         (palindromic pair sum). Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1's
///         TwoFactor and F68's SumCoefficient.</item>
///   <item><b>ArTanhCoefficient = 2 = a_0</b>: in ln(d_fast/d_slow) =
///         2·artanh(...). Same a_0 anchor; the "2" multiplier on artanh comes
///         from the standard identity ln((1+x)/(1−x)) = 2·artanh(x).</item>
///   <item><b>F1 palindrome identity (transitive)</b>: d_fast + d_slow = 2·Σγ
///         is F1's identity at the eigenvalue level, the same as F68's
///         "α_b + α_p = 2γ₀".</item>
/// </list>
///
/// <para><b>Inverse-rate linear coefficient.</b> The derivative of the
/// log-ratio with respect to Δd at Δd = 0 is 1/Σγ. Its unit is inverse rate;
/// the derivation supplies no temperature interpretation.</para>
///
/// <para>Tier1Derived: F44 is Tier 1 proven in D08 (one-line algebraic
/// identity from artanh definition); valid for any palindromic Liouvillian,
/// all N. Pi2-Foundation anchoring is composition through Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F44 (line 879) +
/// <c>docs/proofs/derivations/D08_*.md</c> +
/// <c>experiments/ENTROPY_PRODUCTION.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (palindrome
/// identity at eigenvalue level) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F44AlgebraicPairRateLogIdentityPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No bespoke typed bit_a-twin: F44's rate-ratio identity is a Lindblad eigenvalue-pair
    /// statement, covered by the Hadamard X↔Z duality (see <see cref="BitATwinStatus"/>).</summary>
    public Claim? BitATwin => null;

    /// <summary>Covered by the Hadamard X↔Z duality (Case 2 of docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md):
    /// F44's ln(d_fast/d_slow)=2·artanh(Δd/2Σγ) is a pure algebraic relation between palindromic Liouvillian
    /// eigenvalues, generalizing the already-covered F68; U_H carries
    /// the Z-dephasing eigenvalue pair to the X-dephasing one, so no bespoke twin is owed.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.CoveredByHadamardDuality;
    public Pi2DyadicLadderClaim Ladder { get; }
    public F1Pi2Inheritance F1 { get; }

    /// <summary>The "2" in d_fast + d_slow = 2·Σγ (palindromic pair sum coefficient).
    /// Live from Pi2DyadicLadder a_0. Same anchor as F1's TwoFactor.</summary>
    public double SumCoefficient => Ladder.Term(0);

    /// <summary>The "2" in ln(d_fast/d_slow) = 2·artanh(Δd/(2Σγ)); comes from the
    /// identity ln((1+x)/(1−x)) = 2·artanh(x). Same a_0 anchor as SumCoefficient.</summary>
    public double ArTanhCoefficient => Ladder.Term(0);

    /// <summary>The closed-form rate-ratio logarithm:
    /// <c>ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ))</c>.
    /// Requires finite Σγ &gt; 0 and finite 0 ≤ Δd &lt; 2·Σγ. The ordering
    /// Δd = d_fast − d_slow is part of the public contract.</summary>
    public double LogRatio(double deltaD, double totalGamma)
    {
        ValidateOrderedPairDomain(deltaD, totalGamma);
        double argument = (deltaD / totalGamma) / SumCoefficient;
        return ArTanhCoefficient * Math.Atanh(argument);
    }

    /// <summary>The linear-in-Δd approximation valid for Δd/Σγ ≪ 1:
    /// <c>ln(d_fast/d_slow) ≈ Δd / Σγ</c>. This is the leading-order Taylor
    /// expansion of <see cref="LogRatio"/> around Δd = 0 and retains that
    /// method's finite ordered-pair domain 0 ≤ Δd &lt; 2·Σγ.</summary>
    public double LogRatioLinearApproximation(double deltaD, double totalGamma)
    {
        ValidateOrderedPairDomain(deltaD, totalGamma);
        return deltaD / totalGamma;
    }

    /// <summary>The inverse-rate linear coefficient 1/Σγ. This is the derivative
    /// of the algebraic log-ratio with respect to Δd at the pair centre.</summary>
    public double InverseRateLinearCoefficient(double totalGamma)
    {
        ValidateTotalGamma(totalGamma);
        return 1.0 / totalGamma;
    }

    private void ValidateOrderedPairDomain(double deltaD, double totalGamma)
    {
        ValidateTotalGamma(totalGamma);
        if (!double.IsFinite(deltaD) || deltaD < 0.0 || deltaD / totalGamma >= SumCoefficient)
            throw new ArgumentOutOfRangeException(nameof(deltaD), deltaD,
                "Δd = d_fast − d_slow must be finite and satisfy 0 ≤ Δd < 2·Σγ.");
    }

    private static void ValidateTotalGamma(double totalGamma)
    {
        if (!double.IsFinite(totalGamma) || totalGamma <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(totalGamma), totalGamma, "Σγ must be finite and > 0.");
    }

    /// <summary>Drift check: verifies the palindromic-pair sum d_fast + d_slow = 2·Σγ
    /// holds within tolerance. This is the F1 palindrome identity at the eigenvalue
    /// level; F44 asserts it as the precondition for the artanh closed form.</summary>
    public bool PalindromicSumHolds(double dFast, double dSlow, double totalGamma, double tolerance = 1e-12)
    {
        double sum = dFast + dSlow;
        double expected = SumCoefficient * totalGamma;
        return Math.Abs(sum - expected) < tolerance;
    }

    public F44AlgebraicPairRateLogIdentityPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1)
        : base("F44 algebraic pair-rate log identity: ln(d_fast/d_slow) = 2·artanh(Δd/(2Σγ)) for finite positive palindromic pairs with Δd=d_fast−d_slow and 0≤Δd<2Σγ; no probability/work ensemble or temperature is defined",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F44 + " +
               "docs/proofs/derivations/D08 + " +
               "experiments/ENTROPY_PRODUCTION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F44 algebraic pair-rate log identity as Pi2-Foundation a_0 + F1 inheritance";

    public override string Summary =>
        $"ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ)) for finite positive palindromic pairs, Δd=d_fast−d_slow and 0≤Δd<2Σγ; both 2's = a_0 (= {SumCoefficient}); inverse-rate coefficient 1/Σγ; no probability/work ensemble or temperature is defined ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F44 closed form",
                summary: "ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ)) for a finite ordered palindromic pair (d_fast, d_slow), Δd=d_fast−d_slow, 0≤Δd<2Σγ, and d_fast+d_slow=2Σγ; valid for any palindromic Liouvillian, all N");
            yield return InspectableNode.RealScalar("SumCoefficient (= a_0 = 2)", SumCoefficient);
            yield return InspectableNode.RealScalar("ArTanhCoefficient (= a_0 = 2)", ArTanhCoefficient);
            yield return new InspectableNode("F1 palindrome identity (transitive)",
                summary: $"F44's d_fast + d_slow = 2·Σγ is F1's palindrome identity at the eigenvalue level. F1's TwoFactor (= {F1.TwoFactor}) is the same '2' as F44's SumCoefficient.");
            yield return new InspectableNode("linear approximation",
                summary: "Δd/Σγ ≪ 1: ln(d_fast/d_slow) ≈ Δd/Σγ; inverse-rate linear coefficient = 1/Σγ");
            yield return new InspectableNode("F68 specialization",
                summary: "F68 (α_b + α_p = 2γ₀) is F44 applied to the F67 bonding-mode pair under endpoint Z-dephasing where Σγ = γ₀. F44 generalizes to any palindromic pair.");
            yield return new InspectableNode("scope boundary",
                summary: "Algebraic decay-rate identity only: no process probabilities, work distribution, free-energy change, bath temperature, or dimensionless thermodynamic exponent is defined.");
            yield return new InspectableNode("verified example (F68 N=3, γ=0.05)",
                summary: $"d_fast = 0.075 (F68 partner), d_slow = 0.025 (F68 bonding); Δd = 0.05; sum = 0.1 = 2·Σγ ✓; ln(d_fast/d_slow) = 2·artanh(0.5) = {LogRatio(0.05, 0.05):G6}");
        }
    }
}
