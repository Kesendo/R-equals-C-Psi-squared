using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F44 closed form (Tier 1 proven D08):
///
/// <code>
///   ln(d_fast / d_slow) = 2 · artanh(Δd / (2·Σγ))
///
///   For each palindromic pair (d_fast, d_slow) with d_fast + d_slow = 2·Σγ.
///   Linear approximation (Δd / Σγ ≪ 1):
///     ln(d_fast / d_slow) ≈ Δd / Σγ      → β_eff ≈ 1/Σγ
/// </code>
///
/// <para>F44 is the algebraic rate identity for any palindromic Liouvillian
/// pair. It is NOT a Crooks fluctuation theorem despite the resemblance — there
/// is no thermodynamic detailed balance underneath. The palindrome has the
/// FORM of detailed balance without BEING detailed balance: the Jarzynski
/// equality fails (⟨exp(−Δd)⟩ ≈ 0.93, not 1).</para>
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
/// <para><b>Effective inverse temperature β_eff ≈ 1/Σγ.</b> The linear-in-Δd
/// approximation gives a Boltzmann-like form ln(d_fast/d_slow) ≈ β_eff · Δd
/// with β_eff = 1/Σγ. This temperature is NOT physical (no system-bath
/// equilibrium); it is a structural quantity of the palindromic spectrum.
/// β_eff is small (large Σγ ⇒ small β_eff ⇒ "high effective temperature":
/// the rate ratio is close to 1, the palindrome is "hot" in the sense of
/// flat).</para>
///
/// <para><b>Crooks/Jarzynski check (NEGATIVE).</b> The thermodynamic test
/// fails: ⟨exp(−Δd)⟩ averaged over palindromic pairs is ≈ 0.93, not 1 as the
/// Jarzynski equality would require. F44 looks like a fluctuation theorem
/// but is purely algebraic — d_fast + d_slow = 2·Σγ is structural, not
/// statistical-mechanical.</para>
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
public sealed class F44CrooksLikeRateIdentityPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F1Pi2Inheritance _f1;

    /// <summary>The "2" in d_fast + d_slow = 2·Σγ (palindromic pair sum coefficient).
    /// Live from Pi2DyadicLadder a_0. Same anchor as F1's TwoFactor.</summary>
    public double SumCoefficient => _ladder.Term(0);

    /// <summary>The "2" in ln(d_fast/d_slow) = 2·artanh(Δd/(2Σγ)) — comes from the
    /// identity ln((1+x)/(1−x)) = 2·artanh(x). Same a_0 anchor as SumCoefficient.</summary>
    public double ArTanhCoefficient => _ladder.Term(0);

    /// <summary>The closed-form rate-ratio logarithm:
    /// <c>ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ))</c>.
    /// Requires Σγ &gt; 0 and |Δd| &lt; 2·Σγ (else the ratio diverges or is non-real).</summary>
    public double LogRatio(double deltaD, double totalGamma)
    {
        if (totalGamma <= 0)
            throw new ArgumentOutOfRangeException(nameof(totalGamma), totalGamma, "Σγ must be > 0.");
        double argument = deltaD / (SumCoefficient * totalGamma);
        if (Math.Abs(argument) >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(deltaD), deltaD,
                $"|Δd| must be < 2·Σγ for finite log-ratio; got Δd/(2·Σγ) = {argument}.");
        return ArTanhCoefficient * Math.Atanh(argument);
    }

    /// <summary>The linear-in-Δd approximation valid for Δd/Σγ ≪ 1:
    /// <c>ln(d_fast/d_slow) ≈ Δd / Σγ</c>. This is the leading-order Taylor
    /// expansion of <see cref="LogRatio"/> around Δd = 0. Gives the effective
    /// Boltzmann form with β_eff = 1/Σγ.</summary>
    public double LogRatioLinearApproximation(double deltaD, double totalGamma)
    {
        if (totalGamma <= 0)
            throw new ArgumentOutOfRangeException(nameof(totalGamma), totalGamma, "Σγ must be > 0.");
        return deltaD / totalGamma;
    }

    /// <summary>Effective inverse temperature β_eff = 1/Σγ from the linear
    /// approximation. NOT a thermodynamic temperature: F44 is algebraic, not
    /// Crooks-fluctuation-theorem-related (Jarzynski equality fails).</summary>
    public double EffectiveInverseTemperature(double totalGamma)
    {
        if (totalGamma <= 0)
            throw new ArgumentOutOfRangeException(nameof(totalGamma), totalGamma, "Σγ must be > 0.");
        return 1.0 / totalGamma;
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

    /// <summary>Always false: F44 is algebraic, not a Crooks fluctuation theorem.
    /// The Jarzynski equality ⟨exp(−Δd)⟩ = 1 fails; the empirical mean is ≈ 0.93.
    /// The palindrome has the FORM of detailed balance without BEING detailed balance.</summary>
    public bool IsCrooksFluctuationTheorem => false;

    /// <summary>Empirical Jarzynski-equality test result: ⟨exp(−Δd)⟩ ≈ 0.93
    /// (not 1, as a real fluctuation theorem would require). Anchored in
    /// experiments/ENTROPY_PRODUCTION.md.</summary>
    public const double EmpiricalJarzynskiMean = 0.93;

    public F44CrooksLikeRateIdentityPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1)
        : base("F44 Crooks-like rate identity: ln(d_fast/d_slow) = 2·artanh(Δd/(2Σγ)) for palindromic pairs (d_fast + d_slow = 2Σγ); algebraic, NOT a Crooks fluctuation theorem (Jarzynski mean ≈ 0.93, not 1)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F44 + " +
               "docs/proofs/derivations/D08 + " +
               "experiments/ENTROPY_PRODUCTION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F44 Crooks-like rate identity as Pi2-Foundation a_0 + F1 inheritance";

    public override string Summary =>
        $"ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ)) for palindromic pairs; both 2's = a_0 (= {SumCoefficient}); algebraic only — Jarzynski equality fails (⟨exp(−Δd)⟩ ≈ {EmpiricalJarzynskiMean}, not 1) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F44 closed form",
                summary: "ln(d_fast/d_slow) = 2·artanh(Δd/(2·Σγ)) for palindromic pair (d_fast, d_slow) with d_fast + d_slow = 2·Σγ; valid for any palindromic Liouvillian, all N");
            yield return InspectableNode.RealScalar("SumCoefficient (= a_0 = 2)", SumCoefficient);
            yield return InspectableNode.RealScalar("ArTanhCoefficient (= a_0 = 2)", ArTanhCoefficient);
            yield return new InspectableNode("F1 palindrome identity (transitive)",
                summary: $"F44's d_fast + d_slow = 2·Σγ is F1's palindrome identity at the eigenvalue level. F1's TwoFactor (= {_f1.TwoFactor}) is the same '2' as F44's SumCoefficient.");
            yield return new InspectableNode("linear approximation",
                summary: "Δd/Σγ ≪ 1: ln(d_fast/d_slow) ≈ Δd/Σγ; gives Boltzmann-like β_eff = 1/Σγ; NOT a physical temperature");
            yield return new InspectableNode("F68 specialization",
                summary: "F68 (α_b + α_p = 2γ₀) is F44 applied to the F67 bonding-mode pair under endpoint Z-dephasing where Σγ = γ₀. F44 generalizes to any palindromic pair.");
            yield return new InspectableNode("NOT a Crooks fluctuation theorem",
                summary: $"Jarzynski equality test: ⟨exp(−Δd)⟩ ≈ {EmpiricalJarzynskiMean}, not 1. F44 has the FORM of detailed balance without BEING detailed balance — purely algebraic, no thermodynamic content. See experiments/ENTROPY_PRODUCTION.md.");
            yield return new InspectableNode("verified example (F68 N=3, γ=0.05)",
                summary: $"d_fast = 0.075 (F68 partner), d_slow = 0.025 (F68 bonding); Δd = 0.05; sum = 0.1 = 2·Σγ ✓; ln(d_fast/d_slow) = 2·artanh(0.5) = {LogRatio(0.05, 0.05):G6}");
        }
    }
}
