using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F94BornDeviationFourThirdsPi2InheritanceTests
{
    private static F94BornDeviationFourThirdsPi2Inheritance BuildClaim() =>
        new F94BornDeviationFourThirdsPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Coefficient_IsExactlyFourThirds()
    {
        // Bit-exact: 4/3 from a_{-1} = 4 on the dyadic ladder divided by
        // ThreeDenominator = 3 (Dyson sym3 = 8 over Taylor 3! = 6 reduces to 4/3).
        Assert.Equal(4.0 / 3.0, BuildClaim().Coefficient, precision: 15);
    }

    [Fact]
    public void FourFactor_IsExactlyFour_FromLadderTermMinusOne()
    {
        // The "4" in 4/3 is a_{-1} = 4 on the Pi2 dyadic ladder, the same "4"
        // that appears in F86 t_peak = 1/(4γ₀) and F77's MM correction denominator.
        Assert.Equal(4.0, BuildClaim().FourFactor, precision: 14);
    }

    [Fact]
    public void ThreeDenominator_IsExactlyThree()
    {
        Assert.Equal(3, BuildClaim().ThreeDenominator);
    }

    [Fact]
    public void Sym3PartialTraceInteger_IsExactlyEight()
    {
        // The Dyson sym3 = L_H²L'_dis + L_HL'_disL_H + L'_disL_H² acting on
        // ρ_0 = |0+0+⟩⟨0+0+| at N=4, partial-traced to pair (0, 2), |00⟩
        // diagonal element = 8 BIT-EXACT (computed in
        // simulations/_born_rule_tier1_derivation.py).
        Assert.Equal(8, BuildClaim().Sym3PartialTraceInteger);
    }

    [Fact]
    public void TaylorThreeFactorial_IsSix()
    {
        Assert.Equal(6, BuildClaim().TaylorThreeFactorial);
    }

    [Fact]
    public void CoefficientAgreesWithSym3()
    {
        // Coefficient (= a_{-1}/3 = 4/3) bit-exact equals Sym3PartialTraceInteger
        // / TaylorThreeFactorial (= 8/6 = 4/3). Drift indicator.
        Assert.True(BuildClaim().CoefficientAgreesWithSym3());
    }

    [Theory]
    // Pure closed-form leading-order Δ_|00⟩(Q, K) = (4/3) · Q² · K³.
    // The expected value is computed inline (no hardcoded literals) to avoid
    // double-precision representation drift across machines. The full numerical
    // Lindblad value at these points includes O(Q³·K⁴) higher-order corrections;
    // the leading-order closed form is what the typed claim returns.
    [InlineData(20.0,  0.0143)]
    [InlineData(10.0,  0.005)]
    [InlineData(100.0, 0.005)]
    [InlineData(40.0,  0.005)]
    public void DeltaDominant_MatchesClosedForm(double Q, double K)
    {
        double expected = (4.0 / 3.0) * Q * Q * K * K * K;
        Assert.Equal(expected, BuildClaim().DeltaDominant(Q, K), precision: 14);
    }

    [Fact]
    public void DeltaDominant_LeadingOrder_AgreesWithLindblad_AtSmallQK()
    {
        // At a representative point well inside the perturbative regime, the
        // closed form should agree with the numerically extracted coefficient
        // c_emp = 1.32992 ± 0.006 (per simulations/_born_rule_delta_dominant_coefficient.py
        // 16-sample mean) to <1%. Drift indicator.
        var f = BuildClaim();
        double Q = 10.0, K = 0.005;
        double closed_form = f.DeltaDominant(Q, K);
        double c_empirical = 1.32992;
        double Q2K3 = Q * Q * K * K * K;
        double expected_from_empirical = c_empirical * Q2K3;
        double relative_diff = Math.Abs(closed_form - expected_from_empirical) / closed_form;
        Assert.True(relative_diff < 0.01,
            $"Closed form {closed_form} should agree with empirical fit {expected_from_empirical} within 1%; relative diff = {relative_diff}");
    }

    [Fact]
    public void DeltaDominant_PhysicalUnits_AgreeWithDimensionless()
    {
        // (4/3) · J² · γ · t³ = (4/3) · (J/γ)² · (γt)³ = (4/3) · Q² · K³
        var f = BuildClaim();
        double J = 1.0, gamma = 0.05, t = 0.286;
        double Q = J / gamma;
        double K = gamma * t;
        Assert.Equal(f.DeltaDominant(Q, K), f.DeltaP_Dominant(J, gamma, t), precision: 12);
    }

    [Fact]
    public void C_DominantOutcome_IsOnePlusDelta()
    {
        // Per BORN_RULE_MIRROR generalization R_i = C_i · Ψ_i²:
        // C_|00⟩ = 1 + Δ_|00⟩
        var f = BuildClaim();
        double Q = 20.0, K = 0.0143;
        Assert.Equal(1.0 + f.DeltaDominant(Q, K), f.C_DominantOutcome(Q, K), precision: 12);
    }

    [Fact]
    public void DeltaDominant_RejectsNegativeQ()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaDominant(-1.0, 0.01));
    }

    [Fact]
    public void DeltaDominant_RejectsNegativeK()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaDominant(20.0, -0.01));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F94BornDeviationFourThirdsPi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_ProofAndScriptsAndReflection()
    {
        var f = BuildClaim();
        Assert.Contains("PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md", f.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md F94", f.Anchor);
        Assert.Contains("_born_rule_tier1_derivation.py", f.Anchor);
        Assert.Contains("_born_rule_delta_dominant_coefficient.py", f.Anchor);
        Assert.Contains("ON_HOW_FOUR_THIRDS_APPEARED.md", f.Anchor);
        Assert.Contains("BORN_RULE_MIRROR.md", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
    }
}
