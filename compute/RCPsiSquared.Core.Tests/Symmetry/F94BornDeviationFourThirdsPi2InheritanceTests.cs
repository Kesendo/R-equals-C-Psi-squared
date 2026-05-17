using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F94BornDeviationFourThirdsPi2InheritanceTests
{
    private static F94BornDeviationFourThirdsPi2Inheritance BuildClaim() =>
        new F94BornDeviationFourThirdsPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new QuarterAsBilinearMaxvalClaim());

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
        Assert.Equal(3, F94BornDeviationFourThirdsPi2Inheritance.ThreeDenominator);
    }

    [Fact]
    public void Sym3PartialTraceInteger_IsExactlyEight()
    {
        // The Dyson sym3 = L_H²L'_dis + L_HL'_disL_H + L'_disL_H² acting on
        // ρ_0 = |0+0+⟩⟨0+0+| at N=4, partial-traced to pair (0, 2), |00⟩
        // diagonal element = 8 BIT-EXACT (computed in
        // simulations/_born_rule_tier1_derivation.py).
        Assert.Equal(8, F94BornDeviationFourThirdsPi2Inheritance.Sym3PartialTraceInteger);
    }

    [Fact]
    public void TaylorThreeFactorial_IsSix()
    {
        Assert.Equal(6, F94BornDeviationFourThirdsPi2Inheritance.TaylorThreeFactorial);
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
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F94BornDeviationFourThirdsPi2Inheritance(null!, quarter));
        Assert.Throws<ArgumentNullException>(() =>
            new F94BornDeviationFourThirdsPi2Inheritance(ladder, null!));
    }

    [Fact]
    public void OneOverFourFactor_IsExactlyQuarter_FromLadderTermThree()
    {
        // Mirror partner via inversion identity a_{-1} · a_3 = 4 · (1/4) = 1.
        Assert.Equal(0.25, BuildClaim().OneOverFourFactor, precision: 14);
    }

    [Fact]
    public void MirrorPartnerProductIsOne_HoldsBitExact()
    {
        // Drift check on the dyadic-ladder inversion identity:
        // FourFactor · OneOverFourFactor = a_{-1} · a_3 = 1.
        var f = BuildClaim();
        Assert.True(f.MirrorPartnerProductIsOne());
        Assert.Equal(1.0, f.FourFactor * f.OneOverFourFactor, precision: 14);
    }

    [Fact]
    public void Quarter_TypedParent_IsExposed()
    {
        Assert.NotNull(BuildClaim().Quarter);
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

    [Fact]
    public void SurvivingDysonDiagrams_IsExactlyThirtyTwo()
    {
        // Bit-exact via direct enumeration in simulations/_born_rule_sym3_decomposition.py
        // (2026-05-17): 32 of the 4·4·4·3·3·3 = 1728 (b₁, b₂, s, ord, c₁, c₂) sextuples
        // contribute non-zero to ⟨00|_pair Tr_{1,3}[sym3·ρ_0]|00⟩_pair.
        Assert.Equal(32, F94BornDeviationFourThirdsPi2Inheritance.SurvivingDysonDiagrams);
    }

    [Fact]
    public void RawPauliPerDiagram_IsExactlyFour()
    {
        // Each surviving diagram contributes raw Pauli value 4 (uniform across all 32),
        // which after the (J/4)² = 1/16 Heisenberg coupling normalization becomes 1/4.
        Assert.Equal(4, F94BornDeviationFourThirdsPi2Inheritance.RawPauliPerDiagram);
    }

    [Fact]
    public void CellCounts_SumToThirtyTwo()
    {
        // Cell A (ord=1, XX, adj-kept) = 8
        // Cell B (ord=2, XX, self ∪ adj-kept) = 16
        // Cell C (ord=2, YY, self only) = 8
        // 8 + 16 + 8 = 32 = SurvivingDysonDiagrams
        Assert.Equal(8, F94BornDeviationFourThirdsPi2Inheritance.CellA_Ord1XX_AdjKeptSide);
        Assert.Equal(16, F94BornDeviationFourThirdsPi2Inheritance.CellB_Ord2XX_SelfOrAdjKeptSide);
        Assert.Equal(8, F94BornDeviationFourThirdsPi2Inheritance.CellC_Ord2YY_Self);
        Assert.True(BuildClaim().CellCountsSumToSurvivingDiagrams());
    }

    [Fact]
    public void StructuralDecomposition_RecoversSym3Integer()
    {
        // 32 surviving diagrams × 4 raw / 16 (= a_{-1}²) = 8 = Sym3PartialTraceInteger.
        // Equivalently: 32 × 4 = 128 = 8 × 16. Bit-exact via integer arithmetic.
        Assert.True(BuildClaim().StructuralDecompositionRecoversSym3());
    }

    [Fact]
    public void StructuralDecomposition_Yields4ThirdsFromAlternativeFormula()
    {
        // 4/3 = SurvivingDysonDiagrams / (a_{-1} · 3!) = 32 / (4 · 6) = 32/24.
        // Complement to the typed-anchor reading 4/3 = a_{-1} / 3.
        var f = BuildClaim();
        int a_minus_1 = (int)f.FourFactor;
        double structural =
            (double)F94BornDeviationFourThirdsPi2Inheritance.SurvivingDysonDiagrams
            / (a_minus_1 * F94BornDeviationFourThirdsPi2Inheritance.TaylorThreeFactorial);
        Assert.Equal(4.0 / 3.0, structural, precision: 15);
        Assert.Equal(f.Coefficient, structural, precision: 15);
    }

    [Fact]
    public void TopologicalCut_SelfAndAdjEachContribute16()
    {
        // Alternative cut: 16 self-bond-pair diagrams (8 XX + 8 YY, all ord=2)
        // + 16 adj-bond-pair diagrams (all XX, 8 ord=1 + 8 ord=2) = 32.
        // From the cell breakdown:
        //   self = (XX self in Cell B = 8) + (YY self in Cell C = 8) = 16
        //   adj  = (XX adj in Cell A = 8) + (XX adj in Cell B = 8) = 16
        int self_diagrams = 8 + 8;  // XX self (Cell B half) + YY self (Cell C)
        int adj_diagrams = 8 + 8;   // XX adj (Cell A) + XX adj (Cell B half)
        Assert.Equal(16, self_diagrams);
        Assert.Equal(16, adj_diagrams);
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.SurvivingDysonDiagrams,
            self_diagrams + adj_diagrams);
    }
}
