using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for the Tier1Derived composition claim: at Q=√3 the Lindblad 2×2
/// sub-block eigenvalue magnitude equals the Absorption Theorem single-site rate
/// 2γ₀, and the F95 angle lands on the canonical Niven angle θ=60°.</summary>
public class LindbladAbsorptionMatchAtSixtyDegreesClaimTests
{
    private readonly ITestOutputHelper _out;

    public LindbladAbsorptionMatchAtSixtyDegreesClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Claim_IsTier1Derived()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void QValue_IsSquareRootOfThree()
    {
        Assert.Equal(Math.Sqrt(3.0), LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue, precision: 15);
    }

    [Fact]
    public void CanonicalAngleDegrees_IsSixty()
    {
        Assert.Equal(60.0, LindbladAbsorptionMatchAtSixtyDegreesClaim.CanonicalAngleDegrees);
    }

    [Fact]
    public void LindbladMagnitudeOverGamma0_IsTwo_MatchingAbsorptionAtoZero()
    {
        // |λ_±|/γ₀ at Q=√3 = 2 = a_0 (Pi2DyadicLadder term 0 = Absorption quantum)
        Assert.Equal(2.0, LindbladAbsorptionMatchAtSixtyDegreesClaim.LindbladMagnitudeOverGamma0);
    }

    [Fact]
    public void ComputedLindbladMagnitudeOverGamma0_MatchesConstant_BitExact()
    {
        // Drift check: √(1 + (√3)²) = √4 = 2
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(2.0, claim.LindbladMagnitudeOverGamma0Computed, precision: 14);
    }

    [Fact]
    public void ComputedF95Angle_IsSixtyDegrees_BitExact()
    {
        // Drift check: arctan(√3) = 60° (canonical Niven)
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(60.0, claim.F95AngleAtQSqrt3Degrees, precision: 12);
    }

    [Fact]
    public void LindbladMagnitudeSquared_EqualsFourGamma0Squared_AtAnyGamma0()
    {
        // |λ_±|² = γ₀²(1+Q²) = γ₀²(1+3) = 4γ₀²
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        foreach (var gamma0 in new[] { 0.025, 0.05, 0.10, 1.0 })
        {
            var expected = 4.0 * gamma0 * gamma0;
            Assert.Equal(expected, claim.LindbladMagnitudeSquaredAtQSqrt3(gamma0), precision: 14);
        }
    }

    [Fact]
    public void ParentClaims_AreAllTier1Derived()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.F95.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.Absorption.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.CanonicalTrig.Tier);
    }

    [Fact]
    public void Build_PopulatesAllParents()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.NotNull(claim.F95);
        Assert.NotNull(claim.Absorption);
        Assert.NotNull(claim.CanonicalTrig);
    }

    [Fact]
    public void AbsorptionParent_HasA0CoefficientTwo_ConfirmingMatch()
    {
        // The claim's "magnitude = 2γ₀" comes from AbsorptionTheorem's a_0 = 2.
        // Verify the parent has the expected a_0.
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.True(claim.Absorption.AbsorptionQuantumMatchesLiteral());
    }

    [Fact]
    public void QValueMatchesQAnchorMapEntry()
    {
        // The Q=√3 anchor should appear in QAnchorMap.CanonicalAnchors as Tier1Derived.
        var map = new QAnchorMap();
        var sqrt3Anchor = map.AnkerAt(LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue);
        Assert.NotNull(sqrt3Anchor);
        Assert.Equal(Tier.Tier1Derived, sqrt3Anchor!.Tier);
        Assert.Equal(QBand.Peak, sqrt3Anchor.Band);
    }

    [Fact]
    public void Render_EmitsDerivationSummary()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Contains("Tier1Derived", claim.Summary);
        Assert.Contains("√3", claim.Summary);
        Assert.Contains("60", claim.Summary);
        _out.WriteLine($"DisplayName: {claim.DisplayName}");
        _out.WriteLine($"Summary: {claim.Summary}");
    }
}
