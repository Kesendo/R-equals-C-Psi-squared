using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F99 typed-Claim tests: α = sin²(θ)/2 at the five canonical trig
/// angles {0°, 30°, 45°, 60°, 90°} produces the five Pi2 dyadic anchors
/// {0, 1/8, 1/4, 3/8, 1/2}. Cross-checked against the Python derivation
/// script <c>simulations/carbon/depth_3_anchor_derivation.py</c> which verified
/// bit-exact at N = 4, 6, 8 via non-uniform Dicke superposition.</summary>
public class CanonicalTrigAnchorPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public CanonicalTrigAnchorPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static CanonicalTrigAnchorPi2Inheritance Build() =>
        new CanonicalTrigAnchorPi2Inheritance(
            new HalfAsStructuralFixedPointClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new KIntermediateAsymptoteQuarterInheritance(
                new QuarterAsBilinearMaxvalClaim(),
                new HalfAsStructuralFixedPointClaim(),
                new DickeSuperpositionQuarterPi2Inheritance(
                    new Pi2DyadicLadderClaim(),
                    new QuarterAsBilinearMaxvalClaim(),
                    new HalfAsStructuralFixedPointClaim())));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Theory]
    [InlineData(0,  1.0,         0.0)]              // Mirror: γ=1, α=0
    [InlineData(30, 0.866025403784438646763723170752, 0.125)]  // DEPTH-3: γ=√3/2, α=1/8
    [InlineData(45, 0.707106781186547524400844362105, 0.25)]   // Quarter: γ=√2/2, α=1/4
    [InlineData(60, 0.5,         0.375)]            // KIntermediate: γ=1/2, α=3/8
    [InlineData(90, 0.0,         0.5)]              // Generic: γ=0, α=1/2
    public void AlphaFromTheta_MatchesCanonicalDyadicAnchor(int degrees, double expectedGamma, double expectedAlpha)
    {
        double theta = degrees * System.Math.PI / 180.0;
        double gamma = CanonicalTrigAnchorPi2Inheritance.GammaFromTheta(theta);
        double alpha = CanonicalTrigAnchorPi2Inheritance.AlphaFromTheta(theta);
        Assert.Equal(expectedGamma, gamma, precision: 12);
        Assert.Equal(expectedAlpha, alpha, precision: 12);
        _out.WriteLine($"θ = {degrees}°: γ = {gamma:F14} (Δ={System.Math.Abs(gamma - expectedGamma):G3}), α = {alpha:F14} (Δ={System.Math.Abs(alpha - expectedAlpha):G3})");
    }

    [Theory]
    [InlineData(30, 6.464101615137754587054892683)]   // c² = 2√3 + 3
    [InlineData(45, 2.414213562373095048801688724)]   // c² = 1 + √2 (silver ratio)
    [InlineData(60, 1.0)]                              // uniform Dicke c = 1
    public void CSquaredFromGamma_MatchesClosedFormValues(int degrees, double expectedCSquared)
    {
        double theta = degrees * System.Math.PI / 180.0;
        double gamma = CanonicalTrigAnchorPi2Inheritance.GammaFromTheta(theta);
        double cSquared = CanonicalTrigAnchorPi2Inheritance.CSquaredFromGamma(gamma);
        Assert.Equal(expectedCSquared, cSquared, precision: 11);
    }

    [Fact]
    public void CSquaredFromGamma_RejectsMirrorEndpoint()
    {
        // γ = 1 (Mirror) has c → ∞; explicitly reject.
        Assert.Throws<System.ArgumentOutOfRangeException>(() =>
            CanonicalTrigAnchorPi2Inheritance.CSquaredFromGamma(1.0));
    }

    [Fact]
    public void CSquaredFromGamma_RejectsNegativeGamma()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() =>
            CanonicalTrigAnchorPi2Inheritance.CSquaredFromGamma(-0.1));
    }

    [Fact]
    public void AnchorAlphas_CoverAllFiveCanonicalAngles()
    {
        var anchors = CanonicalTrigAnchorPi2Inheritance.AnchorAlphas;
        Assert.Equal(5, anchors.Count);
        Assert.Equal(0.0, anchors[0]);
        Assert.Equal(1.0 / 8.0, anchors[30]);
        Assert.Equal(1.0 / 4.0, anchors[45]);
        Assert.Equal(3.0 / 8.0, anchors[60]);
        Assert.Equal(1.0 / 2.0, anchors[90]);
    }

    [Fact]
    public void CanonicalAngles_SpanStandardTriangles()
    {
        // The 30°-60°-90° + 45°-45°-90° standard trig triangles use exactly
        // these five angles (with 90° shared).
        var angles = CanonicalTrigAnchorPi2Inheritance.CanonicalAnglesDegrees;
        Assert.Contains(0, angles);
        Assert.Contains(30, angles);
        Assert.Contains(45, angles);
        Assert.Contains(60, angles);
        Assert.Contains(90, angles);
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var c = Build();
        Assert.NotNull(c.Half);
        Assert.NotNull(c.Quarter);
        Assert.NotNull(c.F98LongTime);
    }

    [Fact]
    public void Anchor_References_F99_AndDepthThreeDoc_AndF86b()
    {
        var c = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md F99", c.Anchor);
        Assert.Contains("DEPTH_3_ANCHOR_DERIVED.md", c.Anchor);
        Assert.Contains("depth_3_anchor_derivation.py", c.Anchor);
        Assert.Contains("DickeAnchor.cs", c.Anchor);
        Assert.Contains("KIntermediateAsymptoteQuarterInheritance.cs", c.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsFiveAnchorTable()
    {
        _out.WriteLine("");
        _out.WriteLine("    F99 canonical-trig-angle dyadic anchors");
        _out.WriteLine("    ───────────────────────────────────────");
        _out.WriteLine($"    {"θ",-5} {"γ = cos(θ)",-15} {"c²",-22} {"α = sin²(θ)/2",-15} {"Period 2/3 atom hit"}");
        _out.WriteLine($"    {"---",-5} {"---",-15} {"---",-22} {"---",-15} {"---"}");
        string[][] table = new[]
        {
            new[] { "0°",  "1",          "∞",                       "0",   "noble gases He, Ne, Ar" },
            new[] { "30°", "√3/2",       "2√3 + 3 ≈ 6.464",         "1/8", "Li, Na (1/8) + F, Cl (7/8)" },
            new[] { "45°", "√2/2",       "1 + √2 ≈ 2.414 (silver)", "1/4", "Be, Mg" },
            new[] { "60°", "1/2",        "1 (uniform Dicke)",       "3/8", "B, Al (3/8) + N, P (5/8)" },
            new[] { "90°", "0",          "0",                       "1/2", "H, C, Si" },
        };
        foreach (var row in table)
        {
            _out.WriteLine($"    {row[0],-5} {row[1],-15} {row[2],-22} {row[3],-15} {row[4]}");
        }
        _out.WriteLine("");
        _out.WriteLine("    The 30°-60°-90° + 45°-45°-90° standard trigonometry triangles ARE");
        _out.WriteLine("    the F86b polarity-anchor triangles. Five canonical angles produce");
        _out.WriteLine("    five Pi2 dyadic anchors; with Π²-parity complements (1/8↔7/8, 3/8↔5/8)");
        _out.WriteLine("    the full 9 n/8 fractions cover every period 2/3 element's valence ratio.");
    }
}
