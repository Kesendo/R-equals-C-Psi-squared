using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class NivenRationalityRootClaimTests
{
    [Fact]
    public void Build_WiresBothSeFaceParents_AndIsTier1Derived()
    {
        var claim = NivenRationalityRootClaim.Build();
        // Tier1Derived: pure number theory, both parents are Tier1Derived (no cap).
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.NotNull(claim.BandEdge);
        Assert.NotNull(claim.Rates);
        Assert.Equal(Tier.Tier1Derived, claim.BandEdge.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.Rates.Tier);

        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.BandEdge, kids);
        Assert.Contains(claim.Rates, kids);
    }

    [Fact]
    public void BandEdgeDegree_IsEulerPhi2mOver2()
    {
        // [Q(2cos(π/(N+1))):Q] = φ_euler(2(N+1))/2: rational (deg 1) for N≤2, quadratic (deg 2) for N=3,4,5
        Assert.Equal(1, NivenRationalityRootClaim.BandEdgeDegree(1));  // 2cos(π/2)=0
        Assert.Equal(1, NivenRationalityRootClaim.BandEdgeDegree(2));  // 2cos(π/3)=1
        Assert.Equal(2, NivenRationalityRootClaim.BandEdgeDegree(3));  // √2
        Assert.Equal(2, NivenRationalityRootClaim.BandEdgeDegree(4));  // φ
        Assert.Equal(2, NivenRationalityRootClaim.BandEdgeDegree(5));  // √3
        Assert.Equal(3, NivenRationalityRootClaim.BandEdgeDegree(6));  // first cubic
        // degree at most 2 iff N≤5: rational for N≤2, quadratic surd exactly for N∈{3,4,5}
        Assert.Equal(new[] { 1, 2 }, Enumerable.Range(1, 8).Where(n => NivenRationalityRootClaim.BandEdgeDegree(n) == 1));
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, Enumerable.Range(1, 8).Where(n => NivenRationalityRootClaim.BandEdgeDegree(n) <= 2));
    }

    [Fact]
    public void BandEdge_N2IsRational_NotAQuadraticSurd_OnEveryClaimSurface()
    {
        var claim = NivenRationalityRootClaim.Build();
        var imFace = Assert.Single(
            ((IInspectable)claim).Children,
            child => child.DisplayName.StartsWith("IM-face:"));

        Assert.Equal(1, NivenRationalityRootClaim.BandEdgeDegree(2));
        Assert.Equal(new[] { 3, 4, 5 },
            Enumerable.Range(1, 8).Where(n => NivenRationalityRootClaim.BandEdgeDegree(n) == 2));
        foreach (var surface in new[] { claim.Name, claim.Summary, imFace.Summary })
        {
            Assert.Contains("degree at most 2 iff N≤5", surface);
            Assert.Contains("quadratic surd exactly for N∈{3,4,5}", surface);
        }
    }

    [Fact]
    public void FirstOrderRateCombAllRational_IffNPlus1InCrystallographicSet()
    {
        // First-order F65 comb coefficients are rational iff N+1 ∈ {1,2,3,4,6}, i.e. N ∈ {1,2,3,5} for N≥1.
        Assert.Equal(new[] { 1, 2, 3, 5 }, Enumerable.Range(1, 8).Where(NivenRationalityRootClaim.FirstOrderRateCombAllRational));
        Assert.False(NivenRationalityRootClaim.FirstOrderRateCombAllRational(4));   // N=4 the first irrational (golden)
    }

    [Fact]
    public void LiveClaim_UsesCanonicalPositiveF65FirstOrderEndpointCombScope()
    {
        var claim = NivenRationalityRootClaim.Build();
        var reFace = Assert.Single(
            ((IInspectable)claim).Children,
            child => child.DisplayName.StartsWith("RE-face:"));

        Assert.Equal(0.5, claim.Rates.FirstOrderRateCoefficient(3, 1), precision: 12);
        Assert.Equal(0.5, claim.Rates.FirstOrderRateTerm(3, 1, 1.0), precision: 12);
        Assert.Equal(
            claim.Rates.FirstOrderRateTerm(3, 1, 1.0),
            claim.Rates.SingleExcitationRate(3, 1, 1.0),
            precision: 12);
        Assert.Contains("a_k = (4/(N+1))·sin²(kπ/(N+1))", claim.Name);
        Assert.Contains("α_k^full = γ₀·a_k + O(γ₀³/J²)", claim.Name);
        Assert.Contains("α_k^full/γ₀ = a_k + O((γ₀/J)²)", claim.Name);
        Assert.Contains("Exact Niven rationality belongs to the first-order coefficient comb", claim.Name);
        Assert.Contains("uniform open XX chain with one dephased endpoint", claim.Name);
        Assert.Contains("first-order coefficient comb", claim.Name);
        Assert.Contains("relative full-L rate shift is O((γ₀/J)²)", claim.Name);
        Assert.Contains("absolute shift δα_k = O(γ₀³/J²)", claim.Name);
        Assert.Contains("no exact finite-γ₀/J full-L rationality is claimed", claim.Name);
        Assert.Contains("a_k = (4/(N+1))·sin²(kπ/(N+1))", claim.Summary);
        Assert.Contains("α_k^full = γ₀·a_k + O(γ₀³/J²)", claim.Summary);
        Assert.Contains("no exact finite-γ₀/J full-L rationality", claim.Summary);
        Assert.Contains("uniform open XX chain with one dephased endpoint", reFace.Summary);
        Assert.Contains("exact Niven rationality", reFace.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("a_k = (4/(N+1))·sin²(kπ/(N+1))", reFace.Summary);
        Assert.Contains("α_k^full = γ₀·a_k + O(γ₀³/J²)", reFace.Summary);
        Assert.Contains("relative full-L rate shift is O((γ₀/J)²)", reFace.Summary);
        Assert.Contains("absolute shift δα_k = O(γ₀³/J²)", reFace.Summary);
        Assert.DoesNotContain("−2γ·sin²", claim.Name);
        Assert.DoesNotContain("α_k/γ₀ = (4/(N+1))", claim.Name);
    }

    [Fact]
    public void ThirdFace_IsTheF6QEdgeGain_NotACouplingComplexityIdentity()
    {
        var claim = NivenRationalityRootClaim.Build();
        Assert.Contains("F6 Q-edge gain", claim.Summary);

        var face = Assert.Single(
            ((IInspectable)claim).Children,
            child => child.DisplayName.StartsWith("F6 face: the Q-edge gain"));
        Assert.Contains("1+cos(π/N)", face.Summary);
        Assert.DoesNotContain("coupling-created complexity", face.Summary);
    }
}
