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
        // rational iff N≤2; single quadratic surd iff N≤5
        Assert.Equal(new[] { 1, 2 }, Enumerable.Range(1, 8).Where(n => NivenRationalityRootClaim.BandEdgeDegree(n) == 1));
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, Enumerable.Range(1, 8).Where(n => NivenRationalityRootClaim.BandEdgeDegree(n) <= 2));
    }

    [Fact]
    public void RatesAllRational_IffNPlus1InCrystallographicSet()
    {
        // rates rational iff N+1 ∈ {1,2,3,4,6}, i.e. N ∈ {1,2,3,5} for N≥1
        Assert.Equal(new[] { 1, 2, 3, 5 }, Enumerable.Range(1, 8).Where(NivenRationalityRootClaim.RatesAllRational));
        Assert.False(NivenRationalityRootClaim.RatesAllRational(4));   // N=4 the first irrational (golden)
    }
}
