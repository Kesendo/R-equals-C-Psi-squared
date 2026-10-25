using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="CatTurnedFieldSelfEncounterClaim"/> (F168): it lands in the default
/// registry as Tier1Derived with its two typed parents, F63LCommutesPi2Pi2Inheritance (the parity) and
/// XGlobalEigenstateMirrorPi2Inheritance (the two points); its closed form reads the page's numbers.</summary>
public class CatTurnedFieldSelfEncounterClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<CatTurnedFieldSelfEncounterClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<CatTurnedFieldSelfEncounterClaim>().Tier);

    [Fact]
    public void Ancestors_ContainBothParents()
    {
        var ancestors = Default().AncestorsOf<CatTurnedFieldSelfEncounterClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F63LCommutesPi2Pi2Inheritance), ancestors);
        Assert.Contains(typeof(XGlobalEigenstateMirrorPi2Inheritance), ancestors);
    }

    [Fact]
    public void LeadingOddMoment_ReadsThePagesNumbers()
    {
        // N = 2, gamma = 1 on both bits (sigma = 2), h = 1, sin^2 = 1: (2/3) 5! = 80, the page's 80 gamma h^4 sin^2(2 phi)
        Assert.Equal(new BigRational(80), CatTurnedFieldSelfEncounterClaim.LeadingOddMoment(2, 2, 1, 1));
        // N = 3, sigma = 3 gamma at gamma = 1, h = 1, sin^2 = 1: (3/3) 7! = 5040
        Assert.Equal(new BigRational(5040), CatTurnedFieldSelfEncounterClaim.LeadingOddMoment(3, 3, 1, 1));
        Assert.Equal(7, CatTurnedFieldSelfEncounterClaim.FirstOddOrder(3));
    }
}
