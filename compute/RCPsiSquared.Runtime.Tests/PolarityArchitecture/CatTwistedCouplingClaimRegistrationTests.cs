using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="CatTwistedCouplingClaim"/> (F168 under a twisted coupling): it lands
/// in the default registry as Tier1Derived with F168 as its typed parent.</summary>
public class CatTwistedCouplingClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<CatTwistedCouplingClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<CatTwistedCouplingClaim>().Tier);

    [Fact]
    public void Ancestors_ContainFoneSixtyEight()
    {
        var ancestors = Default().AncestorsOf<CatTwistedCouplingClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(CatTurnedFieldSelfEncounterClaim), ancestors);
    }

    [Fact]
    public void OrderSevenWeight_IsSixOrderingsTimesSixteen() =>
        Assert.Equal(6 * 16, CatTwistedCouplingClaim.OrderSevenWeight);
}
