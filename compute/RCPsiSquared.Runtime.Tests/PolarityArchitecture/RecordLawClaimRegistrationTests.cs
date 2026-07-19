using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Registration gates for the record laws (F135 <see cref="RecordParityLawClaim"/> and
/// F136 <see cref="RecordLetterLawClaim"/>): registered, Tier 1 derived, and carrying the typed
/// spine <see cref="AbsorptionTheoremClaim"/> → F135 → F136 (the pair reduction is the
/// absorption substrate on a pair page; every letter channel is Proposition 1 on a channel
/// class).</summary>
public class RecordLawClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim();

    [Fact]
    public void RegisterRecordParityLaw_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterRecordParityLawClaim()
            .Build();
        Assert.True(registry.Contains<RecordParityLawClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<RecordParityLawClaim>().Tier);
    }

    [Fact]
    public void RegisterRecordLetterLaw_AddsClaim_Tier1Derived_WithParityParent()
    {
        var registry = BuildBaseRegistry()
            .RegisterRecordParityLawClaim()
            .RegisterRecordLetterLawClaim()
            .Build();
        Assert.True(registry.Contains<RecordLetterLawClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<RecordLetterLawClaim>().Tier);
        Assert.Same(registry.Get<RecordParityLawClaim>(), registry.Get<RecordLetterLawClaim>().ParityLaw);
    }

    [Fact]
    public void RecordLetterLaw_AncestorsReachAbsorptionThroughParity()
    {
        var registry = BuildBaseRegistry()
            .RegisterRecordParityLawClaim()
            .RegisterRecordLetterLawClaim()
            .Build();
        var ancestors = registry.AncestorsOf<RecordLetterLawClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(RecordParityLawClaim), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void BuildDefault_ContainsBothRecordLawClaims()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<RecordParityLawClaim>());
        Assert.True(registry.Contains<RecordLetterLawClaim>());
    }
}
