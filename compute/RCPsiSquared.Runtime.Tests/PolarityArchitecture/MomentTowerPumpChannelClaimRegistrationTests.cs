using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MomentTowerPumpChannelClaimRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus its two typed parents and
    /// their own transitive parents (F113's F112 + F108 Part 1/Part 2 + BitA-twin chain,
    /// and F84's ladder + F82 + F81 Pi2-Foundation chain).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            // F84 chain (ladder + Pi2 family + F81 + F82).
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop()
            .RegisterF1PalindromeIdentity()
            .RegisterF81Pi2Inheritance()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            // F113 chain (F108 Parts 1+2, BitA twin, F112).
            .RegisterF108Part2Pi2XEvenAlwaysPalindromic()
            .RegisterF108Part1Pi2EvenAlwaysPalindromic()
            .RegisterLindbladBitAPiBalance()
            .RegisterLindbladBitBPiBalance()
            .RegisterLindbladBitBPiBreakMagnitude()
            .RegisterMomentTowerPumpChannelClaim()
            .Build();

    [Fact]
    public void RegisterMomentTowerPumpChannelClaim_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<MomentTowerPumpChannelClaim>());
    }

    [Fact]
    public void RegisterMomentTowerPumpChannelClaim_TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<MomentTowerPumpChannelClaim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainBothTypedParents()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<MomentTowerPumpChannelClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(LindbladBitBPiBreakMagnitude), ancestors);
        Assert.Contains(typeof(F84ThermalAmplitudeDampingPi2Inheritance), ancestors);
        // Transitive: F113's F112 foundation, and F84's F82 mother claim.
        Assert.Contains(typeof(LindbladBitBPiBalance), ancestors);
        Assert.Contains(typeof(F82T1AmplitudeDampingPi2Inheritance), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var registry = BuildMinimalRegistry();
        var claim = registry.Get<MomentTowerPumpChannelClaim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsMomentTowerPumpChannelClaim()
    {
        // The claim is wired into the production registry via KnowledgeRegistryFactory
        // (registered after the antilinear triangle; both parents, F113 and F84, are
        // registered earlier, 2026-06-11).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<MomentTowerPumpChannelClaim>());
        var claim = registry.Get<MomentTowerPumpChannelClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
