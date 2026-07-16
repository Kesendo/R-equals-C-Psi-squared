using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class DeadSetLawClaimRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus its three typed parents and
    /// their transitive parents (the F131 chain reuses the MirrorOrderSorting minimal
    /// registry's ingredients; ChiralK is a root).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterChiralK()
            .Register<KleinFourCellClaim>(_ => new KleinFourCellClaim())
            .RegisterKleinEightCellClaim()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .RegisterCommutatorDConjugationSign()
            .RegisterMirrorGroupD4Claim()
            .RegisterF108Part2Pi2XEvenAlwaysPalindromic()
            .RegisterF108Part1Pi2EvenAlwaysPalindromic()
            .RegisterLindbladBitAPiBalance()
            .RegisterLindbladBitBPiBalance()
            .RegisterAntilinearTriangleClaim()
            .RegisterSymmetryFamilyInventory()
            .RegisterJointPopcountSectors()
            .RegisterF71MirrorBlockRefinement()
            .RegisterF71AntiPalindromicGammaSpectralInvariance()
            .RegisterF92BondAntiPalindromicJSpectralInvariance()
            .RegisterF93DetuningAntiPalindromicSpectralInvariance()
            .RegisterMirrorOrderSortingClaim()
            .RegisterDeadSetLawClaim()
            .Build();

    [Fact]
    public void RegisterDeadSetLawClaim_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<DeadSetLawClaim>());
    }

    [Fact]
    public void RegisterDeadSetLawClaim_TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<DeadSetLawClaim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainAllThreeTypedParents()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<DeadSetLawClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(ChiralKClaim), ancestors);
        Assert.Contains(typeof(AntilinearTriangleClaim), ancestors);
        Assert.Contains(typeof(MirrorOrderSortingClaim), ancestors);
        // Transitive through F131: the F91-family axes it owns.
        Assert.Contains(typeof(RCPsiSquared.Core.BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var registry = BuildMinimalRegistry();
        var claim = registry.Get<DeadSetLawClaim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsDeadSetLawClaim()
    {
        // The claim is wired into the production registry via KnowledgeRegistryFactory
        // (registered after its last parent MirrorOrderSortingClaim, 2026-07-16).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DeadSetLawClaim>());
        var claim = registry.Get<DeadSetLawClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
