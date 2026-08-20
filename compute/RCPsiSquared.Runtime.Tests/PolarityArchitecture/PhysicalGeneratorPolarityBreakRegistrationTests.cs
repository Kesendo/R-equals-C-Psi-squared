using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Registration of F155 (<see cref="PhysicalGeneratorPolarityBreak"/>) and the shape of
/// the edge to its typed parent F113.</summary>
public class PhysicalGeneratorPolarityBreakRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus F113 and F113's own transitive
    /// parents (F112 and the F108 Part 1 / Part 2 + BitA-twin chain).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF108Part2Pi2XEvenAlwaysPalindromic()
            .RegisterF108Part1Pi2EvenAlwaysPalindromic()
            .RegisterLindbladBitAPiBalance()
            .RegisterLindbladBitBPiBalance()
            .RegisterLindbladBitBPiBreakMagnitude()
            .RegisterPhysicalGeneratorPolarityBreak()
            .Build();

    [Fact]
    public void RegisterPhysicalGeneratorPolarityBreak_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<PhysicalGeneratorPolarityBreak>());
    }

    [Fact]
    public void TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<PhysicalGeneratorPolarityBreak>().Tier);
    }

    [Fact]
    public void TheParentPointerIsTheRegisteredF113Instance()
    {
        var registry = BuildMinimalRegistry();

        Assert.Same(
            registry.Get<LindbladBitBPiBreakMagnitude>(),
            registry.Get<PhysicalGeneratorPolarityBreak>().BreakMagnitude);
    }

    [Fact]
    public void Ancestors_ReachF113AndThroughItF112()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<PhysicalGeneratorPolarityBreak>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(LindbladBitBPiBreakMagnitude), ancestors);
        Assert.Contains(typeof(LindbladBitBPiBalance), ancestors);
    }

    [Fact]
    public void ItSitsOnTheBitBAxis_WithTheHadamardTwinVerdict()
    {
        // The verdict differs from its own parent's, which is the point rather than a slip:
        // F113 is BitBSpecific because sigma+/sigma- single out the Z energy ladder, and F155
        // makes no channel assumption at all. See the claim's BitATwinStatus for the argument
        // and gate S10 for the measurement.
        var claim = BuildMinimalRegistry().Get<PhysicalGeneratorPolarityBreak>();

        Assert.Equal(Z2Axis.BitB, claim.Z2Axis);
        Assert.Equal(BitATwinClassification.CoveredByHadamardDuality, claim.BitATwinStatus);
        Assert.Equal(BitATwinClassification.BitBSpecific, claim.BreakMagnitude.BitATwinStatus);
    }

    [Fact]
    public void BuildDefault_ContainsPhysicalGeneratorPolarityBreak()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        Assert.True(registry.Contains<PhysicalGeneratorPolarityBreak>());
        Assert.Same(
            registry.Get<LindbladBitBPiBreakMagnitude>(),
            registry.Get<PhysicalGeneratorPolarityBreak>().BreakMagnitude);
    }

    [Fact]
    public void BuildDefault_CountsItInThePolarityCubeMap()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        Assert.Contains(cubeMap.BitBClaims, c => c is PhysicalGeneratorPolarityBreak);
    }
}
