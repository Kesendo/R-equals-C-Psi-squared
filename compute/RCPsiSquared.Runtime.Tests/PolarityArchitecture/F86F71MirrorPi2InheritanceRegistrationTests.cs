using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86F71MirrorPi2InheritanceRegistrationTests
{
    // Minimal registry for F86F71MirrorPi2Inheritance (Welle 9 cleanup
    // 2026-05-26): the Claim's direct ctor parents are F71MirrorSymmetry +
    // F86MirrorGeneralisationLink — NOT F1. Previously this test also
    // registered the F1 chain transitively (RegisterF1Family +
    // RegisterF1Pi2Inheritance + F88b + F38/F63/F61 BitA chain) which made
    // the Welle 7 F1→F61 propagation surface here as a false-positive
    // breakage. Keeping BuildBaseRegistry minimal-to-the-claim is the
    // F86F71MirrorPi2InheritanceRegistration pattern; transitive deps only
    // appear when the Claim's own ctor needs them.
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterHalfIntegerMirror(5)
            .RegisterF71Family(N: 5)
            .RegisterF71MirrorSymmetryPi2Inheritance();

    [Fact]
    public void RegisterF86F71Mirror_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F86F71MirrorPi2Inheritance>());
    }

    [Fact]
    public void RegisterF86F71Mirror_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F86F71MirrorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF86F71Mirror_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F86F71MirrorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
        Assert.Contains(typeof(F86MirrorGeneralisationLink), ancestors);
    }

    [Theory]
    [InlineData(5, 0, 3)]
    [InlineData(6, 2, 2)]
    public void RegisterF86F71Mirror_MirrorPartnerBondAcrossRegistry(int N, int b, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F86F71MirrorPi2Inheritance>().MirrorPartnerBond(N, b));
    }

    [Fact]
    public void RegisterF86F71Mirror_WithoutF71Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterHalfIntegerMirror(5)
                .RegisterF71MirrorSymmetryPi2Inheritance()
                // Missing: RegisterF71Family
                .RegisterF86F71MirrorPi2Inheritance()
                .Build());
    }
}
