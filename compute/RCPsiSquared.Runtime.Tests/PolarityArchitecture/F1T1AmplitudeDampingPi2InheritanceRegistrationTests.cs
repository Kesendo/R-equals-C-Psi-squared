using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F1T1AmplitudeDampingPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F1T1AmplitudeDampingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_TierIsTier1Candidate()
    {
        // Closed form is Tier1Candidate (verified empirically N=3..6, analytic
        // derivation open per F1 open-question Item 2).
        var registry = BuildBaseRegistry()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Candidate,
            registry.Get<F1T1AmplitudeDampingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F1T1AmplitudeDampingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 32.0, 16.0)]    // N=3: H-part = 2^5 = 32; T1-part = 4^2 = 16
    [InlineData(4, 64.0, 64.0)]
    [InlineData(5, 128.0, 256.0)]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_PrefactorsAgreeAcrossRegistry(int N, double hPart, double t1Part)
    {
        // Cross-registry verification: the H-part and T1-part prefactors computed via
        // the registered claim match the algebraic 2^(N+2) and 4^(N-1) values.
        var registry = BuildBaseRegistry()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .Build();

        var f = registry.Get<F1T1AmplitudeDampingPi2Inheritance>();
        Assert.Equal(hPart, f.HPartPrefactor(N), precision: 12);
        Assert.Equal(t1Part, f.T1PartPrefactor(N), precision: 12);
    }

    [Fact]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_FourMultiplierIsAMinusOne()
    {
        var registry = BuildBaseRegistry()
            .RegisterF1T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(4.0, registry.Get<F1T1AmplitudeDampingPi2Inheritance>().FourMultiplier, precision: 14);
    }

    [Fact]
    public void RegisterF1T1AmplitudeDampingPi2Inheritance_WithoutF1Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                // Missing: RegisterF1Family
                .RegisterF1T1AmplitudeDampingPi2Inheritance()
                .Build());
    }
}
