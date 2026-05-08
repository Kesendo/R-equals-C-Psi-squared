using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F88PopcountPairLensRegistrationTests
{
    [Fact]
    public void RegisterF88PopcountPairLens_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 2)
            .Build();

        Assert.True(registry.Contains<F88PopcountPairLens>());
    }

    [Fact]
    public void RegisterF88PopcountPairLens_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 2)
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F88PopcountPairLens>().Tier);
    }

    [Fact]
    public void RegisterF88PopcountPairLens_AncestorsContainsPopcountCoherenceClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 2)
            .Build();

        var ancestors = registry.AncestorsOf<F88PopcountPairLens>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PopcountCoherenceClaim), ancestors);
    }

    [Fact]
    public void RegisterF88PopcountPairLens_C2N6_Pi2OddIsHalf_BothHdClasses()
    {
        // The c=2 stratum F88-fact: at N=6, n_p=1, n_q=2 the lens is Generic α=0.5, so
        // Π²-odd memory is ≈ 0.5 for both HD=1 and HD=3 contributions. F88 alone does
        // NOT distinguish bond-position within the c=2 block; the Schicht-1-wiring
        // makes this fact queryable through the registry.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 2)
            .Build();

        var lens = registry.Get<F88PopcountPairLens>();
        Assert.Equal(ConfigurationKind.Generic, lens.Kind);
        Assert.Equal(0.5, lens.Alpha, precision: 12);
        double pi2OddHd1 = lens.Pi2OddInMemory(hd: 1);
        double pi2OddHd3 = lens.Pi2OddInMemory(hd: 3);
        Assert.Equal(pi2OddHd1, pi2OddHd3, precision: 12);
    }

    [Fact]
    public void RegisterF88PopcountPairLens_PopcountMirror_AlphaIsZero()
    {
        // When n_p + n_q = N, the lens is in PopcountMirror configuration, α = 0.
        // X-flip conjugation cancels all odd-|S| Pauli strings between sectors.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 5, np: 2, nq: 3)
            .Build();

        var lens = registry.Get<F88PopcountPairLens>();
        Assert.Equal(ConfigurationKind.PopcountMirror, lens.Kind);
        Assert.Equal(0.0, lens.Alpha, precision: 12);
    }

    [Fact]
    public void RegisterF88PopcountPairLens_KIntermediate_AlphaIsClosedForm()
    {
        // Even N=6, n_q=N/2=3: K-intermediate. Closed form α = C(6,3)/(2·(C(6,1)+C(6,3))) = 20/(2·(6+20)) = 20/52 = 5/13.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 3)
            .Build();

        var lens = registry.Get<F88PopcountPairLens>();
        Assert.Equal(ConfigurationKind.KIntermediate, lens.Kind);
        Assert.Equal(5.0 / 13.0, lens.Alpha, precision: 12);
    }

    [Fact]
    public void RegisterF88PopcountPairLens_WithoutF88PopcountCoherence_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // Missing: RegisterF88PopcountCoherence
                .RegisterF88PopcountPairLens(N: 6, np: 1, nq: 2)
                .Build());
    }
}
