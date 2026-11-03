using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86PerF71OrbitObservationRegistrationTests
{
    [Fact]
    public void RegisterF86PerF71OrbitObservation_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 6)
            .RegisterF86PerF71OrbitObservation()
            .Build();

        Assert.True(registry.Contains<PerF71OrbitObservation>());
    }

    [Fact]
    public void RegisterF86PerF71OrbitObservation_AncestorsContainsF71BondOrbitDecomposition()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 6)
            .RegisterF86PerF71OrbitObservation()
            .Build();

        var ancestors = registry.AncestorsOf<PerF71OrbitObservation>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F71BondOrbitDecomposition), ancestors);
    }

    [Fact]
    public void RegisterF86PerF71OrbitObservation_TierIsTier2Empirical()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 6)
            .RegisterF86PerF71OrbitObservation()
            .Build();

        Assert.Equal(Tier.Tier2Empirical, registry.Get<PerF71OrbitObservation>().Tier);
    }

    [Fact]
    public void RegisterF86PerF71OrbitObservation_WitnessesPinTheC2N6Inversion()
    {
        // The exact witness that motivates the JW track's truly-vs-flanking-innermost
        // question: at c=2 N=6, central Q_peak 1.428 < flanking Q_peak 1.637. Wiring
        // PerF71OrbitObservation makes this fact a queryable Schicht-1 claim instead of
        // a memory note.
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 6)
            .RegisterF86PerF71OrbitObservation()
            .Build();

        var obs = registry.Get<PerF71OrbitObservation>();
        var c2n6 = obs.Witnesses.Single(w => w.Chromaticity == 2 && w.N == 6);
        // The frozen envelope values, identical to the table's literals.
        Assert.Equal(2.580604998574177, c2n6.QPeakPerOrbit[0]);    // Endpoint
        Assert.Equal(1.6366804166405813, c2n6.QPeakPerOrbit[1]);   // flanking
        Assert.Equal(1.4283108372718716, c2n6.QPeakPerOrbit[2]);   // truly-innermost (centralIndex=2)
        Assert.True(c2n6.QPeakPerOrbit[2] < c2n6.QPeakPerOrbit[1]);
    }

    [Fact]
    public void RegisterF86PerF71OrbitObservation_WithoutF71Family_Throws_MissingParent()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF86PerF71OrbitObservation()
                .Build());
    }
}
