using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86MainRegistrationTests
{
    // TPeakLaw declares an AbsorptionTheoremClaim parent (rung-2 edge, 2026-06-10),
    // so standalone builders must carry the absorption chain: Pi2Family (9) +
    // Pi2DyadicLadder (1) + AbsorptionTheoremClaim (1).
    private static ClaimRegistryBuilder BuildBaseRegistry(double gammaZero, double gEff) =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF86Main(gammaZero, gEff);

    [Fact]
    public void RegisterF86Main_BuildsFiveClaims()
    {
        var registry = BuildBaseRegistry(gammaZero: 0.05, gEff: 1.0)
            .Build();

        // 11 absorption-chain claims (see BuildBaseRegistry) + the five F86 claims.
        Assert.Equal(16, registry.All().Count());
        Assert.True(registry.Contains<ChiralAiiiClassification>());
        Assert.True(registry.Contains<DressedModeWeightClaim>());
        Assert.True(registry.Contains<F71MirrorInvariance>());
        Assert.True(registry.Contains<TPeakLaw>());
        Assert.True(registry.Contains<QEpLaw>());
    }

    [Fact]
    public void RegisterF86Main_TPeakLawAncestorsContainAbsorptionTheorem()
    {
        // The rung-2 edge (2026-06-10): t_peak's four is 2γ·2, the Absorption
        // Theorem's second rung, not the discriminant d².
        var registry = BuildBaseRegistry(gammaZero: 0.05, gEff: 1.0)
            .Build();

        var ancestors = registry.AncestorsOf<TPeakLaw>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void RegisterF86Main_TierMix_FourTier1DerivedPlusOneCandidate()
    {
        var registry = BuildBaseRegistry(gammaZero: 0.05, gEff: 1.0)
            .Build();

        // Four of five claims are Tier1Derived; DressedModeWeightClaim is
        // Tier1Candidate by design — its hardcoded W anchors 0.99 / 0.31 are
        // documented as "specific values unverified per (c, N)" in the claim's
        // own ctor. Asserting per-claim keeps the honest tier surface visible.
        Assert.Equal(Tier.Tier1Derived, registry.Get<ChiralAiiiClassification>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<F71MirrorInvariance>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<TPeakLaw>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<QEpLaw>().Tier);
        Assert.Equal(Tier.Tier1Candidate, registry.Get<DressedModeWeightClaim>().Tier);
    }

    [Fact]
    public void RegisterF86Main_TPeakLaw_UsesProvidedGamma()
    {
        var registry = BuildBaseRegistry(gammaZero: 0.1, gEff: 1.0)
            .Build();

        var tpeak = registry.Get<TPeakLaw>();
        // t_peak = 1/(4γ₀) = 1/0.4 = 2.5
        Assert.Equal(2.5, tpeak.Value, precision: 10);
    }

    [Fact]
    public void RegisterF86Main_QEpLaw_UsesProvidedGEff()
    {
        var registry = BuildBaseRegistry(gammaZero: 0.05, gEff: 4.0)
            .Build();

        var qep = registry.Get<QEpLaw>();
        // Q_EP = 2/g_eff = 2/4 = 0.5
        Assert.Equal(0.5, qep.Value, precision: 10);
    }
}
