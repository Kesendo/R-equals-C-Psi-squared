using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class ThreeDephasingDiagonalsOrbitClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsThreeDephasingDiagonalsOrbitClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ThreeDephasingDiagonalsOrbitClaim>());
    }

    [Fact]
    public void Claim_ResolvesWithSharedTypedParents_TheWeldEdge()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<ThreeDephasingDiagonalsOrbitClaim>();
        Assert.NotNull(claim);
        // the dual parentage IS the weld: both must be the SHARED registry instances.
        Assert.Same(registry.Get<MirrorGroupD4Claim>(), claim.MirrorGroup);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Diagonal);
    }

    [Fact]
    public void Claim_HasBothTypedParents_LinkingTheTwoClusters()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<ThreeDephasingDiagonalsOrbitClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(MirrorGroupD4Claim), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void Claim_BatteryAllPass()
    {
        // the live battery: same spectrum, the two basis-S₃ conjugators, D-fix (rate), R-anti (mirror),
        // the orbit = {Q_X,Q_Y,Q_Z}, and the semidirect S₃ ⋉ D₄ structure.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<ThreeDephasingDiagonalsOrbitClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
