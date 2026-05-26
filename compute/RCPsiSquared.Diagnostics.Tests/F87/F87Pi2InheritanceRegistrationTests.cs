using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F87Pi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            // Welle 7 (2026-05-25) wired F1 → F61 (bit_a twin), and the F61
            // chain pulls in F38 → F63 → BitA references → Pi2OperatorSpaceMirror
            // → F88b dyadic anchors. Mirror the F1Pi2InheritanceRegistrationTests
            // base registry exactly so this test runs with the same Pi²
            // inheritance topology as production.
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterF38BitAInvolutionInheritance()
            .RegisterF63BitAReference()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF1Pi2Inheritance()
            .RegisterF87Family();

    [Fact]
    public void RegisterF87Pi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F87Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F87Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        // KleinFourCellClaim was historically asserted "not an ancestor" because
        // F87 entstand vor Klein cells, sie sind konzeptuell unterschiedlich.
        // Post Welle 7 (2026-05-25), Welle 8 (2026-05-26) wired F1 → F61 → F63
        // → F38 → Pi2OperatorSpaceMirror → F88bStaticDyadicAnchor →
        // PopcountCoherence → KleinFourCellClaim, so KleinFour is now a
        // TRANSITIVE ancestor of F87 (not a direct parent). The conceptual
        // distinction holds (F87 doesn't directly use Klein cells); the
        // inheritance-graph relationship is now structural.
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_TransitivelyInheritedTwoFactorIsTwo()
    {
        // The "2" in F87's discriminator M = F1's residual; via F1Pi2Inheritance the "2"
        // is a_0 on the Pi2 ladder. Cross-registry verification.
        var registry = BuildBaseRegistry()
            .RegisterF87Pi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F87Pi2Inheritance>().TransitivelyInheritedTwoFactor, precision: 14);
    }

    [Fact]
    public void RegisterF87Pi2Inheritance_WithoutF1Inheritance_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                .RegisterF88bPopcountCoherence()
                .RegisterF88bStaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                .RegisterF38BitAInvolutionInheritance()
                .RegisterF63BitAReference()
                .RegisterF38Pi2InvolutionPi2Inheritance()
                .RegisterF63LCommutesPi2Pi2Inheritance()
                .RegisterF61BitAParityPi2Inheritance()
                // Missing: RegisterF1Pi2Inheritance — the throw must surface F1
                // itself as the unresolved parent, not an upstream Pi² dep.
                .RegisterF87Family()
                .RegisterF87Pi2Inheritance()
                .Build());
    }
}
