using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F127 (CrossTripleOrthogonalityClaim / CrossTripleOrthogonalityWitness):
/// (a) the live witness passes (on-variety zeros over both wall primes, nonzero controls);
/// (b) registry wiring: the claim resolves from the default registry, Tier1Candidate (the
/// code-trust caveat is the named reason), with its typed parent SeedExistenceCountingClaim;
/// (c) the claim's constants match the F127 registry entry.</summary>
public class CrossTripleOrthogonalityTests
{
    [Fact]
    public void Witness_SlicePasses()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.StartsWith("PASS", w.Summary);
        Assert.Contains("controls nonzero", w.Summary);
    }

    [Fact]
    public void Claim_ResolvesFromDefaultRegistry_WithTypedParent()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.NotNull(claim);
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.SeedExistence);
        Assert.Contains("F127", claim.DisplayName);
    }

    [Fact]
    public void Claim_ConstantsMatchTheRegistryEntry()
    {
        Assert.Equal(527, CrossTripleOrthogonalityClaim.CertifiedTasks);
        Assert.Equal(17, CrossTripleOrthogonalityClaim.PrimeCount);
        Assert.Equal(510.0, CrossTripleOrthogonalityClaim.Log2PrimeProduct);
        Assert.True(CrossTripleOrthogonalityClaim.AssemblyDIsSymbolic);
    }

    [Fact]
    public void Claim_AnchorNamesTheWallAndTheAssembly()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.Contains("grid_proof_sweep.py", claim.Anchor);
        Assert.Contains("assembly_d_symbolic.py", claim.Anchor);
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md", claim.Anchor);
    }
}
