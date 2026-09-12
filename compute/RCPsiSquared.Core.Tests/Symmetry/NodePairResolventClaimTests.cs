using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class NodePairResolventClaimTests
{
    private static NodePairResolventClaim BuildClaim() =>
        new(new SeatCutBlindnessClaim(new F4KernelDimensionByComponentsClaim()));

    [Fact]
    public void Claim_IsTier1Derived_AndAnchoredToTheExactProof()
    {
        var claim = BuildClaim();

        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("docs/proofs/PROOF_NODE_PAIR_RESOLVENT.md", claim.Anchor);
        Assert.Contains("Theorem 1", claim.Anchor);
        Assert.Contains("Corollary B", claim.Anchor);
        Assert.Contains("Corollary C", claim.Anchor);
        Assert.DoesNotContain("F65", $"{claim.Name} {claim.Anchor} {claim.Summary}");
    }

    [Fact]
    public void Claim_StatesSymmetryFreeTheorem_AndScopedUniformCentreIff()
    {
        var claim = BuildClaim();
        string surface = $"{claim.Name} {claim.Summary} {claim.Scope}";

        Assert.Contains("zero-free", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("both nodes", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("sufficient", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("epsilon != -1", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("non-incident bond", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("uniform centre-watched", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("if and only if", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("r != 0,+1,-1", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("off-centre", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("nonuniform", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("bound state in the continuum", surface, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Claim_HoldsItsSingleCoreOwnedParent_AndRejectsNull()
    {
        var parent = new SeatCutBlindnessClaim(new F4KernelDimensionByComponentsClaim());
        var claim = new NodePairResolventClaim(parent);

        Assert.Same(parent, claim.SeatBlindness);
        Assert.Throws<ArgumentNullException>(() => new NodePairResolventClaim(null!));
    }

    [Fact]
    public void Claim_BreadcrumbsTheLiveDiagnosticsRoot_WithoutOwningItsType()
    {
        var claim = BuildClaim();
        string children = string.Join("\n", claim.Children.Select(child =>
            $"{child.DisplayName} {child.Summary}"));

        Assert.Contains("inspect --root nodepair", children, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain(typeof(NodePairResolventClaim).Assembly.GetReferencedAssemblies(),
            reference => string.Equals(reference.Name, "RCPsiSquared.Diagnostics",
                StringComparison.Ordinal));
    }
}
