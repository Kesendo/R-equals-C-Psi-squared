using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

[Trait("Category", "ROUTE_B_A2_N6_REMAINDER")]
public class RouteBN6RemainderBoundClaimRegistrationTests
{
    private static readonly Lazy<RCPsiSquared.Runtime.ObjectManager.ClaimRegistry> Registry =
        new(KnowledgeRegistryFactory.BuildDefault);

    private static Claim RegisteredClaim()
    {
        var claim = Registry.Value.All()
            .SingleOrDefault(c => c.GetType().Name == nameof(RouteBN6RemainderBoundClaim));
        Assert.NotNull(claim);
        return claim;
    }

    [Fact]
    public void ClaimIsRegisteredAsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, RegisteredClaim().Tier);
    }

    /// <summary>Exactly one typed parent. F163 is a premise of the ball computation; the N=4
    /// certificates share the contraction design but are a method precedent, not a premise,
    /// and must not appear as an edge.</summary>
    [Fact]
    public void TheOnlyTypedParentIsTheN6UnfoldingClaim()
    {
        var parents = Registry.Value.AllEdges()
            .Where(e => e.Child == typeof(RouteBN6RemainderBoundClaim))
            .Select(e => e.Parent.Name)
            .ToList();
        Assert.Equal(new[] { nameof(RouteBN6A2UnfoldingClaim) }, parents);
    }

    /// <summary>Every anchored path is a file that exists; a moved producer or artifact breaks
    /// this rather than leaving the claim pointing at nothing. The house anchor audit skips this
    /// claim entirely, having no parameterless constructor to build it with, and checks only .md
    /// tokens where it does run, so this is the claim's only anchor coverage and it covers the
    /// producers and the certificate as well. Tokenisation follows the house audit exactly, so a
    /// later anchor written with its other separators does not fail a correct claim.</summary>
    [Fact]
    public void EveryAnchorPathResolves()
    {
        var repoRoot = FindRepoRoot();
        var paths = RegisteredClaim().Anchor
            .Split(new[] { " + ", " / ", ", " }, StringSplitOptions.RemoveEmptyEntries)
            .Select(raw => raw.Trim())
            .Select(token =>
            {
                int strip = token.IndexOfAny(new[] { ' ', '#', '(', ':' });
                return strip > 0 ? token[..strip] : token;
            })
            .Where(path => !string.IsNullOrEmpty(path))
            .ToList();
        Assert.Equal(4, paths.Count);
        foreach (var relative in paths)
            Assert.True(File.Exists(Path.Combine(repoRoot, relative)), $"anchor missing: {relative}");
    }

    private static string FindRepoRoot()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
            if (Directory.Exists(Path.Combine(directory.FullName, "docs", "proofs")))
                return directory.FullName;
        throw new DirectoryNotFoundException("Cannot locate the repository root.");
    }
}
