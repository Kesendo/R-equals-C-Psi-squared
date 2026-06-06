using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Witnesses for the soft-certifier's structural ceiling (PROOF_F103 §7.12), AFTER the hidden-Q
/// routing strategy receded it. The colouring strategies are all 2-colourings of the basis-state graph, so
/// they reach only SOFT cases whose graph is bipartite. XX+XZ is soft yet its basis-state graph is
/// NON-bipartite, so no colouring can express it; but the non-diagonal routing now CERTIFIES it (the shared
/// uniform family {P4}), so it is no longer the ceiling. The remaining ceiling is the k-body routed-soft
/// frontier (Stufe B): XZX+XZY+YZX is soft yet beyond the 2-body routing table, so NotCertified. These
/// tests pin both down with the spectral authority, and also pin that XY+YX+XZ+ZX is HARD on the
/// chain.</summary>
public class PalindromeSoftCertifierCeilingTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(l => T(l)).ToList();

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void NonBipartiteSoft_XXxz_IsCertifiedByRouting(int n)
    {
        // XX+XZ is soft by the spectral authority, yet its basis-state graph is non-bipartite, so no chiral
        // K = diag(±1) exists and no colouring can certify it. The hidden-Q routing reaches it anyway
        // (the shared uniform family {P4}): the ceiling has receded, this is now certified, not the wall.
        var terms = H("XX", "XZ");
        var chain = new ChainSystem(n, 1.0, 0.05);
        var bc = BipartiteChirality.Classify(chain, terms);
        Assert.Equal(TrichotomyClass.Soft, bc.ActualClass);   // genuinely soft
        Assert.False(bc.IsBipartite);                         // non-bipartite basis-state graph
        // The non-diagonal routing certifies it; the colourings could not (and the certifier never lies).
        var cert = PalindromeSoftCertifier.Certify(terms, n);
        Assert.True(cert.Certified);
        Assert.Equal(PalindromeSoftCertifier.SoftStrategy.Routing, cert.Strategy);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void KBodyRoutedSoft_IsRealAndBeyondTheRoutingTable(int n)
    {
        // XZX+XZY+YZX is a 3-body routed-soft case: soft by the spectral authority, yet the routing family
        // table is 2-body and cannot reach it, so the certifier returns NotCertified. The new ceiling
        // (Stufe B). Checked at N=4,5 (not N=3, where k=3 fills the whole chain, a different regime).
        var terms = H("XZX", "XZY", "YZX");
        var chain = new ChainSystem(n, 1.0, 0.05);
        Assert.Equal(TrichotomyClass.Soft, PauliPairTrichotomy.Classify(chain, terms));   // genuinely soft
        // The certifier is sound: it does not (and cannot) certify this k-body routed-soft case.
        Assert.False(PalindromeSoftCertifier.Certify(terms, n).Certified);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void PairPlusOddMix_IsHardOnTheChain(int n)
    {
        // XY+YX+XZ+ZX is soft only off-chain (e.g. a triangle); on the chain it is HARD. NotCertified here
        // is correct because it is hard, not because of the ceiling. (Corrects an earlier example.)
        var terms = H("XY", "YX", "XZ", "ZX");
        var chain = new ChainSystem(n, 1.0, 0.05);
        Assert.Equal(TrichotomyClass.Hard, PauliPairTrichotomy.Classify(chain, terms));
        Assert.False(PalindromeSoftCertifier.Certify(terms, n).Certified);
    }
}
