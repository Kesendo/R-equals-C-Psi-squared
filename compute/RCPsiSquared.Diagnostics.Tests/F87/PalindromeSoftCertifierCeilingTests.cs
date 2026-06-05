using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Witnesses for the soft-certifier's structural ceiling (PROOF_F103 §7.12). The certifier's
/// strategies are all 2-colourings of the basis-state graph, so they can only ever reach SOFT cases whose
/// graph is bipartite. There are soft Hamiltonians whose basis-state graph is NON-bipartite (XX+XZ): soft
/// by a non-diagonal mechanism no colouring can express. The certifier returns NotCertified for them and
/// must (it is sound, not complete). These tests pin that down with the spectral authority, and also
/// correct an earlier mislabelled example by pinning that XY+YX+XZ+ZX is HARD on the chain.</summary>
public class PalindromeSoftCertifierCeilingTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(l => T(l)).ToList();

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void NonBipartiteSoft_IsRealAndBeyondAnyColouring(int n)
    {
        // XX+XZ is soft by the spectral authority, yet its basis-state graph is non-bipartite, so no chiral
        // K = diag(±1) exists. No colouring (scalable or not) can certify it -- the structural ceiling.
        var terms = H("XX", "XZ");
        var chain = new ChainSystem(n, 1.0, 0.05);
        var bc = BipartiteChirality.Classify(chain, terms);
        Assert.Equal(TrichotomyClass.Soft, bc.ActualClass);   // genuinely soft
        Assert.False(bc.IsBipartite);                         // non-bipartite basis-state graph
        Assert.False(bc.Agrees);                              // so the bipartite criterion mispredicts (Hard)
        // The certifier is sound: it does not (and cannot) certify this soft case.
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
