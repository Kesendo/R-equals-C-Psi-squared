using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;
using Strategy = RCPsiSquared.Diagnostics.F87.PalindromeSoftCertifier.SoftStrategy;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The Liouvillian-free soft-certifier: the σ± pure-pairing test, the two strategies, the
/// orchestration, and the safety property that a certificate is never a false positive.</summary>
public class PalindromeSoftCertifierTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static PauliTerm T(string label, Complex c) => new(PauliLabel.Parse(label), c);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(l => T(l)).ToList();

    [Fact]
    public void IsPurePairing_TrueForTheSymmetricPairTerm()
    {
        Assert.True(PalindromeSoftCertifier.IsPurePairing(H("XY", "YX")));   // -2i(s+s+ - s-s-)
    }

    [Fact]
    public void IsPurePairing_FalseForHoppingMixedOrDiagonal()
    {
        Assert.False(PalindromeSoftCertifier.IsPurePairing(H("XY")));        // pair + hop
        Assert.False(PalindromeSoftCertifier.IsPurePairing(
            new List<PauliTerm> { T("XY"), T("YX", -Complex.One) }));        // XY - YX, pure hopping
        Assert.False(PalindromeSoftCertifier.IsPurePairing(H("XX", "YY")));  // symmetric hopping
        Assert.False(PalindromeSoftCertifier.IsPurePairing(H("ZI")));        // pure diagonal (no flip)
    }

    [Fact]
    public void LinearSiteColoring_TrueForChainBipartiteHopping_FalseForHardOrDiagonal()
    {
        // XY-YX is pure hopping; on the chain its flip-masks are bipartite -> certified.
        Assert.True(PalindromeSoftCertifier.CertifyByLinearSiteColoring(
            new List<PauliTerm> { T("XY"), T("YX", -Complex.One) }, 4));
        // XYI+YIX is the F87 windowed hard pair: non-bipartite flip-masks -> not certified.
        Assert.False(PalindromeSoftCertifier.CertifyByLinearSiteColoring(H("XYI", "YIX"), 4));
        // A pure-diagonal term (ZZ) lifts the diagonal, which the chiral K cannot negate -> not certified
        // even though the flip part (XX) is bipartite. (This guards against a false-positive certificate.)
        Assert.False(PalindromeSoftCertifier.CertifyByLinearSiteColoring(H("XX", "YY", "ZZ"), 4));
    }
}
