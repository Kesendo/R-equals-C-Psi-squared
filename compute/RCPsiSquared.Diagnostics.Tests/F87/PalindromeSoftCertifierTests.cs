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

    [Fact]
    public void Certify_PrefersExcitation_FallsBackToLinear_ElseNotCertified()
    {
        // Pure pairing -> ExcitationPairing (the stronger, topology-independent certificate).
        Assert.Equal(Strategy.ExcitationPairing, PalindromeSoftCertifier.Certify(H("XY", "YX"), 4).Strategy);
        // Pure hopping, chain-bipartite, not a pairing -> LinearSiteColoring.
        Assert.Equal(Strategy.LinearSiteColoring, PalindromeSoftCertifier.Certify(
            new List<PauliTerm> { T("XY"), T("YX", -Complex.One) }, 4).Strategy);
        // F87 windowed hard pair -> neither strategy -> not certified.
        var hard = PalindromeSoftCertifier.Certify(H("XYI", "YIX"), 4);
        Assert.False(hard.Certified);
        Assert.Equal(Strategy.None, hard.Strategy);
    }

    [Fact]
    public void Certify_NeverFalsePositive_AgainstTheSpectralVerdict()
    {
        // A certificate must imply not-hard. Check against the actual trichotomy at N=4.
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var battery = new[]
        {
            H("XY", "YX"),          // pairing, soft
            H("XX", "YY"),          // truly
            H("XY", "YX", "XX"),    // mixed cell-ish; whatever it is, certificate must not lie
            H("XYI", "YIX"),        // hard
            H("XIY", "IXY"),        // hard
            H("XZ", "ZX"),          // soft
        };
        foreach (var terms in battery)
        {
            var cert = PalindromeSoftCertifier.Certify(terms, 4);
            if (cert.Certified)
            {
                var spectral = PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: PauliLetter.Z);
                Assert.NotEqual(TrichotomyClass.Hard, spectral);   // certified => not hard
            }
        }
    }

    [Fact]
    public void Certify_ScalesToLargeN_WithNoLiouvillian()
    {
        // Pure pairing certified at N=40 (excitation strategy is N-independent; no 2^40 matrix).
        var cert = PalindromeSoftCertifier.Certify(H("XY", "YX"), 40);
        Assert.True(cert.Certified);
        Assert.Equal(Strategy.ExcitationPairing, cert.Strategy);
    }
}
