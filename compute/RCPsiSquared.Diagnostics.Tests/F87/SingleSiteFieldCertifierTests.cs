using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The SingleSiteField strategy: a sum of single-site TRANSVERSE fields (every term weight-1 with
/// letter X or Y) is constructively soft-local. H = Σ_i (a_i X_i + b_i Y_i) gives L = Σ_i L_i over commuting
/// single-site Liouvillians, palindromized by Q = ⊗_i M_i (each M_i the per-site crossover map). Sound by
/// derivation: weight-1 transverse is always soft; a weight-1 Z (longitudinal) field is HARD (the partner of
/// the 0 eigenvalue is absent), so Z is excluded. This is the 4→2 correction: the two I-heavy cases are now
/// certified local, no longer the non-local ceiling.</summary>
public class SingleSiteFieldCertifierTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    private static TrichotomyClass Spectral(IReadOnlyList<PauliTerm> terms, int n) =>
        PauliPairTrichotomy.Classify(new ChainSystem(n, 1.0, 0.05), terms, dephaseLetter: PauliLetter.Z);

    public static IEnumerable<object[]> IHeavyCases() => new[]
    {
        new object[] { new[] { "IXI", "IIY", "YII" } },
        new object[] { new[] { "IYI", "IIX", "XII" } },
    };

    [Theory]
    [MemberData(nameof(IHeavyCases))]
    public void Certify_IHeavy_IsSingleSiteField_AndSpectralSoft(string[] labels)
    {
        var terms = H(labels);
        var cert = PalindromeSoftCertifier.Certify(terms, n: 4);
        Assert.True(cert.Certified);
        Assert.Equal(PalindromeSoftCertifier.SoftStrategy.SingleSiteField, cert.Strategy);
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms, 4));
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms, 5));
    }

    [Theory]
    [InlineData(4)]
    [InlineData(1000)]
    [InlineData(1_000_000)]
    public void Certify_IHeavy_IsSingleSiteField_AtAnyN_TheClassifierIsNFree(int n)
    {
        // The whole point of the certifier is the shortcut around the 2^N wall: it is structural and k-local
        // and NEVER builds the 4^N Liouvillian. The SingleSiteField detection reads only the term labels (it
        // does not use n at all); the earlier strategies are scalable; RoutingKBody is N-independent by
        // additivity. So a chain of a million sites certifies as fast as four. (Contrast: brute-forcing the
        // full 4^N palindrome residual stops near N=7; that wall is exactly what the certifier exists to
        // avoid. Certifying at N = 1_000_000 here completes in microseconds.)
        var terms = H("IXI", "IIY", "YII");
        var cert = PalindromeSoftCertifier.Certify(terms, n);
        Assert.True(cert.Certified);
        Assert.Equal(PalindromeSoftCertifier.SoftStrategy.SingleSiteField, cert.Strategy);
    }

    [Fact]
    public void CertifyBySingleSiteField_True_ForWeightOneTransverse()
    {
        Assert.True(PalindromeSoftCertifier.CertifyBySingleSiteField(H("IXI", "IIY", "YII")));
        Assert.True(PalindromeSoftCertifier.CertifyBySingleSiteField(H("XII")));
        Assert.True(PalindromeSoftCertifier.CertifyBySingleSiteField(H("IYI", "XII")));
    }

    [Fact]
    public void CertifyBySingleSiteField_False_ForLongitudinalZField()
    {
        Assert.False(PalindromeSoftCertifier.CertifyBySingleSiteField(H("IZI")));
        Assert.False(PalindromeSoftCertifier.CertifyBySingleSiteField(H("IXI", "IZI")));
        Assert.Equal(TrichotomyClass.Hard, Spectral(H("IZI"), 4));
    }

    [Fact]
    public void CertifyBySingleSiteField_False_ForWeightTwoOrMore()
    {
        Assert.False(PalindromeSoftCertifier.CertifyBySingleSiteField(H("XIX")));
        Assert.False(PalindromeSoftCertifier.CertifyBySingleSiteField(H("XZX", "XZY")));
        Assert.False(PalindromeSoftCertifier.CertifyBySingleSiteField(new List<PauliTerm>()));
    }

    [Fact]
    public void EveryWeightOneTransverseSet_IsNotHard_SoundnessSweep()
    {
        var probes = new[]
        {
            H("XII"), H("IXI"), H("IIX"), H("YII"), H("IYI"),
            H("IXI", "YII"), H("XII", "IIY"), H("IXI", "IIY", "YII"), H("IYI", "IIX", "XII"),
        };
        foreach (var terms in probes)
        {
            Assert.True(PalindromeSoftCertifier.CertifyBySingleSiteField(terms));
            Assert.NotEqual(TrichotomyClass.Hard, Spectral(terms, 4));
        }
    }
}
