using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;
using Strategy = RCPsiSquared.Diagnostics.F87.PalindromeSoftCertifier.SoftStrategy;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The Liouvillian-free soft-certifier: the σ± pure-pairing test, the colouring strategies, the
/// non-diagonal hidden-Q routing strategy, the orchestration, and the safety property that a certificate
/// is never a false positive.</summary>
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
    public void Certify_RejectsMixedCellAndNonPairingWitnesses()
    {
        // Bug 1: XX (bit_b=0) + XY (bit_b=1) is a mixed Klein cell and spectrally Hard -> not certified.
        Assert.False(PalindromeSoftCertifier.CertifyByLinearSiteColoring(H("XX", "XY"), 4));
        Assert.False(PalindromeSoftCertifier.Certify(H("XX", "XY"), 4).Certified);
        // Bug 2: a pairing plus an odd-flip term (|Δn| mixes 2 and 1) is not a pure pairing and is Hard.
        var pairingPlusOddFlip = new List<PauliTerm> { T("XY"), T("YX"), T("YZ") };
        Assert.False(PalindromeSoftCertifier.IsPurePairing(pairingPlusOddFlip));
        Assert.False(PalindromeSoftCertifier.Certify(pairingPlusOddFlip, 4).Certified);
    }

    [Fact]
    public void ExcitationParity_CertifiesAllOddFlip_TopologyIndependent()
    {
        // Every term flips an odd number of sites (here one), so every basis-edge has odd Δn and the
        // excitation parity n mod 2 two-colours the basis graph -> soft on any topology.
        Assert.True(PalindromeSoftCertifier.IsAllOddFlip(H("XZ", "ZX")));
        Assert.Equal(Strategy.ExcitationParity, PalindromeSoftCertifier.Certify(H("XZ", "ZX"), 4).Strategy);
        Assert.Equal(Strategy.ExcitationParity, PalindromeSoftCertifier.Certify(H("XZ", "ZX"), 40).Strategy);
        // Even or mixed k_xy is not the parity strategy.
        Assert.False(PalindromeSoftCertifier.IsAllOddFlip(H("XY", "YX")));               // k_xy = 2 (even)
        Assert.False(PalindromeSoftCertifier.IsAllOddFlip(H("XY", "YX", "XZ", "ZX")));   // mixed even + odd
        Assert.False(PalindromeSoftCertifier.IsAllOddFlip(H("ZZ")));                     // k_xy = 0 (diagonal)
    }

    [Fact]
    public void IsReversalSymmetric_TrueWhenReversingEachLabelFixesTheMultiset()
    {
        // {XY, YX}: reversing each label maps XY -> YX and YX -> XY, the multiset is unchanged.
        Assert.True(PalindromeSoftCertifier.IsReversalSymmetric(H("XY", "YX")));
        // {XX, XY, YX}: XX is its own reversal; XY <-> YX swap; multiset preserved.
        Assert.True(PalindromeSoftCertifier.IsReversalSymmetric(H("XX", "XY", "YX")));
        // {XX, XY}: XY reverses to YX, which is absent -> not symmetric.
        Assert.False(PalindromeSoftCertifier.IsReversalSymmetric(H("XX", "XY")));
        // {XY} alone reverses to {YX}, a different multiset -> not symmetric.
        Assert.False(PalindromeSoftCertifier.IsReversalSymmetric(H("XY")));
    }

    [Fact]
    public void SiteSwapSymmetry_CertifiesMixedBipartiteReversalSymmetric2Body()
    {
        // The Tom case: XX (bit_b=0) + XY,YX (bit_b=1) is bit_b-MIXED, so the three colourings decline it,
        // but it is 2-body, mask-bipartite, and reversal-symmetric -> certified by the site-swap strategy.
        Assert.True(PalindromeSoftCertifier.CertifyBySiteSwapSymmetry(H("XX", "XY", "YX"), 4));
        Assert.Equal(Strategy.SiteSwapSymmetry, PalindromeSoftCertifier.Certify(H("XX", "XY", "YX"), 4).Strategy);
        // YY (bit_b=0) + XY,YX (bit_b=1) is the same shape -> also certified by the site-swap strategy.
        Assert.Equal(Strategy.SiteSwapSymmetry, PalindromeSoftCertifier.Certify(H("YY", "XY", "YX"), 4).Strategy);
    }

    [Fact]
    public void SiteSwapSymmetry_RejectsThe3BodyKiller_NonAdjacent2Body_AndAsymmetricHardPair()
    {
        // The 2-body gate rejects the 3-body killer XXX+XXY+YXX: reversal-symmetric, bit_b-MIXED, and
        // mask-bipartite, yet spectrally HARD at N=4,5. It is 3-body, so it must not be certified.
        Assert.False(PalindromeSoftCertifier.CertifyBySiteSwapSymmetry(H("XXX", "XXY", "YXX"), 4));
        Assert.False(PalindromeSoftCertifier.Certify(H("XXX", "XXY", "YXX"), 4).Certified);
        // The CONTIGUITY gate rejects the non-adjacent 2-body killer {XIIX, XY, YX}: every term has exactly
        // two non-identity letters, it is bit_b-MIXED, mask-bipartite, and reversal-symmetric, yet XIIX's two
        // X's are 3 sites apart (outside the verified adjacent scope) and the set is spectrally HARD at N=5.
        Assert.False(PalindromeSoftCertifier.CertifyBySiteSwapSymmetry(H("XIIX", "XY", "YX"), 5));
        Assert.False(PalindromeSoftCertifier.Certify(H("XIIX", "XY", "YX"), 5).Certified);
        // The asymmetric pair XX+XY is bit_b-MIXED but not reversal-symmetric (XY reverses to absent YX)
        // and is spectrally hard -> not certified by the site-swap strategy either.
        Assert.False(PalindromeSoftCertifier.CertifyBySiteSwapSymmetry(H("XX", "XY"), 4));
        Assert.False(PalindromeSoftCertifier.Certify(H("XX", "XY"), 4).Certified);
    }

    [Fact]
    public void TryGetUniformFamilyMask_RecognizedBilinearsReturnTheirMask_UnknownLabelsReturnFalse()
    {
        // fam(XX) = {P1, P4} = 1 | 2 = 3 (both crossovers route a same-letter bilinear).
        Assert.True(TwoTermPalindromeRouting.TryGetUniformFamilyMask("XX", out int xx));
        Assert.Equal(3, xx);
        // fam(XZ) = {P4} = 2 (only the Y-router survives the X·Z mix).
        Assert.True(TwoTermPalindromeRouting.TryGetUniformFamilyMask("XZ", out int xz));
        Assert.Equal(2, xz);
        // fam(XY) = {} = 0, but XY is RECOGNIZED, so it returns true with an empty family.
        Assert.True(TwoTermPalindromeRouting.TryGetUniformFamilyMask("XY", out int xy));
        Assert.Equal(0, xy);
        // A non-bilinear single-letter-ish label and a 3-body template are not recognized: false, mask 0.
        Assert.False(TwoTermPalindromeRouting.TryGetUniformFamilyMask("XI", out int xi));
        Assert.Equal(0, xi);
        Assert.False(TwoTermPalindromeRouting.TryGetUniformFamilyMask("XZX", out int xzx));
        Assert.Equal(0, xzx);
    }

    [Fact]
    public void CertifyByRouting_CertifiesSharedUniformFamilyAndTwoTermEscapes()
    {
        // XX (fam {P1,P4}) + XZ (fam {P4}) share {P4}: one uniform Q palindromizes both -> certified.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("XX", "XZ")));
        // Multi-term uniform: XX,XZ,ZX all admit P4 (XX={P1,P4}, XZ=ZX={P4}); intersection {P4} -> certified.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("XX", "XZ", "ZX")));
        // YY (fam {P1,P4}) + YZ (fam {P1}) share {P1} -> certified.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("YY", "YZ")));
        // XY + ZZ: no shared uniform family, but it is a two-term alternating escape -> certified via Classify.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("XY", "ZZ")));
    }

    [Fact]
    public void CertifyByRouting_DeclinesNoSharedFamily_AllMother_AndKBodyLabels()
    {
        // XY (fam {}) + XZ (fam {P4}): no shared family AND not a two-term soft escape (spectrally hard).
        Assert.False(PalindromeSoftCertifier.CertifyByRouting(H("XY", "XZ")));
        // XX + YY is all-Mother: TRULY (canonical Π pairs the spectrum), not soft -> must decline.
        Assert.False(PalindromeSoftCertifier.CertifyByRouting(H("XX", "YY")));
        // A 3-body set has non-recognized labels: out of the 2-body scope (Stufe B) -> decline, no throw.
        Assert.False(PalindromeSoftCertifier.CertifyByRouting(H("XZX", "XZY", "YZX")));
    }

    [Fact]
    public void CertifyByRouting_SingleNonMotherBilinearIsInScope_SingleMotherIsExcluded()
    {
        // A single recognized non-Mother bilinear routes via its own uniform family: XZ alone -> {P4}.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("XZ")));
        // A single Mother bilinear is excluded by the all-Mother gate: ZZ alone is truly, not soft.
        Assert.False(PalindromeSoftCertifier.CertifyByRouting(H("ZZ")));
    }

    [Fact]
    public void CertifyByRouting_IsCoefficientIndependent()
    {
        // The routing strategy reads LABELS only: one uniform Q palindromizes any real linear combination.
        // XX (fam {P1,P4}) + XZ (fam {P4}) share {P4}; varying the coefficients does not change the verdict.
        var terms = new List<PauliTerm> { T("XX", new Complex(2, 0)), T("XZ", new Complex(-3, 0)) };
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(terms));
        Assert.Equal(Strategy.Routing, PalindromeSoftCertifier.Certify(terms, 4).Strategy);
    }

    [Fact]
    public void CertifyByRouting_CertifiesP1NarrowingMultiTerm()
    {
        // The existing routing positives all land on {P4}; this 3-term set narrows to the OTHER family bit.
        // YY={P1,P4} ∩ YZ={P1} ∩ ZY={P1} = {P1}: one uniform Q palindromizes all three -> certified.
        Assert.True(PalindromeSoftCertifier.CertifyByRouting(H("YY", "YZ", "ZY")));
        Assert.Equal(Strategy.Routing, PalindromeSoftCertifier.Certify(H("YY", "YZ", "ZY"), 4).Strategy);
    }

    [Fact]
    public void Certify_LabelsUniformRoutingAsRouting()
    {
        // XX+XZ: non-diagonal hidden-Q soft, missed by all four colourings, certified by Routing.
        Assert.Equal(Strategy.Routing, PalindromeSoftCertifier.Certify(H("XX", "XZ"), 4).Strategy);
        // The multi-term uniform case is also labelled Routing.
        Assert.Equal(Strategy.Routing, PalindromeSoftCertifier.Certify(H("XX", "XZ", "ZX"), 4).Strategy);
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
            H("XZ", "ZX"),          // all-odd-flip, soft (parity strategy)
            H("XZ", "ZX", "YZ", "ZY"),  // all-odd-flip but bit_b-MIXED: soft, yet the parity gate declines it (a miss, not a lie)
            H("XX", "XY"),          // mixed Klein cell, hard (Bug 1 witness)
            new List<PauliTerm> { T("XY"), T("YX"), T("YZ") },  // pairing + odd-flip, hard (Bug 2 witness)
            new List<PauliTerm> { T("XY"), T("YX"), T("XZ"), T("ZX") },  // pairing + odd mix, hard on chain
            H("XZ", "ZXZ"),         // all-odd-flip but bit_b-MIXED ({1,0}) and hard (Bug 3 witness): the
                                    // parity strategy needs the bit_b-homogeneity gate or it false-positives
            H("XX", "XY", "YX"),    // bit_b-MIXED, mask-bipartite, reversal-symmetric, soft (site-swap strategy)
            H("YY", "XY", "YX"),    // same shape, soft (site-swap strategy)
            H("XY", "YX", "XZ", "ZX"),  // bit_b-homogeneous but non-mask-bipartite on the chain: hard, declined
            H("XXX", "XXY", "YXX"), // 3-body killer: reversal-symmetric, mask-bipartite, yet spectrally hard
            H("XIIX", "XY", "YX"),  // non-adjacent 2-body killer: every term 2 non-I letters but XIIX spans 4 sites; hard
            H("XX", "XZ"),          // hidden-Q routing: non-bipartite-soft, certified by Routing (shared {P4})
            H("XX", "XZ", "ZX"),    // multi-term uniform routing-soft (intersection {P4})
            H("YY", "YZ", "ZY"),    // multi-term uniform routing-soft narrowing to the OTHER bit (intersection {P1})
            H("XY", "XZ"),          // routing-declined: no shared family, not a two-term escape, spectrally hard
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
