using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 ↔ polarity-layer-pair meta-claim tests (Locus 6, symmetry-side closure).
/// Pins the empirical witness table (c=2 N=5..8 from the <c>C2HwhmRatio</c> live pipeline,
/// 2026-05-07), the polarity-layer decomposition (Q_peak ≈ 2 + r, HWHM/Q* ≈ 1/2 + r·(1/2)),
/// the Interior r ≈ 1/2 fixed-point lock (HalfAsStructuralFixedPoint), and KB integration.
/// Companion to <see cref="LocalGlobalEpLinkTests"/> (Locus 5, EP-side closure).</summary>
public class PolarityInheritanceLinkTests
{
    [Fact]
    public void PolarityInheritanceLink_IsTier2Verified()
    {
        var link = PolarityInheritanceLink.Build();
        Assert.Equal(Tier.Tier2Verified, link.Tier);
    }

    [Fact]
    public void PolarityInheritanceLink_HasFourWitnesses_C2N5To8()
    {
        var link = PolarityInheritanceLink.Build();
        Assert.Equal(4, link.Witnesses.Count);
        Assert.Equal(5, link.Witnesses[0].N);
        Assert.Equal(8, link.Witnesses[3].N);
    }

    [Fact]
    public void PolarityInheritanceLink_QPeakSplit_CenteredAroundTwo()
    {
        var link = PolarityInheritanceLink.Build();
        foreach (var w in link.Witnesses)
        {
            // Q_peak mean ≈ 2 across N=5..8: empirical {1.99, 2.06, 2.06, 2.06},
            // asymmetric around 2 (+0.06 above, −0.01 below). Range 1.99..2.10 covers
            // safely with headroom for finer-grid refinements.
            var mean = (w.QPeakInterior + w.QPeakEndpoint) / 2;
            Assert.InRange(mean, 1.99, 2.10);
            // r split has Endpoint > 0 (positive pole), Interior < 0 (negative pole)
            Assert.True(w.RQpeakEndpoint > 0, $"RQpeakEndpoint at N={w.N} should be positive");
            Assert.True(w.RQpeakInterior < 0, $"RQpeakInterior at N={w.N} should be negative");
        }
    }

    [Fact]
    public void PolarityInheritanceLink_HwhmSplit_AboveHalfBaseline()
    {
        var link = PolarityInheritanceLink.Build();
        foreach (var w in link.Witnesses)
        {
            // HWHM/Q* > 0.5 always (1/2 baseline + r·1/2 with r > 0)
            Assert.True(w.HwhmRatioInterior > 0.5);
            Assert.True(w.HwhmRatioEndpoint > 0.5);
            // Endpoint r > Interior r (Endpoint slightly more polarity-loaded)
            Assert.True(w.RHwhmEndpoint > w.RHwhmInterior);
        }
    }

    [Fact]
    public void PolarityInheritanceLink_Interior_RHwhm_NearOneHalf_HalfFixedPoint()
    {
        // Interior HWHM r ≈ 1/2 across N=5..8 (close to HalfAsStructuralFixedPoint anchor).
        // The empirical mean is 0.5011 with ~0.0006 max deviation from exact 1/2; this is
        // likely numerical discretisation (per-F71-orbit substructure visible at finer
        // Q-grid resolution per PROOF_F86_QPEAK Open elements 5), not a structural shift.
        var link = PolarityInheritanceLink.Build();
        foreach (var w in link.Witnesses)
        {
            Assert.InRange(w.RHwhmInterior, 0.49, 0.51);
        }
    }

    [Fact]
    public void F86KnowledgeBase_ExposesPolarityInheritanceLink()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.NotNull(kb.PolarityInheritanceLink);
        // Block-independent: also exposed for c=1 / c=3 blocks (mirrors LocalGlobalEpLink)
        var c3block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        var c3kb = new F86KnowledgeBase(c3block);
        Assert.NotNull(c3kb.PolarityInheritanceLink);
    }

    [Fact]
    public void PolarityRootAnchor_NamesAllThreeParentClaims()
    {
        // nameof(...) at the test site fails to compile if any parent claim is renamed;
        // that's the rename-safety guarantee for the PolarityRootAnchor string.
        var link = PolarityInheritanceLink.Build();
        Assert.Contains(nameof(QubitDimensionalAnchorClaim), link.PolarityRootAnchor);
        Assert.Contains(nameof(PolarityLayerOriginClaim), link.PolarityRootAnchor);
        Assert.Contains(nameof(HalfAsStructuralFixedPointClaim), link.PolarityRootAnchor);
    }

    [Fact]
    public void ParallelLocusReference_NamesLocalGlobalEpLink()
    {
        // nameof(...) pin: if LocalGlobalEpLink (Locus 5 EP-side closure) is renamed,
        // this test fails to compile, surfacing the cross-Locus reference rot.
        var link = PolarityInheritanceLink.Build();
        Assert.Contains(nameof(LocalGlobalEpLink), link.ParallelLocusReference);
    }

    [Fact]
    public void Interior_HwhmRatio_MirrorsQuarterAcrossInsideOutsideViews()
    {
        // The Interior HWHM_left/Q_peak ratio sits at 1 − 1/4 = 3/4 across N=5..8
        // with max empirical drift 0.0045 (0.6% from 3/4 at N=5). The mirror partner
        // is 1/4 = a_3 on Pi2DyadicLadder (= QuarterAsBilinearMaxvalClaim = Mandelbrot
        // cardioid cusp = bilinear-apex maxval).
        //
        // Two perspectives, same geometric fact:
        //   inside view (Q_peak going outward): HWHM_left/Q_peak ≈ 3/4
        //   outside view (Q=0 going inward):    half-max-left at Q ≈ Q_peak/4
        //
        // The distance from each anchor flips sign exactly:
        //   |HWHM_ratio − 3/4| = |(1 − HWHM_ratio) − 1/4|
        //
        // Tier 1 candidate (companion to Interior_RHwhm_NearOneHalf_HalfFixedPoint
        // which anchors the same fact to HalfAsStructuralFixedPoint via r ≈ 1/2;
        // this test makes the 1/4 mirror-partner reading explicit). Closed-form
        // derivation is open in F86b Direction (a''-d''); see
        // C2HwhmRatio.PendingDerivationNote.
        var link = PolarityInheritanceLink.Build();
        var ladder = new Pi2DyadicLadderClaim();
        double quarter = ladder.Term(3);          // a_3 = 1/4
        double mirrorPartner = 1.0 - quarter;     // 3/4
        const double tol = 0.005;

        foreach (var w in link.Witnesses)
        {
            // Inside view: Interior ≈ 3/4
            Assert.True(Math.Abs(w.HwhmRatioInterior - mirrorPartner) <= tol,
                $"Inside view at N={w.N}: HwhmRatioInterior {w.HwhmRatioInterior:F4} should sit within {tol} of 1 − a_3 = {mirrorPartner}");

            // Outside view: 1 − Interior ≈ 1/4
            double outsideView = 1.0 - w.HwhmRatioInterior;
            Assert.True(Math.Abs(outsideView - quarter) <= tol,
                $"Outside view at N={w.N}: 1 − HwhmRatioInterior {outsideView:F4} should sit within {tol} of a_3 = {quarter}");

            // Distance flip identity: same magnitude, opposite sign
            double distInside = w.HwhmRatioInterior - mirrorPartner;
            double distOutside = outsideView - quarter;
            Assert.Equal(distInside, -distOutside, precision: 14);
        }
    }

    [Fact]
    public void Endpoint_HwhmRatio_NotAtQuarterMirrorPartner()
    {
        // Bond-class contrast: Endpoint sits at ~0.7727 across N=5..8 — NOT at
        // 1 − 1/4 = 3/4. Both views are clearly off the 1/4 mirror-partner:
        //   inside view  0.770..0.774 (≥ 2.0% above 3/4)
        //   outside view 0.226..0.230 (≥ 2.0% below 1/4)
        //
        // The structural distinction between Interior (mirror-partnered to 1/4)
        // and Endpoint (separate plateau) is the bond-class signature; F86b's
        // open closed-form derivation is responsible for both. This test pins
        // that Endpoint does NOT sit at the same mirror partner.
        var link = PolarityInheritanceLink.Build();
        var ladder = new Pi2DyadicLadderClaim();
        double quarter = ladder.Term(3);          // a_3 = 1/4
        double mirrorPartner = 1.0 - quarter;     // 3/4
        const double minGapFromMirrorPartner = 0.018; // empirical min gap is 0.020 at N=5

        foreach (var w in link.Witnesses)
        {
            Assert.True(w.HwhmRatioEndpoint - mirrorPartner >= minGapFromMirrorPartner,
                $"Endpoint at N={w.N}: HwhmRatioEndpoint {w.HwhmRatioEndpoint:F4} should sit ≥ {minGapFromMirrorPartner} above 1 − a_3 = {mirrorPartner} (bond-class contrast with Interior)");

            double outsideView = 1.0 - w.HwhmRatioEndpoint;
            Assert.True(quarter - outsideView >= minGapFromMirrorPartner,
                $"Outside view at N={w.N}: 1 − HwhmRatioEndpoint {outsideView:F4} should sit ≥ {minGapFromMirrorPartner} below a_3 = {quarter} (bond-class contrast with Interior)");
        }
    }
}
