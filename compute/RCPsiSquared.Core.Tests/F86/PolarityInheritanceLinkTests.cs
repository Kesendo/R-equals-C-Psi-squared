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
}
