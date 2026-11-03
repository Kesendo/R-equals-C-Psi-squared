using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 ↔ polarity-layer-pair meta-claim tests (Locus 6, symmetry-side closure).
/// Pins the witness table (c=2 N=5..8 class envelope values, checked against a live
/// <c>ResonanceScan</c>), the polarity-layer reading (Q_peak ≈ 2 + r, HWHM/Q* ≈ 1/2 + r·(1/2)),
/// the Interior r_HWHM drift above 1/2, and KB integration.
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
    public void PolarityInheritanceLink_QPeakSplit_PolesOnEitherSideOfTwo()
    {
        // Endpoint above 2, Interior below; symmetric about 2 at N = 5 (mean 2.0001), the
        // Interior pole drifting toward 2 from N = 6 (mean 2.073..2.079).
        var link = PolarityInheritanceLink.Build();
        foreach (var w in link.Witnesses)
        {
            Assert.True(w.RQpeakEndpoint > 0, $"RQpeakEndpoint at N={w.N} should be positive");
            Assert.True(w.RQpeakInterior < 0, $"RQpeakInterior at N={w.N} should be negative");
        }
        for (int k = 1; k < link.Witnesses.Count; k++)
            Assert.True(link.Witnesses[k].RQpeakInterior > link.Witnesses[k - 1].RQpeakInterior,
                $"Interior pole should move toward 2 with N (N={link.Witnesses[k].N})");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void PolarityInheritanceLink_Table_IsTheLiveEnvelope(int n)
    {
        // The table holds what ResonanceScan computes (gate 1e-10, the per-bond pin law).
        var w = PolarityInheritanceLink.Build().Witnesses.Single(x => x.N == n);
        var curve = new RCPsiSquared.Core.Resonance.ResonanceScan(new CoherenceBlock(N: n, n: 1, gammaZero: 0.05)).ComputeKCurve();
        var pi = curve.Peak(RCPsiSquared.Core.Resonance.BondClass.Interior);
        var pe = curve.Peak(RCPsiSquared.Core.Resonance.BondClass.Endpoint);
        Assert.Equal(w.QPeakInterior, pi.QPeak, 1e-10);
        Assert.Equal(w.QPeakEndpoint, pe.QPeak, 1e-10);
        Assert.Equal(w.HwhmRatioInterior, pi.HwhmLeftOverQPeak!.Value, 1e-10);
        Assert.Equal(w.HwhmRatioEndpoint, pe.HwhmLeftOverQPeak!.Value, 1e-10);
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
    public void PolarityInheritanceLink_Interior_RHwhm_RisesAwayFromOneHalf()
    {
        // Interior r_HWHM sits above the 1/2 fixed point and moves away from it with N
        // (0.505, 0.514, 0.518, 0.520 at N = 5..8): not locked to HalfAsStructuralFixedPoint.
        var link = PolarityInheritanceLink.Build();
        for (int k = 0; k < link.Witnesses.Count; k++)
        {
            Assert.True(link.Witnesses[k].RHwhmInterior > 0.5);
            if (k > 0) Assert.True(link.Witnesses[k].RHwhmInterior > link.Witnesses[k - 1].RHwhmInterior);
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
