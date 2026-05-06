using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.Tests.F86;

public class F86ExtendedClaimsTests
{
    // Reduced grid + small location subset to keep tests under ~30 sec.
    // Range extended to 6.0 so x=+1.0 plateau (Q ≈ 2·Q_peak ≈ 5 for Endpoint) stays in grid.
    private static readonly double[] TestQGrid = ResonanceScan.LinearQGrid(0.20, 6.00, 40);
    private static readonly IReadOnlyList<(int C, int N)> TestLocations = new[] { (2, 5), (2, 6), (2, 7) };

    private static WitnessCache TestCache() => new(TestQGrid);

    [Fact]
    public void TwoLevelEpModel_HigherK_ShiftsCenterAndDecayTime()
    {
        const double γ0 = 0.05, gEff = 2.83, j = 0.025;
        var k1 = new TwoLevelEpModel(γ0, j, gEff, k: 1);
        var k2 = new TwoLevelEpModel(γ0, j, gEff, k: 2);
        var k3 = new TwoLevelEpModel(γ0, j, gEff, k: 3);

        Assert.Equal(-4.0 * γ0 * 1, k1.TraceCentre, 12);
        Assert.Equal(-4.0 * γ0 * 2, k2.TraceCentre, 12);
        Assert.Equal(-4.0 * γ0 * 3, k3.TraceCentre, 12);

        Assert.Equal(1.0 / (4.0 * γ0 * 1), k1.DecayTimeAtEp, 12);
        Assert.Equal(1.0 / (4.0 * γ0 * 2), k2.DecayTimeAtEp, 12);
        Assert.Equal(1.0 / (4.0 * γ0 * 3), k3.DecayTimeAtEp, 12);

        // Discriminant is k-independent — same EP location for every level.
        Assert.Equal(k1.Discriminant, k2.Discriminant, 12);
        Assert.Equal(k1.Discriminant, k3.Discriminant, 12);
    }

    [Fact]
    public void TwoLevelEpModel_RejectsKLessThan1()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new TwoLevelEpModel(0.05, 0.025, 2.83, k: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new TwoLevelEpModel(0.05, 0.025, 2.83, k: -1));
    }

    [Fact]
    public void TwoLevelEpModel_TauSquared_IsZeroAtEp_AndPositiveAfter()
    {
        const double γ0 = 0.05, gEff = 2.0;
        var atEp = new TwoLevelEpModel(γ0, j: 0.05, gEff);  // J·gEff = 2·γ₀
        Assert.Equal(0.0, atEp.TauSquared, 12);
        var post = new TwoLevelEpModel(γ0, j: 0.10, gEff);
        Assert.True(post.TauSquared > 0);
    }

    [Fact]
    public void PerBlockQPeakClaim_StandardList_HasC2_C3_C4_C5()
    {
        var standard = PerBlockQPeakClaim.Standard;
        Assert.Equal(4, standard.Count);
        Assert.Contains(standard, c => c.Chromaticity == 2 && c.Caveat is not null);
        Assert.Equal(1.6, standard.First(c => c.Chromaticity == 3).QPeakValue);
        Assert.Equal(1.8, standard.First(c => c.Chromaticity == 4).QPeakValue);
        Assert.Equal(1.8, standard.First(c => c.Chromaticity == 5).QPeakValue);
        Assert.True(standard.First(c => c.Chromaticity == 4).Saturated);
    }

    [Fact]
    public void PerBlockQPeakClaim_GammaZeroExtraction_GivesGammaFromMeasuredJStar()
    {
        var c4 = PerBlockQPeakClaim.Standard.First(c => c.Chromaticity == 4);
        // If measured J* = 0.09 and Q_peak(c=4) = 1.8, then γ₀ = 0.09 / 1.8 = 0.05.
        Assert.Equal(0.05, c4.ExtractGammaZero(0.09), 12);
    }

    [Fact]
    public void PerBondQPeakWitnessTable_EndpointAndInterior_HaveLocationsList()
    {
        Assert.Equal(10, F86StandardLocations.Full.Count);
        // Building the table creates self-computing witnesses (no scan triggered yet).
        var endpoint = PerBondQPeakWitnessTable.BuildEndpoint();
        var interior = PerBondQPeakWitnessTable.BuildInterior();
        Assert.Equal(F86StandardLocations.Full.Count, endpoint.Witnesses.Count);
        Assert.Equal(F86StandardLocations.Full.Count, interior.Witnesses.Count);
        Assert.All(endpoint.Witnesses, w => Assert.Equal(BondClass.Endpoint, w.BondClass));
        Assert.All(interior.Witnesses, w => Assert.Equal(BondClass.Interior, w.BondClass));
    }

    [Fact]
    public void Sigma0Scaling_TrajectoryCrossingValueIsTwoSqrtTwoCMinusOne_AndComputesSelfFromInterChannelSvd()
    {
        // σ_0 witnesses use small InterChannelSvd computations — fast, no need for reduced grid.
        // SigmaZeroChromaticityScaling subsumes the old c=2 Sigma0AsymptoticClaim by covering all c.
        // Note (2026-05-08): the 2√(2(c−1)) value is a trajectory CROSSING (sweet-spot at c=2 N=7,
        // bit-exact), not an asymptote; see SigmaZeroChromaticityScaling class summary. The test
        // name was renamed accordingly. Default Ns extended from {5,6,7,8} to {5,6,7,8,9} to
        // expose the post-crossing region where ratio > 2.0 at c=2.
        var claim = new SigmaZeroChromaticityScaling(0.05, chromaticities: new[] { 2 });
        Assert.Equal(2.0 * Math.Sqrt(2.0), SigmaZeroChromaticityScaling.Asymptote(c: 2), 12);
        Assert.Equal(5, claim.Witnesses.Count);
        // The σ_0 values are now computed live from InterChannelSvd.Build(block).Sigma0.
        Assert.True(claim.Witnesses[0].Sigma0 > 2.7 && claim.Witnesses[0].Sigma0 < 2.9,
            $"σ_0 at c=2 N=5 = {claim.Witnesses[0].Sigma0:F4} should be near 2√2 ≈ 2.828");
        // c=2 N=7 should hit 2√2 to ~4 decimals (the structural sweet spot crossing).
        Assert.InRange(claim.Witnesses[2].Sigma0,
            SigmaZeroChromaticityScaling.Asymptote(2) - 0.001,
            SigmaZeroChromaticityScaling.Asymptote(2) + 0.001);
        // c=2 N=9 is post-crossing: σ_0 > 2√2.
        Assert.True(claim.Witnesses[4].Sigma0 > 2.0 * Math.Sqrt(2.0),
            $"σ_0 at c=2 N=9 = {claim.Witnesses[4].Sigma0:F4} should be > 2√2 (post-crossing region)");
    }

    [Fact]
    public void ShapeFunctionWitnesses_InteriorAndEndpoint_HavePeakAtX0_AndComputeLive()
    {
        var cache = TestCache();
        var interior = new ShapeFunctionWitnesses(BondClass.Interior, 0.05, TestLocations, xPoints: null, cache);
        var endpoint = new ShapeFunctionWitnesses(BondClass.Endpoint, 0.05, TestLocations, xPoints: null, cache);

        // Peak at x=0: y range collapses to 1.0 (within ~0.01 due to linear interpolation
        // between grid points; definitionally Q*·(1+0) = Q_peak hits the curve maximum but
        // the interpolated value depends on the grid spacing).
        var interiorPeak = interior.Points.First(p => p.X == 0);
        Assert.InRange(interiorPeak.YMin, 0.99, 1.001);
        Assert.InRange(interiorPeak.YMax, 0.99, 1.001);

        // Plateau values computed from reduced (c=2 only) sample with 40-point grid;
        // ranges are wider than full-grid Python to accommodate sampling differences.
        Assert.InRange(interior.PostPeakPlateau, 0.75, 0.92);
        Assert.InRange(endpoint.PostPeakPlateau, 0.83, 0.97);

        // Endpoint plateau at x=+1 is higher than Interior plateau (structural).
        var interiorAtX1 = interior.Points.First(p => p.X == 1.0);
        var endpointAtX1 = endpoint.Points.First(p => p.X == 1.0);
        Assert.True(endpointAtX1.YMin > interiorAtX1.YMax,
            $"Endpoint plateau ({endpointAtX1.YMin:F3}) should exceed Interior plateau ({interiorAtX1.YMax:F3})");
    }

    [Fact]
    public void ShapeFunctionWitnesses_SpreadIsBoundedNearPeak()
    {
        // Near x=0 the shape collapse is tightest. With reduced grid (30 points) the spread
        // is wider than the full-grid 2% but still bounded; structural observation.
        var cache = TestCache();
        var interior = new ShapeFunctionWitnesses(BondClass.Interior, 0.05, TestLocations, xPoints: null, cache);
        var nearPeakPoints = interior.Points.Where(p => Math.Abs(p.X) <= 0.2);
        foreach (var pt in nearPeakPoints)
            Assert.True(pt.SpreadPercent < 5.0,
                $"Interior near-peak spread at x={pt.X} too large: {pt.SpreadPercent}%");
    }

    [Fact]
    public void DressedModeWeightClaim_99PctAtPeak_31PctAtPlateau()
    {
        var claim = new DressedModeWeightClaim();
        Assert.Equal(0.99, claim.WeightAtQPeak);
        Assert.Equal(0.31, claim.WeightAtPlateau);
        Assert.Equal(20.0, claim.PlateauQ);
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void OpenQuestion_StandardList_HasItem1Prime_Item4Prime_Item5()
    {
        var open = F86OpenQuestions.Standard;
        Assert.Equal(3, open.Count);
        Assert.All(open, q => Assert.Equal(Tier.OpenQuestion, q.Tier));
        Assert.Contains(open, q => q.Name.Contains("Item 1'"));
        Assert.Contains(open, q => q.Name.Contains("Item 4'"));
        Assert.Contains(open, q => q.Name.Contains("Item 5"));
    }

    [Fact]
    public void ChiralAiiiClassification_IsTier1Derived_AndNotPT()
    {
        var classification = new ChiralAiiiClassification();
        Assert.Equal(Tier.Tier1Derived, classification.Tier);
        Assert.Contains("AIII", classification.DisplayName);
        Assert.Contains("NOT", classification.DisplayName);
    }

    [Fact]
    public void F86KnowledgeBase_ForC2N5_HasSigma0Scaling_AndPerBlockTable()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);

        Assert.NotNull(kb.Sigma0Scaling);
        Assert.Equal(4, kb.PerBlockQPeaks.Count);
        Assert.Equal(BondClass.Endpoint, kb.EndpointPerBondTable.BondClass);
        Assert.Equal(BondClass.Interior, kb.InteriorPerBondTable.BondClass);
        Assert.NotNull(kb.DressedModeWeight);
        Assert.NotNull(kb.AlgebraicClass);
        Assert.Equal(3, kb.OpenQuestions.Count);
    }

    [Fact]
    public void F86KnowledgeBase_ForC3Block_HasHigherKLevels()
    {
        // c=3 has 2 channel pairs (k=1, k=2); the chromaticity-scaling claim covers all c.
        var block = new CoherenceBlock(N: 7, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        var kb = new F86KnowledgeBase(block);

        // Higher-k hierarchy: 2 levels (k=1, k=2).
        Assert.Equal(2, kb.HigherKLevels.Count);
        Assert.Equal(1, kb.HigherKLevels[0].K);
        Assert.Equal(2, kb.HigherKLevels[1].K);
        // Decay times decrease with k: t_k = 1/(4γ₀·k)
        Assert.True(kb.HigherKLevels[0].DecayTimeAtEp > kb.HigherKLevels[1].DecayTimeAtEp);
    }
}
