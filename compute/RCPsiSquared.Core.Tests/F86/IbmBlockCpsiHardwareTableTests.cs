using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F86;

public class IbmBlockCpsiHardwareTableTests
{
    private static IbmBlockCpsiHardwareTable Table() => new();

    [Fact]
    public void Claim_IsTier2Verified()
    {
        Assert.Equal(Tier.Tier2Verified, Table().Tier);
    }

    [Fact]
    public void Claim_AnchorReferences_ProofBlockCpsiQuarter_AndPipelineScript()
    {
        var t = Table();
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER.md", t.Anchor);
        Assert.Contains("_block_cpsi_lens_ibm_snapshots.py", t.Anchor);
    }

    [Fact]
    public void Witnesses_HasExpected_FourBackends_TimesFourScenarios_TimesTwoBlocks()
    {
        var ws = Table().Witnesses;
        Assert.Equal(IbmBlockCpsiHardwareTable.Backends.Count
                     * IbmBlockCpsiHardwareTable.Scenarios.Count
                     * 2,
                     ws.Count);
        Assert.Equal(32, ws.Count);
    }

    [Fact]
    public void Witnesses_EveryBackendScenarioBlock_PresentExactlyOnce()
    {
        var ws = Table().Witnesses;
        foreach (string backend in IbmBlockCpsiHardwareTable.Backends)
        foreach (string scen in IbmBlockCpsiHardwareTable.Scenarios)
        foreach (int blockN in new[] { 0, 1 })
        {
            var matches = ws.Where(w => w.Backend == backend
                                     && w.Scenario == scen
                                     && w.BlockN == blockN).ToList();
            Assert.Single(matches);
        }
    }

    [Fact]
    public void EveryWitness_SatisfiesTheorem2_CBlockBelowOrAtQuarter()
    {
        // The structural test: every pinned C_block measurement must respect the
        // Tier-1-derived universal ceiling (PROOF_BLOCK_CPSI_QUARTER Theorem 2).
        const double tol = 1e-12;
        foreach (var w in Table().Witnesses)
        {
            Assert.True(w.CBlockMeasured <= BlockCoherenceContent.Quarter + tol,
                $"Theorem 2 violation: {w.Backend} {w.Scenario} block({w.BlockN},{w.BlockN + 1}) " +
                $"= {w.CBlockMeasured:F6} > {BlockCoherenceContent.Quarter}");
            Assert.True(w.CBlockIdealContinuous <= BlockCoherenceContent.Quarter + tol,
                $"Theorem 2 violation in ideal column: {w.Scenario} block({w.BlockN},{w.BlockN + 1}) " +
                $"= {w.CBlockIdealContinuous:F6} > {BlockCoherenceContent.Quarter}");
        }
    }

    [Fact]
    public void IdealAsymmetryRatios_AreSymmetric_ForPi2EvenScenarios()
    {
        // heisenberg (XX+YY+ZZ) and truly (XX+YY) are Π²-even: block(0,1) = block(1,2)
        // at the continuous-time noiseless level → ratio = 1.000.
        var t = Table();
        Assert.Equal(1.0, t.IdealAsymmetryRatios["heisenberg"], precision: 4);
        Assert.Equal(1.0, t.IdealAsymmetryRatios["truly"], precision: 4);
    }

    [Fact]
    public void IdealAsymmetryRatios_AreAboveOne_ForPi2OddScenarios()
    {
        // soft (XY+YX, Π²-odd) and hard (XX+XY, mixed) break popcount-block symmetry:
        // block(0,1) > block(1,2) → ratio > 1 strictly. soft is more asymmetric than hard.
        var t = Table();
        double soft = t.IdealAsymmetryRatios["soft"];
        double hard = t.IdealAsymmetryRatios["hard"];
        Assert.True(soft > 1.0, $"soft asymmetry should be > 1; got {soft:F3}");
        Assert.True(hard > 1.0, $"hard asymmetry should be > 1; got {hard:F3}");
        Assert.True(soft > hard,
            $"soft asymmetry should exceed hard (Π²-odd more strongly than mixed); got soft={soft:F3} hard={hard:F3}");
    }

    [Fact]
    public void IdealAsymmetryRatios_PinnedValues_MatchTheRecordedRatios()
    {
        // Pinned values from simulations/results/_block_cpsi_lens_ibm_snapshots.txt:
        // ideal heisenberg block(0,1)=0.0642 / block(1,2)=0.0642 ≈ 1.000
        // ideal truly      block(0,1)=0.1064 / block(1,2)=0.1064 ≈ 1.000
        // ideal soft       block(0,1)=0.1838 / block(1,2)=0.0499 ≈ 3.683
        // ideal hard       block(0,1)=0.1707 / block(1,2)=0.0785 ≈ 2.175
        var r = Table().IdealAsymmetryRatios;
        Assert.Equal(1.000, r["heisenberg"], precision: 3);
        Assert.Equal(1.000, r["truly"],      precision: 3);
        Assert.Equal(3.683, r["soft"],       precision: 3);
        Assert.Equal(2.175, r["hard"],       precision: 3);
    }

    [Fact]
    public void HardwareResolution_RangeAcrossIbmWitnesses_IsOrderOfMagnitude()
    {
        // Tom's resolution worry refuted: hardware values do NOT all sit "near 1/4".
        // The 24 ibm_* witnesses should span at least one order of magnitude in
        // FractionOfQuarter.
        var ibmWits = Table().Witnesses
            .Where(w => w.Backend.StartsWith("ibm_"))
            .ToList();
        Assert.Equal(24, ibmWits.Count);
        double minFrac = ibmWits.Min(w => w.FractionOfQuarter);
        double maxFrac = ibmWits.Max(w => w.FractionOfQuarter);
        Assert.True(maxFrac / minFrac > 10.0,
            $"hardware resolution should be > 10× across the lens; got {maxFrac / minFrac:F1}× (min {minFrac:F3}, max {maxFrac:F3})");
    }

    [Fact]
    public void Witness_ExposesFractionOfQuarter_AndDeltaFromIdeal()
    {
        var w = Table().Witnesses.First(x => x.Backend == "ibm_marrakesh"
                                          && x.Scenario == "hard"
                                          && x.BlockN == 0);
        Assert.Equal(w.CBlockMeasured / 0.25, w.FractionOfQuarter, precision: 12);
        Assert.Equal(w.CBlockMeasured - w.CBlockIdealContinuous, w.DeltaFromIdeal, precision: 12);
    }

    [Fact]
    public void Children_IncludeWitnessesGroup_AndOpenQuestionAndAsymmetryNodes()
    {
        IInspectable c = Table();
        var labels = c.Children.Select(ch => ch.DisplayName).ToList();
        Assert.Contains("source pipeline", labels);
        Assert.Contains("Theorem 2 ceiling", labels);
        Assert.Contains("hardware resolution", labels);
        Assert.Contains("OpenQuestion: Π²-odd asymmetry derivation", labels);
        Assert.Contains("witnesses (32 rows: 4 backends × 4 scenarios × 2 blocks)", labels);
    }

    [Fact]
    public void F86KnowledgeBase_ExposesIbmBlockCpsiHardwareTable_AsTier2Empirical()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.NotNull(kb.IbmBlockCpsiHardwareTable);
        Assert.Equal(Tier.Tier2Verified, kb.IbmBlockCpsiHardwareTable.Tier);
    }
}
