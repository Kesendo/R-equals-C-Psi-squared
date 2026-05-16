using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.F86.JordanWigner;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

/// <summary>Diagnostic probe (Phase B test for Γ_pair derivation, 2026-05-16):
/// dump per-bond cluster-pair Q_EP distributions at N=5 to see whether the JW-
/// cluster-pair data physically support the 2-pair superposition lift hypothesis
/// from `simulations/_f86_lift_via_superposition_probe.py`.
///
/// <para>The Python superposition probe found that 2-pair sums at specific
/// (Q_EP_1, Q_EP_2) tuples reproduce Interior 0.7506 / Endpoint 0.7728:
/// Interior target: Q1≈1.01, Q2≈2.86 (ratio 2.84). Endpoint target: Q1≈0.58,
/// Q2≈1.86 (ratio 3.21). If the actual JW cluster-pair distribution per bond
/// does NOT have similar Q_EP spreads, the superposition hypothesis is
/// structurally wrong, not just numerically off.</para></summary>
public class F86ClusterPairQEpDistributionProbe
{
    private readonly ITestOutputHelper _out;
    public F86ClusterPairQEpDistributionProbe(ITestOutputHelper output) => _out = output;

    [Fact]
    public void EmitClusterPairQEpDistributionN5()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var unified = JwBondQPeakUnified.Build(block);
        var empirical = C2HwhmRatio.Build(block);

        _out.WriteLine($"N=5 c=2 cluster-pair Q_EP distribution per bond.");
        _out.WriteLine($"Empirical anchors: Interior 0.7506, Endpoint 0.7728 (HWHM_left/Q_peak).");
        _out.WriteLine($"Phase 4 superposition search: Interior needs Q1≈1.01, Q2≈2.86;");
        _out.WriteLine($"                              Endpoint needs Q1≈0.58, Q2≈1.86.");
        _out.WriteLine("");

        for (int b = 0; b < unified.Bonds.Count; b++)
        {
            var bond = unified.Bonds[b];
            var emp = empirical.Witnesses[b];
            _out.WriteLine($"=== Bond {b} ({bond.BondClass}) ===");
            _out.WriteLine($"  Empirical: Q_peak={emp.QPeak:F4}, HWHM_left={emp.HwhmLeft:F4}, ratio={emp.HwhmLeftOverQPeak:F4}");
            _out.WriteLine($"  Predicted Q_peak (current Tier2): {bond.QPeakPredicted:F4} (regime {bond.Regime})");
            _out.WriteLine($"  Top1/Top2 frob² ratio: {bond.Top1Top2FrobRatio:F3}");
            _out.WriteLine($"  Cluster pairs with non-NaN Q_EP (sorted by Q_EP):");
            _out.WriteLine($"    (this output relies on internal GetPairQEp — emit raw cluster-pair sequence)");
            _out.WriteLine($"  See full RankedPairs structure via JwBondClusterPairAffinity inspect.");
            _out.WriteLine("");
        }

        // For each bond emit the cluster-pair affinity ranking
        var affinity = JwBondClusterPairAffinity.Build(block);
        _out.WriteLine("=== JwBondClusterPairAffinity raw ranking (top 5 per bond) ===");
        _out.WriteLine("");
        foreach (var bondAff in affinity.Bonds)
        {
            _out.WriteLine($"Bond {bondAff.Bond} ({bondAff.BondClass}):");
            int shown = 0;
            foreach (var pair in bondAff.RankedPairs)
            {
                if (shown >= 5) break;
                _out.WriteLine($"  pair ({pair.Cluster1Index}, {pair.Cluster2Index}): " +
                    $"|Δδ|={pair.AbsoluteDeltaδ:F3}, sizes=({pair.Cluster1Size},{pair.Cluster2Size}), " +
                    $"frob²={pair.FrobeniusSquared:F4}");
                shown++;
            }
            _out.WriteLine("");
        }
    }
}
