using System;
using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="SectorBraidWitness"/>, the live face of the N-dependent multi-sector monodromy
/// verdict. At N=4 the census node must read CONFINED (no dense diagonal core carries the braid) and the F89d
/// cross-fold self-fold must be machine zero, with the dilute (1,1) control reading node (not braid). The N=5
/// SPREAD verdict and the diagonal-spectator byte-identity are pinned by <see cref="MultiSectorMonodromyCensus"/>
/// tests and shown live via <c>--root sectorbraid --N 5</c>.</summary>
public class SectorBraidWitnessTests
{
    [Fact]
    public void Witness_RejectsBelowAnchor()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SectorBraidWitness(3));
    }

    [Fact(DisplayName = "SectorBraidWitness N=4: CONFINED verdict + cross-fold gate")]
    [Trait("Category", "SLOW_MSM")]
    public void Witness_N4_ReadsConfined_AndCrossFoldGate()
    {
        var witness = new SectorBraidWitness(4);
        var nodes = witness.Children.ToList();

        // Node 1 — the census verdict is CONFINED at N=4 (the dense diagonal core is braid-free), never SPREAD.
        Assert.Contains(nodes, n => n.DisplayName.Contains("CONFINED"));
        Assert.DoesNotContain(nodes, n => n.DisplayName.Contains("SPREAD"));

        // Node 2 — the additivity + cross-fold gate: the dilute (1,1) control is a node (not a braid), and at N=4
        // the diagonal-spectator byte-identity is (correctly) reported as an N≥5 phenomenon, not asserted here.
        var gate = nodes.First(n => n.DisplayName.Contains("cross-fold gate"));
        Assert.Contains("HasDefective=False", gate.Summary);
        Assert.Contains("N≥5 phenomenon", gate.Summary);
    }
}
