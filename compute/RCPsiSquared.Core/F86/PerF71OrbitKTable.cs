using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F71-orbit-grouped K-resonance witness table for c=2 K-resonance data.
/// Tier2Verified. <see cref="Build"/> consumes
/// <see cref="F71BondOrbitDecomposition"/> (F71, Tier1Derived) and
/// <see cref="C2HwhmRatio.Build"/>.Witnesses (F86 per-bond), groups by orbit, and
/// averages 2-bond orbits with F71 mirror invariance guarded at 1e-6 max-deviation.
/// Live-computed counterpart to <see cref="PerF71OrbitObservation"/> (which is
/// frozen 9-case sweep across (c, N)). Anchor: PROOF_F86_QPEAK Statement 2 +
/// Statement 3 (F71 mirror invariance).</summary>
public sealed class PerF71OrbitKTable : Claim
{
    public CoherenceBlock Block { get; }
    public IReadOnlyList<OrbitKWitness> OrbitWitnesses { get; }

    public PerF71OrbitKTable(
        CoherenceBlock block,
        IReadOnlyList<OrbitKWitness> orbitWitnesses)
        : base("c=2 F71-orbit K-resonance witness table",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Statement 2 + Statement 3")
    {
        Block = block;
        OrbitWitnesses = orbitWitnesses;
    }

    public static PerF71OrbitKTable Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"PerF71OrbitKTable applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var orbits = new F71BondOrbitDecomposition(block.N).Orbits;
        var hwhmRatio = C2HwhmRatio.Build(block);
        var bondWitnesses = hwhmRatio.Witnesses;

        var orbitWitnesses = new List<OrbitKWitness>(orbits.Count);
        foreach (var orbit in orbits)
        {
            if (orbit.IsSelfPaired)
            {
                var w = bondWitnesses[orbit.BondA];
                orbitWitnesses.Add(new OrbitKWitness(
                    orbit,
                    QPeak: w.QPeak,
                    HwhmLeft: w.HwhmLeft,
                    HwhmLeftOverQPeak: w.HwhmLeftOverQPeak,
                    KMax: w.KMax));
            }
            else
            {
                var wA = bondWitnesses[orbit.BondA];
                var wB = bondWitnesses[orbit.BondB!.Value];
                double maxDev = Math.Max(
                    Math.Max(Math.Abs(wA.QPeak - wB.QPeak),
                             Math.Abs(wA.HwhmLeft - wB.HwhmLeft)),
                    Math.Abs(wA.KMax - wB.KMax));
                if (maxDev > 1e-6)
                    throw new InvalidOperationException(
                        $"F71 mirror invariance violated at orbit {{b={orbit.BondA}, b={orbit.BondB}}}: " +
                        $"max-Δ = {maxDev:E2} exceeds 1e-6");
                orbitWitnesses.Add(new OrbitKWitness(
                    orbit,
                    QPeak: 0.5 * (wA.QPeak + wB.QPeak),
                    HwhmLeft: 0.5 * (wA.HwhmLeft + wB.HwhmLeft),
                    HwhmLeftOverQPeak: 0.5 * (wA.HwhmLeftOverQPeak + wB.HwhmLeftOverQPeak),
                    KMax: 0.5 * (wA.KMax + wB.KMax)));
            }
        }
        return new PerF71OrbitKTable(block, orbitWitnesses);
    }

    public override string DisplayName =>
        $"F71-orbit K-resonance table (N={Block.N}, {OrbitWitnesses.Count} orbits)";

    public override string Summary =>
        $"{OrbitWitnesses.Count} F71 orbits at c=2 N={Block.N}; live counterpart to PerF71OrbitObservation ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("orbit count", summary: OrbitWitnesses.Count.ToString());
            yield return InspectableNode.Group("orbits",
                OrbitWitnesses.Cast<IInspectable>().ToArray());
        }
    }
}
