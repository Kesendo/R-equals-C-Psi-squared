using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F71-orbit-grouped K-resonance witnesses for one c=2 block. Tier2Verified;
/// live-computed counterpart to <see cref="PerF71OrbitObservation"/> (frozen 9-case sweep).</summary>
public sealed class PerF71OrbitKTable : Claim
{
    public CoherenceBlock Block { get; }
    public IReadOnlyList<OrbitKWitness> OrbitWitnesses { get; }

    public PerF71OrbitKTable(
        CoherenceBlock block,
        IReadOnlyList<OrbitKWitness> orbitWitnesses)
        : base("c=2 F71-orbit K-resonance witness table",
               Tier.Tier2Verified,
               F86Anchors.Statement2Plus3)
    {
        Block = block;
        OrbitWitnesses = orbitWitnesses;
    }

    /// <summary>Builds the table for one block: walks
    /// <see cref="F71BondOrbitDecomposition"/> orbits and joins each with the
    /// <see cref="C2HwhmRatio"/> per-bond witnesses. 2-bond orbits average their two
    /// bond witnesses; the F71-mirror invariance guard rejects max-Δ &gt; 1e-6 (which
    /// would indicate a numerical regression, since R-paired bonds are bit-identical
    /// by construction). Pass <paramref name="cache"/> to share the underlying
    /// <see cref="C2HwhmRatio"/> with sibling claims (e.g.
    /// <see cref="C2UniversalShapeDerivation"/>) and avoid a duplicate Q-scan.
    ///
    /// <para><paramref name="throwOnGridEdgeSnap"/> propagates to
    /// <see cref="Item1Derivation.C2HwhmRatio.Build"/>. When using a shared
    /// <paramref name="cache"/>, the cached <see cref="Item1Derivation.C2HwhmRatio"/>
    /// retains its own throw-on-snap setting from when it was first computed.</para></summary>
    /// <seealso cref="Item1Derivation.GridEdgeEscapeException"/>
    public static PerF71OrbitKTable Build(
        CoherenceBlock block,
        WitnessCache? cache = null,
        bool throwOnGridEdgeSnap = false)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"PerF71OrbitKTable applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var orbits = new F71BondOrbitDecomposition(block.N).Orbits;
        var hwhmRatio = cache?.GetOrComputeC2HwhmRatio(block) ?? C2HwhmRatio.Build(block, throwOnGridEdgeSnap: throwOnGridEdgeSnap);
        var bondWitnesses = hwhmRatio.Witnesses;

        var orbitWitnesses = new List<OrbitKWitness>(orbits.Count);
        foreach (var orbit in orbits)
        {
            orbitWitnesses.Add(orbit.IsSelfPaired
                ? FromSingleBond(orbit, bondWitnesses[orbit.BondA])
                : FromMirrorPair(orbit, bondWitnesses[orbit.BondA], bondWitnesses[orbit.BondB!.Value]));
        }
        return new PerF71OrbitKTable(block, orbitWitnesses);
    }

    private static OrbitKWitness FromSingleBond(F71BondOrbit orbit, HwhmRatioWitness w) =>
        new(orbit, QPeak: w.QPeak, HwhmLeft: w.HwhmLeft, KMax: w.KMax);

    private static OrbitKWitness FromMirrorPair(F71BondOrbit orbit, HwhmRatioWitness wA, HwhmRatioWitness wB)
    {
        double maxDev = Math.Max(
            Math.Max(Math.Abs(wA.QPeak - wB.QPeak), Math.Abs(wA.HwhmLeft - wB.HwhmLeft)),
            Math.Abs(wA.KMax - wB.KMax));
        if (maxDev > 1e-6)
            throw new InvalidOperationException(
                $"F71 mirror invariance violated at orbit {{b={orbit.BondA}, b={orbit.BondB}}}: " +
                $"max-Δ = {maxDev:E2} exceeds 1e-6");
        return new OrbitKWitness(
            orbit,
            QPeak: 0.5 * (wA.QPeak + wB.QPeak),
            HwhmLeft: 0.5 * (wA.HwhmLeft + wB.HwhmLeft),
            KMax: 0.5 * (wA.KMax + wB.KMax));
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
