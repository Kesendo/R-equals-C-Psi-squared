using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 surface for the five canonical F87 trichotomy witnesses
/// (XX+YY truly, Heisenberg truly, YZ+ZY soft / EQ-030 Marrakesh-confirmed,
/// XX+XY hard, XY+YX soft) as a single registry-queryable claim. Wraps the
/// existing <see cref="F87CanonicalWitness.StandardSet"/> with a typed parent
/// edge to <see cref="F87TrichotomyClassification"/>.
///
/// <para>The witnesses inside still classify lazily via
/// <see cref="PauliPairTrichotomy.Classify"/>; this wrapper does NOT trigger
/// classification at construction. Consumers that touch <c>Witnesses[i].ActualClass</c>
/// pay the L-build cost at first access only.</para>
///
/// <para>Tier: Tier2Empirical. The Pauli-pair classification rule itself is
/// Tier1Derived (in <see cref="F87TrichotomyClassification"/>), but the named
/// canonical examples + the EQ-030 Marrakesh hardware match are empirical
/// pin-downs — the registered set is the empirical surface of the law.</para>
///
/// <para>Anchor: <c>docs/ANALYTICAL_FORMULAS.md</c> F87 + EQ-030 anchor in
/// <c>review/EMERGING_QUESTIONS.md</c>; per-witness anchors live on each
/// inner <see cref="F87CanonicalWitness.Anchor"/>.</para></summary>
public sealed class F87StandardWitnessSet : Claim
{
    public ChainSystem Chain { get; }
    public IReadOnlyList<F87CanonicalWitness> Witnesses { get; }

    public F87StandardWitnessSet(ChainSystem chain)
        : base("F87 standard canonical witness set (XX+YY, Heisenberg, EQ-030 soft, XX+XY hard, XY+YX soft)",
               Tier.Tier2Empirical,
               "docs/ANALYTICAL_FORMULAS.md F87 + review/EMERGING_QUESTIONS.md (EQ-030)")
    {
        Chain = chain;
        Witnesses = F87CanonicalWitness.StandardSet(chain);
    }

    public override string DisplayName =>
        $"F87 standard witness set (N={Chain.N}, {Witnesses.Count} witnesses)";

    public override string Summary =>
        $"{Witnesses.Count} canonical Pauli-pair witnesses spanning truly/soft/hard ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("chain", summary: $"N={Chain.N}, {Chain.HType}, {Chain.Topology}");
            foreach (var w in Witnesses) yield return w;
        }
    }
}
