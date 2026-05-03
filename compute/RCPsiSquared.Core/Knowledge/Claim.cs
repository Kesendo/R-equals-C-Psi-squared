using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Abstract base for any typed knowledge claim across the F-theorem family
/// (F86, F71, F1, F87, etc.). Carries a name, a <see cref="Tier"/>, and a pointer to the
/// canonical written anchor (<c>docs/ANALYTICAL_FORMULAS.md</c> entry, <c>docs/proofs/</c>
/// file, or memory-graph entry).
///
/// <para>Each concrete subclass IS its own computation: the F-theorem-specific claim
/// classes derive their values from chain or block primitives lazily. The abstraction here
/// is the "what we know" frame: every claim places itself in tier × anchor space and
/// exposes its computed content as <see cref="IInspectable.Children"/>.</para>
///
/// <para>Subclasses should override <see cref="ExtraChildren"/> to contribute their
/// F-theorem-specific children; the base prepends tier + anchor metadata automatically.</para>
/// </summary>
public abstract class Claim : IInspectable
{
    public string Name { get; }
    public Tier Tier { get; }
    public string Anchor { get; }

    protected Claim(string name, Tier tier, string anchor)
    {
        Name = name;
        Tier = tier;
        Anchor = anchor;
    }

    public abstract string DisplayName { get; }
    public abstract string Summary { get; }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("tier", summary: Tier.Label());
            yield return new InspectableNode("anchor", summary: Anchor);
            foreach (var c in ExtraChildren) yield return c;
        }
    }

    /// <summary>Subclass-specific children, appended after the tier/anchor metadata. Default
    /// is empty.</summary>
    protected virtual IEnumerable<IInspectable> ExtraChildren => Array.Empty<IInspectable>();

    public virtual InspectablePayload Payload => InspectablePayload.Empty;
}
