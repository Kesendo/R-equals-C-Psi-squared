using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>Abstract base for any typed F86 knowledge object. Carries name, tier, and a
/// pointer to the canonical written anchor (<c>docs/ANALYTICAL_FORMULAS.md</c> entry,
/// <c>docs/proofs/</c> file, or memory-graph entry).
///
/// <para>Each concrete subclass IS its own computation: <see cref="TPeakLaw"/> takes γ₀ and
/// computes 1/(4γ₀); <see cref="TwoLevelEpModel"/> takes (γ₀, J, g_eff) and computes the
/// pre/at/post-EP eigenvalue state; etc. The abstraction here is the "what we know" frame —
/// each claim places itself in tier × anchor space and exposes its computed content as
/// inspectable children.</para>
/// </summary>
public abstract class F86Claim : IInspectable
{
    public string Name { get; }
    public Tier Tier { get; }
    public string Anchor { get; }

    protected F86Claim(string name, Tier tier, string anchor)
    {
        Name = name;
        Tier = tier;
        Anchor = anchor;
    }

    public abstract string DisplayName { get; }
    public abstract string Summary { get; }

    /// <summary>The base prepends tier + anchor metadata before each subclass's
    /// <see cref="ExtraChildren"/>. Override <see cref="ExtraChildren"/> instead of
    /// <see cref="Children"/> to keep the prefix consistent.</summary>
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
