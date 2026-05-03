using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>A substantive theoretical item that has not yet been resolved: the OOP form
/// of "what's missing for full Tier 1 promotion". Encoded as a typed <see cref="Claim"/> at
/// <see cref="Tier.OpenQuestion"/> so the OM tree surfaces them automatically.
///
/// <para>Generic across F-theorems: each F-family knowledge base defines its own canonical
/// list (e.g. <c>F86OpenQuestions.Standard</c>, <c>F71OpenItems.Standard</c>) anchoring at
/// the relevant proof / docs entry.</para>
/// </summary>
public sealed class OpenQuestion : Claim
{
    public string Description { get; }
    public string Approach { get; }

    public OpenQuestion(string name, string description, string approach, string anchor)
        : base(name, Tier.OpenQuestion, anchor)
    {
        Description = description;
        Approach = approach;
    }

    public override string DisplayName => $"[OPEN] {Name}";
    public override string Summary => Description;

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("description", summary: Description);
            yield return new InspectableNode("approach", summary: Approach);
        }
    }
}
