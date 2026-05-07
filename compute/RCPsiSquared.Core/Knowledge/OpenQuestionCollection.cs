using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Generic wrapper Claim that holds a list of <see cref="OpenQuestion"/> items
/// for one F-theorem family. Distinguished from other family collections by a marker type
/// parameter <typeparamref name="TFamilyMarker"/>: each F-family declares its own empty
/// marker class (e.g. <c>F1Marker</c>, <c>F86Marker</c>) and registers
/// <c>OpenQuestionCollection&lt;F1Marker&gt;</c> as a separate Type in the Runtime registry.
///
/// <para>Solves the multi-instance constraint: <see cref="ObjectManager.ClaimRegistryBuilder"/>
/// keys factories by Type, so multiple bare <see cref="OpenQuestion"/> instances cannot
/// all register. The generic marker makes each F-family collection a distinct Type, while
/// the items inside stay typed as <see cref="OpenQuestion"/>.</para>
///
/// <para>Tier is fixed at <see cref="Tier.OpenQuestion"/> because the collection itself
/// represents what is currently unresolved for the family. Individual items may lift to
/// stronger tiers as their resolutions land; the collection's tier reflects the open
/// surface as a whole.</para></summary>
public sealed class OpenQuestionCollection<TFamilyMarker> : Claim
{
    public IReadOnlyList<OpenQuestion> Items { get; }

    public OpenQuestionCollection(string familyName, IReadOnlyList<OpenQuestion> items, string anchor)
        : base($"{familyName} open questions", Tier.OpenQuestion, anchor)
    {
        Items = items;
    }

    public override string DisplayName =>
        $"{Name} ({Items.Count} item{(Items.Count == 1 ? "" : "s")})";

    public override string Summary =>
        Items.Count == 0
            ? $"{Name}: no open questions"
            : $"{Name}: {Items.Count} open question{(Items.Count == 1 ? "" : "s")} pending";

    protected override IEnumerable<IInspectable> ExtraChildren =>
        Items;
}
