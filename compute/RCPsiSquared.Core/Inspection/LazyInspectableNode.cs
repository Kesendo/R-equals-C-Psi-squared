namespace RCPsiSquared.Core.Inspection;

/// <summary>An <see cref="IInspectable"/> whose subtree is produced on demand by a factory,
/// not at construction time. <see cref="DisplayName"/> and <see cref="Summary"/> are cheap and
/// always available; the heavy work hides behind <see cref="Children"/>, which only runs the
/// factory the first time it is enumerated.
///
/// <para>This is the lazy seam the world root sits on: a catalog entry can carry its name and
/// one-line description without ever building its (potentially OOM-prone) root object, so a
/// shallow render that never enumerates this node's children fires no computation.</para>
/// </summary>
public sealed class LazyInspectableNode : IInspectable
{
    private readonly Func<IInspectable> _factory;
    private IInspectable? _built;

    public LazyInspectableNode(string displayName, string summary, Func<IInspectable> factory)
    {
        DisplayName = displayName;
        Summary = summary;
        _factory = factory;
    }

    public string DisplayName { get; }
    public string Summary { get; }

    public IEnumerable<IInspectable> Children
    {
        get { yield return _built ??= _factory(); }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
