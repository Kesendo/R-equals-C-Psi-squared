namespace MirrorWorld;

// The base. An object produces outputs. Each output is exactly one of two things, nothing more:
//   left  = own       , the object produces it itself.
//   right = inherited  , the object produces it because it inherits from its parent.
public abstract class GameObject
{
    public GameObject? Parent { get; }

    protected GameObject(GameObject? parent) => Parent = parent;

    // left
    public abstract IReadOnlyList<string> Own { get; }

    // right , everything the parents own, walked up the chain.
    public IReadOnlyList<string> Inherited =>
        Parent is null
            ? Array.Empty<string>()
            : Parent.Own.Concat(Parent.Inherited).ToList();
}
