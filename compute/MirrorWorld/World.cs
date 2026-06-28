namespace MirrorWorld;

// The root. No parent, inference-free. Everything any object holds in its "right" bucket comes from
// here. Its own outputs are the frame read at an object: x, y, z. Nothing more.
public sealed class World : GameObject
{
    public World() : base(null) { }

    public override IReadOnlyList<string> Own => new[] { "x", "y", "z" };
}
