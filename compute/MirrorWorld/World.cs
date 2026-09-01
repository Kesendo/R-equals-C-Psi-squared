namespace MirrorWorld;

// The root. No parent, inference-free. Everything any object holds in its "right" bucket comes from
// here, with one exception since 2026-09-01: Crack hangs on Cyclotomy, whose parent is deliberately
// open, so its right bucket carries the two combs and no frame. The World's own outputs are the frame
// read at an object: x, y, z. Nothing more.
public sealed class World : GameObject
{
    public World() : base(null) { }

    public override IReadOnlyList<string> Own => new[] { "x", "y", "z" };
}
