using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Tests.Inspection;

public class LazyInspectableNodeTests
{
    [Fact]
    public void NameAndSummary_AreAvailableWithoutFiringFactory()
    {
        int calls = 0;
        var node = new LazyInspectableNode("lbl", "sum", () => { calls++; return Leaf(); });

        Assert.Equal("lbl", node.DisplayName);
        Assert.Equal("sum", node.Summary);
        Assert.Equal(0, calls);
    }

    [Fact]
    public void Factory_FiresOnlyWhenChildrenEnumerated()
    {
        int calls = 0;
        var node = new LazyInspectableNode("lbl", "sum", () => { calls++; return Leaf(); });

        Assert.Equal(0, calls);
        _ = node.Children.ToList();
        Assert.Equal(1, calls);
    }

    [Fact]
    public void Factory_RunsAtMostOnceAcrossRepeatedEnumeration()
    {
        int calls = 0;
        var node = new LazyInspectableNode("lbl", "sum", () => { calls++; return Leaf(); });

        _ = node.Children.ToList();
        _ = node.Children.ToList();
        Assert.Equal(1, calls);
    }

    [Fact]
    public void DepthOneWalk_DoesNotFireFactory()
    {
        // A walk that stops at the lazy node itself (depth 0 relative to it) never touches
        // its children, so the heavy factory behind it stays cold. This is the property the
        // world root leans on: a shallow render never enumerates a catalog node's children.
        int calls = 0;
        var node = new LazyInspectableNode("lbl", "sum", () => { calls++; return Leaf(); });

        _ = node.Walk(maxDepth: 0).ToList();
        Assert.Equal(0, calls);
    }

    private static IInspectable Leaf() => new InspectableNode("leaf", "computed");
}
