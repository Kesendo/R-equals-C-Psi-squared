using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class EdgeTests
{
    private sealed class A { }
    private sealed class B { }

    [Fact]
    public void Edge_RecordEquality()
    {
        var e1 = new Edge(typeof(A), typeof(B), "reason");
        var e2 = new Edge(typeof(A), typeof(B), "reason");
        Assert.Equal(e1, e2);
    }

    [Fact]
    public void Edge_DifferentReason_NotEqual()
    {
        var e1 = new Edge(typeof(A), typeof(B), "r1");
        var e2 = new Edge(typeof(A), typeof(B), "r2");
        Assert.NotEqual(e1, e2);
    }
}
