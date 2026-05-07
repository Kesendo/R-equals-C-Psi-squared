using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests;

public class SmokeTests
{
    [Fact]
    public void Runtime_AssemblyLoads()
    {
        var asm = typeof(ClaimRegistry).Assembly;
        Assert.NotNull(asm);
        Assert.Equal("RCPsiSquared.Runtime", asm.GetName().Name);
    }
}
