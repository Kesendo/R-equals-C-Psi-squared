namespace RCPsiSquared.Orchestration.Tests;

public class SmokeTests
{
    [Fact]
    public void Orchestration_AssemblyLoads()
    {
        var asm = typeof(Placeholder).Assembly;
        Assert.NotNull(asm);
        Assert.Equal("RCPsiSquared.Orchestration", asm.GetName().Name);
    }
}
