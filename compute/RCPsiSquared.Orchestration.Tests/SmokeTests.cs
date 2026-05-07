using RCPsiSquared.Orchestration.Cli;

namespace RCPsiSquared.Orchestration.Tests;

public class SmokeTests
{
    [Fact]
    public void Orchestration_AssemblyLoads()
    {
        var asm = typeof(KnowledgeCli).Assembly;
        Assert.NotNull(asm);
        Assert.Equal("RCPsiSquared.Orchestration", asm.GetName().Name);
    }
}
