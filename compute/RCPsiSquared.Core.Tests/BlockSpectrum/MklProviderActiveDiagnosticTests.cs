using MathNet.Numerics;
using MathNet.Numerics.Providers.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Diagnostic witness that the RCPsiSquared.Core ModuleInitializer correctly
/// activates MathNet's native MKL provider for the test process. If this test ever fails,
/// the F1 N=7/N=8 block-spectrum tests (and any other Evd-touching test) will silently
/// fall back to the managed single-threaded provider, with 100-400× perf degradation on
/// large blocks. 2026-05-18 incident anchor: see the LiouvillianBlockSpectrum.cs Contract
/// paragraph + the F1 general topology dogfood test family.</summary>
public class MklProviderActiveDiagnosticTests
{
    private readonly ITestOutputHelper _out;
    public MklProviderActiveDiagnosticTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void LinearAlgebraProvider_IsMklAfterCoreLoad()
    {
        // Touch any Core type to ensure the ModuleInitializer has run. The static class
        // itself works without instantiation.
        _ = MathNetSetup.MklActive;

        var provider = LinearAlgebraControl.Provider;
        _out.WriteLine($"Provider: {provider}");
        _out.WriteLine($"MaxDegreeOfParallelism: {Control.MaxDegreeOfParallelism}");
        _out.WriteLine($"MathNetSetup.MklActive: {MathNetSetup.MklActive}");

        var typeName = provider.GetType().Name;
        Assert.True(
            typeName.Contains("Mkl", System.StringComparison.OrdinalIgnoreCase),
            $"Expected an MKL-backed LinearAlgebraProvider, got '{typeName}'. " +
            "ModuleInitializer in CoreModuleInitializer.cs failed to activate MKL, or " +
            "MKL native DLLs are missing from the test output (check runtimes/win-x64/native/).");
    }
}
