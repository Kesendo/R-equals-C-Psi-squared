using System.Runtime.CompilerServices;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Module-load initializer for RCPsiSquared.Core. Ensures
/// <see cref="MathNetSetup.EnsureInitialized"/> runs exactly once when the assembly is
/// loaded by any consumer (CLI, Compute, test runners, future hosts). This fixes the
/// 2026-05-18 N=8 perf incident where test processes silently fell back to MathNet's
/// managed single-threaded Evd because no caller in the test assembly invoked the
/// initialization. Idempotent and thread-safe via <see cref="MathNetSetup"/>'s existing
/// <see cref="LazyThreadSafetyMode.ExecutionAndPublication"/> guard.</summary>
internal static class CoreModuleInitializer
{
    // CA2255 (ModuleInitializer in library code) is the intended use here: every consumer of
    // the Core assembly, including xUnit test runners that never explicitly initialize MathNet,
    // must run native-MKL activation on assembly load. There is no application entry point to
    // hook for the test process.
#pragma warning disable CA2255
    [ModuleInitializer]
#pragma warning restore CA2255
    internal static void Initialize() => MathNetSetup.EnsureInitialized();
}
