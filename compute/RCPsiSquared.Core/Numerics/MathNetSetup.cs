using MathNet.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>One-time MathNet.Numerics initialisation: native MKL provider + threading.
///
/// Without this, MathNet falls back to managed implementations that are single-threaded
/// for matrix operations — eigendecomposition and inverse stay on a single CPU core.
/// With MKL active, all linear algebra calls saturate the available cores.
/// </summary>
public static class MathNetSetup
{
    private static readonly Lazy<bool> _init = new(() =>
    {
        bool mkl;
        try
        {
            Control.UseNativeMKL();
            mkl = true;
        }
        catch
        {
            mkl = false;
        }
        Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
        return mkl;
    }, LazyThreadSafetyMode.ExecutionAndPublication);

    public static bool MklActive => _init.Value;

    /// <summary>Idempotent and thread-safe. Returns true if MKL was successfully activated.</summary>
    public static bool EnsureInitialized() => _init.Value;
}
