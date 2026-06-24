using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Committed F_d(λ,2) oracle literals (lowest-first, monic, over Z[i]), reproducible via
/// `python simulations/f89_pathk_galois.py emit k`. These are NOT trusted blindly: the semi-live
/// path-k witnesses verify each one divides the LIVE block characteristic polynomial with a coprime
/// degree-AtDegree(k) complement (the validation triple), proving it IS the block's H_B-mixed factor
/// before reading its Galois group. path-3's octic is recomputed fully live (F89Path3OcticBlock), so
/// it is not stored here. (path-5/6 literals exceed Int64 and await a committed loader.)</summary>
public static class F89PathKFdOracle
{
    /// <summary>Whether an F_d(λ,2) literal is currently wired for path-k.</summary>
    public static bool Has(int k) => k == 4;

    public static int Degree(int k) => k switch
    {
        4 => 18, 5 => 32, 6 => 53,
        _ => throw new ArgumentOutOfRangeException(nameof(k)),
    };

    public static int SymDim(int k) => k switch
    {
        4 => 26, 5 => 45, 6 => 75,
        _ => throw new ArgumentOutOfRangeException(nameof(k)),
    };

    /// <summary>Degree of the AT-locked complement = symmetric dimension − deg F_d (8 / 13 / 22).</summary>
    public static int AtDegree(int k) => SymDim(k) - Degree(k);

    /// <summary>F_d with roots ×2 (matching the cleared block from <see cref="F89PathKSeDeBlock"/>),
    /// over Z[i].</summary>
    public static GaussianInteger[] FdScaled(int k)
    {
        var (re, im) = Fd(k);
        var scaled = new GaussianInteger[re.Length];
        for (int j = 0; j < re.Length; j++)
        {
            var pow = BigInteger.Pow(2, re.Length - 1 - j);
            scaled[j] = new GaussianInteger(re[j] * pow, im[j] * pow);
        }
        return scaled;
    }

    /// <summary>F_d(λ,2) as (re, im) BigInteger coefficients, lowest-first, monic.</summary>
    public static (BigInteger[] Re, BigInteger[] Im) Fd(int k)
    {
        if (k != 4) throw new NotImplementedException($"path-{k} oracle literal not yet wired.");
        var re = Path4Re.Select(BigInteger.Parse).ToArray();
        return (re, new BigInteger[re.Length]);   // path-4 F_18(λ,2) is real (im = 0)
    }

    // path-4 F_18(λ,2): degree 18, monic, real (≈1.4·10^16). Reproduce: `emit 4`.
    private static readonly string[] Path4Re =
    {
        "14035467334057984", "23312947804372992", "19130559561793536", "10187231015206912",
        "3924428179963904", "1158171077967872", "271111943831552", "51492033363968", "8056754182144",
        "1048846436352", "114242698752", "10427209728", "795425024", "50329856", "2605120", "107648",
        "3412", "76", "1",
    };
}
