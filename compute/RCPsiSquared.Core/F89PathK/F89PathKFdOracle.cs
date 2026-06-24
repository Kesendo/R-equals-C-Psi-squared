using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Committed F_d(λ,2) oracle literals for path-4/5/6 (lowest-first, monic, over Z[i]),
/// generated reproducibly by `python simulations/f89_pathk_galois.py gen-cs` into
/// <see cref="F89PathKFdLiterals"/>. These are NOT trusted blindly: the semi-live path-k witnesses
/// verify each one divides the LIVE block characteristic polynomial with a coprime degree-AtDegree(k)
/// complement (the validation triple), proving it IS the block's H_B-mixed factor before reading its
/// Galois group. path-3's octic is recomputed fully live (F89Path3OcticBlock), so it is not here.</summary>
public static class F89PathKFdOracle
{
    /// <summary>Whether an F_d(λ,2) literal is wired for path-k (k = 4,5,6).</summary>
    public static bool Has(int k) => k is 4 or 5 or 6;

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
        var (re, im) = k switch
        {
            4 => (F89PathKFdLiterals.Path4Re, F89PathKFdLiterals.Path4Im),
            5 => (F89PathKFdLiterals.Path5Re, F89PathKFdLiterals.Path5Im),
            6 => (F89PathKFdLiterals.Path6Re, F89PathKFdLiterals.Path6Im),
            _ => throw new ArgumentOutOfRangeException(nameof(k)),
        };
        return (re.Select(BigInteger.Parse).ToArray(), im.Select(BigInteger.Parse).ToArray());
    }
}
