using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Semi-live D for the higher paths: the committed F_d(λ,2) oracle literal (reproducible:
/// python simulations/f89_pathk_galois.py emit 4) is verified to be the LIVE block's H_B-mixed
/// factor (it divides the live characteristic polynomial with a coprime degree-AT complement), and
/// its Frobenius cycle types certify S_d. (path-3 is the fully-live anchor; path-4 shown here, the
/// path-5/6 literals exceed Int64 and use the BigInteger path identically.)</summary>
public class F89PathKOracleTests
{
    // path-4 F_18(λ,2): degree 18, monic, real integer coefficients (≈1.4·10^16, fits long).
    private static readonly long[] Path4Re =
    {
        14035467334057984, 23312947804372992, 19130559561793536, 10187231015206912,
        3924428179963904, 1158171077967872, 271111943831552, 51492033363968, 8056754182144,
        1048846436352, 114242698752, 10427209728, 795425024, 50329856, 2605120, 107648, 3412, 76, 1
    };

    [Fact]
    public void Path4_FdOracle_CertifiesS18()
    {
        var im = new long[Path4Re.Length];
        var types = new List<int[]>();
        GaloisGroupCertificate cert = default;
        foreach (int p in SplitPrimesUpTo(1500))
        {
            var ct = OcticGaloisCertificate.CycleType(Path4Re, im, p);
            if (ct is not null) types.Add(ct);
            cert = OcticGaloisCertificate.JordanVerdict(types, 18);
            if (cert.IsFullSymmetric) break;            // early exit once S_18 is certified
        }
        Assert.True(cert.IsFullSymmetric);              // Gal(F_18 / Q(i)(q)) = S_18
    }

    [Fact]
    public void Path4_OracleFd_IsTheLiveBlocksHbMixedFactor()
    {
        // Build the path-4 block LIVE and confirm the oracle F_18 (×2-scaled to match the cleared
        // block) divides its characteristic polynomial exactly, with a coprime degree-8 AT
        // complement — the triple proving F_18 IS the live block's H_B-mixed factor (semi-live D).
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 5);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);   // degree 26, live

        var fdScaled = new GaussianInteger[Path4Re.Length];
        for (int j = 0; j < Path4Re.Length; j++)                        // roots ×2: coeff·2^(18−j)
            fdScaled[j] = new GaussianInteger((BigInteger)Path4Re[j] * BigInteger.Pow(2, Path4Re.Length - 1 - j), 0);

        var at = F89HbMixedIsolation.Isolate(charpoly, fdScaled, expectedDegree: 8);
        Assert.Equal(8, GaussianPolynomial.Degree(at));                 // F_18 + coprime degree-8 AT
    }

    private static IEnumerable<int> SplitPrimesUpTo(int hi) =>
        Enumerable.Range(5, System.Math.Max(0, hi - 5)).Where(n => n % 4 == 1 && IsPrime(n));

    private static bool IsPrime(int n)
    {
        if (n < 2) return false;
        for (int d = 2; (long)d * d <= n; d++) if (n % d == 0) return false;
        return true;
    }
}
