using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The committed F_d(λ,2) oracle literal (reproducible: python simulations/f89_pathk_galois.py
/// emit 4) must read as S_18 through the C# Frobenius engine — validating the oracle and the
/// generalised Jordan verdict end-to-end for a higher path (path-3 is the live anchor; the live-D
/// isolation for path-4/5/6 will be cross-checked against these same literals).</summary>
public class F89PathKOracleTests
{
    [Fact]
    public void Path4_FdOracle_CertifiesS18()
    {
        // path-4 F_18(λ,2): degree 18, monic, real integer coefficients (≈1.4·10^16, fits long).
        var re = new long[]
        {
            14035467334057984, 23312947804372992, 19130559561793536, 10187231015206912,
            3924428179963904, 1158171077967872, 271111943831552, 51492033363968, 8056754182144,
            1048846436352, 114242698752, 10427209728, 795425024, 50329856, 2605120, 107648, 3412, 76, 1
        };
        var im = new long[re.Length];

        var types = new List<int[]>();
        GaloisGroupCertificate cert = default;
        foreach (int p in SplitPrimesUpTo(1500))
        {
            var ct = OcticGaloisCertificate.CycleType(re, im, p);
            if (ct is not null) types.Add(ct);
            cert = OcticGaloisCertificate.JordanVerdict(types, 18);
            if (cert.IsFullSymmetric) break;            // early exit once S_18 is certified
        }
        Assert.True(cert.IsFullSymmetric);              // Gal(F_18 / Q(i)(q)) = S_18
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
