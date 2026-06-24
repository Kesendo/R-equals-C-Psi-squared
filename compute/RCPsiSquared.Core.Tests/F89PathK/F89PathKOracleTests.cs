using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Semi-live D for the higher paths (k = 4,5,6): the committed F_d(λ,2) oracle literal
/// (reproducible: python simulations/f89_pathk_galois.py gen-cs) is verified to be the LIVE block's
/// H_B-mixed factor — it divides the live characteristic polynomial with a coprime degree-AT
/// complement (the triple) — and its Frobenius cycle types certify S_d. (path-3 is the fully-live
/// anchor, tested via the block builder + isolation.)</summary>
public class F89PathKOracleTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void OracleFd_IsTheLiveBlocksHbMixedFactor(int k)
    {
        // Build the path-k block LIVE; the oracle F_d (×2-scaled to match the cleared block) must
        // divide its Berkowitz charpoly exactly, with a coprime degree-AtDegree(k) AT complement.
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: k + 1);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);
        var at = F89HbMixedIsolation.Isolate(charpoly, F89PathKFdOracle.FdScaled(k), F89PathKFdOracle.AtDegree(k));
        Assert.Equal(F89PathKFdOracle.AtDegree(k), GaussianPolynomial.Degree(at));
    }

    [Theory]
    [InlineData(4, 1500)]
    [InlineData(5, 2500)]
    [InlineData(6, 4000)]
    public void OracleFd_CertifiesFullSymmetric(int k, int primeHi)
    {
        var (re, im) = F89PathKFdOracle.Fd(k);
        int d = F89PathKFdOracle.Degree(k);
        var types = new List<int[]>();
        GaloisGroupCertificate cert = default;
        foreach (int p in SplitPrimesUpTo(primeHi))
        {
            var ct = OcticGaloisCertificate.CycleType(re, im, p);
            if (ct is null) continue;
            types.Add(ct);
            cert = OcticGaloisCertificate.JordanVerdict(types, d);
            if (cert.IsFullSymmetric) break;            // early exit once S_d is certified
        }
        Assert.True(cert.IsFullSymmetric);              // Gal(F_d / Q(i)(q)) = S_d
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
