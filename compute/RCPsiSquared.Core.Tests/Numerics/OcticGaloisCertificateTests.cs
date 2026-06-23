using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Known-answer gates for the finite-field Frobenius-cycle-type engine that
/// certifies the F89 octic Galois group. Mirrors the gate-first Python verifier
/// (simulations/f89_path3_octic_galois.py): the engine must read S_8 polynomials as
/// "has a 5-cycle" and solvable/imprimitive ones as "no 5-cycle".</summary>
public class OcticGaloisCertificateTests
{
    // ≡1 (mod 4) primes (split in Z[i]); the only ones the engine can reduce over.
    private static readonly int[] SplitPrimes =
        { 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113 };

    [Fact]
    public void F89OcticAtQ2_AtSplitPrime5_FactorsAs_5_2_1()
    {
        // THE CERTIFICATE: 𝔭|5 (F_5, i↦2) factors F_8(·,2) to cycle type (5,2,1).
        // (5,2,1) ⟹ a 5-cycle (square) ⟹ primitive ⟹ ⊇A_8; odd ⟹ ⊄A_8 ⟹ S_8.
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        Assert.Equal(new[] { 5, 2, 1 }, OcticGaloisCertificate.CycleType(re, im, 5));
    }

    [Fact]
    public void F89OcticAtQ2_AtSplitPrime41_IsAnEightCycle()
    {
        // An 8-cycle on all 8 roots certifies the group is TRANSITIVE.
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        Assert.Equal(new[] { 8 }, OcticGaloisCertificate.CycleType(re, im, 41));
    }

    [Fact]
    public void SelmerTrinomial_S8_HasAFiveCycleOverSomeSplitPrime()
    {
        // x^8 - x - 1 has Galois group S_8 over Q — the engine must find a 5-cycle.
        var re = new long[] { -1, -1, 0, 0, 0, 0, 0, 0, 1 };  // lowest-first: −1 −x + x^8
        var im = new long[9];
        bool hasFive = SplitPrimes
            .Select(p => OcticGaloisCertificate.CycleType(re, im, p))
            .Where(ct => ct is not null)
            .Any(ct => ct!.Contains(5));
        Assert.True(hasFive);
    }

    [Fact]
    public void XEighthMinusTwo_Solvable_NeverHasAFiveCycle()
    {
        // x^8 - 2 has the solvable group QD_16 (no element of order 5) — no 5-cycle, ever.
        var re = new long[] { -2, 0, 0, 0, 0, 0, 0, 0, 1 };
        var im = new long[9];
        bool hasFive = SplitPrimes
            .Select(p => OcticGaloisCertificate.CycleType(re, im, p))
            .Where(ct => ct is not null)
            .Any(ct => ct!.Contains(5));
        Assert.False(hasFive);
    }

    [Fact]
    public void NonSplitPrime_ReturnsNull()
    {
        // p ≡ 3 (mod 4) does NOT split in Z[i]: the engine cannot reduce over F_p.
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        Assert.Null(OcticGaloisCertificate.CycleType(re, im, 7));
    }

    [Fact]
    public void CycleType_DegreesSumToPolynomialDegree()
    {
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        foreach (int p in SplitPrimes)
        {
            var ct = OcticGaloisCertificate.CycleType(re, im, p);
            if (ct is not null) Assert.Equal(8, ct.Sum());
        }
    }

    // ---- JordanVerdict: the reusable transitive + ⊇A_d + S_d reading (path-3..6) ----

    [Fact]
    public void JordanVerdict_OcticFiveTwoOne_ReadsAsFullSymmetric()
    {
        // (5,2,1) on 8 points: 5 is prime in (4, 5] ⟹ ⊇A_8; (8 − 3) is odd ⟹ ⊄A_8 ⟹ S_8.
        var cert = OcticGaloisCertificate.JordanVerdict(new[] { new[] { 5, 2, 1 } }, 8);
        Assert.Equal(5, cert.JordanPrime);
        Assert.True(cert.HasOddCycleType);
        Assert.True(cert.ContainsAlternating);
        Assert.True(cert.IsNonSolvable);
        Assert.True(cert.IsFullSymmetric);
    }

    [Fact]
    public void JordanVerdict_OnlyAnEightCycle_IsTransitiveButNotAlternating()
    {
        // A single 8-cycle certifies transitivity, but 8 is not prime ⟹ no Jordan prime.
        var cert = OcticGaloisCertificate.JordanVerdict(new[] { new[] { 8 } }, 8);
        Assert.True(cert.Transitive);
        Assert.Null(cert.JordanPrime);
        Assert.False(cert.ContainsAlternating);
        Assert.False(cert.IsFullSymmetric);
    }

    [Fact]
    public void JordanVerdict_SelmerDegree18_ReadsAsFullSymmetric()
    {
        // x^18 − x − 1 has Gal/Q = S_18 (Osada). Gathering Frobenius cycle types over split
        // primes, the generalised Jordan verdict (a prime cycle in (9,15] + an odd type) reads S_18.
        var re = new long[19]; re[0] = -1; re[1] = -1; re[18] = 1;   // −1 − x + x^18
        var im = new long[19];
        var cycleTypes = SplitPrimesUpTo(600)
            .Select(p => OcticGaloisCertificate.CycleType(re, im, p))
            .Where(ct => ct is not null).Select(ct => ct!).ToList();
        var cert = OcticGaloisCertificate.JordanVerdict(cycleTypes, 18);
        Assert.True(cert.IsFullSymmetric);
    }

    [Fact]
    public void JordanVerdict_XEighteenMinusTwo_IsNotAlternating()
    {
        // x^18 − 2 is solvable (group order 108 = 2²·3³, no element of order 11 or 13),
        // so no prime cycle in (9, 15] ever appears ⟹ the verdict is not ⊇A_18.
        var re = new long[19]; re[0] = -2; re[18] = 1;
        var im = new long[19];
        var cycleTypes = SplitPrimesUpTo(600)
            .Select(p => OcticGaloisCertificate.CycleType(re, im, p))
            .Where(ct => ct is not null).Select(ct => ct!).ToList();
        var cert = OcticGaloisCertificate.JordanVerdict(cycleTypes, 18);
        Assert.Equal(18, cert.Degree);              // verdict is computed (fails the default stub)
        Assert.True(cert.Transitive);               // x^18 − 2 is irreducible ⟹ an 18-cycle appears
        Assert.False(cert.ContainsAlternating);     // ... but no prime cycle in (9, 15]
    }

    // ---- CycleTypeOfFpPoly: DDF directly on an already-reduced F_p polynomial ----

    [Fact]
    public void CycleTypeOfFpPoly_IrreducibleQuadratic_IsASingleTwoCycle()
    {
        // x² + 2 over F_5: −2 = 3 is a non-square mod 5 ⟹ irreducible ⟹ cycle type (2).
        Assert.Equal(new[] { 2 }, OcticGaloisCertificate.CycleTypeOfFpPoly(new[] { 2, 0, 1 }, 5));
    }

    [Fact]
    public void CycleTypeOfFpPoly_SplitQuadratic_IsTwoOneCycles()
    {
        // x² + 1 = (x − 2)(x − 3) over F_5 (2² = −1) ⟹ cycle type (1,1).
        Assert.Equal(new[] { 1, 1 }, OcticGaloisCertificate.CycleTypeOfFpPoly(new[] { 1, 0, 1 }, 5));
    }

    [Fact]
    public void CycleTypeOfFpPoly_NotSquarefree_ReturnsNull()
    {
        // (x − 1)² = x² − 2x + 1 over F_5 = [1, 3, 1]; a repeated root ⟹ Dedekind n/a ⟹ null.
        Assert.Null(OcticGaloisCertificate.CycleTypeOfFpPoly(new[] { 1, 3, 1 }, 5));
    }

    // ---- BigInteger overload: the path-5/6 regime (coefficients exceed Int64) ----

    [Fact]
    public void CycleType_BigInteger_MatchesLongOverload_OnTheOctic()
    {
        // The BigInteger path must read the same Frobenius cycle types as the long[] path.
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        var bigRe = re.Select(x => (BigInteger)x).ToArray();
        var bigIm = im.Select(x => (BigInteger)x).ToArray();
        foreach (int p in SplitPrimes)
            Assert.Equal(OcticGaloisCertificate.CycleType(re, im, p),
                         OcticGaloisCertificate.CycleType(bigRe, bigIm, p));
    }

    [Fact]
    public void CycleType_BigInteger_HandlesCoefficientsBeyondLong()
    {
        // (x − 10^25)(x − (10^25+1)) has coefficients up to ~10^50 (far beyond Int64); mod 5 the
        // roots reduce to 0 and 1 ⟹ cycle type (1,1). The long[] path would overflow.
        var a = BigInteger.Pow(10, 25);
        var b = a + 1;
        var re = new[] { a * b, -(a + b), BigInteger.One };   // x² − (a+b)x + ab, lowest-first
        var im = new[] { BigInteger.Zero, BigInteger.Zero, BigInteger.Zero };
        Assert.Equal(new[] { 1, 1 }, OcticGaloisCertificate.CycleType(re, im, 5));
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
