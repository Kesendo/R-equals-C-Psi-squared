using System.Linq;
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
}
