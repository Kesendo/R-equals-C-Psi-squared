using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>CrossFormCertificate (the F127 GF(p) slice engine): (a) Tonelli-Shanks known
/// answers; (b) the CROSS-LANGUAGE known-answer pin: the same fixed sample evaluated by the
/// committed Python transcription (`_cross_form_generic`, cross_triple_orthogonality.py) at the
/// wall's first prime, on-variety zeros and the exact nonzero control value; (c) the slice
/// certificate passes at the witness's parameters with a discriminating control.</summary>
public class CrossFormCertificateTests
{
    private const long P = 1073741833;   // core_grid.gen_primes(17)[0], = 1 mod 4

    [Fact]
    public void Tonelli_KnownAnswers()
    {
        Assert.Equal(4, CrossFormCertificate.PowMod(2, 2, 61));
        long r = CrossFormCertificate.Tonelli(10, 13);          // 10 = 6^2 mod 13
        Assert.True(r == 6 || r == 7);
        Assert.Equal(-1, CrossFormCertificate.Tonelli(5, 13));  // 5 is a non-residue mod 13
        long i = CrossFormCertificate.SqrtMinusOne(P);
        Assert.Equal(715816966, i);                             // Python _ModPField(P).i
        Assert.Equal(P - 1, i * i % P);
    }

    [Fact]
    public void BothRoots_ProductIsOne_MatchesPython()
    {
        var rz = CrossFormCertificate.BothRoots(5, 7, P);
        Assert.NotNull(rz);
        // Python both_roots(5, 7) at this prime
        Assert.Equal(438303469, rz!.Value.R1);
        Assert.Equal(758151704, rz.Value.R2);
        Assert.Equal(1, rz.Value.R1 * rz.Value.R2 % P);
    }

    [Fact]
    public void Evaluate_CrossLanguagePin_OnVarietyZero_ControlNonZero()
    {
        long i = CrossFormCertificate.SqrtMinusOne(P);
        var rz = CrossFormCertificate.BothRoots(5, 7, P)!.Value;
        var rw = CrossFormCertificate.BothRoots(11, 13, P)!.Value;
        foreach (long z3 in new[] { rz.R1, rz.R2 })
            foreach (long w3 in new[] { rw.R1, rw.R2 })
                Assert.Equal(0, CrossFormCertificate.Evaluate(P, i, new[] { 5L, 7, z3, 11, 13, w3 }));
        // the off-variety control, pinned to the exact Python value (w3 = 17)
        Assert.Equal(251518001,
            CrossFormCertificate.Evaluate(P, i, new[] { 5L, 7, rz.R1, 11, 13, 17 }));
    }

    [Fact]
    public void CertifySlice_PassesWithDiscriminatingControl()
    {
        foreach (long p in new[] { 1073741833L, 1073741857L })   // the witness's two wall primes
        {
            var (onVariety, bad, controls, samples) = CrossFormCertificate.CertifySlice(p, 6);
            Assert.Equal(0, bad);
            Assert.Equal(4 * samples, onVariety);
            Assert.True(controls >= (int)(0.9 * samples),
                $"control too weak mod {p}: {controls}/{samples}");
        }
    }
}
