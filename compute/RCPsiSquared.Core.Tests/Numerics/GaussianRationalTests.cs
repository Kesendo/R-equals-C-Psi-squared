using System;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>ℚ(i), the field <see cref="GaussianInteger"/> deliberately is not: Z[i] stays division-free,
/// and a Taylor recurrence divides by k+1 at every order. The operation to gate is therefore division,
/// and it is gated on inputs that can break it (thirds and elevenths, never a dyadic), against the
/// float route that does break it.</summary>
public class GaussianRationalTests
{
    private static GaussianRational Q(int nRe, int dRe, int nIm, int dIm)
        => new(new BigRational(nRe, dRe), new BigRational(nIm, dIm));

    [Fact]
    public void Multiply_ISquared_IsMinusOne()
    {
        Assert.Equal(new GaussianRational(-1, 0), GaussianRational.I * GaussianRational.I);
    }

    [Fact]
    public void Divide_ThenMultiply_ReturnsTheNumeratorExactly()
    {
        // Non-representable in binary floating point on BOTH components of BOTH operands, so a
        // rounding route cannot pass this by accident: 1/3 + 7/11 i over 13/17 − 5/7 i.
        var a = Q(1, 3, 7, 11);
        var b = Q(13, 17, -5, 7);

        var roundTrip = a / b * b;

        Assert.True(roundTrip == a, $"(a/b)·b = {roundTrip}, expected exactly {a}.");
    }

    [Fact]
    public void Divide_ThenMultiply_HoldsOverASweep_WhereTheFloatRouteDoesNot()
    {
        // The single pair above round-trips in double as well, so on its own it does not show that
        // the exact assertion is measuring anything. This is the control that does: the SAME
        // identity over 400 pseudo-random Gaussian rationals holds exactly EVERY time in Q(i) and
        // fails for a large fraction of them in double. One passing input proves nothing about
        // exactness; the two counts side by side do.
        var rng = new Random(20260905);
        int floatFailures = 0;

        for (int trial = 0; trial < 400; trial++)
        {
            var a = Q(rng.Next(-99, 100), rng.Next(1, 100), rng.Next(-99, 100), rng.Next(1, 100));
            var b = Q(rng.Next(-99, 100), rng.Next(1, 100), rng.Next(-99, 100), rng.Next(1, 100));
            if (b.IsZero) continue;

            Assert.True(a / b * b == a, $"trial {trial}: ({a})/({b})·({b}) != {a}");

            var (aRe, aIm) = ((double)a.Re.Numerator / (double)a.Re.Denominator,
                              (double)a.Im.Numerator / (double)a.Im.Denominator);
            var (bRe, bIm) = ((double)b.Re.Numerator / (double)b.Re.Denominator,
                              (double)b.Im.Numerator / (double)b.Im.Denominator);
            var norm = bRe * bRe + bIm * bIm;
            var (qRe, qIm) = ((aRe * bRe + aIm * bIm) / norm, (aIm * bRe - aRe * bIm) / norm);
            if (qRe * bRe - qIm * bIm != aRe || qRe * bIm + qIm * bRe != aIm) floatFailures++;
        }

        Assert.True(floatFailures > 100,
            $"only {floatFailures} of 400 float round-trips lost the value; the exact gate above " +
            "would then be testing an identity the float route also satisfies, so it would show nothing.");
    }

    [Fact]
    public void Divide_ByZero_Throws()
    {
        Assert.Throws<DivideByZeroException>(() => Q(1, 3, 7, 11) / GaussianRational.Zero);
    }

    [Fact]
    public void NormSquared_IsTheExactRationalModulusSquared()
    {
        // |3/5 + 4/5 i|² = 9/25 + 16/25 = 1, exactly.
        Assert.Equal(BigRational.One, Q(3, 5, 4, 5).NormSquared);
    }

    [Fact]
    public void Conjugate_TimesItself_IsTheNormSquared()
    {
        var z = Q(1, 3, 7, 11);
        Assert.Equal(new GaussianRational(z.NormSquared, BigRational.Zero), z * z.Conjugate);
    }

    [Fact]
    public void LiftFromGaussianInteger_KeepsBothComponents()
    {
        GaussianRational lifted = new GaussianInteger(new BigInteger(-4), new BigInteger(9));
        Assert.Equal(new GaussianRational(new BigRational(-4), new BigRational(9)), lifted);
    }

    [Fact]
    public void DivisionByTheStepIndex_StaysExactOverManyOrders()
    {
        // The Taylor recurrence's own shape: divide by 1, 2, …, 20 in turn, then multiply back by
        // 20!. A field that rounds anywhere in the chain cannot return the seed.
        var seed = Q(1, 3, 7, 11);
        var running = seed;
        for (int k = 1; k <= 20; k++) running /= new GaussianRational(new BigRational(k), BigRational.Zero);

        BigInteger factorial = BigInteger.One;
        for (int k = 1; k <= 20; k++) factorial *= k;
        var restored = running * new GaussianRational(new BigRational(factorial), BigRational.Zero);

        Assert.True(restored == seed, $"after 20 divisions and one multiply back: {restored}, expected {seed}.");
    }
}
