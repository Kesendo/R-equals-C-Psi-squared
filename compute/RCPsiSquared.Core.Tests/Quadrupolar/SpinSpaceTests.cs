using System;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The spin space: index 0 is m = I; the gauge makes I₊ ones on the superdiagonal and
/// puts a(m)² on I₋'s subdiagonal, so the commutation relations hold exactly at every I; the float
/// physical operators and tensors are the gauge's image, checked entry by entry through σ.</summary>
public class SpinSpaceTests
{
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(7)] [InlineData(12)]
    public void GaugedOperators_SatisfyTheSpinAlgebraExactly(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var iz = s.GaugedIz();
        var ip = s.GaugedIplus();
        var im = s.GaugedIminus();
        // [I₊, I₋] = 2 I_z and [I_z, I₊] = I₊, exactly
        Assert.True(M.AreEqual(M.Subtract(M.Multiply(ip, im), M.Multiply(im, ip)), M.Scale(iz, 2)));
        Assert.True(M.AreEqual(M.Subtract(M.Multiply(iz, ip), M.Multiply(ip, iz)), ip));
        // the Casimir I_z² + (I₊I₋ + I₋I₊)/2 = I(I + 1)·1
        var casimir = M.Add(M.Multiply(iz, iz),
                            M.Scale(M.Add(M.Multiply(ip, im), M.Multiply(im, ip)), new BigRational(1, 2)));
        Assert.True(M.AreEqual(casimir, M.Scale(M.Identity(s.Levels), s.SpinSquared)));
    }

    [Fact]
    public void LevelsAndRaiseSquares_AtSpinThreeHalves()
    {
        var s = SpinSpace.Of(3);
        Assert.Equal(4, s.Levels);
        Assert.Equal(new BigRational(3, 2), s.M(0));
        Assert.Equal(new BigRational(-3, 2), s.M(3));
        Assert.Equal((BigRational)3, s.RaiseSquared(1));   // ⟨3/2|I₊|1/2⟩² = 15/4 − 3/4
        Assert.Equal((BigRational)4, s.RaiseSquared(2));
        Assert.Equal((BigRational)3, s.RaiseSquared(3));
        Assert.Equal((BigRational)36, s.GaugeWeightSquared(3));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpinSpace.Of(0));
    }

    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(5)] [InlineData(12)]
    public void FloatPhysical_IsTheGaugesImage(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var ipF = s.Iplus();
        var imF = s.Iminus();
        var ipG = s.GaugedIplus();
        var imG = s.GaugedIminus();
        double eps = Math.Pow(2, -52);
        for (int i = 0; i < s.Levels; i++)
            for (int j = 0; j < s.Levels; j++)
            {
                // physical = σ_i⁻¹·gauged·σ_j, σ_k = √GaugeWeightSquared(k)
                double f = Math.Sqrt(SpinSpace.ToDouble(s.GaugeWeightSquared(j)) / SpinSpace.ToDouble(s.GaugeWeightSquared(i)));
                double ip = SpinSpace.ToDouble(ipG[i, j].Re) * f, im = SpinSpace.ToDouble(imG[i, j].Re) * f;
                Assert.True(Math.Abs(ipF[i, j].Real - ip) <= 64 * eps * Math.Max(1, Math.Abs(ip)));
                Assert.True(Math.Abs(imF[i, j].Real - im) <= 64 * eps * Math.Max(1, Math.Abs(im)));
            }
    }

    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(9)]
    public void Tensors_GaugedAreThePhysicalImage_AndTheAdjointRelationHolds(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        double eps = Math.Pow(2, -52);
        foreach (int q in new[] { 1, -1, 2, -2 })
        {
            var g = QuadrupolarTensors.Gauged(s, q);
            var p = QuadrupolarTensors.Physical(s, q);
            var pm = QuadrupolarTensors.Physical(s, -q);
            for (int i = 0; i < s.Levels; i++)
                for (int j = 0; j < s.Levels; j++)
                {
                    double f = Math.Sqrt(SpinSpace.ToDouble(s.GaugeWeightSquared(j)) / SpinSpace.ToDouble(s.GaugeWeightSquared(i)));
                    double expected = SpinSpace.ToDouble(g[i, j].Re) * f;
                    Assert.True(Math.Abs(p[i, j].Real - expected) <= 64 * eps * Math.Max(1, Math.Abs(expected)));
                    // T₂q† = (−1)^q T₂,−q (real entries, so † is the transpose)
                    double sign = (q % 2 == 0) ? 1 : -1;
                    Assert.True(Math.Abs(p[j, i].Real - sign * pm[i, j].Real) <= 64 * eps * Math.Max(1, Math.Abs(p[j, i].Real)));
                }
        }
        // √6·T₂₀ = 3I_z² − I(I + 1), and the float T₂₀ is it over √6
        var a = QuadrupolarTensors.Root6T20(s);
        var t0 = QuadrupolarTensors.Physical(s, 0);
        for (int k = 0; k < s.Levels; k++)
            Assert.True(Math.Abs(t0[k, k].Real - SpinSpace.ToDouble(a[k, k].Re) / Math.Sqrt(6)) <= 64 * eps * 10);
    }

    [Fact]
    public void LiteratureConversion_IsAFactorThree()
    {
        var r = QuadrupolarRates.FromLiterature(3, 6, 9);
        Assert.Equal(new QuadrupolarRates(1, 2, 3), r);
        Assert.Equal((BigRational)2, r.ByOrder(-1));
        Assert.Equal((BigRational)3, r.ByOrder(2));
    }
}
