using System;
using System.Linq;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The generator: trace preserved and coherence order conserved exactly, the mean rate of the
/// claim at 3/2, the population rates of the sodium page, one rate per rank at extreme narrowing at every
/// I; the physical view exact where rational (D symmetric, the Hamiltonian part skew), refused at spin 2;
/// the float view the gauged one taken back through σ, within 64·eps of the norm.</summary>
public class QuadrupolarGeneratorTests
{
    private static QuadrupolarGenerator Gen(int twoI, long j0, long j1, long j2, long q = 0, long w = 0) =>
        new(SpinSpace.Of(twoI), new QuadrupolarRates(j0, j1, j2), SecularHamiltonian.Of(q, w));

    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(7)] [InlineData(12)]
    public void TracePreserved_AndCoherenceOrderConserved_Exactly(int twoI)
    {
        var g = Gen(twoI, 1, 2, 3, 5, 7);
        var l = g.Gauged();
        int d = g.Space.Levels;
        for (int c = 0; c < g.Dim; c++)
        {
            var sum = GaussianRational.Zero;
            for (int k = 0; k < d; k++) sum += l[d * k + k, c];
            Assert.True(sum.IsZero);
        }
        for (int r = 0; r < g.Dim; r++)
            for (int c = 0; c < g.Dim; c++)
                if (QuadrupolarGenerator.OrderOf(d, r) != QuadrupolarGenerator.OrderOf(d, c))
                    Assert.True(l[r, c].IsZero);
    }

    [Theory]
    [InlineData(1, 2, 3)] [InlineData(0, 2, 3)] [InlineData(4, 0, 1)] [InlineData(1, 1, 1)]
    public void MeanRate_AtThreeHalves_IsTheClaimsGlobalCentre(long j0, long j1, long j2)
    {
        var g = Gen(3, j0, j1, j2, 2, 3);
        Assert.Equal((BigRational)Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(j0, j1, j2), 2 * g.MeanRate());
    }

    /// <summary>The spin-1 rates the claim's cone is read from, block by block at J = (5, 2, 3): populations 0, 3J₁/2,
    /// (J₁ + 4J₂)/2; single quantum (3J₀ + J₁ + 2J₂)/4 and (3J₀ + 5J₁ + 2J₂)/4; double quantum (J₁ + 2J₂)/2.</summary>
    [Fact]
    public void SpinOneRates_BlockByBlock()
    {
        var g = Gen(2, 5, 2, 3);
        Assert.Equal(new BigRational[] { 0, 21, 10, 1 }, BigRationalMatrixCharpoly.Characteristic(g.DissipatorBlock(0)));          // x(x + 3)(x + 7)
        Assert.Equal(new BigRational[] { new BigRational(713, 16), new BigRational(27, 2), 1 }, BigRationalMatrixCharpoly.Characteristic(g.DissipatorBlock(1)));  // (x + 23/4)(x + 31/4)
        Assert.Equal(new BigRational[] { 4, 1 }, BigRationalMatrixCharpoly.Characteristic(g.DissipatorBlock(2)));                  // x + 4
    }

    [Fact]
    public void PopulationBlock_AtThreeHalves_HasRates0_6J1_6J2_6J1Plus6J2()
    {
        var g = Gen(3, 5, 2, 3);
        var block = g.DissipatorBlock(0);
        var charpoly = BigRationalMatrixCharpoly.Characteristic(block);   // lowest first, monic
        // x(x + 12)(x + 18)(x + 30) = x⁴ + 60x³ + 1116x² + 6480x
        Assert.Equal(new BigRational[] { 0, 6480, 1116, 60, 1 }, charpoly);
    }

    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)] [InlineData(8)]
    [InlineData(9)] [InlineData(10)] [InlineData(11)] [InlineData(12)]
    public void ExtremeNarrowing_OneRatePerRank_InEveryOrderBlock(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var g = new QuadrupolarGenerator(s, new QuadrupolarRates(4, 4, 4), SecularHamiltonian.None);
        var rates = Enumerable.Range(1, twoI).Select(k => NumericalRangeBound.ExtremeNarrowingRankRate(s, k, 4)).Distinct().ToArray();
        // two ranks share a rate exactly when k(k + 1) + k′(k′ + 1) = 4I(I + 1) − 1: among 2I = 2 … 12 only at 3/2, 7/2 and 9/2
        Assert.Equal(twoI is 3 or 7 or 9, rates.Length < twoI);
        BigRational R(int k) => NumericalRangeBound.ExtremeNarrowingRankRate(s, k, 4);
        if (twoI == 3) Assert.Equal(R(1), R(3));
        if (twoI == 7) { Assert.Equal(R(2), R(7)); Assert.Equal(R(4), R(6)); }
        if (twoI == 9) Assert.Equal(R(6), R(7));
        for (int p = -twoI; p <= twoI; p++)
        {
            var b = g.DissipatorBlock(p);
            int n = b.GetLength(0);
            // B·Π_k (B + R_k·1) = 0: every eigenvalue of the block is 0 or −R_k for a rank k
            var prod = b;
            foreach (var r in rates)
            {
                var shifted = (BigRational[,])b.Clone();
                for (int i = 0; i < n; i++) shifted[i, i] += r;
                prod = BigRationalLinearAlgebra.Multiply(prod, shifted);
            }
            foreach (var x in prod) Assert.True(x.IsZero);
        }
        // every rank's rate is present: the populations block holds T_k0 for every k = 0 … 2I, so W + R_k·1 is singular
        var w0 = g.DissipatorBlock(0);
        for (int k = 1; k <= twoI; k++)
        {
            var shifted = (BigRational[,])w0.Clone();
            var rk = NumericalRangeBound.ExtremeNarrowingRankRate(s, k, 4);
            for (int i = 0; i < w0.GetLength(0); i++) shifted[i, i] += rk;
            Assert.NotEmpty(BigRationalLinearAlgebra.Nullspace(shifted));
        }
        // the trace: −Σ_k (2k + 1)·R_k
        var trace = BigRational.Zero;
        var dg = g.GaugedDissipator();
        for (int i = 0; i < g.Dim; i++) trace += dg[i, i].Re;
        var expected = BigRational.Zero;
        for (int k = 1; k <= twoI; k++) expected += (2 * k + 1) * NumericalRangeBound.ExtremeNarrowingRankRate(s, k, 4);
        Assert.Equal(-expected, trace);
    }

    [Theory]
    [InlineData(2)] [InlineData(3)]
    public void Physical_IsExactWhereRational_DissipatorSymmetric_HamiltonianPartSkew(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var g = new QuadrupolarGenerator(s, new QuadrupolarRates(1, 2, 3), SecularHamiltonian.Of(5, 7).WithLambShift(1, 3));
        var gd = new QuadrupolarGenerator(s, new QuadrupolarRates(1, 2, 3), SecularHamiltonian.None);
        var dPhys = gd.Physical();
        Assert.True(M.AreEqual(dPhys, M.Transpose(dPhys)));                 // D real symmetric = self-adjoint
        var lh = M.Subtract(g.Physical(), dPhys);
        Assert.True(M.AreEqual(M.Adjoint(lh), M.Scale(lh, -1)));             // −i[H, ·] skew-adjoint
    }

    [Fact]
    public void Physical_RefusesSpinTwo_WhereTheEntriesAreIrrational()
    {
        var g = Gen(4, 1, 2, 3);
        Assert.Throws<NotSupportedException>(() => g.Physical());
        var lock32 = new QuadrupolarGenerator(SpinSpace.Of(3), new QuadrupolarRates(1, 2, 3), SecularHamiltonian.Of(1, 1).WithLock(2));
        Assert.Throws<NotSupportedException>(() => lock32.Physical());     // ν·I_x carries √3 at 3/2
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)]
    public void Physical_IsTheGaugedViewItself_WhereItIsRational(int twoI)
    {
        foreach (var (j, q, w, k1, k2) in new[] { ((1L, 2L, 3L), 5L, 7L, 1L, 3L), ((0L, 4L, 1L), -2L, 3L, 0L, 0L), ((2L, 0L, 5L), 1L, 0L, -1L, 2L) })
        {
            var g = new QuadrupolarGenerator(SpinSpace.Of(twoI), new QuadrupolarRates(j.Item1, j.Item2, j.Item3),
                                             SecularHamiltonian.Of(q, w).WithLambShift(k1, k2));
            Assert.True(M.AreEqual(g.Physical(), g.Gauged()));
        }
    }

    /// <summary>Error model: each float entry is a few flops of terms no larger than max|L|, so the deviation
    /// from the exact gauged view taken back through the shared factor <see cref="QuadrupolarGenerator.GaugeFactorSquared"/>
    /// stays within 64·eps·max|L|; the three settings put the Hamiltonian and the dissipator three decades apart
    /// in both directions, so entries of very different size are read in one matrix. From 2I = 4 on the factor
    /// is not 1, so this is also the gate on its direction.</summary>
    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(6)] [InlineData(9)] [InlineData(12)]
    public void Float_IsTheGaugedViewBackTransformed_Within64EpsOfTheNorm(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        foreach (var (rs, hs) in new[] { (1L, 1L), (1L, 1000L), (1000L, 1L) })
        {
            var h = SecularHamiltonian.Of(5 * hs, 7 * hs).WithLambShift(hs, 3 * hs).WithLock(2 * hs);
            var g = new QuadrupolarGenerator(s, new QuadrupolarRates(rs, 2 * rs, 3 * rs), h);
            var lg = g.Gauged();
            var lf = g.Float();
            double norm = 0, dev = 0;
            for (int r = 0; r < g.Dim; r++)
                for (int c = 0; c < g.Dim; c++)
                {
                    var expected = System.Numerics.Complex.Zero;
                    if (!lg[r, c].IsZero)
                    {
                        double f = Math.Sqrt(SpinSpace.ToDouble(QuadrupolarGenerator.GaugeFactorSquared(s, r, c)));
                        expected = new System.Numerics.Complex(SpinSpace.ToDouble(lg[r, c].Re) * f, SpinSpace.ToDouble(lg[r, c].Im) * f);
                    }
                    norm = Math.Max(norm, expected.Magnitude);
                    dev = Math.Max(dev, (lf[r, c] - expected).Magnitude);
                }
            Assert.True(dev <= 64 * Math.Pow(2, -52) * norm, $"2I = {twoI}, rates ×{rs}, H ×{hs}: deviation {dev:e3} against norm {norm:e3}");
        }
    }
}
