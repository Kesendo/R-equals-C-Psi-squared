using System.Linq;
using RCPsiSquared.Core.Numerics;
using Xunit;
using K = RCPsiSquared.Core.Numerics.Cyclotomic24;
using KM = RCPsiSquared.Core.Numerics.Cyclotomic24Matrix;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>ℚ(ζ₂₄) on <see cref="CyclotomicRing"/>'s basis: the identities that define its named numbers, conjugation
/// as ζ ↦ ζ⁻¹, and ranks through the regular representation against matrices of known rank.</summary>
public class Cyclotomic24Tests
{
    private static K Q(long n, long d = 1) => K.FromRational(new BigRational(n, d));

    [Fact]
    public void Zeta_HasOrder24_AndHalfTurnMinusOne()
    {
        Assert.Equal(K.One, K.RootPower(24));
        Assert.Equal(-K.One, K.RootPower(12));
        Assert.Equal(K.RootPower(5), K.RootPower(-19));
        Assert.Equal(K.RootPower(7) * K.RootPower(11), K.RootPower(18));
        Assert.All(Enumerable.Range(1, 23), k => Assert.NotEqual(K.One, K.RootPower(k)));
        // ζ multiplied by itself through the fold by Φ₂₄, not read off the table
        var power = K.One;
        for (int k = 1; k <= 24; k++)
        {
            power *= K.RootPower(1);
            if (k == 12) Assert.Equal(-K.One, power);
        }
        Assert.Equal(K.One, power);
    }

    [Fact]
    public void NamedNumbers_SquareToWhatTheyName()
    {
        Assert.Equal(-K.One, K.I * K.I);
        Assert.Equal(Q(2), K.Sqrt2 * K.Sqrt2);
        Assert.Equal(Q(3), K.Sqrt3 * K.Sqrt3);
        Assert.NotEqual(K.Sqrt2, K.Sqrt3);
        // e^{iπ/4} = (1 + i)/√2 and e^{iπ/6} = (√3 + i)/2
        Assert.Equal(K.RootPower(3) * K.Sqrt2, K.One + K.I);
        Assert.Equal(K.RootPower(2) * new BigRational(2), K.Sqrt3 + K.I);
    }

    [Fact]
    public void Conjugation_IsZetaToItsInverse_AndKeepsTheRealOnes()
    {
        Assert.Equal(-K.I, K.I.Conjugate);
        Assert.Equal(K.Sqrt2, K.Sqrt2.Conjugate);
        Assert.Equal(K.Sqrt3, K.Sqrt3.Conjugate);
        Assert.All(Enumerable.Range(0, 24), k =>
        {
            Assert.Equal(K.RootPower(-k), K.RootPower(k).Conjugate);
            Assert.Equal(K.One, K.RootPower(k) * K.RootPower(k).Conjugate);
        });
        var z = K.Sqrt2 * Q(3, 7) + K.I * K.Sqrt3 - K.RootPower(5);
        var w = K.RootPower(9) * Q(-2) + K.Sqrt3;
        Assert.Equal(z.Conjugate * w.Conjugate, (z * w).Conjugate);
    }

    [Fact]
    public void CosAndSin_AtEveryFifteenDegrees()
    {
        Assert.All(Enumerable.Range(0, 24), k =>
        {
            Assert.Equal(K.One, K.Cos(k) * K.Cos(k) + K.Sin(k) * K.Sin(k));
            Assert.Equal(K.RootPower(k), K.Cos(k) + K.I * K.Sin(k));
        });
        Assert.Equal(K.Sqrt3 * new BigRational(1, 2), K.Cos(2));      // cos 30°
        Assert.Equal(K.Sqrt2 * new BigRational(1, 2), K.Sin(3));      // sin 45°
        Assert.Equal(K.One, K.Sin(6));                               // sin 90°
    }

    [Fact]
    public void FromGaussian_AndTheRegularRepresentation()
    {
        Assert.Equal(Q(3, 4) - K.I * new BigRational(2), K.FromGaussian(new GaussianRational(new BigRational(3, 4), -2)));
        // multiplication by i on {1, ζ, …, ζ⁷}: ζʲ ↦ ζʲ⁺⁶, folded by Φ₂₄; the matrix squares to −1 and has the trace of i's
        // conjugates, 0
        var m = K.I.RegularRepresentation();
        var sq = BigRationalLinearAlgebra.Multiply(m, m);
        Assert.All(Enumerable.Range(0, 8), r => Assert.All(Enumerable.Range(0, 8), c =>
            Assert.Equal(r == c ? new BigRational(-1) : BigRational.Zero, sq[r, c])));
        Assert.Equal(BigRational.Zero, Enumerable.Range(0, 8).Aggregate(BigRational.Zero, (s, j) => s + m[j, j]));
        var col = K.I * K.RootPower(3);
        Assert.All(Enumerable.Range(0, 8), r => Assert.Equal(col[r], m[r, 3]));
    }

    [Fact]
    public void TryRational_OnlyOnTheConstants()
    {
        Assert.True((K.Sqrt2 * K.Sqrt2).TryRational(out var two));
        Assert.Equal(new BigRational(2), two);
        Assert.False(K.Sqrt2.TryRational(out _));
        Assert.False(K.I.TryRational(out _));
    }

    [Fact]
    public void Rank_ReadsTheFieldNotTheEntries()
    {
        K[,] M(K a, K b, K c, K d) => new[,] { { a, b }, { c, d } };
        Assert.Equal(1, KM.Rank(M(K.One, K.Sqrt2, K.Sqrt2, Q(2))));
        Assert.Equal(1, KM.Rank(M(K.Sqrt3, Q(3), K.One, K.Sqrt3)));
        Assert.Equal(1, KM.Rank(M(K.One, K.RootPower(1), K.RootPower(1), K.RootPower(2))));
        Assert.Equal(2, KM.Rank(M(K.One, K.Sqrt2, K.Sqrt2, Q(3))));
        Assert.Equal(0, KM.Rank(KM.Zero(3)));
        Assert.Equal(4, KM.Rank(KM.Identity(4)));
        // a rectangular one: two rows, the second √2·i times the first
        var r = new[,] { { K.One, K.Sqrt3, K.RootPower(7) }, { K.Sqrt2 * K.I, K.Sqrt2 * K.I * K.Sqrt3, K.Sqrt2 * K.I * K.RootPower(7) } };
        Assert.Equal(1, KM.Rank(r));
    }

    [Fact]
    public void Matrices_AdjointReversesProducts_KronMultipliesRanks()
    {
        var a = new[,] { { K.One, K.Sqrt2 }, { K.I, K.RootPower(5) } };
        var b = new[,] { { K.Sqrt3, Q(0) }, { K.RootPower(-3), Q(1, 2) } };
        Assert.True(KM.AreEqual(KM.Adjoint(KM.Multiply(a, b)), KM.Multiply(KM.Adjoint(b), KM.Adjoint(a))));
        Assert.Equal(a[1, 0].Conjugate, KM.Adjoint(a)[0, 1]);                       // conjugated, not only transposed
        Assert.NotEqual(KM.Transpose(a)[0, 1], KM.Adjoint(a)[0, 1]);
        // the Kronecker product with A's index the high digit, rectangular factors included
        var row = new[,] { { K.One, K.I } };
        var col = new[,] { { K.Sqrt2 }, { K.Sqrt3 } };
        var rc = KM.Kron(row, col);
        Assert.Equal(2, rc.GetLength(0));
        Assert.Equal(2, rc.GetLength(1));
        Assert.Equal(K.I * K.Sqrt3, rc[1, 1]);
        Assert.Equal(K.Sqrt3, rc[1, 0]);
        Assert.Equal(K.I * K.Sqrt2, rc[0, 1]);
        Assert.Equal(b[0, 0] * a[1, 0], KM.Kron(a, b)[2, 0]);
        var rank1 = new[,] { { K.One, K.Sqrt2 }, { K.Sqrt2, Q(2) } };
        Assert.Equal(2, KM.Rank(KM.Kron(rank1, KM.Identity(2))));
        Assert.Equal(KM.Trace(a) * KM.Trace(b), KM.Trace(KM.Kron(a, b)));
    }
}
