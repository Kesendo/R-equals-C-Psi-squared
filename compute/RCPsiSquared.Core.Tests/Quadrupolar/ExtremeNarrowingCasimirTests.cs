using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>At extreme narrowing the dissipator is a polynomial in the rotations' Casimir C = Σ_a ad(I_a)²,
/// D = −(J/4)·C·(4I(I + 1) − 1 − C), here as an exact matrix identity at every 2I from 2 to 12 (at 2I = 1 both
/// sides vanish; the sodium page proves it for every I). C is k(k + 1) on rank k, each rank once (its eigenspaces
/// counted exactly), so R_k = (J/4)·k(k + 1)·(4I(I + 1) − k(k + 1) − 1); two ranks share a rate at the half-integer
/// spins listed. Off extreme narrowing the dissipator keeps the field's axis without its direction: it commutes
/// with ad(I_z) and with the π turn about y, and with ad(I₊) for no rate triple off the line J₀ = J₁ = J₂ (the
/// commutator is linear in the rates; its kernel is read exactly), so there it is no polynomial in C.</summary>
public class ExtremeNarrowingCasimirTests
{
    private static GaussianRational[,] Dissipator(SpinSpace s, long j0, long j1, long j2) =>
        new QuadrupolarGenerator(s, new QuadrupolarRates(j0, j1, j2), SecularHamiltonian.None).GaugedDissipator();

    private static GaussianRational[,] Ad(GaussianRational[,] a)
    {
        var one = M.Identity(a.GetLength(0));
        return M.Subtract(M.Sandwich(a, one), M.Sandwich(one, a));
    }

    private static bool Commute(GaussianRational[,] a, GaussianRational[,] b) =>
        M.AreEqual(M.Multiply(a, b), M.Multiply(b, a));

    private static BigInteger Factorial(int n)
    {
        BigInteger f = BigInteger.One;
        for (int i = 2; i <= n; i++) f *= i;
        return f;
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    [InlineData(8)] [InlineData(9)] [InlineData(10)] [InlineData(11)] [InlineData(12)]
    public void Dissipator_IsAQuadraticInTheRotationCasimir(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        Assert.True(M.AreEqual(NumericalRangeBound.ExtremeNarrowingFromCasimir(s, 4), Dissipator(s, 4, 4, 4)));
    }

    /// <summary>C − k(k + 1)·1 has nullity 2k + 1 for every k = 0 … 2I; since Σ_k (2k + 1) = d², C is diagonalizable
    /// with exactly these eigenspaces, the ranks.</summary>
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(7)] [InlineData(12)]
    public void Casimir_IsKTimesKPlusOne_OnEachRank_EachRankOnce(int twoI)
    {
        Assert.Equal(Enumerable.Range(0, twoI + 1).Select(k => 2 * k + 1), SpinSpace.Of(twoI).RotationCasimirNullities());
    }

    /// <summary>Two ranks share a rate exactly when k(k + 1) + k′(k′ + 1) = 4I(I + 1) − 1: from R_k, among 2I = 2 … 30 at
    /// 2I = 3, 7, 9, 13, 15, 19, 21 and 27 (never at an integer spin), the pairs listed.</summary>
    [Fact]
    public void SharedRankRates_UpTo2IThirty()
    {
        var found = new System.Collections.Generic.List<string>();
        for (int twoI = 2; twoI <= 30; twoI++)
        {
            var s = SpinSpace.Of(twoI);
            for (int k = 1; k <= twoI; k++)
                for (int k2 = k + 1; k2 <= twoI; k2++)
                    if (NumericalRangeBound.ExtremeNarrowingRankRate(s, k, 4) == NumericalRangeBound.ExtremeNarrowingRankRate(s, k2, 4))
                        found.Add($"{twoI}:{k},{k2}");
        }
        Assert.Equal(new[] { "3:1,3", "7:2,7", "7:4,6", "9:6,7", "13:3,13", "15:8,13", "19:7,18", "21:4,21", "21:14,16", "27:11,25", "27:13,24" },
                     found);
    }

    /// <summary>At extreme narrowing D commutes with ad(I_z) and ad(I₊); being self-adjoint in the physical basis
    /// (its jumps come in adjoint pairs at equal rates), it then commutes with ad(I₋) too, so with every rotation. Off it: the
    /// dissipator is linear in the rates, D = J₀·D₀ + J₁·D₁ + J₂·D₂, and the three commutators [D_i, ad(I₊)] have a
    /// one-dimensional kernel in (J₀, J₁, J₂), the line (1, 1, 1): off extreme narrowing D commutes with ad(I₊) for no
    /// rate triple. Every D_i commutes with ad(I_z) (the coherence order) and with the π turn about y, which reverses
    /// the field (it carries ad(I_z) to −ad(I_z) and ad(I₊) to −ad(I₋)); in the gauge that turn is rational,
    /// σ_{2I−k}/σ_k = (2I − k)!/k!.</summary>
    [Theory]
    [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(7)] [InlineData(12)]
    public void ExtremeNarrowing_KeepsEveryRotation_OffIt_TheFieldsAxisWithoutItsDirection(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        int d = s.Levels;
        var adZ = Ad(s.GaugedIz());
        var adPlus = Ad(s.GaugedIplus());
        var adMinus = Ad(s.GaugedIminus());
        var en = Dissipator(s, 4, 4, 4);
        Assert.True(Commute(en, adZ));
        Assert.True(Commute(en, adPlus));

        // the π turn about y, |m⟩ ↦ (−1)^(I − m)|−m⟩, and its inverse, in the gauge
        var r = M.Zero(d);
        var rInv = M.Zero(d);
        for (int k = 0; k < d; k++)
        {
            int sign = k % 2 == 0 ? 1 : -1;
            r[d - 1 - k, k] = new BigRational(sign * Factorial(twoI - k), Factorial(k));
            rInv[k, d - 1 - k] = new BigRational(sign * Factorial(k), Factorial(twoI - k));
        }
        Assert.True(M.AreEqual(M.Multiply(r, rInv), M.Identity(d)));
        var flip = M.Sandwich(r, rInv);
        Assert.True(M.AreEqual(M.Multiply(flip, adZ), M.Scale(M.Multiply(adZ, flip), -1)));
        Assert.True(M.AreEqual(M.Multiply(flip, adPlus), M.Scale(M.Multiply(adMinus, flip), -1)));

        var parts = new[] { Dissipator(s, 1, 0, 0), Dissipator(s, 0, 1, 0), Dissipator(s, 0, 0, 1) };
        foreach (var part in parts)
        {
            Assert.True(Commute(part, adZ));
            Assert.True(Commute(part, flip));
        }
        // the kernel of (J₀, J₁, J₂) ↦ Σ J_i·[D_i, ad(I₊)], the commutators flattened into three columns
        var commutators = parts.Select(p => M.Subtract(M.Multiply(p, adPlus), M.Multiply(adPlus, p))).ToArray();
        int n = d * d;
        var columns = new BigRational[n * n, 3];
        for (int i = 0; i < 3; i++)
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                {
                    Assert.True(commutators[i][a, b].Im.IsZero);
                    columns[n * a + b, i] = commutators[i][a, b].Re;
                }
        var kernel = BigRationalLinearAlgebra.Nullspace(columns);
        Assert.Single(kernel);
        Assert.Equal(kernel[0][0], kernel[0][1]);
        Assert.Equal(kernel[0][0], kernel[0][2]);
        Assert.False(kernel[0][0].IsZero);
    }
}
