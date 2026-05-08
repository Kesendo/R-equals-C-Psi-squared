using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class XyJordanWignerModesTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void SineModes_AreOrthonormal(int N)
    {
        var modes = XyJordanWignerModes.Build(N);
        var M = modes.SineModeMatrix;
        var product = M * M.Transpose();
        var identity = Matrix<double>.Build.DenseIdentity(N);
        var residual = (product - identity).FrobeniusNorm();
        Assert.True(residual < 1e-12,
            $"sine-mode matrix must be row-orthonormal; got Frobenius residual {residual:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    [InlineData(10)]
    public void Dispersion_MatchesXyEigenvalues_AtJEqualsOne(int N)
    {
        var modes = XyJordanWignerModes.Build(N, J: 1.0);
        var hop = Matrix<double>.Build.Dense(N, N);
        for (int i = 0; i < N - 1; i++) { hop[i, i + 1] = 1.0; hop[i + 1, i] = 1.0; }
        var eigs = hop.Evd().EigenValues.Select(z => z.Real).OrderBy(x => x).ToArray();
        var theirs = modes.Dispersion.OrderBy(x => x).ToArray();
        for (int i = 0; i < N; i++)
            Assert.Equal(eigs[i], theirs[i], precision: 10);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void Dispersion_ScalesLinearlyWithJ(int N)
    {
        var modesA = XyJordanWignerModes.Build(N, J: 1.0);
        var modesB = XyJordanWignerModes.Build(N, J: 2.5);
        for (int k = 0; k < N; k++)
            Assert.Equal(2.5 * modesA.Dispersion[k], modesB.Dispersion[k], precision: 12);
    }

    [Theory]
    [InlineData(5, 1, 0)]
    [InlineData(5, 3, 2)]
    [InlineData(8, 5, 4)]
    public void SineMode_MatchesAnalyticalFormula(int N, int k, int j)
    {
        var modes = XyJordanWignerModes.Build(N);
        double expected = Math.Sqrt(2.0 / (N + 1)) * Math.Sin(Math.PI * k * (j + 1) / (N + 1));
        Assert.Equal(expected, modes.SineMode(k, j), precision: 12);
    }

    [Fact]
    public void Build_RejectsZeroOrNegativeN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => XyJordanWignerModes.Build(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => XyJordanWignerModes.Build(-1));
    }

    [Fact]
    public void SineMode_RejectsOutOfRangeIndices()
    {
        var modes = XyJordanWignerModes.Build(5);
        Assert.Throws<ArgumentOutOfRangeException>(() => modes.SineMode(0, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => modes.SineMode(6, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => modes.SineMode(1, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => modes.SineMode(1, 5));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var modes = XyJordanWignerModes.Build(5);
        Assert.Equal(Tier.Tier1Derived, modes.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK()
    {
        var modes = XyJordanWignerModes.Build(5);
        Assert.Contains("PROOF_F86_QPEAK", modes.Anchor);
    }
}
