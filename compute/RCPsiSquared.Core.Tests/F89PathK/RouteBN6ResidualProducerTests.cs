using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Tests.F89PathK;

[Trait("Category", "ROUTE_B_A2_N6_PRODUCER")]
public sealed class RouteBN6ResidualProducerTests
{
    private static readonly Lazy<RouteBN6ExactPencil> Even = new(() => FoldResultantCertificate.ExportRouteBN6ExactPencil(false));
    private static readonly Lazy<RouteBN6ExactPencil> Odd = new(() => FoldResultantCertificate.ExportRouteBN6ExactPencil(true));

    [Fact]
    public void IndependentParityExportsHaveClearedDegreesAndExactTransport()
    {
        Assert.False(Even.Value.ROdd);
        Assert.True(Odd.Value.ROdd);
        foreach (var pencil in new[] { Even.Value, Odd.Value })
        {
            Assert.Equal(45, pencil.SectorDimension);
            Assert.Equal(32, pencil.ResidualLambdaDegree);
            Assert.Equal(13, pencil.AtFactorInT.Length - 1);
            Assert.Equal(new[] { BigInteger.One }, pencil.ResidualInT[^1]);
            Assert.Equal(new[] { BigInteger.One }, pencil.AtFactorInT[^1]);
        }
        Assert.True(FoldResultantCertificate.RouteBN6TransportMatches(Even.Value.ResidualInT, Odd.Value.ResidualInT));
    }

    [Fact]
    public void OraclePinsQuarterTurnSignAndClearedLambdaOrientation()
    {
        var oracle = F89PathKFdOracle.FdScaled(5);
        Assert.Equal(oracle, Evaluate(Even.Value.ResidualInT, 2 * GaussianInteger.I));
        Assert.False(oracle.SequenceEqual(Evaluate(Even.Value.ResidualInT, -2 * GaussianInteger.I)));
        Assert.False(oracle.SequenceEqual(Evaluate(Even.Value.ResidualInT, 2)));
        Assert.False(oracle.SequenceEqual(Evaluate(Even.Value.ResidualInT, -2)));
    }

    [Fact]
    public void TransportRejectsCoefficientMutation()
    {
        var mutated = Odd.Value.ResidualInT.Select(row => row.ToArray()).ToArray();
        mutated[0][0] += 1;
        Assert.False(FoldResultantCertificate.RouteBN6TransportMatches(Even.Value.ResidualInT, mutated));
    }

    [Fact]
    public void ActualHoppingSlopeMutationFailsExactAtDivision()
    {
        var b0 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(0, 6);
        var b1 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(1, 6);
        var b2 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(2, 6);
        bool changed = false;
        for (int i = 0; i < 45 && !changed; i++)
            for (int j = i + 1; j < 45 && !changed; j++)
            {
                var slope = b1[i, j] - b0[i, j];
                if (slope == GaussianInteger.Zero || slope != b1[j, i] - b0[j, i]) continue;
                b1[i, j] -= 2 * slope;
                b1[j, i] -= 2 * slope;
                b2[i, j] -= 4 * slope;
                b2[j, i] -= 4 * slope;
                changed = true;
            }
        Assert.True(changed);
        var error = Assert.Throws<InvalidOperationException>(() =>
            FoldResultantCertificate.ExportRouteBN6ExactPencilFromSamples(false, b0, b1, b2));
        Assert.Contains("AT division left a remainder", error.Message);
    }

    [Fact]
    public void NonlinearSampleIsRejectedBeforePolynomialConstruction()
    {
        var b0 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(0, 6);
        var b1 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(1, 6);
        var b2 = F89PathKSeDeBlock.BuildTwoTimesSymBlock(2, 6);
        b2[0, 0] += 1;
        var error = Assert.Throws<InvalidOperationException>(() =>
            FoldResultantCertificate.ExportRouteBN6ExactPencilFromSamples(false, b0, b1, b2));
        Assert.Contains("not linear", error.Message);
    }

    [Fact]
    public void N7ControlHasUnequalSectorsAndPublicExportIsN6Only()
    {
        Assert.Equal(75, F89PathKSeDeBlock.BuildTwoTimesSymBlock(0, 7).GetLength(0));
        Assert.Equal(72, F89PathKSeDeBlock.BuildTwoTimesROddBlock(0, 7).GetLength(0));
        var method = Assert.Single(typeof(FoldResultantCertificate).GetMethods(),
            m => m.Name == nameof(FoldResultantCertificate.ExportRouteBN6ExactPencil));
        Assert.Equal(typeof(bool), Assert.Single(method.GetParameters()).ParameterType);
    }

    private static GaussianInteger[] Evaluate(BigInteger[][] polynomial, GaussianInteger t)
        => polynomial.Select(row => row.Reverse().Aggregate(GaussianInteger.Zero,
            (acc, coefficient) => acc * t + coefficient)).ToArray();
}
