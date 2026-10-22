using RCPsiSquared.Core.F86;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

public class TwoLevelEpModelRegimeTests
{
    private const double Eps = 2.220446049250313e-16;

    [Fact]
    public void FarFromTheEp_IsNeverAtEp_AtAnyScaleOfGamma()
    {
        // an absolute threshold on the discriminant 4γ₀² − J²g² called this AtEp: at γ₀ = 10⁻⁶, J = 0 the discriminant
        // is 4·10⁻¹², below any fixed 10⁻¹⁰
        foreach (double gamma in new[] { 1e-9, 1e-6, 1e-3, 1.0, 1e3, 1e6 })
        {
            Assert.Equal(EpRegime.PreEp, new TwoLevelEpModel(gamma, 0.0, 1.0).Regime);
            Assert.Equal(EpRegime.PostEp, new TwoLevelEpModel(gamma, 4.0 * gamma, 1.0).Regime);
        }
    }

    [Fact]
    public void OnTheEp_BuiltAsAtQ_IsAtEp_AndTheResidualStaysBelowTheThreeRoundingBoundInEveryDecade()
    {
        // the error model: three roundings, |x − 1| ≤ 1.5ε; the ratio to ε is read per decade and must stay below that bound
        var rng = new Random(7);
        for (int e = -9; e <= 9; e++)
        {
            double worst = 0;
            for (int i = 0; i < 2000; i++)
            {
                double g = 0.1 + 9.9 * rng.NextDouble();
                double gamma = (1 + 9 * rng.NextDouble()) * Math.Pow(10, e);
                var m = TwoLevelEpModel.AtQ(gamma, 2.0 / g, g);
                Assert.Equal(EpRegime.AtEp, m.Regime);
                worst = Math.Max(worst, m.RelativeDistanceToEp / Eps);
            }
            Assert.True(worst <= 1.5, $"decade 1e{e}: worst |x − 1|/ε = {worst}, above the three-rounding bound 1.5");
        }
    }

    [Fact]
    public void OneStepOutsideTheBand_IsNotAtEp()
    {
        // x = 1 ± 10⁻¹² is far outside 2ε (≈ 4.4·10⁻¹⁶) and on the correct side
        Assert.Equal(EpRegime.PostEp, new TwoLevelEpModel(1.0, 2.0 * (1 + 1e-12), 1.0).Regime);
        Assert.Equal(EpRegime.PreEp, new TwoLevelEpModel(1.0, 2.0 * (1 - 1e-12), 1.0).Regime);
        // a point passed directly with J·g = 2γ₀ exact lands exactly; an AtQ point need not
        Assert.Equal(0.0, new TwoLevelEpModel(1.0, 2.0, 1.0).RelativeDistanceToEp);
        Assert.NotEqual(0.0, TwoLevelEpModel.AtQ(1.0, 2.0 / 49, 49).RelativeDistanceToEp);
        Assert.Equal(EpRegime.AtEp, TwoLevelEpModel.AtQ(1.0, 2.0 / 49, 49).Regime);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-1.0)]
    [InlineData(double.NaN)]
    public void NonPositiveGamma_IsRejected(double gamma)
    {
        // a negative 2γ₀ would make the relative distance negative and read as AtEp
        Assert.Throws<ArgumentOutOfRangeException>(() => new TwoLevelEpModel(gamma, 1.0, 1.0));
    }
}
