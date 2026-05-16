using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F83;
using RCPsiSquared.Diagnostics.F87;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Asserting tests for <see cref="F89F87BreakPredictionFromF83"/>: the
/// closed-form prediction ‖M(H_F89 + δH)‖²_F = 2^(N+2)·(HOdd(δH)² + 2·HEvenNonTruly(δH)²)
/// matches the actual L palindrome residual bit-exactly across the five Pi2Class regimes,
/// and Pitfall 1 detection flags the YZ+ZY family correctly.</summary>
public class F89F87BreakPredictionFromF83Tests
{
    private readonly ITestOutputHelper _out;
    public F89F87BreakPredictionFromF83Tests(ITestOutputHelper output) => _out = output;

    private const double BitExactTolerance = 1e-9;

    private static readonly PauliPairBondTerm Xx = new(PauliLetter.X, PauliLetter.X);
    private static readonly PauliPairBondTerm Yy = new(PauliLetter.Y, PauliLetter.Y);
    private static readonly PauliPairBondTerm Xy = new(PauliLetter.X, PauliLetter.Y);
    private static readonly PauliPairBondTerm Yx = new(PauliLetter.Y, PauliLetter.X);
    private static readonly PauliPairBondTerm Yz = new(PauliLetter.Y, PauliLetter.Z);
    private static readonly PauliPairBondTerm Zy = new(PauliLetter.Z, PauliLetter.Y);

    private static readonly IReadOnlyList<PauliPairBondTerm> F89Only = new[] { Xx, Yy };
    private static readonly IReadOnlyList<PauliPairBondTerm> F89PlusOddPure = new[] { Xx, Yy, Xy, Yx };
    private static readonly IReadOnlyList<PauliPairBondTerm> F89PlusEvenNonTruly = new[] { Xx, Yy, Yz, Zy };
    private static readonly IReadOnlyList<PauliPairBondTerm> PureOddPure = new[] { Xy, Yx };
    private static readonly IReadOnlyList<PauliPairBondTerm> PureEvenNonTruly = new[] { Yz, Zy };

    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    private static double ActualBreakNorm(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        var L = PauliDephasingDissipator.Build(H, gammaList, PauliLetter.Z);
        var M = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, PauliLetter.Z);
        return M.FrobeniusNorm();
    }

    [Fact]
    public void PredictBreakNorm_F89Only_IsZero()
    {
        var chain = MakeChain(N: 4);
        var forecast = PiDecompositionPrediction.Predict(chain, F89Only);
        double predicted = F89F87BreakPredictionFromF83.PredictBreakNorm(chain.N, forecast);
        Assert.Equal(0.0, predicted, precision: 12);
        // Cross-check: actual L palindrome residual is also numerical zero.
        Assert.True(ActualBreakNorm(chain, F89Only) < BitExactTolerance);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void PredictBreakNorm_AcrossAllFiveRegimes_MatchesActual(int N)
    {
        var chain = MakeChain(N);
        var cases = new (string Label, IReadOnlyList<PauliPairBondTerm> Terms)[]
        {
            ("F89 only", F89Only),
            ("F89 + (XY+YX)", F89PlusOddPure),
            ("F89 + (YZ+ZY)", F89PlusEvenNonTruly),
            ("Pure XY+YX", PureOddPure),
            ("Pure YZ+ZY", PureEvenNonTruly),
        };

        foreach (var (label, terms) in cases)
        {
            var forecast = PiDecompositionPrediction.Predict(chain, terms);
            double predicted = F89F87BreakPredictionFromF83.PredictBreakNorm(N, forecast);
            double actual = ActualBreakNorm(chain, terms);

            double residual = Math.Abs(predicted - actual);
            _out.WriteLine($"N={N} {label}: predicted={predicted:G10}, actual={actual:G10}, residual={residual:E2}");

            Assert.True(residual < BitExactTolerance,
                $"N={N} {label}: |predicted − actual| = {residual:E2} exceeds tolerance {BitExactTolerance:E2}");
        }
    }

    [Fact]
    public void WouldAntiFractionAloneMissBreak_OnYzZyFamily_True()
    {
        var chain = MakeChain(N: 4);
        var forecastYz = PiDecompositionPrediction.Predict(chain, PureEvenNonTruly);
        Assert.Equal(0.0, forecastYz.AntiFraction);
        Assert.True(forecastYz.HEvenNonTrulySquared > 0);
        Assert.True(F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak(forecastYz));
    }

    [Fact]
    public void WouldAntiFractionAloneMissBreak_OnPureF89_False()
    {
        var chain = MakeChain(N: 4);
        var forecast = PiDecompositionPrediction.Predict(chain, F89Only);
        Assert.Equal(0.0, forecast.AntiFraction);
        Assert.Equal(0.0, forecast.HEvenNonTrulySquared);
        Assert.False(F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak(forecast));
    }

    [Fact]
    public void WouldAntiFractionAloneMissBreak_OnOddPure_False()
    {
        // Pi2OddPure has HOdd² > 0, AntiFraction = 0.5 > 0, so the "alone" miss test is false.
        var chain = MakeChain(N: 4);
        var forecast = PiDecompositionPrediction.Predict(chain, PureOddPure);
        Assert.True(forecast.AntiFraction > 0);
        Assert.False(F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak(forecast));
    }

    [Fact]
    public void PredictBreakNormSquared_NullForecast_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            F89F87BreakPredictionFromF83.PredictBreakNormSquared(N: 4, perturbationForecast: null!));
    }

    [Fact]
    public void PredictBreakNormSquared_NTooSmall_Throws()
    {
        var forecast = new PiDecompositionForecast(0, 0, 0, 0, 0, 0, 0);
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89F87BreakPredictionFromF83.PredictBreakNormSquared(N: 1, perturbationForecast: forecast));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var trichotomy = new F87TrichotomyClassification();
        var trulyBridge = new F89F87TrulyInheritance(trichotomy);
        var bridge = new F89F87BreakPredictionFromF83(trulyBridge);
        Assert.Equal(Tier.Tier1Derived, bridge.Tier);
    }
}
