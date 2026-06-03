using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class SlowManifoldPauliContentTests
{
    private readonly ITestOutputHelper _out;
    public SlowManifoldPauliContentTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void CrossoverCore_Is_IZ_On_Lit_Sites_Rotating_Carries_XY()
    {
        // The fan resolved the N = 3 crossover slow manifold into a 4-dimensional invariant core plus
        // a rotating remainder. The reading to nail: the core is {I, Z} on the LIT sites (0, 1, the
        // X/Y carriers the turn rotates), and the rotating directions carry an X or a Y on a lit site
        // (PTF's far bank). The Pauli projection confirms it: the core's mass sits on {I, Z}-on-lit
        // strings (lit-XY weight ≈ 0), the rotating mass carries lit-XY content.
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);

        var reading = SlowManifoldPauliContent.Compute(
            sweep.SlowBasis[0], sweep.SlowBasis[^1], axis.N, axis.LitSites);

        _out.WriteLine($"core dim {reading.CoreDim}, lit-XY weight {reading.CoreLitXYWeight:E3}");
        foreach (var s in reading.CoreTop)
            _out.WriteLine($"  core   {s.Label}  {s.Weight:0.000}{(s.LitXY ? "  [lit-XY]" : "")}");
        _out.WriteLine($"rotating dim {reading.RotatingDim}, lit-XY weight {reading.RotatingLitXYWeight:E3}");
        foreach (var s in reading.RotatingTop)
            _out.WriteLine($"  rot    {s.Label}  {s.Weight:0.000}{(s.LitXY ? "  [lit-XY]" : "")}");

        Assert.Equal(4, reading.CoreDim);
        Assert.Equal(12, reading.RotatingDim);

        // The core is {I, Z} on the lit sites: essentially no XY-on-lit content (to the Evd floor).
        Assert.True(reading.CoreLitXYWeight < 1e-6,
            $"core should be {{I,Z}} on lit sites; lit-XY weight {reading.CoreLitXYWeight:E3}");
        // The rotating part carries XY on the lit sites: a substantial fraction.
        Assert.True(reading.RotatingLitXYWeight > 0.1,
            $"rotating part should carry XY on lit sites; lit-XY weight {reading.RotatingLitXYWeight:E3}");
    }
}
