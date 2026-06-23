using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>(D) the felt_time closure functional, graduated to a live witness
/// (<c>inspect --root gradient</c>): the C# landing of <c>simulations/felt_time_amplitude_law.py</c>.
///
/// <para>The survivor's first-order bond rate shift dRe(b) = (density-mode gradient at bond b)^2 -- the
/// diffusion Rayleigh quotient ("amplitude^2"). The slow survivor is a DENSITY/diffusion mode; a delta-J
/// bond defect perturbs the local diffusion coefficient, so dRe(b) = dlambda/dD_b ~ (n(j)-n(j+1))^2.
/// Gate-first (2026-06-19): dRe/grad^2 is bond-INDEPENDENT and the log-log slope dRe vs |grad| ~ 2.</para></summary>
public class SurvivorDiffusionGradientWitnessTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void RateShift_is_amplitude_squared_in_the_density_gradient(int n)
    {
        var w = new SurvivorDiffusionGradientWitness(n);
        Assert.InRange(w.PowerSlope, 1.6, 2.4);     // dRe ~ |grad|^2 (the diffusion Rayleigh quotient)
        Assert.True(w.RatioCv < 0.15, $"dRe/grad^2 CV = {w.RatioCv} (want bond-independent < 0.15)");
        Assert.True(w.LawHolds);
    }

    [Fact]
    public void Survivor_is_a_density_mode_in_the_strong_dephasing_limit()
    {
        // The diffusion-Rayleigh derivation's premise: at Q -> 0 (gamma >> J) the survivor's mode
        // operator is (nearly) diagonal, a pure density mode (off-diagonal weight -> 0). At finite Q it
        // picks up coherence dressing, so the off-diagonal weight grows with Q (and the law degrades).
        var strong = new SurvivorDiffusionGradientWitness(5, 0.3);
        var dressed = new SurvivorDiffusionGradientWitness(5, 1.5);
        Assert.True(strong.OffDiagonalWeight < 0.1,
            $"strong-dephasing off-diag = {strong.OffDiagonalWeight} (want < 0.1: a clean density mode)");
        Assert.True(dressed.OffDiagonalWeight > strong.OffDiagonalWeight,
            "off-diagonal coherence dressing should grow with Q");
    }

    [Fact]
    public void Density_gradient_is_quiet_at_the_no_flux_ends()
    {
        var w = new SurvivorDiffusionGradientWitness(5);
        var b = w.Bonds;
        double endMax = System.Math.Max(b[0].GradSq, b[^1].GradSq);
        double interiorMax = b.Skip(1).Take(b.Count - 2).Max(x => x.GradSq);
        Assert.True(endMax < 0.5 * interiorMax,     // reflecting (no-flux) boundary: gradient -> 0 at the ends
            $"end grad^2 max {endMax} should be < 0.5 * interior max {interiorMax}");
    }
}
