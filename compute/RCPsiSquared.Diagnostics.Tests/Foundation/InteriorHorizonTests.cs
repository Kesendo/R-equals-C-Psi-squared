using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class InteriorHorizonTests
{
    private readonly ITestOutputHelper _out;
    public InteriorHorizonTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Heading_HasCanonicalAnchors()
    {
        // F95: theta = arctan(sqrt(4*CPsi - 1)). At the cusp 1/4 theta = 0 (degenerate root);
        // at 1/3 theta = 30 deg; at the anchor 1/2 theta = 45 deg.
        Assert.Equal(0.0, InteriorHorizon.HeadingDegrees(0.25), 6);
        Assert.Equal(30.0, InteriorHorizon.HeadingDegrees(1.0 / 3.0), 4);
        Assert.Equal(45.0, InteriorHorizon.HeadingDegrees(0.5), 6);
        // Below the cusp the interior heading is not defined; clamp to 0.
        Assert.Equal(0.0, InteriorHorizon.HeadingDegrees(0.10), 12);
    }

    [Fact]
    public void Regime_FollowsDiscriminantSign()
    {
        Assert.Equal("quantum", InteriorHorizon.Regime(0.30));   // D = 1 - 4*0.30 = -0.2 < 0
        Assert.Equal("classical", InteriorHorizon.Regime(0.20)); // D = 1 - 4*0.20 = +0.2 > 0
        Assert.Equal("cusp", InteriorHorizon.Regime(0.25));      // D = 0
    }

    [Fact]
    public void Recursion_LiveK_MatchesClosedForm()
    {
        // Live Mandelbrot iteration count K = n*sqrt(eps) must match the closed form
        // K(eps,tol) = (1/2)*ln(4*eps/tol) + alpha(tol)*sqrt(eps) (CRITICAL_SLOWING_AT_THE_CUSP.md).
        const double tol = 1e-12;
        foreach (double eps in new[] { 1e-4, 1e-5 })
        {
            double cpsi = 0.25 - eps;
            int n = InteriorHorizon.RecursionIterations(cpsi, tol);
            double liveK = n * System.Math.Sqrt(eps);
            double closedK = InteriorHorizon.RecursionKClosedForm(cpsi, tol);
            _out.WriteLine($"eps={eps:E0} n={n} liveK={liveK:F3} closedK={closedK:F3}");
            Assert.True(System.Math.Abs(liveK - closedK) < 0.3,
                $"liveK {liveK:F3} vs closedK {closedK:F3} at eps={eps:E0}");
        }
    }

    [Fact]
    public void Recursion_IterationCount_DivergesAtTheHorizon()
    {
        // The raw iteration count diverges as CPsi -> 1/4 from below (the recursion crawls; time stops).
        const double tol = 1e-12;
        int nFar = InteriorHorizon.RecursionIterations(0.25 - 1e-2, tol);
        int nNear = InteriorHorizon.RecursionIterations(0.25 - 1e-4, tol);
        _out.WriteLine($"nFar(eps=1e-2)={nFar}  nNear(eps=1e-4)={nNear}");
        Assert.True(nNear > 5 * nFar, $"recursion should crawl near the horizon: nNear {nNear} vs nFar {nFar}");
    }

    [Fact]
    public void Recursion_RelativeStop_GivesConstantK()
    {
        // With a relative stop tol = k*eps, the rescaled K is a pure constant (1/2)*ln(4/k):
        // the critical slowing belonged to the stop criterion, not the cusp.
        const double k = 1e-3;
        double expected = 0.5 * System.Math.Log(4.0 / k); // = 4.1447...
        foreach (double eps in new[] { 1e-6, 1e-7 })
        {
            double cpsi = 0.25 - eps;
            int n = InteriorHorizon.RecursionIterationsRelative(cpsi, k);
            double kRel = n * System.Math.Sqrt(eps);
            _out.WriteLine($"eps={eps:E0} kRel={kRel:F3} expected={expected:F3}");
            Assert.True(System.Math.Abs(kRel - expected) < 0.05,
                $"relative-stop K {kRel:F3} should equal {expected:F3} at eps={eps:E0}");
        }
    }

    [Fact]
    public void DwellK_IsGammaInvariant()
    {
        // K_dwell = gamma * t_dwell = prefactor * delta, independent of gamma (F57).
        const double delta = 1e-3;
        double kRef = InteriorHorizon.DwellK(delta);
        foreach (double gamma in new[] { 0.1, 1.0, 10.0 })
        {
            double kFromTime = gamma * InteriorHorizon.DwellTime(gamma, delta);
            Assert.Equal(kRef, kFromTime, 12);
        }
        Assert.Equal(1.080088 * delta, kRef, 9); // Bell+ pure-Z prefactor
    }

    [Fact]
    public void BellPlusCpsi_StartsAtThird_CrossesQuarter()
    {
        // F25: CPsi(t) = f*(1+f^2)/6, f = exp(-4*gamma*t). At t=0 (f=1) CPsi = 1/3; it decays
        // monotonically through the cusp 1/4 to 0.
        Assert.Equal(1.0 / 3.0, InteriorHorizon.BellPlusCpsi(gamma: 0.5, t: 0.0), 9);
        Assert.True(InteriorHorizon.BellPlusCpsi(0.5, 5.0) < 0.25, "CPsi should be below 1/4 at large t");
    }
}
