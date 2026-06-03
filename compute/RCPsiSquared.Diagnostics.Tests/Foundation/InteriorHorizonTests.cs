using System.Linq;
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

    [Fact]
    public void Field_Surfaces_Marks_Heading_Recursion_Seam_Dwell()
    {
        var field = new InteriorHorizonField(epsLo: 1e-4, epsHi: 0.25, epsPoints: 11,
            tol: 1e-12, relK: 1e-3, gamma: 0.5);
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)field).Children.ToList();
        var labels = children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("marks"));
        Assert.Contains(labels, l => l.Contains("heading"));
        Assert.Contains(labels, l => l.Contains("recursion"));
        Assert.Contains(labels, l => l.Contains("ours"));
        Assert.Contains(labels, l => l.Contains("dwell"));

        // The heading curve falls to ~0 at the horizon (its interior end, closest to 1/4).
        var heading = children.First(c => c.DisplayName.Contains("heading"));
        var hc = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(heading.Payload);
        Assert.True(hc.Y.Min() < 1.0, "heading should reach ~0 deg at the horizon");
        Assert.True(hc.Y.Max() > 30.0, "heading should rise toward 45 deg at the anchor");

        // The recursion curve diverges toward the horizon (its max count is large).
        var recursion = children.First(c => c.DisplayName.Contains("recursion"));
        var rc = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(recursion.Payload);
        Assert.True(rc.Y.Max() > 100.0, "the recursion should crawl (large iteration count) near the horizon");

        // The seam (the field's most distinctive reading): near the horizon the relative-stop K is
        // ~flat at the constant 1/2*ln(4/k), while the absolute-tol K drifts above it. The slowing
        // belonged to the stop criterion, not the cusp.
        var ours = children.First(c => c.DisplayName.Contains("ours"));
        var seam = ours.Children.ToList();
        var relCurve = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(
            seam.First(c => c.DisplayName.Contains("relative")).Payload);
        var absCurve = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(
            seam.First(c => c.DisplayName.Contains("absolute")).Payload);
        double kConst = 0.5 * System.Math.Log(4.0 / 1e-3); // = 4.147, the relative-stop constant
        Assert.True(System.Math.Abs(relCurve.Y[0] - kConst) < 0.2,
            $"relative-stop K near the horizon should be ~{kConst:F2}; got {relCurve.Y[0]:F2}");
        Assert.True(absCurve.Y[0] > relCurve.Y[0] + 2.0,
            $"absolute-tol K should drift above the relative-stop constant; abs {absCurve.Y[0]:F2} rel {relCurve.Y[0]:F2}");

        // Smoke: renders to JSON without throwing and carries the horizon story.
        var json = RCPsiSquared.Core.Inspection.InspectionJsonExporter.ToJson(field);
        Assert.Contains("horizon", json);
        Assert.Contains("recursion", json);
    }
}
