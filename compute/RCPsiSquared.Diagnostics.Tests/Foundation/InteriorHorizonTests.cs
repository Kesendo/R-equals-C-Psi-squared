using System;
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
        Assert.Equal("cusp", InteriorHorizon.Regime(0.25));      // D = 0, exactly (1/4 is dyadic)

        // Just off the cusp the three readings must agree: no tolerance band may call "cusp" a
        // point where the discriminant is negative and the heading is a real angle.
        Assert.Equal("quantum", InteriorHorizon.Regime(0.25 + 1e-13));
        Assert.True(InteriorHorizon.Discriminant(0.25 + 1e-13) < 0.0);
        Assert.True(InteriorHorizon.Heading(0.25 + 1e-13) > 0.0);
        Assert.Equal("classical", InteriorHorizon.Regime(0.25 - 1e-13));
        Assert.Equal(0.0, InteriorHorizon.Heading(0.25 - 1e-13));
    }

    [Fact]
    public void Recursion_LiveK_ApproachesAsymptoticFormula()
    {
        // Live Mandelbrot iteration count K = n*sqrt(eps) approaches the asymptotic expansion
        // K(eps,tol) = (1/2)*ln(4*eps/tol) + alpha(tol)*sqrt(eps) (CRITICAL_SLOWING_AT_THE_CUSP.md).
        const double tol = 1e-12;
        foreach (double eps in new[] { 1e-4, 1e-5 })
        {
            double cpsi = 0.25 - eps;
            int n = InteriorHorizon.RecursionIterations(cpsi, tol);
            double liveK = n * System.Math.Sqrt(eps);
            double asymptoticK = InteriorHorizon.RecursionKAsymptotic(cpsi, tol);
            _out.WriteLine($"eps={eps:E0} n={n} liveK={liveK:F3} asymptoticK={asymptoticK:F3}");
            Assert.True(System.Math.Abs(liveK - asymptoticK) < 0.3,
                $"liveK {liveK:F3} vs asymptoticK {asymptoticK:F3} at eps={eps:E0}");
        }
    }

    [Fact]
    public void Recursion_AsymptoticFormula_IsNotPresentedAsExactAtFiniteEpsilon()
    {
        const double eps = 1e-4, tol = 1e-12;
        int n = InteriorHorizon.RecursionIterations(0.25 - eps, tol);
        double residual = n * System.Math.Sqrt(eps)
                        - InteriorHorizon.RecursionKAsymptotic(0.25 - eps, tol);

        Assert.True(System.Math.Abs(residual) > 1e-3, $"finite-epsilon residual unexpectedly vanished: {residual:E3}");
        Assert.True(System.Math.Abs(residual) < 0.3, $"asymptotic residual too large in its stated regime: {residual:E3}");
    }

    [Theory]
    [InlineData(1e-3)]
    [InlineData(1e-4)]
    [InlineData(1e-5)]
    public void Recursion_AsymptoticAlphaTerm_EarnsItsPlace(double eps)
    {
        // F56's content beyond the leading log is alpha(tol) = -4 + (1/2)*ln(16*tol). The <0.3 gate
        // above passes with alpha dropped entirely, so it certifies only the logarithm. This one
        // compares the shipped form against its two truncations on the same live count: alpha = 0
        // (no correction at all) and alpha = -4 (transient only, no Modified-Equation term). The
        // shipped form is closer by more than 5x at every eps, and halving alpha also fails that.
        const double tol = 1e-12;
        double cpsi = 0.25 - eps;
        double liveK = InteriorHorizon.RecursionIterations(cpsi, tol) * System.Math.Sqrt(eps);
        double leadingLog = 0.5 * System.Math.Log(4.0 * eps / tol);

        double resFull = System.Math.Abs(liveK - InteriorHorizon.RecursionKAsymptotic(cpsi, tol));
        double resNoAlpha = System.Math.Abs(liveK - leadingLog);
        double resTransientOnly = System.Math.Abs(liveK - (leadingLog - 4.0 * System.Math.Sqrt(eps)));

        _out.WriteLine($"eps={eps:E0} full={resFull:F4} alpha=0 {resNoAlpha:F4} alpha=-4 {resTransientOnly:F4}");
        Assert.True(resFull * 5.0 < resNoAlpha, $"dropping alpha should hurt: {resFull:F4} vs {resNoAlpha:F4}");
        Assert.True(resFull * 5.0 < resTransientOnly, $"the ln(16 tol) half of alpha should earn its place: {resFull:F4} vs {resTransientOnly:F4}");
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

    [Theory]
    [InlineData(0.0)]        // the cusp itself
    [InlineData(1e-13)]      // above it, and inside the tolerance
    [InlineData(0.99e-9)]    // just inside
    public void Recursion_AboveTheCusp_StillReturnsAFiniteCount_WhileDeltaIsBelowTol(double delta)
    {
        // The -1 branch is not the CPsi >= 1/4 boundary. Writing c = 1/4 + delta and u = 1/2 + v
        // gives step = v^2 + delta, whose minimum over the crawl is delta; so the increment still
        // falls under tol above the cusp, and the count is finite exactly while delta < tol.
        const double tol = 1e-9;
        int n = InteriorHorizon.RecursionIterations(0.25 + delta, tol);
        _out.WriteLine($"delta={delta:E2} n={n}");
        Assert.True(n > 0, $"expected a finite count at delta={delta:E2} (< tol), got {n}");
    }

    [Theory]
    [InlineData(1.01e-9)]    // just outside
    [InlineData(1e-6)]
    [InlineData(1e-2)]
    public void Recursion_AboveTheCusp_ReturnsMinusOne_OnceDeltaExceedsTol(double delta)
    {
        // The other side of the same threshold: the crawl never gets slower than tol, so the orbit
        // escapes and the reading is -1. The two theories bracket delta = tol from both sides at 1 %,
        // which is what makes the pair a gate on the law rather than on two convenient points.
        const double tol = 1e-9;
        Assert.Equal(-1, InteriorHorizon.RecursionIterations(0.25 + delta, tol));
    }

    [Fact]
    public void Recursion_TheEscapeThreshold_MovesWithTolAndNotWithTheCusp()
    {
        // Same CPsi, two tolerances, opposite answers: the threshold belongs to the stop criterion.
        // If -1 marked "CPsi >= 1/4" these two would have to agree.
        const double cpsi = 0.25 + 1e-10;
        Assert.True(InteriorHorizon.RecursionIterations(cpsi, tol: 1e-9) > 0);
        Assert.Equal(-1, InteriorHorizon.RecursionIterations(cpsi, tol: 1e-11));
    }

    [Theory]
    [InlineData(0.25)]
    [InlineData(0.3)]
    public void Recursion_RelativeStop_RefusesTheCuspAndBeyond(double cpsi)
    {
        // eps <= 0 makes the stop tolerance non-positive, which no increment can undercut: without
        // the guard the call runs silently to the 10-million cap and returns it as a count.
        Assert.Throws<ArgumentOutOfRangeException>(() => InteriorHorizon.RecursionIterationsRelative(cpsi, 1e-3));
    }

    [Fact]
    public void Recursion_RelativeStop_RefusesNonPositiveK()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => InteriorHorizon.RecursionIterationsRelative(0.24, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => InteriorHorizon.RecursionIterationsRelative(0.24, -1e-3));
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
        var field = new InteriorHorizonField();   // the shipped defaults, inside F56's committed range
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
        // The first point IS the cusp, where the closed form clamps to exactly 0, so a "< 1 deg"
        // gate there would hold for any heading formula whatsoever. Read it exactly instead, and
        // gate the approach on the interior rungs, which are the ones that carry the shape.
        Assert.Equal(0.0, hc.Y[0]);
        Assert.True(hc.Y[1] > 0.0, "the first interior rung must be off the cusp");
        for (int i = 2; i < hc.Y.Count; i++)
            Assert.True(hc.Y[i] > hc.Y[i - 1], $"the heading must rise away from the horizon; fell at rung {i}");

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
        Assert.True(System.Math.Abs(relCurve.Y[0] - kConst) < 0.01,
            $"relative-stop K near the horizon should be ~{kConst:F3}; got {relCurve.Y[0]:F3}");

        // The content is that one curve is flat where the other is not, and "drifts above" is not
        // it: which of the two is larger depends on eps against tol/k, and at the near end the
        // absolute K is the SMALLER one. Gate the spread over the same rungs instead.
        int rungs = System.Math.Min(6, relCurve.Y.Count);
        double relSpread = relCurve.Y.Take(rungs).Max() - relCurve.Y.Take(rungs).Min();
        double absSpread = absCurve.Y.Take(rungs).Max() - absCurve.Y.Take(rungs).Min();
        _out.WriteLine($"over {rungs} rungs: rel spread {relSpread:F4}, abs spread {absSpread:F4}");
        Assert.True(relSpread < 0.02, $"the relative-stop K should be flat near the horizon; spread {relSpread:F4}");
        Assert.True(absSpread > 20.0 * relSpread, $"the absolute-tol K should drift; abs {absSpread:F4} vs rel {relSpread:F4}");

        // Smoke: renders to JSON without throwing and carries the horizon story.
        var json = RCPsiSquared.Core.Inspection.InspectionJsonExporter.ToJson(field);
        Assert.Contains("horizon", json);
        Assert.Contains("recursion", json);
    }

    [Fact]
    public void Field_Ctor_RejectsInvalidInputs()
    {
        // Each input flows straight into the closed forms (the ε-ladder, log/sqrt/divide), so the ctor
        // guards them. One bad value per case, the rest left at their valid defaults. ParamName is
        // asserted too: a guard that blames the wrong argument sends the reader to the wrong knob,
        // and only naming the parameter can catch that.
        AssertRejects("epsPoints", () => new InteriorHorizonField(epsPoints: 1));                  // < 2
        AssertRejects("epsLo", () => new InteriorHorizonField(epsLo: 0.0));                        // ≤ 0
        AssertRejects("epsHi", () => new InteriorHorizonField(epsHi: 0.30));                       // > 1/4
        AssertRejects("epsHi", () => new InteriorHorizonField(epsLo: 0.2, epsHi: 0.1));            // hi ≤ lo
        AssertRejects("tol", () => new InteriorHorizonField(tol: 0.0));                            // ≤ 0
        AssertRejects("relK", () => new InteriorHorizonField(relK: -1.0));                         // ≤ 0
        AssertRejects("gamma", () => new InteriorHorizonField(gamma: 0.0));                        // ≤ 0
    }

    private static void AssertRejects(string paramName, Func<InteriorHorizonField> build)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => build());
        Assert.Equal(paramName, ex.ParamName);
    }
}
