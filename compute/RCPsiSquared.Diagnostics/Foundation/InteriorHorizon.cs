namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The cusp mathematics of the ¼-to-½ interior, the third navigation axis. Closed forms
/// (the F95 heading, the F56 critical-slowing iteration count, the F57 dwell, the F25 Bell+ geodesic)
/// plus the one live computation: running the Mandelbrot recursion u → u² + c and counting the steps,
/// which diverge as CΨ → ¼ (the recursion crawling at its own cusp, the horizon where time stops).
///
/// <para>The horizon CΨ = ¼ is the cardioid cusp / saddle-node fold of the R = C(Ψ + R)² recursion,
/// a structural / dephasing threshold (hardware-observed on IBM Kingston), NOT a gravitational
/// horizon. The critical slowing here is closed-form and hardware-confirmed; see
/// experiments/CRITICAL_SLOWING_AT_THE_CUSP.md and the Confirmations entries f25_cusp_trajectory,
/// f57_kdwell_gamma_invariance.</para></summary>
public static class InteriorHorizon
{
    /// <summary>The cusp / horizon: CΨ = 1/4 (the saddle-node fold, the quantum-classical door).</summary>
    public const double Cusp = 0.25;

    /// <summary>The dimension anchor: CΨ = 1/2 = 1/d at d = 2.</summary>
    public const double Anchor = 0.5;

    /// <summary>The Bell+ pure-Z dwell prefactor (F57): K_dwell = γ·t_dwell = 1.080088·δ.</summary>
    public const double BellPlusDwellPrefactor = 1.080088;

    /// <summary>The discriminant D = 1 − 4·CΨ of the recursion's quadratic. D &gt; 0 classical (two
    /// real fixed points), D = 0 the cusp, D &lt; 0 quantum (complex-conjugate roots).</summary>
    public static double Discriminant(double cpsi) => 1.0 - 4.0 * cpsi;

    /// <summary>"classical" / "cusp" / "quantum" from the discriminant sign.</summary>
    public static string Regime(double cpsi)
    {
        double d = Discriminant(cpsi);
        if (Math.Abs(d) < 1e-12) return "cusp";
        return d > 0 ? "classical" : "quantum";
    }

    /// <summary>The F95 heading θ = arctan(√(4·CΨ − 1)) in radians, the compass on the interior side.
    /// θ = 0 at the cusp ¼, 45° at the anchor ½. Clamped to 0 below the cusp (no interior heading).</summary>
    public static double Heading(double cpsi)
    {
        double arg = 4.0 * cpsi - 1.0;
        if (arg <= 0.0) return 0.0;
        return Math.Atan(Math.Sqrt(arg));
    }

    /// <summary>The heading in degrees.</summary>
    public static double HeadingDegrees(double cpsi) => Heading(cpsi) * 180.0 / Math.PI;

    /// <summary>The live Mandelbrot iteration count: u_{n+1} = u_n² + c with c = CΨ, u_0 = c, stopping
    /// when |u_{n+1} − u_n| &lt; <paramref name="tol"/>. For CΨ &lt; ¼ (the classical side) the iteration
    /// converges and the count diverges as CΨ → ¼ (critical slowing). Returns −1 if the iteration
    /// diverges (CΨ ≥ ¼, the quantum side has no fixed point) and <paramref name="maxIter"/> if it does
    /// not converge within the cap. The count is 1-based (the reference tables in
    /// experiments/CRITICAL_SLOWING_AT_THE_CUSP.md are 0-based, so this returns one more); the rescaled
    /// K = n·√ε is unaffected at the relevant ε.</summary>
    public static int RecursionIterations(double cpsi, double tol, int maxIter = 10_000_000)
    {
        double c = cpsi;
        double u = c;
        for (int n = 1; n <= maxIter; n++)
        {
            double uNext = u * u + c;
            if (Math.Abs(uNext - u) < tol) return n;
            u = uNext;
            if (double.IsInfinity(u) || u > 1e6) return -1; // diverged: CΨ ≥ ¼, no convergence
        }
        return maxIter;
    }

    /// <summary>The recursion count with a relative stop criterion tol = k·ε (ε = ¼ − CΨ). The rescaled
    /// K = n·√ε then converges to the constant ½·ln(4/k): the slowing was the stop criterion's.</summary>
    public static int RecursionIterationsRelative(double cpsi, double k, int maxIter = 10_000_000)
    {
        double eps = Cusp - cpsi;
        return RecursionIterations(cpsi, k * eps, maxIter);
    }

    /// <summary>The closed-form rescaled iteration count (F56), zero fit parameters:
    /// K(ε, tol) = ½·ln(4ε/tol) + α(tol)·√ε, α(tol) = −4 + ½·ln(16·tol), ε = ¼ − CΨ.</summary>
    public static double RecursionKClosedForm(double cpsi, double tol)
    {
        double eps = Cusp - cpsi;
        double alpha = -4.0 + 0.5 * Math.Log(16.0 * tol);
        return 0.5 * Math.Log(4.0 * eps / tol) + alpha * Math.Sqrt(eps);
    }

    /// <summary>F25 Bell+ coherence under Z-dephasing: CΨ(t) = f·(1 + f²)/6, f = exp(−4·γ·t). Starts
    /// at 1/3 (t = 0), decays monotonically through the cusp ¼ to 0. The geodesic the hardware
    /// confirmed crossing CΨ = ¼ on IBM Kingston.</summary>
    public static double BellPlusCpsi(double gamma, double t)
    {
        double f = Math.Exp(-4.0 * gamma * t);
        return f * (1.0 + f * f) / 6.0;
    }

    /// <summary>The dwell time near the cusp (F57): t_dwell = prefactor·δ/γ (Bell+, pure Z). Scales as
    /// 1/γ; faster dephasing compresses the window.</summary>
    public static double DwellTime(double gamma, double delta) => BellPlusDwellPrefactor * delta / gamma;

    /// <summary>The γ-invariant dwell in K-units: K_dwell = γ·t_dwell = prefactor·δ (F57), the fixed
    /// dose carrying a Bell+ state through the fold regardless of how bright γ is.</summary>
    public static double DwellK(double delta) => BellPlusDwellPrefactor * delta;
}
