namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Closed-form readings around the scalar value CΨ = ¼. The recurrence
/// u → u²+c has discriminant 1−4c and a double fixed point at c=¼; the F95 angle is a
/// separate positive-b quadratic coordinate, and the Bell+/pure-Z curve is a named finite
/// trajectory whose scalar readout happens to pass through ¼. This class keeps those objects
/// separate while retaining the F56 iteration-count asymptotic and the F57 dwell calculation.</summary>
public static class InteriorHorizon
{
    /// <summary>The recurrence's double-root value c = 1/4.</summary>
    public const double Cusp = 0.25;

    /// <summary>The dimension anchor: CΨ = 1/2 = 1/d at d = 2.</summary>
    public const double Anchor = 0.5;

    /// <summary>The Bell+ pure-Z dwell prefactor (F57): K_dwell = γ·t_dwell = 1.080088·δ.</summary>
    public const double BellPlusDwellPrefactor = 1.080088;

    /// <summary>The discriminant D = 1 − 4·CΨ of the recurrence fixed-point quadratic.
    /// D &gt; 0 gives two real roots, D = 0 a double root, and D &lt; 0 a complex-conjugate pair.</summary>
    public static double Discriminant(double cpsi) => 1.0 - 4.0 * cpsi;

    /// <summary>The root type from the discriminant sign, read exactly. ¼ is
    /// dyadic, so <c>1 − 4·CΨ</c> is exactly 0.0 there and no tolerance band is needed; a band would
    /// only put this reading at odds with <see cref="Heading"/>, whose own clamp is the exact
    /// <c>4·CΨ − 1 ≤ 0</c>, in whatever window it drew.</summary>
    public static string Regime(double cpsi)
    {
        double d = Discriminant(cpsi);
        if (d == 0.0) return "double-root";
        return d > 0 ? "two-real-roots" : "complex-root-pair";
    }

    /// <summary>The positive-b quadratic coordinate θ = arctan(√(4·CΨ − 1)) in radians.
    /// θ = 0 at the boundary ¼ and 45° at the anchor ½. Clamped to 0 below the boundary,
    /// where this real-angle coordinate is not defined.</summary>
    public static double Heading(double cpsi)
    {
        double arg = 4.0 * cpsi - 1.0;
        if (arg <= 0.0) return 0.0;
        return Math.Atan(Math.Sqrt(arg));
    }

    /// <summary>The heading in degrees.</summary>
    public static double HeadingDegrees(double cpsi) => Heading(cpsi) * 180.0 / Math.PI;

    /// <summary>The live Mandelbrot iteration count: u_{n+1} = u_n² + c with c = CΨ, u_0 = c, stopping
    /// when |u_{n+1} − u_n| &lt; <paramref name="tol"/>. For CΨ &lt; ¼ (the two-real-root side) the iteration
    /// converges and the count diverges as CΨ → ¼ (critical slowing). Returns <paramref name="maxIter"/>
    /// if the stop criterion is not met within the cap, and −1 if the orbit escapes first.
    ///
    /// <para>The −1 branch is a reading of the stop criterion, not of whether a fixed point exists,
    /// and the two do not coincide. Writing c = ¼ + δ and u = ½ + v turns the step into
    /// Δu = v² + δ, whose minimum over the crawl is δ itself; so above the cusp the increment still
    /// falls below tol, and this returns a finite count, exactly while δ &lt; tol. Escape and −1
    /// begin at δ &gt; tol. At tol = 10⁻⁹ that puts CΨ = ¼ (31611) and CΨ = ¼ + 10⁻¹³ (31612) on the
    /// finite side and CΨ = ¼ + 10⁻⁶ on the −1 side. Read the count as "how long the orbit crawled
    /// slower than tol", and read the regime off <see cref="Discriminant"/> instead.</para>
    ///
    /// <para>The count is 1-based (the reference tables in
    /// experiments/CRITICAL_SLOWING_AT_THE_CUSP.md are 0-based, so this returns one more); the rescaled
    /// K = n·√ε is unaffected at the relevant ε.</para></summary>
    public static int RecursionIterations(double cpsi, double tol, int maxIter = 10_000_000)
    {
        double c = cpsi;
        double u = c;
        for (int n = 1; n <= maxIter; n++)
        {
            double uNext = u * u + c;
            if (Math.Abs(uNext - u) < tol) return n;
            u = uNext;
            if (double.IsInfinity(u) || u > 1e6) return -1; // escaped before the increment fell below tol
        }
        return maxIter;
    }

    /// <summary>The recursion count with a relative stop criterion tol = k·ε (ε = ¼ − CΨ). The rescaled
    /// K = n·√ε then converges to the constant ½·ln(4/k) as ε → 0: the slowing was the stop criterion's.
    /// Requires CΨ &lt; ¼ and k &gt; 0; at or above the cusp ε ≤ 0 makes the stop tolerance
    /// non-positive, which no increment can undercut, so the loop would silently run to the cap
    /// instead of reporting anything.</summary>
    public static int RecursionIterationsRelative(double cpsi, double k, int maxIter = 10_000_000)
    {
        if (cpsi >= Cusp)
            throw new ArgumentOutOfRangeException(nameof(cpsi), cpsi, "the relative stop needs ε = ¼ − CΨ > 0; at or above the cusp there is no relative scale to stop against.");
        if (k <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(k), k, "k must be positive.");
        double eps = Cusp - cpsi;
        return RecursionIterations(cpsi, k * eps, maxIter);
    }

    /// <summary>The zero-fit asymptotic expansion for the rescaled iteration count (F56):
    /// K(ε, tol) = ½·ln(4ε/tol) + α(tol)·√ε, α(tol) = −4 + ½·ln(16·tol), ε = ¼ − CΨ.</summary>
    public static double RecursionKAsymptotic(double cpsi, double tol)
    {
        double eps = Cusp - cpsi;
        double alpha = -4.0 + 0.5 * Math.Log(16.0 * tol);
        return 0.5 * Math.Log(4.0 * eps / tol) + alpha * Math.Sqrt(eps);
    }

    /// <summary>F25's named Bell+/pure-Z scalar trajectory:
    /// CΨ(t) = f·(1 + f²)/6, f = exp(−4·γ·t). It starts at 1/3 and decreases through ¼.
    /// The formula is setup-specific; the shared scalar value does not identify it with the recurrence.</summary>
    public static double BellPlusCpsi(double gamma, double t)
    {
        double f = Math.Exp(-4.0 * gamma * t);
        return f * (1.0 + f * f) / 6.0;
    }

    /// <summary>The dwell time near the selected scalar readout (F57): t_dwell = prefactor·δ/γ
    /// for the Bell+/pure-Z setup. It scales as
    /// 1/γ; faster dephasing compresses the window.</summary>
    public static double DwellTime(double gamma, double delta) => BellPlusDwellPrefactor * delta / gamma;

    /// <summary>The γ-invariant dwell in K-units: K_dwell = γ·t_dwell = prefactor·δ (F57), the fixed
    /// dose carrying the named Bell+ trajectory through that readout window.</summary>
    public static double DwellK(double delta) => BellPlusDwellPrefactor * delta;
}
