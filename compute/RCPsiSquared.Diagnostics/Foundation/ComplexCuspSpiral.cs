namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>A named Bell+/pure-Z complex-coherence trajectory. With a common Z-drift Ω,
/// CΨ_com = C·Ψ_com winds inward while its magnitude follows F25. The set |CΨ_com|=¼ is a
/// radial readout ring. It is not the period-one cardioid, whose boundary is parameterized by
/// c=e^{iφ}/2−e^{2iφ}/4 and has only the real point c=+¼ as its cusp.
///
/// <para>Logarithmic is the right word and it was measured, not assumed: over the observed run
/// from |CΨ_com| = 1/3 to the ¼ radial ring, ln|CΨ_com| is linear in the phase to 0.0018 out of a
/// total change of 0.288, i.e. 0.6%. It stays that good further out (2% over thirty e-folds of
/// f), because the only departure from linearity is the slowly varying factor
/// <see cref="MagnitudeEFoldRate"/>/4γ = (1+3f²)/(1+f²), which moves from 2 to 1.85 across the
/// visited range. That factor matters for <see cref="WindingRate"/>, which divides by the wrong
/// e-fold if it is ignored; it does not bend the spiral appreciably.</para>
///
/// <para>The radial magnitude is unchanged and Ω-independent: |CΨ_com|(t) = f(1+f²)/6, f = e^{−4γt}
/// (F25, the same law <see cref="InteriorHorizon.BellPlusCpsi"/> reads on the line). So every spiral
/// reaches the same radial ring at the same time; only the argument φ₀ − Ω·t_cross is free. Saved
/// Kingston rows show that argument changing under an injected Z rotation. They do not realize
/// F95 roots or trace the F97 cardioid.</para>
///
/// <para>The same scalar number may label the recurrence double root and this radial crossing, but
/// the equations and objects remain distinct. Closed-form and N-free within the named setup.</para></summary>
public static class ComplexCuspSpiral
{
    /// <summary>The selected radial readout radius |CΨ_com| = ¼.</summary>
    public const double CircleRadius = 0.25;

    /// <summary>The maximum reachable magnitude, |CΨ_com|(0) = 1/3 (Bell+ at t = 0). No circle of
    /// larger radius is ever crossed.</summary>
    public const double InitialMagnitude = 1.0 / 3.0;

    /// <summary>The radial magnitude |CΨ_com|(t) = f(1+f²)/6, f = e^{−4γt} (F25), Ω-independent. The
    /// same law <see cref="InteriorHorizon.BellPlusCpsi"/> reads on the real line.</summary>
    public static double Magnitude(double gamma, double t) => InteriorHorizon.BellPlusCpsi(gamma, t);

    /// <summary>The winding angle arg(CΨ_com)(t) = φ₀ − Ω·t (radians), linear in t at rate −Ω. The sign
    /// matches the physical evolution ρ₀₃(t) = ρ₀₃(0)·e^{−iΩt} (clockwise for Ω &gt; 0).</summary>
    public static double Argument(double omega, double phi0, double t) => phi0 - omega * t;

    /// <summary>The real part of CΨ_com(t) = |CΨ_com|·cos(arg): the spiral's x in the c-plane.</summary>
    public static double Re(double gamma, double omega, double phi0, double t) =>
        Magnitude(gamma, t) * Math.Cos(Argument(omega, phi0, t));

    /// <summary>The imaginary part of CΨ_com(t) = |CΨ_com|·sin(arg): the spiral's y in the c-plane.</summary>
    public static double Im(double gamma, double omega, double phi0, double t) =>
        Magnitude(gamma, t) * Math.Sin(Argument(omega, phi0, t));

    /// <summary>The winding rate Ω/(4γ): radians of phase winding per e-fold of the COHERENCE
    /// FACTOR f = e^{−4γt}, which is what decays at 4γ. The spiral's tightness; 0 at Ω = 0 (no
    /// winding, the real-axis line).
    ///
    /// <para>Not per e-fold of the magnitude, and the difference is a factor of about two over the
    /// range this class is about. |CΨ_com| = f(1+f²)/6 shrinks at
    /// <see cref="MagnitudeEFoldRate"/> = 4γ·(1+3f²)/(1+f²), which is 8γ at t = 0 and reaches 4γ
    /// only as f → 0. At the ¼ crossing it is 4γ·1.8517012, and that factor is F25's crossing
    /// derivative |dCΨ/dt|/γ itself, not a number that resembles it: the rate is |dCΨ/dt|/|CΨ| and
    /// |CΨ| is ¼ there, so the two agree exactly whenever f(1+f²) = 3/2. Averaged across the run
    /// from the start 1/3 to the ¼ radial ring the true figure is Ω·t_cross / ln(4/3)
    /// = Ω/(4γ) · (−ln f*)/ln(4/3), which is Ω/(4γ) divided by 1.9255760, at every γ and
    /// Ω.</para></summary>
    public static double WindingRate(double gamma, double omega) => omega / (4.0 * gamma);

    /// <summary>The magnitude's own e-folding rate |d ln|CΨ_com|/dt| = 4γ·(1+3f²)/(1+f²) at time t,
    /// the quantity <see cref="WindingRate"/> is NOT divided by. Runs from 8γ at t = 0 down to 4γ
    /// as t → ∞; the two rates coincide only in that limit.</summary>
    public static double MagnitudeEFoldRate(double gamma, double t)
    {
        double f = Math.Exp(-4.0 * gamma * t);
        return 4.0 * gamma * (1.0 + 3.0 * f * f) / (1.0 + f * f);
    }

    /// <summary>Radians of phase swept per e-fold of the MAGNITUDE, averaged over the run from the
    /// start |CΨ_com| = 1/3 to the circle of the given radius: Ω·t_cross / ln((1/3)/radius). This is
    /// the figure <see cref="WindingRate"/> is often mistaken for; at the quarter-radius ring it is
    /// Ω/(4γ) divided by 1.9255760, independently of γ and Ω. NaN where there is no
    /// crossing.</summary>
    public static double PhasePerMagnitudeEFold(double gamma, double omega, double radius = CircleRadius)
    {
        double tc = CrossingTime(gamma, radius);
        if (double.IsNaN(tc) || tc <= 0.0) return double.NaN;
        return omega * tc / Math.Log(InitialMagnitude / radius);
    }

    /// <summary>The number of full turns |Ω|·tMax/(2π) the spiral makes over [0, tMax].</summary>
    public static double WindingNumber(double omega, double tMax) => Math.Abs(omega) * tMax / (2.0 * Math.PI);

    /// <summary>The time the magnitude reaches a radial ring of the given radius (default ¼),
    /// Ω-independent: solves f³ + f = 6·radius for f ∈ (0, 1] (Newton from f = 1), then t = −ln(f)/(4γ).
    /// The radial crossing every spiral shares. NaN if radius is not in (0, 1/3] (the magnitude starts at
    /// 1/3 and only decreases, so no larger circle is crossed). 0 exactly at radius = 1/3 (the start).</summary>
    public static double CrossingTime(double gamma, double radius = CircleRadius)
    {
        if (gamma <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be > 0: without dephasing the magnitude never leaves 1/3 and there is no crossing.");
        if (radius <= 0.0 || radius > InitialMagnitude) return double.NaN;
        double target = 6.0 * radius;             // f³ + f = 6r
        double f = 1.0;
        for (int i = 0; i < 64; i++)
        {
            double g = f * f * f + f - target;
            double gp = 3.0 * f * f + 1.0;
            double step = g / gp;
            f -= step;
            if (Math.Abs(step) < 1e-15) break;
        }
        if (f >= 1.0) return 0.0;                  // radius ≥ 1/3: at or above the start
        return -Math.Log(f) / (4.0 * gamma);
    }

    /// <summary>The argument at which the spiral reaches the radial ring of the given radius:
    /// arg at <see cref="CrossingTime"/> = φ₀ − Ω·t_cross. Saved Kingston rows report this
    /// setup-dependent argument under injected Z rotation. NaN if there is no
    /// crossing.</summary>
    public static double CrossingArgument(double gamma, double omega, double phi0, double radius = CircleRadius)
    {
        double tc = CrossingTime(gamma, radius);
        if (double.IsNaN(tc)) return double.NaN;
        return Argument(omega, phi0, tc);
    }
}
