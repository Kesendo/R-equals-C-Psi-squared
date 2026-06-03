namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The family of approach shapes to the cusp ¼: how the approach depends on the start. For the
/// partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ (Bell+ is α=π/4) under Z-dephasing,
/// the coherence CΨ(t) is a two-exponential at the odd-harmonic rates 4γ and 12γ with closed-form,
/// entanglement-dependent weights (verified bit-exact against the Lindblad evolution):
///
/// <para>CΨ(s,t) = w₀·e^(−4γt) + w₁·e^(−12γt),  w₀ = s(1−s²/2)/3,  w₁ = s³/6,  s = sin2α.</para>
///
/// <para>The scaling laws: the start CΨ(0) = s/3 (the start height is the entanglement); it crosses ¼
/// only if s &gt; 3/4 (a threshold; s = 3/4 starts exactly on ¼); the fast mode (12γ) carries a fraction
/// s²/2 of the start, growing quadratically (Bell+, s = 1, is the 50/50 member); and every member shares
/// the slowest carrier rate 4γ (the eigenvalue −γ₀) and collapses onto it at late time. The shape is the
/// early harmonic transient; the carrier is universal. This is the "the slowing is ours" reading across
/// the whole family. Tier-1 closed form; the Bell+ member reproduces F25
/// (<see cref="InteriorHorizon.BellPlusCpsi"/>) exactly.</para>
///
/// <para>Wired into the typed-knowledge graph (Core) as <c>ApproachFamilyCarrierClaim</c>: the shared
/// carrier 4γ as the Universal Carrier, the c=2 doubled-PTF kinship to <c>C2BareDoubledPtfClosedForm</c>
/// (decay vs susceptibility, a viewpoint not an identity), the algebra/dynamics two readings, and the
/// Bell+ member as F25. Render it with <c>inspect --claim ApproachFamilyCarrierClaim</c>.</para></summary>
public static class OddHarmonicApproach
{
    /// <summary>The cusp / horizon CΨ = ¼.</summary>
    public const double Cusp = 0.25;

    /// <summary>The entanglement threshold s = 3/4: the approach crosses ¼ iff s &gt; 3/4 (CΨ(0)=s/3 &gt; ¼).</summary>
    public const double CrossingThreshold = 0.75;

    /// <summary>The Bell+ member: s = sin2α = 1 (α = π/4), maximally entangled, 50/50 weights.</summary>
    public const double BellPlusS = 1.0;

    /// <summary>The carrier (slowest) rate 4γ: the HD=1 mode, the eigenvalue −γ₀ scaled, shared by every
    /// member of the family.</summary>
    public static double CarrierRate(double gamma) => 4.0 * gamma;

    /// <summary>The harmonic (fast) rate 12γ: the second odd harmonic, the early transient that
    /// distinguishes the shapes.</summary>
    public static double HarmonicRate(double gamma) => 12.0 * gamma;

    /// <summary>The closed-form weights (w₀ carrier, w₁ harmonic) for s = sin2α:
    /// w₀ = s(1−s²/2)/3, w₁ = s³/6.</summary>
    public static (double Carrier, double Harmonic) Weights(double s) =>
        (s * (1.0 - 0.5 * s * s) / 3.0, s * s * s / 6.0);

    /// <summary>CΨ(s,t) = w₀·e^(−4γt) + w₁·e^(−12γt).</summary>
    public static double Cpsi(double s, double gamma, double t)
    {
        var (w0, w1) = Weights(s);
        return w0 * Math.Exp(-CarrierRate(gamma) * t) + w1 * Math.Exp(-HarmonicRate(gamma) * t);
    }

    /// <summary>The start CΨ(0) = w₀ + w₁ = s/3: the start height is the entanglement.</summary>
    public static double InitialCpsi(double s) => s / 3.0;

    /// <summary>Whether the approach crosses ¼ from above: CΨ(0) = s/3 &gt; ¼, i.e. s &gt; 3/4.</summary>
    public static bool Crosses(double s) => InitialCpsi(s) > Cusp;

    /// <summary>The harmonic fraction at t = 0: w₁/(w₀+w₁) = s²/2 (the fast-mode content of the start;
    /// ½ for Bell+).</summary>
    public static double HarmonicFraction(double s) => 0.5 * s * s;

    /// <summary>The time CΨ(s,·) crosses ¼, by bisection on the monotone-decreasing two-exponential
    /// (both terms positive and decreasing). NaN if the approach never reaches ¼ (s ≤ 3/4: the start is
    /// at or below the cusp).</summary>
    public static double CrossingTime(double s, double gamma)
    {
        if (!Crosses(s)) return double.NaN;
        double lo = 0.0, hi = 1.0;
        while (Cpsi(s, gamma, hi) > Cusp) hi *= 2.0;
        for (int i = 0; i < 200; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (Cpsi(s, gamma, mid) > Cusp) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }
}
