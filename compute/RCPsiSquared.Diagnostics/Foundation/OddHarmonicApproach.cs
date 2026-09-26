namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>A named family of two-qubit Z-dephasing readout curves: how the scalar-quarter crossing
/// depends on the start. For the
/// partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ (Bell+ is α=π/4) under Z-dephasing,
/// the coherence CΨ(t) is a two-exponential at the odd-harmonic rates 4γ and 12γ with closed-form,
/// entanglement-dependent weights:
///
/// <para>CΨ(s,t) = w₀·e^(−4γt) + w₁·e^(−12γt),  w₀ = s(1−s²/2)/3,  w₁ = s³/6,  s = sin2α.</para>
///
/// <para>The two exponentials are not fitted. Z-dephasing leaves the populations alone and takes the
/// single coherence ρ₀₃ = cosα·sinα·f, f = e^{−4γt}; with s = sin2α that gives purity
/// 1 − s²/2 + s²f²/2 and ℓ₁-coherence s·f, so CΨ = purity·ℓ₁/(d−1) = (s/3)(1−s²/2)·f + (s³/6)·f³,
/// and f³ is the 12γ term. The 4γ and 12γ exponents arise from the linear and cubic
/// powers of that one coherence factor in this nonlinear observable; they do not identify independent spectral objects.</para>
///
/// <para>Here s is the pure-state concurrence of the initial state. The separate CΨ readout starts at
/// CΨ(0)=s/3, exactly one third of that concurrence. A genuine temporal downward crossing of ¼ occurs
/// iff γ&gt;0 and s&gt;3/4; s=3/4 merely touches ¼ at t=0, and at γ=0 the curve is constant. For every s&gt;0 the cubic 12γ term
/// has relative weight s²/2, growing quadratically (Bell+, s = 1, is the 50/50 member); and every nonzero member
/// has the same late-time 4γ exponential term. The slowing is ours: the log-rate −d ln CΨ/dt is
/// 4γ(w₀ + 3w₁)/(w₀ + w₁) = 4γ(1 + s²) at t = 0 and falls to 4γ, while the coherence factor f decays at a
/// steady 4γ; what slows is the observable, not the carrier. This is a property of the named free two-qubit setup,
/// not an identification with a selected many-body eigenvector. Tier-1 closed form; the Bell+ member reproduces F25
/// (<see cref="InteriorHorizon.BellPlusCpsi"/>) exactly.</para>
///
/// <para>Wired into the typed-knowledge graph (Core) as <c>ApproachFamilyCarrierClaim</c> with two
/// mathematical edges: the Absorption Theorem at n_diff=2 and the exact Bell+ member F25. The C2 3:1
/// resemblance and algebra/dynamics pairing remain prose comparisons rather than typed ancestry.
/// Render it with <c>inspect --claim ApproachFamilyCarrierClaim</c>.</para></summary>
public static class OddHarmonicApproach
{
    /// <summary>The scalar quarter readout CΨ = ¼. The identifier is retained for API stability.</summary>
    public const double Cusp = 0.25;

    /// <summary>The concurrence threshold s = 3/4. For γ&gt;0, a temporal downward crossing occurs
    /// iff s&gt;3/4; equality is only a t=0 touch.</summary>
    public const double CrossingThreshold = 0.75;

    /// <summary>The Bell+ member: s = sin2α = 1 (α = π/4), maximally entangled, 50/50 weights.</summary>
    public const double BellPlusS = 1.0;

    /// <summary>The late-time exponent of every nonzero member, 4γ. At s=0 its coefficient is zero.</summary>
    public static double CarrierRate(double gamma) => 4.0 * gamma;

    /// <summary>The cubic exponent 12γ: f³ in the nonlinear CΨ expression.</summary>
    public static double HarmonicRate(double gamma) => 12.0 * gamma;

    /// <summary>The closed-form weights (w₀ carrier, w₁ harmonic) for s = sin2α:
    /// w₀ = s(1−s²/2)/3, w₁ = s³/6. Thus w₀=w₁=0 at s=0.</summary>
    public static (double Carrier, double Harmonic) Weights(double s) =>
        (s * (1.0 - 0.5 * s * s) / 3.0, s * s * s / 6.0);

    /// <summary>CΨ(s,t) = w₀·e^(−4γt) + w₁·e^(−12γt).</summary>
    public static double Cpsi(double s, double gamma, double t)
    {
        var (w0, w1) = Weights(s);
        return w0 * Math.Exp(-CarrierRate(gamma) * t) + w1 * Math.Exp(-HarmonicRate(gamma) * t);
    }

    /// <summary>The start CΨ(0)=w₀+w₁=s/3: exactly one third of the initial
    /// pure-state concurrence s, and a distinct linear CΨ readout.</summary>
    public static double InitialCpsi(double s) => s / 3.0;

    /// <summary>Whether the curve has a genuine temporal downward crossing of ¼.
    /// This holds iff γ &gt; 0 and CΨ(0)=s/3 &gt; ¼, i.e. iff γ &gt; 0 and s &gt; 3/4.
    /// At γ=0 the curve is constant; at s=3/4 it starts on ¼ but does not cross.</summary>
    public static bool HasDownwardCrossing(double s, double gamma) =>
        gamma > 0.0 && InitialCpsi(s) > Cusp;

    /// <summary>For s &gt; 0, the cubic f³-term fraction at t = 0 is w₁/(w₀+w₁) = s²/2
    /// (the cubic-term share of the start; ½ for Bell+). At s=0 the ratio is 0/0; this method returns
    /// its continuous extension, 0.</summary>
    public static double HarmonicFraction(double s) => 0.5 * s * s;

    /// <summary>The time CΨ(s,·) crosses ¼, by bisection on the monotone-decreasing two-exponential
    /// (both terms positive and decreasing when γ&gt;0). NaN if there is no temporal downward
    /// crossing, including γ=0 or s≤3/4. Throws for γ&lt;0, which is outside the dephasing family.</summary>
    public static double CrossingTime(double s, double gamma)
    {
        if (gamma < 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0 for the dephasing family.");
        if (!HasDownwardCrossing(s, gamma)) return double.NaN;
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
