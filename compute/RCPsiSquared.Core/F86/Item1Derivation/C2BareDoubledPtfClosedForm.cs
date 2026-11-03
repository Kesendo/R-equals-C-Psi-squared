using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 bare-doubled-PTF closed form for the K_b observable (Tier-1-Derived
/// 2026-05-16). Provides the explicit analytical expression that underlies the universal
/// constants <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/> (= 2.196910329331, rounded as
/// 2.196910 in <see cref="C2HwhmRatio"/>) and
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> (= 0.671535517861, rounded as
/// 0.671535 in <see cref="C2HwhmRatio"/>).
///
/// <para><b>Two observables.</b> Those two constants are readings of the FIXED-time slice
/// τ = t/t_peak = 1. The F86 curve of the full block is the envelope max_t |K_b| (see
/// <see cref="Resonance.EnvelopeRefinement"/>), and the bare model's envelope is a different
/// curve: <see cref="EvaluateKb(double, double)"/> gives K_b at any τ in closed form,
/// <see cref="EnvelopePeak"/> reads its peak, x_peak = <see cref="EnvelopeXPeak"/> = 1.7141131714
/// (τ* ≈ 1.30) and HWHM_left/x_peak = <see cref="EnvelopeHwhmLeftOverXPeak"/> = 0.7264892771.
/// A comparison with full-block values takes the envelope pair; a formula defined on the
/// fixed-τ pair keeps it, under that name.</para>
///
/// <para>Setting (γ₀ = 1 dimensionless, t_peak = 1/4 = 1/(4γ₀)):</para>
/// <code>
///   L_2(x) = [[ -2, +2ix ],
///             [ +2ix, -6  ]]               x = Q/Q_EP = Q·g_eff/2
///   ρ_0    = (1, 0)        S = I (probe-block identity)
///   K_b(x) = 2 Re⟨ρ(t_peak)| dρ/dJ⟩
/// </code>
///
/// <para>Eigenvalues split at the EP x=1 into post-EP (complex conjugate pair) and
/// pre-EP (two real). Eigenvector geometry: |v_±|² = 2 universal. Duhamel kernel at
/// t = 1/(4γ₀) reduces via direct sympy algebra (see
/// <c>simulations/f86_doubled_ptf_bare_floor_derivation.py</c>) to two regime-specific
/// closed forms:</para>
///
/// <code>
///   Post-EP (x > 1, ξ = √(x²-1)):
///     K_b(x) = e⁻² · x · [(ξ² + 2)·cos(ξ) − 2] / ξ⁴
///
///   Pre-EP (x &lt; 1, μ = √(1-x²)):
///     K_b(x) = e⁻² · x · [(2 − μ²)·cosh(μ) − 2] / μ⁴
///
///   EP limit (x = 1):
///     K_b(1) = −5·e⁻² / 12
/// </code>
///
/// <para>The post-EP and pre-EP forms are analytic continuations of each other under
/// ξ ↔ iμ (cos → cosh, ξ² → −μ²). Both have the same removable singularity at the EP
/// where (ξ² + 2)·cos(ξ) − 2 ~ −5ξ⁴/12 + 7ξ⁶/180 + O(ξ⁸) — the numerator vanishes at
/// the same order as the denominator. Stable evaluation near the EP uses the Taylor
/// expansion to avoid cancellation.</para>
///
/// <para><b>x_peak and HWHM_left as implicit transcendental equations:</b></para>
///
/// <code>
///   x_peak:  ξ ≥ 1 satisfying dK_b/dξ = 0 (post-EP regime).
///            Numerical solution: ξ_peak = 1.956122438683, x_peak = √(ξ²+1) = 2.196910329331.
///   x_half:  |K_b(x_half)| = |K_b(x_peak)| / 2, with x_half &lt; x_peak.
///            Numerical solution: x_half = 0.721607013629 (PRE-EP regime).
///   ratio:   HWHM_left/x_peak = (x_peak − x_half) / x_peak = 0.671535517861.
/// </code>
///
/// <para><b>The HWHM_left line crosses the EP</b>: x_half &lt; 1 &lt; x_peak. This is
/// the structural reason naïve post-EP-only analyses miss the 0.671535 value — both
/// the post-EP closed form (for K_b at x_peak) AND the pre-EP closed form (for K_b at
/// x_half) are required. Implemented in <see cref="EvaluateKb"/> which dispatches by
/// regime.</para>
///
/// <para><b>No clean algebraic closed form found for the ratio 0.671535517861</b>;
/// candidates tested against include (π−1)/π, 1 − 1/(π + 1/π), e⁻¹ + (1 − e⁻¹)/2,
/// 2/√(2π/e), √(1 − e⁻²), erf(1); closest miss is (π−1)/π = 0.6817 (off by 10⁻²).
/// The transcendental value follows from the implicit equations above and is bit-exact
/// reproducible via the explicit K_b formulas; the rounded <c>0.671535</c> in
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> is the 6-place truncation of this
/// constant.</para>
///
/// <para>Anchors: <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/>,
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/>,
/// <c>docs/superpowers/syntheses/2026-05-06-direction-b-attempt.md</c>,
/// <c>simulations/f86_doubled_ptf_bare_floor_derivation.py</c> (sympy derivation +
/// brute Duhamel verification: closed-form vs eigendecomposition residual ≤ 7·10⁻¹⁷ at
/// 7 sample x values spanning both regimes).</para>
///
/// <para><b>Out of scope</b>: the per-bond-class lift from the bare floor to the full-block
/// values. On the envelope the floor is 0.7264892771 and the c = 2 class values run from 0.7523
/// to 0.7621 (Interior) and from 0.7766 to 0.7815 (Endpoint) over N = 5..12, a lift of +0.026 to
/// +0.036 and +0.050 to +0.055. No closed form for it is held; <see cref="F86HwhmClosedFormClaim"/>
/// tracks the fitted (α, β) description.</para></summary>
public sealed class C2BareDoubledPtfClosedForm : Claim
{
    /// <summary>x_peak to 12 digits of precision: ξ_peak² + 1 with ξ_peak = 1.956122438683
    /// satisfying dK_b/dξ = 0. The 6-place rounding stored in
    /// <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/> is the truncation of this value.</summary>
    public const double XPeakPrecise = 2.196910329331;

    /// <summary>HWHM_left/x_peak to 12 digits of precision. The 6-place rounding stored
    /// in <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> is the truncation of this
    /// value. No clean algebraic closed form among standard π/e/√n candidates matches.</summary>
    public const double HwhmLeftOverXPeakPrecise = 0.671535517861;

    /// <summary>x_peak of the ENVELOPE max_τ |K_b(x, τ)|, the observable the full-block F86 curve
    /// reads (<see cref="EnvelopePeak"/> recomputes it; resolved to about √ε relative, since the
    /// maximum is flat to second order). It sits below <see cref="XPeakPrecise"/>, the fixed-τ = 1
    /// slice's peak, and is reached at τ* ≈ 1.30.</summary>
    public const double EnvelopeXPeak = 1.7141131714;

    /// <summary>HWHM_left/x_peak of the envelope (x_half = 0.468828332730, pre-EP like the
    /// fixed-τ x_half): the bare floor on the observable the full-block curve reads.</summary>
    public const double EnvelopeHwhmLeftOverXPeak = 0.7264892771;

    /// <summary>EP limit value: K_b(x=1) = −5·e⁻²/12. The removable-singularity limit of
    /// both regime closed forms at the EP.</summary>
    public static readonly double KbAtEp = -5.0 * Math.Exp(-2.0) / 12.0;

    /// <summary>Evaluate K_b(x) via the appropriate regime closed form. Stable across
    /// the EP via Taylor expansion when |ξ| or |μ| &lt; 0.05.</summary>
    public static double EvaluateKb(double x)
    {
        if (x <= 0) throw new ArgumentOutOfRangeException(nameof(x), x, "x must be > 0.");
        const double eMinus2 = 0.1353352832366127; // = e^(-2)
        const double small = 0.05;

        if (x > 1.0)
        {
            double xi = Math.Sqrt(x * x - 1.0);
            double bracket;
            if (xi < small)
            {
                double xi2 = xi * xi;
                double xi4 = xi2 * xi2;
                double xi6 = xi4 * xi2;
                bracket = -5.0 * xi4 / 12.0 + 7.0 * xi6 / 180.0;
            }
            else
            {
                bracket = (xi * xi + 2.0) * Math.Cos(xi) - 2.0;
            }
            double xi_4 = xi * xi * xi * xi;
            return eMinus2 * x * bracket / xi_4;
        }

        if (x < 1.0)
        {
            double mu = Math.Sqrt(1.0 - x * x);
            double bracket;
            if (mu < small)
            {
                double mu2 = mu * mu;
                double mu4 = mu2 * mu2;
                double mu6 = mu4 * mu2;
                bracket = -5.0 * mu4 / 12.0 - 7.0 * mu6 / 180.0;
            }
            else
            {
                bracket = (2.0 - mu * mu) * Math.Cosh(mu) - 2.0;
            }
            double mu_4 = mu * mu * mu * mu;
            return eMinus2 * x * bracket / mu_4;
        }

        return KbAtEp; // x == 1
    }

    /// <summary>K_b(x, τ) at any time τ = t/t_peak (τ = 1 is <see cref="EvaluateKb(double)"/>):
    /// <code>
    ///   post-EP: K_b = x·e^(−2τ)·[(τξ² + 2)·cos(τξ) + (τ − 1)·ξ·sin(τξ) − 2] / ξ⁴,   ξ = √(x² − 1)
    ///   pre-EP:  the same with ξ → iμ, μ = √(1 − x²)
    /// </code>
    /// The bracket is one entire series in s = x² − 1 = ξ², K_b = x·e^(−2τ)·Σ_(m≥2) c_m(τ)·s^(m−2) with
    /// c_m = (−1)^(m−1)·[τ^(2m−1)/(2m−2)! + (τ − 1)·τ^(2m−1)/(2m−1)! − 2τ^(2m)/(2m)!]
    /// (c_2 = −τ³/3 − τ⁴/12, which is −5/12 at τ = 1); the series is used where |s|·τ² &lt; 1,
    /// where the closed form cancels, and the closed form elsewhere.</summary>
    public static double EvaluateKb(double x, double tau)
    {
        if (x <= 0) throw new ArgumentOutOfRangeException(nameof(x), x, "x must be > 0.");
        if (tau <= 0) throw new ArgumentOutOfRangeException(nameof(tau), tau, "tau must be > 0.");
        double s = x * x - 1.0;
        double pref = x * Math.Exp(-2.0 * tau);

        if (Math.Abs(s) * tau * tau < 1.0)
        {
            // Σ_(m≥2) c_m s^(m−2), terms falling like (s τ²)^m / (2m)!.
            double sum = 0.0, sPow = 1.0;
            double tPow = tau * tau * tau;          // τ^(2m−1) at m = 2
            double f2m2 = 2.0, f2m1 = 6.0, f2m = 24.0; // (2m−2)!, (2m−1)!, (2m)! at m = 2
            double sign = -1.0;                     // (−1)^(m−1) at m = 2
            for (int m = 2; m <= 40; m++)
            {
                double cm = sign * (tPow / f2m2 + (tau - 1.0) * tPow / f2m1 - 2.0 * tPow * tau / f2m);
                double term = cm * sPow;
                sum += term;
                if (m > 4 && Math.Abs(term) <= 1e-17 * Math.Abs(sum)) break;
                sPow *= s;
                tPow *= tau * tau;
                f2m2 = f2m * 1.0; f2m1 = f2m * (2 * m + 1); f2m = f2m1 * (2 * m + 2);
                sign = -sign;
            }
            return pref * sum;
        }

        if (s > 0)
        {
            double xi = Math.Sqrt(s);
            double bracket = (tau * s + 2.0) * Math.Cos(tau * xi) + (tau - 1.0) * xi * Math.Sin(tau * xi) - 2.0;
            return pref * bracket / (s * s);
        }
        else
        {
            double mu = Math.Sqrt(-s);
            double bracket = (2.0 + tau * s) * Math.Cosh(tau * mu) - (tau - 1.0) * mu * Math.Sinh(tau * mu) - 2.0;
            return pref * bracket / (s * s);
        }
    }

    /// <summary>Seed grid in τ = t/t_peak for the bare envelope: 160 points on [0.05, 8]. The
    /// maximiser τ*(x) runs from about 2 at small x down to about 0.6 at x = 5 (1.30 at the
    /// envelope peak), well inside.</summary>
    private static readonly double[] TauSeedGrid =
        Enumerable.Range(0, 160).Select(i => 0.05 + (8.0 - 0.05) * i / 159.0).ToArray();

    /// <summary>The envelope of the bare model, max over τ of |K_b(x, τ)|, and its maximiser τ*:
    /// the same observable as the full-block F86 curve (<see cref="Resonance.EnvelopeRefinement"/>),
    /// where <see cref="EvaluateKb(double)"/> is the fixed-τ = 1 slice.</summary>
    public static (double KMax, double TauStar) EnvelopeKb(double x)
    {
        int iWin = 0; double kWin = -1.0;
        for (int i = 0; i < TauSeedGrid.Length; i++)
        {
            double k = Math.Abs(EvaluateKb(x, TauSeedGrid[i]));
            if (k > kWin) { kWin = k; iWin = i; }
        }
        var (kMax, tauStar, _) = Resonance.EnvelopeRefinement.Refine(TauSeedGrid, iWin, kWin,
            tau => Math.Abs(EvaluateKb(x, tau)));
        return (kMax, tauStar);
    }

    /// <summary>Peak and left half-maximum of the bare envelope in x: x_peak by golden section
    /// on <see cref="EnvelopeKb"/> seeded on [0.05, 6], x_half by bisection below it. x_peak sits
    /// on a maximum that is flat to second order, so it is resolved only to about √ε relative
    /// (the value KMax to ε); x_half sits on a flank and is resolved to the bisection's end.</summary>
    public static (double XPeak, double KMax, double XHalf, double HwhmLeftOverXPeak) EnvelopePeak()
    {
        const int seeds = 120;
        double lo = 0.05, hi = 6.0;
        int iBest = 0; double kBest = -1.0;
        var xs = new double[seeds];
        for (int i = 0; i < seeds; i++)
        {
            xs[i] = lo + (hi - lo) * i / (seeds - 1);
            double k = EnvelopeKb(xs[i]).KMax;
            if (k > kBest) { kBest = k; iBest = i; }
        }
        double a = xs[Math.Max(iBest - 1, 0)], b = xs[Math.Min(iBest + 1, seeds - 1)];
        double invPhi = (Math.Sqrt(5.0) - 1.0) / 2.0;
        double c = b - invPhi * (b - a), d = a + invPhi * (b - a);
        double fc = EnvelopeKb(c).KMax, fd = EnvelopeKb(d).KMax;
        while (b - a > 1e-12)
        {
            if (fc >= fd) { b = d; d = c; fd = fc; c = b - invPhi * (b - a); fc = EnvelopeKb(c).KMax; }
            else { a = c; c = d; fc = fd; d = a + invPhi * (b - a); fd = EnvelopeKb(d).KMax; }
        }
        double xPeak = fc >= fd ? c : d, kMax = Math.Max(fc, fd);

        double half = kMax / 2.0, xl = lo, xr = xPeak;
        if (EnvelopeKb(xl).KMax >= half)
            throw new InvalidOperationException("bare envelope does not fall to half maximum above x = 0.05");
        while (xr - xl > 1e-14)
        {
            double xm = 0.5 * (xl + xr);
            if (EnvelopeKb(xm).KMax < half) xl = xm; else xr = xm;
        }
        double xHalf = 0.5 * (xl + xr);
        return (xPeak, kMax, xHalf, (xPeak - xHalf) / xPeak);
    }

    public C2BareDoubledPtfClosedForm()
        : base("F86 c=2 bare-doubled-PTF closed form: K_b(x) explicit analytical expression in two regimes (post-EP cos, pre-EP cosh, EP limit −5·e⁻²/12); underlies the universal constants x_peak = 2.196910329331 (= C2HwhmRatio.BareDoubledPtfXPeak) and HWHM_left/x_peak = 0.671535517861 (= C2HwhmRatio.BareDoubledPtfHwhmRatio rounded to 6 places); x_half = 0.721607 sits in pre-EP regime so both closed forms are required for the HWHM evaluation",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (universal constants), docs/superpowers/syntheses/2026-05-06-direction-b-attempt.md (model setup), simulations/f86_doubled_ptf_bare_floor_derivation.py (sympy derivation + brute Duhamel verification; residual ≤ 7e-17 at 7 sample x; transcendental constants reproducible via the explicit closed forms; these are the fixed-τ = 1 slice, the envelope pair 1.7141131714 / 0.7264892771 is computed here by EnvelopePeak), docs/proofs/PROOF_F86B_OBSTRUCTION.md (no closed form for the per-bond-class lift)")
    { }

    public override string DisplayName => "F86 c=2 bare-doubled-PTF K_b closed form (post-EP cos + pre-EP cosh)";

    public override string Summary =>
        $"K_b(x) = e⁻²·x·[(ξ²+2)·cos(ξ) − 2]/ξ⁴ post-EP (ξ=√(x²−1)); = e⁻²·x·[(2−μ²)·cosh(μ) − 2]/μ⁴ pre-EP (μ=√(1−x²)); EP: K_b(1)=−5e⁻²/12. x_peak={XPeakPrecise}, HWHM/x_peak={HwhmLeftOverXPeakPrecise} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Closed form post-EP",
                summary: "K_b(x) = e⁻² · x · [(ξ² + 2)·cos(ξ) − 2] / ξ⁴ with ξ = √(x²−1), valid for x > 1");
            yield return new InspectableNode("Closed form pre-EP",
                summary: "K_b(x) = e⁻² · x · [(2 − μ²)·cosh(μ) − 2] / μ⁴ with μ = √(1−x²), valid for x ∈ (0, 1)");
            yield return new InspectableNode("EP limit",
                summary: "K_b(1) = −5·e⁻² / 12 ≈ −0.0563897 (removable singularity, common limit of both regimes)");
            yield return new InspectableNode("x_peak",
                summary: $"{XPeakPrecise:F12} (ξ_peak = 1.956122438683 from dK/dξ = 0 transcendental)");
            yield return new InspectableNode("x_half (HWHM crossing)",
                summary: "0.721607013629 (pre-EP! HWHM_left line crosses the EP)");
            yield return new InspectableNode("HWHM_left / x_peak",
                summary: $"{HwhmLeftOverXPeakPrecise:F12} (transcendental, no clean algebraic match)");
            yield return new InspectableNode("Envelope (the full-block observable)",
                summary: $"x_peak = {EnvelopeXPeak} at τ* ≈ 1.30, HWHM_left/x_peak = {EnvelopeHwhmLeftOverXPeak}; the fixed-τ = 1 values above are the slice");
            yield return new InspectableNode("Lift out of scope",
                summary: "Envelope floor 0.7264892771 to the c = 2 class values (Interior 0.7523..0.7621, Endpoint 0.7766..0.7815, N = 5..12): +0.026..+0.036 and +0.050..+0.055, no closed form; F86HwhmClosedFormClaim holds the fitted (α, β) description");
            yield return new InspectableNode("Verification",
                summary: "Bit-exact match (residual ≤ 7e-17) against brute Duhamel evaluator at x ∈ {0.5, 0.722, 1.0, 1.5, 2.197, 3.0, 5.0} spanning both regimes; see C2BareDoubledPtfClosedFormTests");
        }
    }
}
