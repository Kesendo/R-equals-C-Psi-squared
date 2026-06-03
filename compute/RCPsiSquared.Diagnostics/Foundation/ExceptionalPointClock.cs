using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The exceptional point read through the clock: the slow-pair effective Liouvillian
/// L_eff = [[−2γ₀, iJ·g_eff],[iJ·g_eff, −6γ₀]] (k=1, J = Q·γ₀) has eigenvalues λ_± = −4γ₀ ± √(4γ₀² −
/// J²g_eff²), two real below the EP, a defective double root −4γ₀ at the EP (Q_EP = 2/g_eff), a complex
/// pair above it. Read as the clock we built (Takt = the radial decay −Re λ, Rotation = the angle
/// arctan(|Im|/|Re|), the F95 angle): the EP is exactly where the Rotation hand lifts off the Takt axis,
/// the rotation born. Closed-form, reusing <see cref="EpAlgebra"/> (the eigenvalues, Q_EP) and
/// <see cref="C2BareDoubledPtfClosedForm"/> (the K_b resonance). N-free (the EP lives in the effective
/// 2-level model, set by g_eff and Q).
///
/// <para>The eigenvector overlap |⟨v₊|v₋⟩| = min(x, 1/x), x = Q/Q_EP: exactly 1 at the EP (the two
/// eigenvectors coalesce, defective, a Jordan block, the Petermann sensitivity spiking), falling off as
/// the modes re-separate on either side. This is the defective signature that makes it an exceptional
/// point and not a mere degeneracy. The cusp (interior axis) is the mirror event at the same F95 zero:
/// there the rotation stills (θ → 0 descending to ¼), here it is born (θ lifts off ascending past Q_EP).
/// The post-EP regime this opens is the BirthCanal of <see cref="PostEpFlowField"/>; hardware-anchored by
/// the IBM Kingston EP-onset run.</para></summary>
public static class ExceptionalPointClock
{
    /// <summary>The resonance peak in x = Q/Q_EP units (the K_b maximum, past the EP).</summary>
    public const double XPeak = C2BareDoubledPtfClosedForm.XPeakPrecise;

    /// <summary>The exceptional point Q_EP = 2/g_eff (where the two slow modes coalesce defectively).</summary>
    public static double QEp(double gEff) => EpAlgebra.QEp(gEff);

    /// <summary>The resonance peak Q_peak = x_peak·Q_EP (where K_b peaks, in the post-EP regime).</summary>
    public static double QPeak(double gEff) => XPeak * QEp(gEff);

    /// <summary>x = Q/Q_EP = Q·g_eff/2, the dimensionless distance past the EP (x=1 is the EP).</summary>
    public static double X(double q, double gEff) => q * gEff / 2.0;

    /// <summary>The slow mode's decay −Re(λ_slow), via <see cref="EpAlgebra"/>: rises from 2γ₀ to 4γ₀ as
    /// Q → Q_EP, then pins at 4γ₀ above the EP (the Takt hand stops moving).</summary>
    public static double Decay(double gamma0, double q, double gEff)
    {
        var (lamP, lamM) = EpAlgebra.SlowestPairEigenvaluesComplex(gamma0, q * gamma0, gEff);
        return Math.Min(-lamP.Real, -lamM.Real);
    }

    /// <summary>The slow mode's oscillation ω = |Im(λ_slow)|: 0 below the EP, growing as 2γ₀√(x²−1)
    /// above it (the rotation born).</summary>
    public static double Omega(double gamma0, double q, double gEff)
    {
        var (lamP, lamM) = EpAlgebra.SlowestPairEigenvaluesComplex(gamma0, q * gamma0, gEff);
        // the slow mode is the one with the smaller decay; both share |Im| above the EP
        return -lamP.Real <= -lamM.Real ? Math.Abs(lamP.Imaginary) : Math.Abs(lamM.Imaginary);
    }

    /// <summary>The Rotation hand angle θ = arctan(ω / decay) in degrees (the F95 angle): 0 below the EP
    /// (pure Takt), lifting off above it.</summary>
    public static double RotationAngleDegrees(double gamma0, double q, double gEff)
    {
        double decay = Decay(gamma0, q, gEff);
        double omega = Omega(gamma0, q, gEff);
        if (decay <= 0.0) return omega > 0.0 ? 90.0 : 0.0;
        return Math.Atan2(omega, decay) * 180.0 / Math.PI;
    }

    /// <summary>The normalized eigenvector overlap |⟨v₊|v₋⟩| = min(x, 1/x), x = Q/Q_EP: 1 at the EP (the
    /// eigenvectors coalesce, defective), falling off on either side. The defective signature.</summary>
    public static double EigenvectorOverlap(double q, double gEff)
    {
        double x = X(q, gEff);
        return x <= 1.0 ? x : 1.0 / x;
    }

    /// <summary>The K_b susceptibility resonance at Q (via <see cref="C2BareDoubledPtfClosedForm"/>),
    /// peaking at Q_peak past the EP.</summary>
    public static double Kb(double q, double gEff) => C2BareDoubledPtfClosedForm.EvaluateKb(X(q, gEff));
}
