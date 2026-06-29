using System;
using System.Collections.Generic;
using System.Linq;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>The honest verdict on a per-sub-class HWHM-lift slope alpha (f86b2_robust_extraction
/// arc, landed 2026-06-29). The slope is a least-squares fit of lift = HWHM_ratio - 0.671535
/// against g_eff over the bond sub-class's anchors; whether it is a real structural slope or an
/// extraction-noise artefact is the question this verdict answers.</summary>
public enum AlphaVerdict
{
    /// <summary>Multi-anchor, real g_eff lever, the slope explains the lift trend and is
    /// jackknife-tight. Marginal at the 0.005 grid floor; promotion to Tier-1-derived still
    /// awaits a grid-convergence study.</summary>
    ResolvedSlope,

    /// <summary>Multi-anchor but a flat lift over a microscopic g_eff lever, so the slope fits
    /// grid noise. Collapses to a per-class constant lift.</summary>
    ConstantNoise,

    /// <summary>One anchor only: no slope is computable. Collapses to a constant lift.</summary>
    SinglePoint,

    /// <summary>Anchors past the default [0.20, 4.00] Q-grid bound (Q_peak > 4), flagged
    /// "grid-snap artefact, not physical" by C2HwhmRatio. Not a structural slope; its fit is
    /// kept only to reproduce the (non-physical) anchors.</summary>
    EscapeArtefact,
}

/// <summary>One sub-class's live extraction diagnostics (recomputed from the anchors).</summary>
public readonly record struct SubClassExtraction(
    BondSubClass SubClass,
    int N,
    int NEscape,
    double GSpan,
    double FittedAlpha,
    double FittedBeta,
    double SigmaAlpha,
    double LiftMean,
    AlphaVerdict Verdict);

/// <summary>The live, honest f86b2 alpha-extraction: recomputes the per-BondSubClass slope alpha
/// and its grid-noise uncertainty from the 22 F90-bridge anchors, and applies the Phase-1 honest
/// reshape (the f86b2_robust_extraction arc, "honest and phased" stance, 2026-06-29).
///
/// <para>Of the six sub-classes only two (Endpoint, Flanking) rest on a real g_eff lever with
/// several physical anchors; their slopes are jackknife-tight but marginal at the 0.005 grid
/// floor. The other four are not defensible: Mid fits a flat lift over a microscopic lever
/// (g_span ~ 0.011), two classes have a single anchor, and Orbit2Escape sits entirely on
/// grid-edge anchors C2HwhmRatio flags non-physical. <see cref="ReshapedParams"/> keeps the two
/// resolved slopes (and the flagged escape fit, for prediction fidelity) and collapses the
/// noise/single-point classes to a per-class constant lift (alpha = 0).</para>
///
/// <para>This replaces the hand-transcribed constants in <see cref="F86HwhmClosedFormClaim"/>
/// with a live recompute: the same anchors, the same least-squares fit, recomputed on every call,
/// so the values cannot drift out of sync with their source.</para></summary>
public static class F86HwhmAlphaExtraction
{
    public const double BareFloor = 0.671535;
    public const double BareDoubledPtfXPeak = 4.39382;
    public const double EscapeQ = 4.0;        // Q_peak past the default [0.20, 4.00] grid bound
    public const double GridSigmaY = 0.005;   // documented Q-grid wobble on HWHM_ratio (inventory B3)

    /// <summary>The 22 F90-bridge anchors (N, bond, Q_peak in F86-J, HWHM_left/Q_peak), N=5..8.
    /// Same table as simulations/f86_hwhm_closed_form_verification.py EMPIRICAL_ANCHORS, sourced
    /// from the PROOF_F90_F86C2_BRIDGE.md verification table.</summary>
    public static readonly (int N, int B, double QPeak, double Hwhm)[] Anchors =
    {
        (5, 0, 2.40, 0.7700), (5, 1, 1.50, 0.7454), (5, 2, 1.50, 0.7454), (5, 3, 2.40, 0.7700),
        (6, 0, 2.52, 0.7737), (6, 1, 1.65, 0.7503), (6, 2, 1.44, 0.7449),
        (6, 3, 1.65, 0.7503), (6, 4, 2.52, 0.7737),
        (7, 0, 2.53, 0.7738), (7, 1, 7.27, 0.9162), (7, 2, 1.54, 0.7469),
        (7, 3, 1.54, 0.7469), (7, 4, 7.27, 0.9162), (7, 5, 2.53, 0.7738),
        (8, 0, 2.53, 0.7734), (8, 1, 8.07, 0.8899), (8, 2, 1.51, 0.7475),
        (8, 3, 16.79, 0.5778), (8, 4, 1.51, 0.7475), (8, 5, 8.07, 0.8899), (8, 6, 2.53, 0.7734),
    };

    // The Phase-1 honest verdict per sub-class. It is a documented judgement, not a magic
    // threshold: the live diagnostics (SigmaAlpha, GSpan, N, NEscape) recomputed in Extract()
    // support it, and the witness surfaces both so a reader can see the verdict match the data.
    private static readonly Dictionary<BondSubClass, AlphaVerdict> _verdict = new()
    {
        [BondSubClass.Endpoint] = AlphaVerdict.ResolvedSlope,
        [BondSubClass.Flanking] = AlphaVerdict.ResolvedSlope,
        [BondSubClass.Mid] = AlphaVerdict.ConstantNoise,
        [BondSubClass.CentralSelfPaired] = AlphaVerdict.SinglePoint,
        [BondSubClass.Orbit2Escape] = AlphaVerdict.EscapeArtefact,
        [BondSubClass.CentralEscapeOrbit3] = AlphaVerdict.SinglePoint,
    };

    public static double GEff(double qPeak) => BareDoubledPtfXPeak / (qPeak + 2.0);

    /// <summary>Recompute every sub-class's slope diagnostics live from the anchors.</summary>
    public static IReadOnlyDictionary<BondSubClass, SubClassExtraction> Extract()
    {
        var groups = new Dictionary<BondSubClass, List<(double G, double Lift, double Q)>>();
        foreach (var (n, b, q, hwhm) in Anchors)
        {
            var sub = BondSubClassExtensions.Classify(n, b, q);
            if (!groups.TryGetValue(sub, out var list)) groups[sub] = list = new();
            list.Add((GEff(q), hwhm - BareFloor, q));
        }

        var result = new Dictionary<BondSubClass, SubClassExtraction>();
        foreach (var (sub, rows) in groups)
        {
            int n = rows.Count;
            int nEsc = rows.Count(r => r.Q > EscapeQ);
            double[] gs = rows.Select(r => r.G).ToArray();
            double[] lifts = rows.Select(r => r.Lift).ToArray();
            (double alpha, double beta) = n >= 2
                ? LeastSquares(gs, lifts)
                : (gs[0] != 0 ? lifts[0] / gs[0] : 0.0, 0.0);
            result[sub] = new SubClassExtraction(
                sub, n, nEsc, gs.Max() - gs.Min(), alpha, beta, SlopeSigma(gs), lifts.Average(),
                _verdict[sub]);
        }
        return result;
    }

    /// <summary>The honest reshaped (alpha, beta) the F86b2 claim uses. A resolved slope keeps
    /// its fit; an escape-artefact keeps its (flagged) fit for fidelity to its two lift levels;
    /// a noise/single-point class collapses to a per-class constant lift (alpha = 0, beta = mean
    /// lift), which still reproduces its anchors within the 0.005 grid floor.</summary>
    public static (double Alpha, double Beta) ReshapedParams(BondSubClass sub)
    {
        var ex = Extract()[sub];
        return ex.Verdict switch
        {
            AlphaVerdict.ResolvedSlope => (ex.FittedAlpha, ex.FittedBeta),
            AlphaVerdict.EscapeArtefact => (ex.FittedAlpha, ex.FittedBeta),
            _ => (0.0, ex.LiftMean),
        };
    }

    private static (double Alpha, double Beta) LeastSquares(double[] x, double[] y)
    {
        double xbar = x.Average(), ybar = y.Average();
        double sxx = 0, sxy = 0;
        for (int i = 0; i < x.Length; i++)
        {
            sxx += (x[i] - xbar) * (x[i] - xbar);
            sxy += (x[i] - xbar) * (y[i] - ybar);
        }
        double alpha = sxx == 0 ? 0.0 : sxy / sxx;
        return (alpha, ybar - alpha * xbar);
    }

    private static double SlopeSigma(double[] gs)
    {
        if (gs.Length < 2) return double.PositiveInfinity;
        double gbar = gs.Average();
        double sxx = gs.Sum(g => (g - gbar) * (g - gbar));
        return sxx == 0 ? double.PositiveInfinity : GridSigmaY / Math.Sqrt(sxx);
    }
}
