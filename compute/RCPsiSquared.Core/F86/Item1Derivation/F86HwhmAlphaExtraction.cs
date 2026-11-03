using System;
using System.Collections.Generic;
using System.Linq;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>What the anchors of one bond sub-class can say about a slope alpha in the form
/// lift = HWHM_ratio − floor = alpha · g_eff + beta: what decides is how many distinct g_eff
/// values a class has.</summary>
public enum AlphaVerdict
{
    /// <summary>Three or more distinct g_eff values: the line is a least-squares fit and its
    /// residual is a measurement of how far the class is from linear in g_eff.</summary>
    Overdetermined,

    /// <summary>Exactly two distinct g_eff values (two chain lengths): the line passes through
    /// both, so the slope is exact and untested.</summary>
    TwoPoint,

    /// <summary>One g_eff value: no slope; the class carries a constant lift.</summary>
    SinglePoint,
}

/// <summary>One sub-class's live extraction diagnostics (recomputed from the anchors).
/// <see cref="N"/> counts bonds, mirror partners included; <see cref="DistinctG"/> counts the
/// distinct g_eff values the fit actually rests on.</summary>
public readonly record struct SubClassExtraction(
    BondSubClass SubClass,
    int N,
    int DistinctG,
    int NEscape,
    double GSpan,
    double FittedAlpha,
    double FittedBeta,
    double MaxResidual,
    double LiftMean,
    AlphaVerdict Verdict);

/// <summary>The live f86b2 alpha extraction: per BondSubClass, the line lift = alpha · g_eff +
/// beta through the c=2 anchors N=5..8, with lift = HWHM_left/Q_peak minus the bare floor.
///
/// <para>Both sides are envelope readings: the anchors are the canonical peaks of
/// max_t |K_b(Q, t)| per bond (the first local maximum, the <see cref="PerF71OrbitObservation"/>
/// values, one per F71 orbit, expanded to both mirror bonds), and the floor is the bare
/// doubled-PTF envelope ratio <see cref="C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak"/>.
/// g_eff is this fit's own coordinate, 4.39382 / (Q_peak + 2), built on the fixed-τ bare peak
/// (2 · 2.196910), so alpha is in fixed-τ units; the envelope's 2 · 1.71411 would rescale alpha
/// and leave every verdict below unchanged.</para>
///
/// <para>What the anchors say. Endpoint and Flanking have four chain lengths each (N=5..8) and
/// neither lies on a line: the largest residuals are 1.28·10⁻³ and 4.2·10⁻⁴. The anchors are
/// grid readings, off the continuous envelope peak by at most 3.5·10⁻⁵ in the ratio; with the
/// continuous values (simulations/f86_envelope_second_route.py --continuous N) the residuals are 1.278·10⁻³ and 4.03·10⁻⁴, so they stand 12 to 36 times
/// above that offset. With four points and two parameters this refutes an exact linear law; it is
/// not a statement about fit quality. Mid (N=7, 8) and CentralSelfPaired (N=6, 8) have two points,
/// a line through them that tests nothing. On the canonical reading no bond peaks past Q = 4, so
/// Orbit2Escape and CentralEscapeOrbit3, which <see cref="BondSubClassExtensions.Classify"/> gives
/// a far maximum (the global reading of <see cref="PerF71OrbitObservation.Witnesses"/>' plateaus),
/// have no anchors here. The linear-in-g_eff form is refuted as an exact law in both classes that
/// can test it.</para></summary>
public static class F86HwhmAlphaExtraction
{
    /// <summary>The bare doubled-PTF envelope floor, the full-block lift's reference.</summary>
    public const double BareFloor = C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak;

    /// <summary>2 · the fixed-τ = 1 bare peak 2.196910: the predictor-defined g_eff scale.</summary>
    public const double BareDoubledPtfXPeak = 4.39382;

    public const double EscapeQ = 4.0;        // Q_peak past the default [0.20, 4.00] grid bound

    /// <summary>The c=2 anchors, one per F71 orbit (N, lower bond of the orbit, canonical Q_peak in
    /// F86-J, HWHM_left/Q_peak), envelope readings at γ₀ = 0.05 on the default Q grid. Second route:
    /// simulations/f86_envelope_second_route.py.</summary>
    public static readonly (int N, int B, double QPeak, double Hwhm)[] OrbitAnchors =
    {
        (5, 0, 2.517071370614458, 0.7765595000527978), (5, 1, 1.4831232668717471, 0.7523433331317328),
        (6, 0, 2.580604998574177, 0.7815203313095082), (6, 1, 1.6366804166405813, 0.7596076842042161),
        (6, 2, 1.4283108372718716, 0.7504884589345521),
        (7, 0, 2.5600734278105572, 0.7813218454972876), (7, 1, 1.698919918615899, 0.7631649971739308),
        (7, 2, 1.5010592087425514, 0.754635079880418),
        (8, 0, 2.543127121124538, 0.7807879816312349), (8, 1, 1.727187504810217, 0.7644478770686235),
        (8, 2, 1.5207477539505165, 0.7562576392988031), (8, 3, 1.5698435864972127, 0.7584673713775594),
    };


    /// <summary>Every bond's anchor: each orbit row expanded to both bonds of its mirror pair.</summary>
    public static IReadOnlyList<(int N, int B, double QPeak, double Hwhm)> Anchors { get; } =
        OrbitAnchors
            .SelectMany(a => a.B == a.N - 2 - a.B
                ? new[] { a }
                : new[] { a, (a.N, a.N - 2 - a.B, a.QPeak, a.Hwhm) })
            .ToArray();

    public static double GEff(double qPeak) => BareDoubledPtfXPeak / (qPeak + 2.0);

    /// <summary>Recompute every sub-class's line and its diagnostics live from the anchors.</summary>
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
            double[] gs = rows.Select(r => r.G).ToArray();
            double[] lifts = rows.Select(r => r.Lift).ToArray();
            // Mirror partners share one literal, so their g_eff values are equal exactly.
            int distinct = gs.Distinct().Count();
            (double alpha, double beta) = distinct >= 2 ? LeastSquares(gs, lifts) : (0.0, lifts.Average());
            double maxRes = rows.Max(r => Math.Abs(r.Lift - (alpha * r.G + beta)));
            var verdict = distinct switch
            {
                >= 3 => AlphaVerdict.Overdetermined,
                2 => AlphaVerdict.TwoPoint,
                _ => AlphaVerdict.SinglePoint,
            };
            result[sub] = new SubClassExtraction(
                sub, rows.Count, distinct, rows.Count(r => r.Q > EscapeQ), gs.Max() - gs.Min(),
                alpha, beta, maxRes, lifts.Average(), verdict);
        }
        return result;
    }

    /// <summary>The (alpha, beta) the F86b2 claim uses: the class's line (the least-squares fit,
    /// exact through two points), or a constant lift for a single point. Throws for a sub-class
    /// with no anchors (the two far-maximum classes on the canonical reading).</summary>
    public static (double Alpha, double Beta) FitParams(BondSubClass sub)
    {
        if (!Extract().TryGetValue(sub, out var ex))
            throw new InvalidOperationException(
                $"{sub} has no anchors: no canonical peak lies past Q = 4, so no anchor classifies as {sub}.");
        return (ex.FittedAlpha, ex.FittedBeta);
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
}
