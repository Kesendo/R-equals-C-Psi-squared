using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 2 universal-shape function values f_class(x) at fixed x points,
/// computed live from the shared <see cref="WitnessCache"/>. For each x ∈
/// {−0.6, −0.4, −0.2, 0, +0.2, +0.4, +1.0} we evaluate y(x) = K(Q*(1+x))/|K|max across the
/// configured (c, N) sample, returning the min/max range.
///
/// <para>The shape collapse manifests as small spreads near the peak (under ~2 % at |x| ≤ 0.4).
/// At x = +1.0 the envelope gives, on a [0.2, 6] grid of step 0.025, c=2 N=5..7 Interior
/// 0.837 / 0.843 / 0.884 and Endpoint
/// 0.864 / 0.874 / 0.876: Endpoint lies above Interior at N=5 and 6 and below it at N=7, where
/// the Interior class mean carries the inner non-central bonds rising toward their plateau
/// (see <see cref="PerF71OrbitObservation"/>), so the x = +1.0 value does not separate the
/// classes at every N. A target Q*(1+x) beyond the curve's Q grid reads as 0: on the default
/// grid [0.2, 4] the Endpoint point x = +1.0 (Q ≈ 5) lies outside it.</para>
/// </summary>
public sealed class ShapeFunctionPoint
{
    public double X { get; }
    public double YMin { get; }
    public double YMax { get; }

    public ShapeFunctionPoint(double x, double yMin, double yMax)
    {
        X = x;
        YMin = yMin;
        YMax = yMax;
    }

    public double YMean => 0.5 * (YMin + YMax);
    public double Spread => YMax - YMin;
    public double SpreadPercent => YMean > 0 ? 100.0 * Spread / YMean : 0;
}

public sealed class ShapeFunctionWitnesses : Claim
{
    public BondClass BondClass { get; }
    public double GammaZero { get; }
    public IReadOnlyList<(int C, int N)> Locations { get; }
    public IReadOnlyList<double> XPoints { get; }
    public WitnessCache Cache { get; }

    private readonly Lazy<IReadOnlyList<ShapeFunctionPoint>> _points;

    public IReadOnlyList<ShapeFunctionPoint> Points => _points.Value;

    /// <summary>The y(x=+1.0) class-specific plateau — the Interior / Endpoint discriminator
    /// in the post-peak tail. Computed.</summary>
    public double PostPeakPlateau => Points.First(p => Math.Abs(p.X - 1.0) < 1e-9).YMean;

    public ShapeFunctionWitnesses(BondClass bondClass, double gammaZero,
        IReadOnlyList<(int C, int N)> locations, IReadOnlyList<double>? xPoints = null,
        WitnessCache? cache = null)
        : base($"f_class(x) shape values ({bondClass})",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F86 universal-shape table + docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md Statement 2")
    {
        BondClass = bondClass;
        GammaZero = gammaZero;
        Locations = locations;
        XPoints = xPoints ?? new[] { -0.60, -0.40, -0.20, 0.00, +0.20, +0.40, +1.00 };
        Cache = cache ?? WitnessCache.Default;
        _points = new Lazy<IReadOnlyList<ShapeFunctionPoint>>(ComputePoints);
    }

    public static ShapeFunctionWitnesses BuildInterior(double gammaZero = 0.05, WitnessCache? cache = null) =>
        new(BondClass.Interior, gammaZero, F86StandardLocations.Full, xPoints: null, cache);

    public static ShapeFunctionWitnesses BuildEndpoint(double gammaZero = 0.05, WitnessCache? cache = null) =>
        new(BondClass.Endpoint, gammaZero, F86StandardLocations.Full, xPoints: null, cache);

    private IReadOnlyList<ShapeFunctionPoint> ComputePoints()
    {
        // For each x point, gather y(x) across all locations; reduce to (min, max).
        var result = new ShapeFunctionPoint[XPoints.Count];
        for (int xi = 0; xi < XPoints.Count; xi++)
        {
            double x = XPoints[xi];
            double yMin = double.PositiveInfinity, yMax = double.NegativeInfinity;
            foreach (var (c, N) in Locations)
            {
                var curve = Cache.GetOrCompute(c, N, GammaZero);
                var peak = curve.Peak(BondClass);
                if (peak.KMax <= 0) continue;
                double qTarget = (1 + x) * peak.QPeak;
                double y = InterpolateK(curve, BondClass, qTarget) / peak.KMax;
                if (y < yMin) yMin = y;
                if (y > yMax) yMax = y;
            }
            result[xi] = new ShapeFunctionPoint(x, yMin, yMax);
        }
        return result;
    }

    /// <summary>Linear interpolation of K(Q) at qTarget on the curve's Q grid for the given
    /// bond class average. Returns 0 if qTarget is outside the grid.</summary>
    private static double InterpolateK(KCurve curve, BondClass cls, double qTarget)
    {
        var classCurve = curve.BondClassAverage(cls);
        var grid = curve.QGrid;
        if (qTarget < grid[0] || qTarget > grid[^1]) return 0;
        // Find bracket
        int lo = 0, hi = grid.Count - 1;
        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (grid[mid] <= qTarget) lo = mid; else hi = mid;
        }
        double t = (qTarget - grid[lo]) / (grid[hi] - grid[lo]);
        return classCurve[lo] * (1 - t) + classCurve[hi] * t;
    }

    public override string DisplayName => $"f_{BondClass}(x) shape function witnesses";

    public override string Summary =>
        $"plateau y(x=+1.0) ≈ {PostPeakPlateau:F2}, max spread = {Points.Max(p => p.SpreadPercent):F1}% ({Tier.Label()}, computed live)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var pt in Points)
                yield return new InspectableNode(
                    $"x = {Formatting.SignedDelta(pt.X, "F2")}",
                    summary: $"y ∈ [{pt.YMin:F3}, {pt.YMax:F3}] ({pt.SpreadPercent:F1}% spread)",
                    payload: new InspectablePayload.Real("y_mean", pt.YMean, "F4"));
            yield return InspectableNode.RealScalar("post-peak plateau y(x=+1.0)", PostPeakPlateau, "F2");
        }
    }
}
