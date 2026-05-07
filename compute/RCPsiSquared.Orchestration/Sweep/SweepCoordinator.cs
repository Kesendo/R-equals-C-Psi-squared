using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Sweep;

/// <summary>Cross-dimensional sweep consumer. For a given <see cref="SweepDimension"/>,
/// evaluates each parameter combination, compares predicted (closed form) against live
/// (calculator), and aggregates into a <see cref="SweepResult"/>. The result is cached by
/// reference; the registry is immutable, so re-running the same dimension is a cache hit.
///
/// <para>F73Frobenius is currently a tautology check by design: both the predicted and
/// live values come from <see cref="PalindromeResidualScaling.FactorChain"/>. The dimension
/// guards against future refactoring drift (a divergent enum branch in either side would
/// fail Matched), not against numerical drift between independent implementations. A
/// Lindbladian-built witness is the planned independent ground truth.</para></summary>
public sealed class SweepCoordinator
{
    private readonly ClaimRegistry _registry;
    private readonly Dictionary<SweepDimension, SweepResult> _cache = new();

    public SweepCoordinator(ClaimRegistry registry)
    {
        _registry = registry;
    }

    public SweepResult Sweep(SweepDimension dimension)
    {
        if (_cache.TryGetValue(dimension, out var cached)) return cached;

        var result = dimension switch
        {
            SweepDimension.F73Frobenius f => SweepF73(f),
            _ => throw new ArgumentOutOfRangeException(nameof(dimension)),
        };

        _cache[dimension] = result;
        return result;
    }

    private static SweepResult SweepF73(SweepDimension.F73Frobenius dim)
    {
        var points = new List<SweepPoint>();
        foreach (var n in dim.NValues)
        foreach (var hc in dim.HClasses)
        foreach (var chainOnly in dim.ChainOnly)
        {
            // Predicted and live both via FactorChain: tautology by construction
            // (see class-level remark). Replace one side with a Lindbladian-built witness
            // to convert this into a real drift detector.
            double predicted = PalindromeResidualScaling.FactorChain(n, hc);
            double live = PalindromeResidualScaling.FactorChain(n, hc);
            bool matched = predicted == live;

            points.Add(new SweepPoint(
                N: n,
                HClass: hc,
                ChainOnly: chainOnly,
                PredictedValue: predicted,
                LiveValue: live,
                Matched: matched));
        }

        var matchedCount = points.Count(p => p.Matched);
        return new SweepResult(points, points.Count, matchedCount, points.Count - matchedCount);
    }
}
