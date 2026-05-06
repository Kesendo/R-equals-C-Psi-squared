using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Sweep;

/// <summary>Layer 3 consumer 3 (sweep). For a given <see cref="SweepDimension"/>, evaluates
/// each parameter combination, compares predicted (closed form) against live (calculator),
/// and aggregates into a <see cref="SweepResult"/>. The result is cached by reference; the
/// registry is immutable, so re-running the same dimension is a cache hit.</summary>
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

    private SweepResult SweepF73(SweepDimension.F73Frobenius dim)
    {
        var points = new List<SweepPoint>();
        foreach (var n in dim.NValues)
        foreach (var hc in dim.HClasses)
        foreach (var chainOnly in dim.ChainOnly)
        {
            // Closed-form prediction.
            double pow = Math.Pow(4, n - 2);
            double predicted = hc switch
            {
                HamiltonianClass.Main => (n - 1) * pow,
                HamiltonianClass.SingleBody => (2 * n - 3) * pow,
                _ => double.NaN,
            };

            // "Live" via the calculator (FactorChain).
            double live = PalindromeResidualScaling.FactorChain(n, hc);

            bool matched = predicted == live;
            points.Add(new SweepPoint(
                Parameters: new Dictionary<string, object>
                {
                    ["N"] = n,
                    ["HClass"] = hc,
                    ["ChainOnly"] = chainOnly,
                },
                PredictedValue: predicted,
                LiveValue: live,
                Matched: matched));
        }

        var matchedCount = points.Count(p => p.Matched);
        return new SweepResult(points, points.Count, matchedCount, points.Count - matchedCount);
    }
}
