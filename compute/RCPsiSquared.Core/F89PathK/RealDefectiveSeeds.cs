using System.Collections.Generic;
using System.Linq;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>A real defective seed EP of the (1,2) coherence block: the containment corollary's one per-N
/// input (docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §7) and the step-3 shell census's probe locus.
/// RParity = the residual sector whose real-root count jumped (+1 R-even, −1 R-odd). LambdaA is the
/// recorded 4-decimal value; census runs REFINE both q* and λ_A from the (1,2) block spectrum before
/// probing (the EP's √|q−q*| sensitivity makes the recorded precision a ~1e-3 member-noise floor;
/// refinement pushes it to ~1e-5).</summary>
public readonly record struct RealSeed(int N, double QStar, double LambdaA, int RParity, string Origin);

/// <summary>The recorded real defective seeds, N = 5..11 (source: the PT-break count-change census,
/// PathKMonodromyScout.FindRealDefectiveByCountChange, gates RealSeedCensusTests / SLOW_SEEDCENSUS;
/// the table lives in experiments/F89_PATH_K_DIABOLIC.md). Census window q ∈ [0.2, 3]. The count-change
/// instrument's blind spots (grazes, real-real crossings: no NET count change) mean these lists are
/// lower bounds; entries found by other instruments carry their own Origin.</summary>
public static class RealDefectiveSeeds
{
    const string Census = "count-change census 2026-07-02";

    public static IReadOnlyList<RealSeed> All { get; } = new RealSeed[]
    {
        new(5, 0.620878, -4.6189, +1, Census),
        new(5, 1.077615, -3.7917, +1, Census),
        new(5, 0.643037, -3.8196, -1, Census),
        new(5, 2.804888, -4.4882, -1, Census),
        new(7, 0.538696, -3.9816, +1, Census),
        new(7, 0.924107, -4.6615, +1, Census),
        new(7, 1.514833, -4.8846, +1, Census),
        new(7, 0.600379, -4.9228, -1, Census),
        new(7, 0.723553, -4.8609, -1, Census),
        new(7, 1.590852, -3.8998, -1, Census),
        new(9, 0.591760, -5.1206, +1, Census),
        new(9, 0.659767, -5.1101, +1, Census),
        new(9, 0.849011, -4.8415, +1, Census),
        new(9, 2.137549, -4.0264, +1, Census),
        new(9, 0.511958, -4.1581, -1, Census),
        new(9, 0.675374, -4.4978, -1, Census),
        new(9, 2.761985, -5.0640, -1, Census),
        new(11, 0.502123, -4.3111, +1, Census),
        new(11, 0.598253, -4.5660, +1, Census),
        new(11, 0.812457, -4.8639, +1, Census),
        new(11, 1.902993, -5.2355, +1, Census),
        new(11, 0.587252, -5.2588, -1, Census),
        new(11, 0.631587, -5.2563, -1, Census),
        new(11, 0.727968, -5.2098, -1, Census),
        new(11, 1.010462, -4.7397, -1, Census),
        new(11, 2.700712, -4.1353, -1, Census),
    };

    public static IEnumerable<RealSeed> ForN(int n) => All.Where(s => s.N == n);
}
