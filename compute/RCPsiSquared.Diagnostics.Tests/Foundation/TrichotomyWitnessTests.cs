using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class TrichotomyWitnessTests
{
    private static IReadOnlyList<double> Uniform(int n, double gamma) =>
        Enumerable.Repeat(gamma, n).ToList();

    [Fact]
    public void Ctor_RejectsOutOfRangeN() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => new TrichotomyWitness(n: 9));

    [Fact]
    public void SurvivorRate_MatchesPostEpFlowField_AtN5_VacBlock()
    {
        // The convention pin (R2): the (0,1) survivor rate on the absolute SectorSlowest scale equals
        // the full-L PostEpFlowField slowest rate to 9 digits. Same (q, profile) convention as the
        // existing SectorReductionWitnessTests pin (VacBlock_SlowestRate_EqualsPostEpFlowField_AtN5):
        // q in {1.5, 1000.0}, tauGrid {0.0, 1.0}, ReadAssembly(q).SlowestRate as ground truth.
        const int n = 5; const double q = 1.5;
        var profile = Uniform(n, 0.5);
        var (pc, pr, rate) = TrichotomyWitness.SurvivorSector(TopologyKind.Chain, n, q, profile);
        var expected = new PostEpFlowField(n, new[] { 1.5, 1000.0 }, new[] { 0.0, 1.0 }, profile)
            .ReadAssembly(q).SlowestRate;
        Assert.Equal(expected, rate, 9);
    }

    [Fact]
    public void Survivor_Star5_IsTheOneOneCommutant()
    {
        var (pc, pr, _) = TrichotomyWitness.SurvivorSector(TopologyKind.Star, 5, 1.5, Uniform(5, 0.5));
        Assert.Equal((1, 1), (pc, pr));
    }
}
