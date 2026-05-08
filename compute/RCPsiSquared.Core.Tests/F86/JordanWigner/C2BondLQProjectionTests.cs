using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class C2BondLQProjectionTests
{
    private readonly ITestOutputHelper _out;

    public C2BondLQProjectionTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F71MirrorInvariance_PerOrbit_BitExact_AtQEqual2(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = C2BondLQProjection.BuildAtQ(block, Q: 2.0);

        Assert.True(proj.MaxF71MirrorDeviation < 1e-10,
            $"N={N}: F71 mirror deviation {proj.MaxF71MirrorDeviation:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointXbNorm_LessThan_InnermostXbNorm_AtQEqual2(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var proj = C2BondLQProjection.BuildAtQ(block, Q: 2.0);

        var endpoint = proj.Bonds.First(b => b.BondClass == BondClass.Endpoint);
        var innermost = ClosestToCenter(proj.Bonds, N);
        if (innermost is null || innermost.Bond == endpoint.Bond) return;

        Assert.True(endpoint.XbFrobeniusNorm < innermost.XbFrobeniusNorm,
            $"N={N}: Endpoint ‖xB‖={endpoint.XbFrobeniusNorm:F4} should be < Innermost ‖xB‖={innermost.XbFrobeniusNorm:F4}");
    }

    private static BondLQProjectionWitness? ClosestToCenter(IReadOnlyList<BondLQProjectionWitness> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        BondLQProjectionWitness? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void AtQEqualsZero_XbFrobeniusNorm_MatchesMhPerBondFrobeniusNorm(int N)
    {
        // Sanity: at Q=0 the L matrix is purely D (diagonal), so R is identity-like
        // (eigenvectors are computational basis vectors, up to permutation). xB(0) =
        // MhPerBond[b] in the same basis, hence ‖xB‖_F = ‖MhPerBond[b]‖_F = bond-uniform
        // per T6. This anchors the Q-dependence.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var projAtZero = C2BondLQProjection.BuildAtQ(block, Q: 0.0);

        double refNorm = projAtZero.Bonds[0].XbFrobeniusNorm;
        for (int b = 1; b < projAtZero.Bonds.Count; b++)
        {
            Assert.Equal(refNorm, projAtZero.Bonds[b].XbFrobeniusNorm, precision: 8);
        }
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => C2BondLQProjection.BuildAtQ(block, Q: 2.0));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var proj = C2BondLQProjection.BuildAtQ(block, Q: 2.0);
        Assert.Equal(Tier.Tier2Verified, proj.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var proj = C2BondLQProjection.BuildAtQ(block, Q: 2.0);
        Assert.Contains("PROOF_F86_QPEAK", proj.Anchor);
        Assert.Contains("Direction (b'')", proj.Anchor);
    }

    // Q-curve exploration: does ‖xB(Q)‖_F have its own peak structure, or does it
    // monotonically grow/shrink with Q? Determines whether Q_peak extraction is meaningful
    // on this observable or whether we need a richer K_b(Q, t) computation.
    [Fact]
    public void Reconnaissance_EmitsXbNorm_QSweep_AtN8()
    {
        int N = 8;
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        double[] qGrid = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0 };

        _out.WriteLine($"=== N={N}: ‖xB(Q)‖_F per bond, Q-sweep ===");
        _out.WriteLine("    b\\Q | " + string.Join(" | ", qGrid.Select(q => $"Q={q,4:F1}")));
        _out.WriteLine("    ----|" + string.Join("|", qGrid.Select(_ => "---------")));

        int numBonds = block.Decomposition.NumBonds;
        var rows = new double[numBonds, qGrid.Length];
        for (int iQ = 0; iQ < qGrid.Length; iQ++)
        {
            var proj = C2BondLQProjection.BuildAtQ(block, qGrid[iQ]);
            for (int b = 0; b < numBonds; b++) rows[b, iQ] = proj.Bonds[b].XbFrobeniusNorm;
        }
        for (int b = 0; b < numBonds; b++)
        {
            string bondClass = (b == 0 || b == numBonds - 1) ? "End" : "Int";
            var line = new System.Text.StringBuilder($"  {b}({bondClass}) |");
            for (int iQ = 0; iQ < qGrid.Length; iQ++) line.Append($" {rows[b, iQ],7:F2}  |");
            _out.WriteLine(line.ToString());
        }
    }

    [Fact]
    public void Reconnaissance_EmitsPerNXbNorm_AtQ2_AcrossN5To8()
    {
        _out.WriteLine("  N | b | class    | ‖xB(Q=2)‖_F | class");
        _out.WriteLine("  --|---|----------|-------------|--------");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var proj = C2BondLQProjection.BuildAtQ(block, Q: 2.0);
            foreach (var b in proj.Bonds)
            {
                _out.WriteLine($"  {N} | {b.Bond} | {b.BondClass,-8} | {b.XbFrobeniusNorm,11:F4}");
            }
            _out.WriteLine($"  -- N={N}: max F71 mirror deviation = {proj.MaxF71MirrorDeviation:G3}");
            _out.WriteLine("");
        }
    }
}
