using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class SlowManifoldPauliContentTests
{
    private readonly ITestOutputHelper _out;
    public SlowManifoldPauliContentTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void CrossoverCore_Is_IZ_On_Lit_Sites_Rotating_Carries_XY()
    {
        // The fan resolved the N = 3 crossover slow manifold into a 4-dimensional invariant core plus
        // a rotating remainder. The reading to nail: the core is {I, Z} on the LIT sites (0, 1, the
        // X/Y carriers the turn rotates), and the rotating directions carry an X or a Y on a lit site
        // (PTF's far bank). The Pauli projection confirms it: the core's mass sits on {I, Z}-on-lit
        // strings (lit-XY weight ≈ 0), the rotating mass carries lit-XY content.
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);

        var reading = SlowManifoldPauliContent.Compute(
            sweep.SlowBasis[0], sweep.SlowBasis[^1], axis.N, axis.LitSites);

        _out.WriteLine($"core dim {reading.CoreDim}, lit-XY weight {reading.CoreLitXYWeight:E3}");
        foreach (var s in reading.CoreTop)
            _out.WriteLine($"  core   {s.Label}  {s.Weight:0.000}{(s.LitXY ? "  [lit-XY]" : "")}");
        _out.WriteLine($"rotating dim {reading.RotatingDim}, lit-XY weight {reading.RotatingLitXYWeight:E3}");
        foreach (var s in reading.RotatingTop)
            _out.WriteLine($"  rot    {s.Label}  {s.Weight:0.000}{(s.LitXY ? "  [lit-XY]" : "")}");

        Assert.Equal(4, reading.CoreDim);
        Assert.Equal(12, reading.RotatingDim);

        // The core is {I, Z} on the lit sites: essentially no XY-on-lit content (to the Evd floor).
        Assert.True(reading.CoreLitXYWeight < 1e-6,
            $"core should be {{I,Z}} on lit sites; lit-XY weight {reading.CoreLitXYWeight:E3}");
        // The rotating part carries XY on the lit sites: a substantial fraction.
        Assert.True(reading.RotatingLitXYWeight > 0.1,
            $"rotating part should carry XY on lit sites; lit-XY weight {reading.RotatingLitXYWeight:E3}");
    }

    [Fact]
    public void IntersectionCore_Retires_The_SlowCount_Window()
    {
        // Edge 2: the threshold-free core. The angle-core inflated 4 → 18 → 27 at slowCount
        // 16 → 24 → 32 (the lens-window caveat: the θ₀ and θ windows re-include each other's rotated
        // images). The intersection core, rank(Π_C·P_slow·Π_C) with Π_C the diagonal projector onto
        // the fixed cell (Δ_l = 0 on every lit site, where Ad_{R_z} is exactly the identity), counts
        // only genuinely fixed directions, needs no angle threshold and no θ-pair, and on the
        // cluster-closed window is window-free in the meaningful sense: every requested slowCount
        // that closes to the same rate manifold reads the same rank.
        //   requested 8, 16, 24  → all close to the 28-dim manifold → rank 10, at θ₀ and θ_final
        //   requested 32         → closes to the 36-dim manifold    → rank 12 (real fixed content
        //                          that became slow; honest membership, not a lens artifact)
        //   requested 6          → already cluster-aligned (6-dim)  → rank 4, the historical core
        // (Pinned 2026-06-10, NumPy twin: ranks θ-stable across all 13 θ on closed windows; the
        // near-1 singular cluster is separated from the next σ by ≥ 0.05, so 1e-6 is a rank cut,
        // not a lens.)
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);

        var sweep6 = DimensionSweep.Compute(axis, slowCount: 6);
        var sweep16 = DimensionSweep.Compute(axis, slowCount: 16);
        var sweep24 = DimensionSweep.Compute(axis, slowCount: 24);
        var sweep32 = DimensionSweep.Compute(axis, slowCount: 32);

        int Rank(DimensionSweepResult s, int p) =>
            SlowManifoldPauliContent.IntersectionCoreRank(s.ClusterClosedBasis[p], axis.N, axis.LitSites);

        // The historical 4 is the closed 6-dim manifold's core (the diagonal {I,Z} shadow).
        Assert.Equal(4, Rank(sweep6, 0));
        Assert.Equal(4, Rank(sweep6, sweep6.Theta.Count - 1));

        // The window-freedom gate: 16 and 24 close to the same 28-dim manifold and read the same 10,
        // where the angle-core read 4 and 18. Stable at θ₀ and at the final θ.
        Assert.Equal(28, sweep16.ClusterClosedBasis[0].ColumnCount);
        Assert.Equal(28, sweep24.ClusterClosedBasis[0].ColumnCount);
        Assert.Equal(10, Rank(sweep16, 0));
        Assert.Equal(10, Rank(sweep16, sweep16.Theta.Count - 1));
        Assert.Equal(10, Rank(sweep24, 0));
        Assert.Equal(10, Rank(sweep24, sweep24.Theta.Count - 1));

        // Requested 32 closes to the next manifold (36-dim) and honestly reads 12: the +2 is real
        // fixed-cell content entering the slow window, not the angle-core's 27.
        Assert.Equal(36, sweep32.ClusterClosedBasis[0].ColumnCount);
        Assert.Equal(12, Rank(sweep32, 0));

        // The raw 16-window (an arbitrary 10-of-22 slice of the Re=−2γ cluster) carries a membership
        // gauge: the NumPy twin reads 4 at twelve of thirteen θ and 5 at θ = π/4 (the slice happens
        // to include a fifth in-cell direction there). The solver-independent guarantees are the
        // bounds: ≥ 4 (the closed 6-dim manifold and its 4 fixed directions are always inside) and
        // ≤ 10 (the raw slice is a subspace of the closed 28-dim manifold). Assert the guarantees,
        // not the gauge.
        int raw0 = SlowManifoldPauliContent.IntersectionCoreRank(sweep16.SlowBasis[0], axis.N, axis.LitSites);
        Assert.InRange(raw0, 4, 10);

        _out.WriteLine($"closed ranks: k6→{Rank(sweep6, 0)}, k16→{Rank(sweep16, 0)}, k24→{Rank(sweep24, 0)}, k32→{Rank(sweep32, 0)}; raw16 θ₀→{raw0}");
    }
}
