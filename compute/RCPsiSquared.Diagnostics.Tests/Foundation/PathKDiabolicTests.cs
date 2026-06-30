using System;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The path-k diabolic investigation (the N=4→N=5 forward edge of the zeros_connecting_structure
/// arc). The N=4 diabolic is the +2.349 σ_T twin pair merging onto the fold at λ_EP=−4+1.318i; these tests
/// build and validate the generalized path-k diabolic tooling against that KNOWN path-3 case before it is
/// trusted at path-4. See docs/superpowers/plans/2026-06-27-f89-path-k-diabolic.md.</summary>
public class PathKDiabolicTests
{
    // Task 1 — ResidualRootsAt(k,q) is the path-k analogue of GaloisMonodromyWitness.OcticRootsAt, valid at
    // q≈2 (R-7 scope). The trusted cross-check: at q=2 the path-3 residual MUST be the path-3 octic exactly.
    [Fact]
    public void ResidualRootsAt_Path3_MatchesOcticWitness()
    {
        var pk = PathKMonodromyScout.ResidualRootsAt(3, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var oc = GaloisMonodromyWitness.OcticRootsAt(new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        Assert.Equal(8, pk.Length);
        for (int i = 0; i < 8; i++)
            Assert.True((pk[i] - oc[i]).Magnitude < 1e-9, $"strand {i}: {pk[i]} vs {oc[i]}");
    }

    // Task 2 — the diabolic classifier MUST re-find the known path-3 diabolic q_EP≈0.659 at λ=−4+iJ·2
    // before it is trusted at path-4. Region [0.5,0.85]×[±0.05] in q contains q_EP (real) and excludes the
    // q≈0.857 defective EP, so exactly one semisimple coalescence is expected.
    [Fact]
    public void DiabolicScan_Path3_FindsKnownDiabolicAtMinus4()
    {
        var found = PathKMonodromyScout.FindDiabolics(k: 3, reLo: 0.5, reHi: 0.85, imLo: -0.05, imHi: 0.05, cell: 0.01);
        var diab = found.Where(d => d.IsSemisimple).ToList();
        Assert.Single(diab);
        var d = diab[0];
        Assert.True(Math.Abs(d.MergeLambda.Real - (-4.0)) < 1e-2, $"Re(λ_d)={d.MergeLambda.Real}");
        Assert.True(Math.Abs(d.QValue.Real - GaloisMonodromyWitness.QEp) < 1e-2, $"q_d={d.QValue.Real}, q_EP={GaloisMonodromyWitness.QEp}");
        Assert.True(d.PairIsResidual, "the coalescing pair must be H_B-mixed residual, not AT-locked");
        Assert.True(d.LoopIsIdentity, "the diabolic loop must be the identity (semisimple, no braid)");
        // gap closes LINEARLY at a diabolic (two sheets crossing) — exponent ≈ 1, not ½ (√-branch).
        Assert.True(d.GapScalingExponent > 0.7 && d.GapScalingExponent < 1.3,
            $"gap-scaling exponent {d.GapScalingExponent} should be ≈1 (linear) for a diabolic");
    }

    // Task 3 — the standalone character discriminant. At the path-3 diabolic (q_EP, λ_EP=−4+2iJ) the
    // residual block is SEMISIMPLE (geo=alg=2, departure≈0, no Jordan block), reproducing F89Path3OcticEpClaim;
    // a generic isolated eigenvalue is NOT a semisimple degeneracy. This is the R-3 load-bearing test.
    [Fact]
    public void CharacterizeAt_Path3_DiabolicAtQEp_NotAtGenericPoint()
    {
        double qEp = GaloisMonodromyWitness.QEp;
        var lamEp = new Complex(-4, 2 * qEp);                       // −4 + 2iJ, J = q_EP·γ (γ=1)
        var rd = PathKMonodromyScout.CharacterizeAt(3, new Complex(qEp, 0), lamEp, radius: 0.1);
        Assert.Equal(EpCharacter.EpKind.Diabolic, rd.Kind);
        Assert.Equal(2, rd.Algebraic);
        Assert.Equal(2, rd.Geometric);
        Assert.True(rd.Departure < 1e-6, $"departure {rd.Departure} should be ≈0 for a diabolic (no Jordan block)");
        Assert.True(PathKMonodromyScout.IsSemisimpleAt(3, new Complex(qEp, 0), lamEp));

        // negative control: a generic isolated octic eigenvalue at q=2 is not a semisimple degeneracy.
        var generic = GaloisMonodromyWitness.OcticRootsAt(new Complex(2, 0))[0];
        Assert.NotEqual(EpCharacter.EpKind.Diabolic, PathKMonodromyScout.CharacterizeAt(3, new Complex(2, 0), generic, radius: 0.05).Kind);
        Assert.False(PathKMonodromyScout.IsSemisimpleAt(3, new Complex(2, 0), generic, radius: 0.05));
    }

    // Residual-only tracking (the N>=6 AT-flood fix): ResidualRootsTracked follows the residual SET by
    // continuity from q0=2, so it returns the F_d residual strands at ANY q (the q-general analogue of
    // ResidualRootsAt, which is q≈2-only). AT-free by construction: the AT-locked strands are excluded, so
    // the AT-AT exact degeneracies that flood the full-block scan at N>=6 never enter this root set.
    [Fact]
    public void ResidualRootsTracked_Path3_MatchesDirectResidual_AtBase()
    {
        var tracked = PathKMonodromyScout.ResidualRootsTracked(3, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var direct = PathKMonodromyScout.ResidualRootsAt(3, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        Assert.Equal(8, tracked.Length);
        for (int i = 0; i < 8; i++)
            Assert.True((tracked[i] - direct[i]).Magnitude < 1e-9, $"strand {i}: {tracked[i]} vs {direct[i]}");
    }

    // The whole point of the fix: at a generic q away from the base, tracking returns EXACTLY the residual
    // strands (F_d degree 32 for path-5), not the full S₂-sym block (45 = 32 residual + 13 AT-locked). The 13
    // AT strands — whose same-⟨n_XY⟩ exact crossings flood the full-block scan at N=6 — are excluded.
    [Fact]
    public void ResidualRootsTracked_Path5_ReturnsResidualStrandsOnly_AwayFromBase()
    {
        var resid = PathKMonodromyScout.ResidualRootsTracked(5, new Complex(1.3, -0.4));
        Assert.Equal(32, resid.Length);
    }

    // Local residual GapRefine (the second perf fix): descend the residual min-gap field from a seed WITHOUT
    // re-tracking from q0=2 at every probe (the nested cost that made the broad N=7 scan ~1h). Track q0→seed
    // once, then identify the residual subset at each probe by local nearest-neighbour continuity. Must still
    // converge a seed to the residual coalescence (gap → ~0), same as the global GapRefine.
    [Fact]
    public void GapRefineResidualLocal_Path5_DrivesResidualGapToZeroAtDiabolic()
    {
        // seed near the known path-5 residual diabolic q≈0.7581+0.260i (broad scan, gap 1.1e-8).
        var qd = PathKMonodromyScout.GapRefineResidualLocalAt(5, new Complex(0.76, 0.27), cell: 0.02);
        double gapAtQd = PathKMonodromyScout.MinGap(PathKMonodromyScout.ResidualRootsTracked(5, qd));
        Assert.True(gapAtQd < 1e-5, $"residual gap at refined q*={qd} is {gapAtQd:E2}, expected ~0 (converged to the diabolic)");
    }

    // The local residual monodromy loop (the perf fix): instead of re-tracking from q0=2 at every one of the
    // 240 loop points (nested O(loopSteps × trackSteps), ~1h at N=7), track q0→loop-start ONCE then walk the
    // tiny circle by LOCAL continuity over the full roots, following only the residual strands. Must read the
    // SAME verdict as the slow global loop: identity at a diabolic (no braid), transposition at a defective EP.
    [Fact]
    public void ResidualLoopIsIdentity_Path5_TrueAtDiabolic_FalseAtDefectiveEp()
    {
        // from the broad k=5 scan: q≈0.709−0.219i is a clean residual diabolic (gap-exp 0.98, loop identity);
        // q≈1.416−0.217i is a residual defective EP (gap-exp 0.50, loop transposition).
        Assert.True(PathKMonodromyScout.ResidualLoopIsIdentityAt(5, new Complex(0.7090, -0.219)),
            "the small residual loop around a diabolic must be the identity (no braid)");
        Assert.False(PathKMonodromyScout.ResidualLoopIsIdentityAt(5, new Complex(1.4159, -0.217)),
            "the small residual loop around a defective EP must braid (a transposition)");
    }

    // The N=6 fix at scan level: the full-block scan FLOODS at path-5 (the DE sector's same-⟨n_XY⟩ AT strands
    // coincide exactly on dense curves, ~800 residual=False gap=0 coalescences, ZERO residual isolated).
    // residualOnly tracks the residual SET from q0=2, so the AT crossings never enter the gap field and every
    // reported coalescence is a genuine residual (H_B-mixed) one.
    [Fact]
    public void FindDiabolics_Path5_ResidualOnly_ReportsOnlyResidualCoalescences()
    {
        var found = PathKMonodromyScout.FindDiabolics(
            k: 5, reLo: 0.5, reHi: 1.0, imLo: 0.05, imHi: 0.45, cell: 0.02, residualOnly: true);
        Assert.NotEmpty(found);
        Assert.All(found, d => Assert.True(d.PairIsResidual,
            $"residualOnly must report only residual coalescences; got an AT-locked one at q={d.QValue}, gap={d.Gap:E2}"));
    }

    // Item-1 regression: the path-4 near-axis diabolic at q≈0.6118+0.012i sits close to a defective EP that a
    // path-3-tuned fixed 0.02 loop would catch (reading a false transposition). The small intrinsic loop radius
    // must classify it correctly: a TRUE diabolic, loop-identity. (The radius-sweep evidence: identity at
    // r<=0.008, the neighbour EP entering only at r>=0.012.)
    [Fact]
    public void Path4_NearAxisDiabolic_IsLoopIdentity_NotContaminatedByNeighbourEp()
    {
        var found = PathKMonodromyScout.FindDiabolics(4, 0.60, 0.625, -0.05, 0.05, 0.005);
        var near = found.Where(d => d.IsSemisimple && Math.Abs(d.QValue.Real - 0.6118) < 0.01).ToList();
        Assert.NotEmpty(near);
        Assert.All(near, d =>
        {
            Assert.True(d.LoopIsIdentity, $"q={d.QValue} should read loop-identity at the small intrinsic radius");
            Assert.True(d.PairIsResidual);
            Assert.True(Math.Abs(d.QValue.Imaginary) > 0.005, $"q={d.QValue} is off the real axis (complex-q diabolic)");
        });
    }

    // ---- the EXACT residual roots (the path-6/N=7 fix) ----

    // The core regression that PROVES the fix. ResidualRootsExact reads the F_d residual roots as the
    // eigenvalues of M(q) compressed onto the orthogonal complement of the AT invariant subspace — AT-free by
    // construction, no nearest-match partition. It returns EXACTLY the F_d count (18/32/53) and the roots are
    // DISTINCT at the base q0=2, where the tracked path at k=6 instead returned a set riddled with exact
    // duplicates and AT-rate (Re=-6) strands (the 293-false-diabolic flood). path-3 (k≤... ) uses the octic
    // witness; the exact split is wired for k=4,5,6.
    [Theory]
    [InlineData(4, 18)]
    [InlineData(5, 32)]
    [InlineData(6, 53)]
    public void ResidualRootsExact_ReturnsFdCount_AndDistinct_AtBase(int k, int fdDegree)
    {
        var res = PathKMonodromyScout.ResidualRootsExact(k, new Complex(2, 0));
        Assert.Equal(fdDegree, res.Length);
        Assert.True(PathKMonodromyScout.MinGap(res) > 1e-6,
            $"path-{k} residual roots must be DISTINCT (the broken tracked path floods with gap-0 duplicates); min gap = {PathKMonodromyScout.MinGap(res):E2}");
    }

    // q-direct correctness away from the base: at a generic complex q the exact reading still returns exactly
    // the 32 residual strands (path-5), distinct — no tracking, so no continuity collision.
    [Fact]
    public void ResidualRootsExact_Path5_AwayFromBase_Returns32Distinct()
    {
        var res = PathKMonodromyScout.ResidualRootsExact(5, new Complex(1.3, -0.4));
        Assert.Equal(32, res.Length);
        Assert.True(PathKMonodromyScout.MinGap(res) > 1e-6, $"min gap = {PathKMonodromyScout.MinGap(res):E2}");
    }

    // The split is sound: AT ⊎ residual reproduces the FULL block spectrum (the block-triangular decomposition),
    // the AT part sits exactly on the rate lines Re λ ∈ {−2,−6}, and the residual part is everything else.
    // This is what ties eig(M22) to the oracle F_d: AT = ForPathK's validated subspace (full-D test), so the
    // complement eigenvalues are the genuine residual roots.
    [Fact]
    public void ExactSplit_Path6_AtPlusResidual_ReproducesFullSpectrum()
    {
        var q = new Complex(2, 0);
        var (a, c) = PathKMonodromyScout.BuildLinear(7);                 // nBlock = k+1 = 7
        var full = PathKMonodromyScout.AllRootsAt(a, c, q);
        var at = PathKMonodromyScout.AtRootsExact(6, q);
        var res = PathKMonodromyScout.ResidualRootsExact(6, q);

        Assert.Equal(75, full.Length);
        Assert.Equal(22, at.Length);
        Assert.Equal(53, res.Length);
        Assert.All(at, z => Assert.True(Math.Abs(z.Real + 2) < 1e-6 || Math.Abs(z.Real + 6) < 1e-6,
            $"AT root {z} must sit on a rate line Re∈{{-2,-6}}"));

        // multiset equality up to tolerance (greedy match; the spectrum has near-degenerate Re=-6 strands
        // whose float-noise Real parts make a sort-and-pair comparison flip pairings — match by distance).
        var remaining = at.Concat(res).ToList();
        foreach (var f in full)
        {
            int idx = remaining.FindIndex(z => (z - f).Magnitude < 1e-6);
            Assert.True(idx >= 0, $"full spectrum root {f} has no match in AT⊎residual");
            remaining.RemoveAt(idx);
        }
        Assert.Empty(remaining);
    }

    // The exact pipeline reproduces the science: on the path-4 near-axis box it finds the same complex-q
    // diabolic (q≈0.6118±0.012i) the tracked path does, semisimple and off-axis — so the new path is gated on
    // a known result before it is trusted to count path-6 (N=7). The full N=5→11 / N=6→16 counts are the CLI
    // re-gate (documented in experiments/F89_PATH_K_DIABOLIC.md, too slow for a unit test).
    [Fact]
    public void FindDiabolicsExact_Path4_FindsNearAxisDiabolic_Semisimple()
    {
        var found = PathKMonodromyScout.FindDiabolicsExact(4, 0.60, 0.625, -0.05, 0.05, 0.005);
        var near = found.Where(d => d.IsSemisimple && Math.Abs(d.QValue.Real - 0.6118) < 0.01).ToList();
        Assert.NotEmpty(near);
        Assert.All(near, d =>
        {
            Assert.True(d.PairIsResidual);
            Assert.True(Math.Abs(d.QValue.Imaginary) > 0.005, $"q={d.QValue} is off the real axis (complex-q diabolic)");
            Assert.True(d.GapScalingExponent > 0.7 && d.GapScalingExponent < 1.3,
                $"gap-scaling exponent {d.GapScalingExponent} should be ≈1 (linear) for a diabolic");
        });
    }
}
