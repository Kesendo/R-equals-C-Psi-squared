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
}
