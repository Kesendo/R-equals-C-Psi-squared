using System;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The q-parametric monodromy tracker: as q traces a loop in the complex plane, follow the
/// d roots by continuity and read the permutation they undergo. This is the spectral face of the
/// Galois group (the monodromy = the Galois group over the function field). The two textbook anchors
/// ARE the physics distinction we apply to the F89 octic: λ²=q is a √-branch point (defective EP,
/// simple discriminant zero) whose loop SWAPS the two roots; λ²=q² is a transversal crossing of two
/// analytic sheets (diabolic, double discriminant zero) whose loop returns IDENTITY.</summary>
public class MonodromyTests
{
    [Fact]
    public void LoopAroundSquareRootBranchPoint_SwapsTheTwoRoots()
    {
        // λ² = q: the two roots ±√q; encircling the branch point q=0 swaps them (the defective signature).
        Func<Complex, Complex[]> roots = q => { var s = Complex.Sqrt(q); return new[] { s, -s }; };
        var perm = Monodromy.Permutation(roots, center: Complex.Zero, radius: 0.5, steps: 400);
        Assert.Equal(new[] { 1, 0 }, perm);
    }

    [Fact]
    public void LoopAroundTransversalCrossing_ReturnsIdentity()
    {
        // λ² = q²: roots ±q coincide at q=0 but as two analytic sheets (double discriminant zero);
        // encircling does NOT swap them — the DIABOLIC signature, exactly the F89 q_EP claim.
        Func<Complex, Complex[]> roots = q => new[] { q, -q };
        var perm = Monodromy.Permutation(roots, center: Complex.Zero, radius: 0.5, steps: 400);
        Assert.Equal(new[] { 0, 1 }, perm);
    }

    [Fact]
    public void LoopEnclosingNoBranchPoint_ReturnsIdentity()
    {
        // a √-system, but the loop is far from the branch point at 0 → nothing braids.
        Func<Complex, Complex[]> roots = q => { var s = Complex.Sqrt(q); return new[] { s, -s }; };
        var perm = Monodromy.Permutation(roots, center: new Complex(5, 0), radius: 0.5, steps: 400);
        Assert.Equal(new[] { 0, 1 }, perm);
    }

    [Fact]
    public void LassoFromADistantBase_AroundTheBranchPoint_SwapsInTheBaseLabelling()
    {
        // λ²=q: from a base q0=1 (roots ±1), lasso out to a small circle around the branch point q=0
        // and back. The net braiding, expressed in the q0 labelling, is the transposition (the basis of
        // assembling many branch points' transpositions into one monodromy = Galois group).
        Func<Complex, Complex[]> roots = q => { var s = Complex.Sqrt(q); return new[] { s, -s }; };
        var lasso = Monodromy.Lasso(new Complex(1, 0), Complex.Zero, radius: 0.05);
        var perm = Monodromy.PermutationAlongPath(roots, lasso);
        Assert.Equal(new[] { 1, 0 }, perm);
    }

    [Fact]
    public void LassoToAPointWithNoBranch_IsIdentity()
    {
        // a lasso around a regular point (no branch enclosed) leaves the labelling unchanged.
        Func<Complex, Complex[]> roots = q => { var s = Complex.Sqrt(q); return new[] { s, -s }; };
        var lasso = Monodromy.Lasso(new Complex(1, 0), new Complex(2, 0), radius: 0.05);
        var perm = Monodromy.PermutationAlongPath(roots, lasso);
        Assert.Equal(new[] { 0, 1 }, perm);
    }

    [Fact]
    public void CubeRootBranchPoint_GivesAThreeCycle()
    {
        // λ³ = q: the three cube roots; encircling q=0 once cycles them (a 3-cycle, not a transposition).
        Func<Complex, Complex[]> roots = q =>
        {
            var r = Complex.Pow(q, 1.0 / 3.0);
            var w = new Complex(Math.Cos(2 * Math.PI / 3), Math.Sin(2 * Math.PI / 3));
            return new[] { r, r * w, r * w * w };
        };
        var perm = Monodromy.Permutation(roots, center: Complex.Zero, radius: 0.5, steps: 600);
        // a single 3-cycle: no fixed points, and applying it three times is the identity.
        Assert.True(perm[0] != 0 && perm[1] != 1 && perm[2] != 2, "a 3-cycle has no fixed points");
        Assert.Equal(0, perm[perm[perm[0]]]);
    }
}
