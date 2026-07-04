using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the parameter-side Klein V4 (F91 + F92 + F93, adopted 2026-07-04 from
// PROOF_F91_GAMMA_NINETY_DEGREES.md / PROOF_F92_BOND_ANTI_PALINDROMIC_J.md /
// PROOF_F93_DETUNING_ANTI_PALINDROMIC.md via the registry entries): on a parameter axis (gamma per
// site, J per bond, h per site) the palindromic mirror F71 and the anti-palindromic reshuffle R90
// are two commuting involutions generating V4 = Z2 x Z2, and the proofs' sharper law is entry-wise
// (MirrorWorld genre, no eigensolver): the F71-refined DIAGONAL-block matrix elements of L depend
// only on the F71 pair-sums, so any two profiles in the anti-palindromic orbit (all pair-sums =
// 2*avg) give IDENTICAL diagonal blocks, entry for entry; the breaking lives in the cross-blocks.
public class ParameterKleinTests
{
    static readonly World W = new();

    static double[] Uniform(int len, double v) => Enumerable.Repeat(v, len).ToArray();

    static void AssertClose(double[] expected, double[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(Math.Abs(expected[i] - actual[i]) < 1e-12, $"element {i}: {expected[i]} vs {actual[i]}");
    }

    // the V4 on a parameter axis: F71 (reverse) and the mean-reflection are commuting involutions,
    // R90 = their product; all four distinct on a generic profile. And the anti-palindromic class
    // is exactly the FIXED-POINT set of R90 (x_l + x_{L-1-l} = 2 avg  <=>  R90 x = x).
    [Fact]
    public void The_Parameter_V4_Is_Two_Commuting_Involutions()
    {
        var x = new[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };
        AssertClose(x, ParameterKlein.Mirror(ParameterKlein.Mirror(x)));
        AssertClose(x, ParameterKlein.Reshuffle(ParameterKlein.Reshuffle(x)));
        AssertClose(x, ParameterKlein.MeanReflection(ParameterKlein.MeanReflection(x)));
        // closure: R90 = MeanReflection o Mirror = Mirror o MeanReflection
        AssertClose(ParameterKlein.Reshuffle(x), ParameterKlein.MeanReflection(ParameterKlein.Mirror(x)));
        AssertClose(ParameterKlein.Reshuffle(x), ParameterKlein.Mirror(ParameterKlein.MeanReflection(x)));
        // all four distinct on a generic profile
        var images = new[] { x, ParameterKlein.Mirror(x), ParameterKlein.Reshuffle(x), ParameterKlein.MeanReflection(x) };
        for (int a = 0; a < 4; a++)
            for (int b = a + 1; b < 4; b++)
                Assert.True(images[a].Zip(images[b]).Any(p => Math.Abs(p.First - p.Second) > 1e-9),
                    $"images {a} and {b} coincide");
        // the anti-palindromic witness profile is fixed pointwise by R90 (and NOT by the mirror)
        var anti = new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 };                   // all pair-sums 0.9
        AssertClose(anti, ParameterKlein.Reshuffle(anti));
        Assert.False(anti.SequenceEqual(ParameterKlein.Mirror(anti)));
    }

    // F91 (gamma axis), the registry's N=6 witness: uniform 0.45, the monotonic ladder and the
    // non-monotonic anti-palindromic profile (all pair-sums 0.9) give bit-identical F71-refined
    // diagonal blocks on every (p,q); the permuted and concentrated profiles are rejected at O(1).
    [Fact]
    public void F91_Gamma_Orbit_Is_EntryWise_Identical_And_The_Rest_Rejected()
    {
        var klein = new ParameterKlein(W, 6);
        var j = Uniform(5, 1.0);
        var h = Uniform(6, 0.0);
        var uniform = Uniform(6, 0.45);
        var monotonic = new[] { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
        var nonMonotonic = new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 };
        var permuted = new[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };               // pair-sums {1.1, 0.8, 0.8}
        var concentrated = new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 2.2 };

        Assert.True(klein.DiagonalBlocksResidual(uniform, j, h, monotonic, j, h) < 1e-12);
        Assert.True(klein.DiagonalBlocksResidual(uniform, j, h, nonMonotonic, j, h) < 1e-12);
        Assert.True(klein.DiagonalBlocksResidual(uniform, j, h, permuted, j, h) > 0.1);
        Assert.True(klein.DiagonalBlocksResidual(uniform, j, h, concentrated, j, h) > 0.1);
    }

    // F91 at N=5 (odd N: the middle site must sit at gamma_avg).
    [Fact]
    public void F91_Gamma_Orbit_Holds_At_Odd_N()
    {
        var klein = new ParameterKlein(W, 5);
        var j = Uniform(4, 1.0);
        var h = Uniform(5, 0.0);
        var anti = new[] { 0.2, 0.4, 0.45, 0.5, 0.7 };                        // pair-sums 0.9, middle 0.45
        Assert.True(klein.DiagonalBlocksResidual(Uniform(5, 0.45), j, h, anti, j, h) < 1e-12);
    }

    // F92 (J axis): bond profiles with constant pair-sum J_b + J_{N-2-b} = 2 J_avg are entry-wise
    // identical on the refined diagonals; a same-sum permutation with unequal pair-sums is not.
    [Fact]
    public void F92_Bond_Orbit_Is_EntryWise_Identical_And_The_Rest_Rejected()
    {
        var klein = new ParameterKlein(W, 5);
        var g = Uniform(5, 0.45);
        var h = Uniform(5, 0.0);
        var anti = new[] { 0.8, 1.0, 1.0, 1.2 };                              // pair-sums 2.0, 2.0
        var broken = new[] { 1.2, 0.8, 1.0, 1.0 };                            // pair-sums 2.2, 1.8
        Assert.True(klein.DiagonalBlocksResidual(g, Uniform(4, 1.0), h, g, anti, h) < 1e-12);
        Assert.True(klein.DiagonalBlocksResidual(g, Uniform(4, 1.0), h, g, broken, h) > 0.01);
    }

    // F93 (h axis, longitudinal detuning only): same law on the third axis.
    [Fact]
    public void F93_Detuning_Orbit_Is_EntryWise_Identical_And_The_Rest_Rejected()
    {
        var klein = new ParameterKlein(W, 5);
        var g = Uniform(5, 0.45);
        var j = Uniform(4, 1.0);
        var anti = new[] { 0.1, 0.25, 0.3, 0.35, 0.5 };                       // pair-sums 0.6, middle 0.3
        var broken = new[] { 0.5, 0.1, 0.3, 0.35, 0.25 };                     // pair-sums 0.75, 0.45
        Assert.True(klein.DiagonalBlocksResidual(g, j, Uniform(5, 0.3), g, j, anti) < 1e-12);
        Assert.True(klein.DiagonalBlocksResidual(g, j, Uniform(5, 0.3), g, j, broken) > 0.01);
    }

    // the breaking lives in the eigenvectors only: for the uniform profile the F71 cross-blocks
    // vanish identically (F71 is an exact L-symmetry there); for the non-monotonic anti-palindromic
    // profile they are O(1) -- while the diagonal blocks match bit for bit (the F91 statement).
    [Fact]
    public void The_Breaking_Lives_In_The_Cross_Blocks()
    {
        var klein = new ParameterKlein(W, 6);
        var j = Uniform(5, 1.0);
        var h = Uniform(6, 0.0);
        Assert.True(klein.CrossBlockNorm(Uniform(6, 0.45), j, h) < 1e-12,
            "uniform gamma must keep F71 exact (no cross-block)");
        Assert.True(klein.CrossBlockNorm(new[] { 0.3, 0.5, 0.4, 0.5, 0.4, 0.6 }, j, h) > 0.05,
            "the anti-palindromic profile must break F71 into the cross-blocks");
    }
}
