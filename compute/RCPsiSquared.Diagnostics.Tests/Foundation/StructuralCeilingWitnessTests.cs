using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first: the live commutant computation reproduces the structural-ceiling closed forms
/// (the C# twin of simulations/topology_ceiling_rep_derivation.py). Small sector matrices, fast.
///
/// <para><b>OPEN, and worth stating here because this is the file that says "verified".</b> What these
/// tests certify is that four small matrices have the eigenvalues 4/N, 4/(N−1), 1 and 2 − 2/√3. What
/// they do NOT certify is that those eigenvalues ARE the gap rate. That leg is the Python verifier's
/// STAGE 0b, which compares the J-independent sector prediction against the full 4^N Liouvillian at two
/// couplings a decade apart and gates the truncation as an error LAW (order 2.00, since ad_H is
/// anti-Hermitian and the first-order correction to the real part cancels). Nothing on the C# side
/// builds that oracle, so the derivation half of "derived, not fit" lives only in Python.
///
/// The machinery is next door and unused here: <c>SecondClockRegimeWitness.G2Full</c> builds exactly
/// that full Liouvillian, but is called only at Q = 2 and Q = 8 for the EP classification, and its own
/// docstring declines the high-Q comparison ("avoids the Q=1000 full-Liouvillian eigenvalue-precision
/// risk"). Closing this would mean either lifting G2Full to a shared helper and adding a high-Q leg
/// with the same two-Q order gate, or saying plainly in the witness that the oracle comparison is not
/// reproduced here. Two further gaps of the same kind: the Python asserts that the winning sector SET
/// contains (1,1) and its X^N partner (an exactly degenerate pair, so an argmin over it reports
/// rounding), and it derives the ring-4 co-occupier value 2√2 from the anti-periodic Jordan-Wigner
/// sector; neither globality nor that derivation has a C# counterpart. Recorded 2026-09-01, from the
/// verifier review round; see docs/CAUGHT_ERRORS.md for that round.</para></summary>
public class StructuralCeilingWitnessTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Complete_OneOneSector_Is4OverN(int n)
    {
        double? g = StructuralCeilingWitness.CommutantDarkest("complete", n, 1, 1);
        Assert.True(g.HasValue);
        Assert.Equal(4.0 / n, g!.Value, 9);
    }

    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Star_OneOneSector_Is4OverNMinus1(int n)
    {
        double? g = StructuralCeilingWitness.CommutantDarkest("star", n, 1, 1);
        Assert.True(g.HasValue);
        Assert.Equal(4.0 / (n - 1), g!.Value, 9);
    }

    [Fact]
    public void N4Outlier_LivesInTheTwoTwoHalfFillingSector()
    {
        // the (1,1) ladder hits 1.0 at N=4 (= 4/4 = the band edge), so it makes no ceiling there
        double? k4Ladder = StructuralCeilingWitness.CommutantDarkest("complete", 4, 1, 1);
        Assert.True(k4Ladder.HasValue);
        Assert.Equal(1.0, k4Ladder!.Value, 9);

        // the K_4 ceiling is the (2,2) half-filling sector = 2 − 2/√3, below the floor
        double? k4Ceiling = StructuralCeilingWitness.CommutantDarkest("complete", 4, 2, 2);
        Assert.True(k4Ceiling.HasValue);
        Assert.Equal(2.0 - 2.0 / Math.Sqrt(3.0), k4Ceiling!.Value, 7);
        Assert.True(k4Ceiling.Value < 1.0);

        // ring-4 is special in the SAME (2,2) sector, but co-occupies the floor (= 1.0)
        double? ring4Ceiling = StructuralCeilingWitness.CommutantDarkest("ring", 4, 2, 2);
        Assert.True(ring4Ceiling.HasValue);
        Assert.Equal(1.0, ring4Ceiling!.Value, 9);
    }

    [Fact]
    public void Chain_NeverCeilings()
    {
        // the chain (no adjacency degeneracy) keeps its (1,1) ladder above 1: the band edge protects
        for (int n = 4; n <= 6; n++)
        {
            double? g = StructuralCeilingWitness.CommutantDarkest("chain", n, 1, 1);
            if (g.HasValue) Assert.True(g.Value > 1.0 - 1e-9, $"chain N={n} should not ceiling, got {g.Value}");
        }
    }
}
