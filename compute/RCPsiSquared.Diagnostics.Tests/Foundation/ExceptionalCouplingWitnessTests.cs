using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live F50 exceptional-coupling witness. Every number pinned here is in the PAULI book (H = J·Σ(XX+YY+ZZ),
/// x = γ/J), the book of <c>simulations/f50_exceptional_couplings.py</c>; the witness reaches it on the repo's own
/// block builders and exact integer arithmetic rather than the gate's hand-built blocks and sympy, so the two meeting
/// is two routes, not one re-run. The exact parts are compared exactly. The eigensolver counts carry the gate's 10⁻⁶
/// window, and the window is a law only where the nullity of L + 2γ equals the count (the eigenvalue is semisimple to
/// the solver, so it moves by O(ε‖L‖) under the backward error and by O(δ) from the 2^−41 reading, far inside the
/// window); every pinned window count AT A POINT OF E below is asserted together with that equality. The generic count
/// carries no nullity: off E the theorem itself says every −2γ mode is a commutant mode with no Jordan chain, so there
/// the window is a law by the proof. At the defective points (count > nullity, both at γ/J = 2, rational, the readings
/// taken at the rational itself) the window count is a reading, so it is bounded rather than pinned: count > nullity
/// is asserted as the regression guard against evaluating off the point (the offset reads 4 < 5 and 8 < 12), and the
/// nullity (geometric) and the exact route over ℚ(i) carry the numbers.</summary>
public class ExceptionalCouplingWitnessTests
{
    private static BigInteger[] Poly(params long[] ascending) => ascending.Select(c => new BigInteger(c)).ToArray();

    [Fact]
    public void Children_AreAllLive()
    {
        var w = new ExceptionalCouplingWitness(2, "chain");
        Assert.All(((IInspectable)w).Children, c => Assert.Equal(NodeProvenance.Live, c.Provenance));
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ExceptionalCouplingWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ExceptionalCouplingWitness(ExceptionalCouplingWitness.MaxN + 1));
        Assert.Throws<ArgumentException>(() => new ExceptionalCouplingWitness(3, "tree"));
    }

    // ---------------------------------------------------------------- the exact polynomial side

    [Theory]
    [InlineData(2, "chain")]
    [InlineData(3, "chain")]
    [InlineData(3, "complete")]
    [InlineData(4, "chain")]
    [InlineData(4, "ring")]
    [InlineData(4, "star")]
    [InlineData(4, "complete")]
    public void EveryEvenBlock_LeadingCoefficient_IsTheHammingProduct_Exactly(int n, string topology)
    {
        var w = new ExceptionalCouplingWitness(n, topology);
        Assert.All(w.DiagonalBlocks.Concat(w.OffDiagonalEvenBlocks), b =>
        {
            Assert.True(b.LeadingMatches, $"block ({b.PKet},{b.QBra}) leading {b.RealPart[^1]} vs Π(2−2h) = {b.PredictedLeading}");
            Assert.NotEqual(BigInteger.Zero, b.PredictedLeading);
        });
        // (a) of the proof needs no computation: an even off-diagonal block has distance ≥ 2 everywhere; the
        // polynomial route agrees, no positive root on any of them
        Assert.All(w.OffDiagonalEvenBlocks, b => Assert.Equal(0, b.DistinctPositiveRoots));
    }

    [Theory]
    [InlineData(2, "chain")]
    [InlineData(3, "chain")]
    [InlineData(3, "complete")]
    [InlineData(4, "chain")]
    [InlineData(4, "ring")]
    [InlineData(4, "star")]
    [InlineData(4, "complete")]
    public void PinnedMinimalPolynomials_AreCertifiedByExactDivision_AndExhaustTheSturmCount(int n, string topology)
    {
        var w = new ExceptionalCouplingWitness(n, topology);
        Assert.True(w.HasPinnedTable);
        Assert.NotEmpty(w.Pinned);
        Assert.All(w.Pinned, c => Assert.True(c.Holds, $"p={c.P} [{string.Join(",", c.MinimalPolynomial)}]^{c.Multiplicity}: divides {c.PowerDivides}, next fails {c.NextPowerFails}"));
        Assert.True(w.PinnedCoversAllRoots);
    }

    [Fact]
    public void N2Chain_OnePoint_AtGammaOverJ2_TheDefectiveDoubleRoot()
    {
        var w = new ExceptionalCouplingWitness(2, "chain");
        var b11 = Assert.Single(w.DiagonalBlocks.Where(b => b.PKet == 1));
        Assert.Equal(4, b11.Dim);
        // Q_1 is (x − 2)·(something with no positive root): one distinct positive root of multiplicity 1 in γ, while
        // the mode at it has algebraic multiplicity 2 and geometric 1 (defective): multiplicity in γ and mode count
        // are two different numbers, read exactly below
        Assert.Equal(1, b11.DistinctPositiveRoots);
        Assert.Equal(1, b11.PositiveRootsWithMultiplicity);
        Assert.Equal(new[] { 2.0 }, w.ExceptionalSet.Select(r => Math.Round(r, 9)));
        Assert.Equal(4, w.GenericCount);
        // the point is rational, so the multiplicity is exact: algebraic 2, geometric 1 (the (1,1) pair's double root
        // is defective); the window count of the eigensolver is no law here (an EP2 splits as √ε), so it is bounded
        // below by the nullity rather than pinned, the bound being the guard against evaluating off the point
        var e = Assert.Single(w.ExactAtRationalPoints);
        Assert.Equal((new BigInteger(2), BigInteger.One), (e.Numerator, e.Denominator));
        Assert.Equal((2, 1), (e.Algebraic, e.Geometric));
        Assert.Equal(5, w.NullitiesAtPoints[0]);                            // 4 generic + 1 geometric, at γ/J = 2 exactly
        Assert.True(w.CountsAtPoints[0] > 5, $"defective: the window count {w.CountsAtPoints[0]} should exceed the nullity 5");
        Assert.Equal(0.5, w.HandoverQ!.Value, 9);
        Assert.True(w.GapSidesHold);
    }

    [Fact]
    public void N3Chain_TwoPoints_EightAtEach_AndTheHandoverIsQGap3()
    {
        var w = new ExceptionalCouplingWitness(3, "chain");
        Assert.Equal(2, w.ExceptionalSet.Count);
        Assert.Equal(Math.Sqrt((Math.Sqrt(17) - 1) / 2), w.ExceptionalSet[0], 10);
        Assert.Equal(Math.Sqrt(3), w.ExceptionalSet[1], 10);
        foreach (var b in w.DiagonalBlocks.Where(b => b.PKet is 1 or 2))
        {
            Assert.Equal(2, b.DistinctPositiveRoots);                  // C(3,p) − 1 = 2
            Assert.Equal(2, b.PositiveRootsWithMultiplicity);
        }
        Assert.Equal(6, w.GenericCount);
        Assert.Equal(new[] { 8, 8 }, w.CountsAtPoints);
        Assert.Equal(w.CountsAtPoints, w.NullitiesAtPoints);                 // semisimple: the window is a law here
        // Q*_gap(3) = √((1+√17)/8) = 0.800243, the six decimals absorption_ladder_regimes.py bisects
        Assert.Equal(Math.Sqrt((1 + Math.Sqrt(17)) / 8), w.HandoverQ!.Value, 10);
        Assert.InRange(Math.Abs(w.HandoverQ.Value - 0.800243), 0, 5e-7);
        Assert.True(w.GapSidesHold);
    }

    [Fact]
    public void K3_OnePoint_AtSqrt3_WithMultiplicityTwoPerBlock_AndCountTwelve()
    {
        var w = new ExceptionalCouplingWitness(3, "complete");
        Assert.Equal(new[] { Math.Round(Math.Sqrt(3), 10) }, w.ExceptionalSet.Select(r => Math.Round(r, 10)));
        foreach (var b in w.DiagonalBlocks.Where(b => b.PKet is 1 or 2))
        {
            Assert.Equal(1, b.DistinctPositiveRoots);
            Assert.Equal(2, b.PositiveRootsWithMultiplicity);
        }
        Assert.Equal(8, w.GenericCount);                                // K₃: 2N + 2
        Assert.Equal(new[] { 12 }, w.CountsAtPoints);
        Assert.Equal(new[] { 12 }, w.NullitiesAtPoints);                // (x²−3)² in γ, yet semisimple in λ (the gate's G4: 4, 4)
        Assert.Null(w.GapSides);                                        // not the chain
    }

    [Fact]
    public void N4Chain_EightPoints_PerBlockBinomialMinusOne_TheTabulatedCounts_AndQGap4()
    {
        var w = new ExceptionalCouplingWitness(4, "chain");
        var perBlock = w.DiagonalBlocks.Where(b => b.PKet is >= 1 and <= 3).ToDictionary(b => b.PKet, b => b.DistinctPositiveRoots);
        Assert.Equal(new Dictionary<int, int> { [1] = 3, [2] = 5, [3] = 3 }, perBlock);   // C(4,p) − 1
        Assert.Equal(8, w.ExceptionalSet.Count);
        var expected = new[] { 0.745022, 0.745439, 0.910180, 1.545936, 1.572303, 1.867978, 2.088800, 2.197368 };
        for (int i = 0; i < 8; i++) Assert.Equal(expected[i], w.ExceptionalSet[i], 6);
        Assert.Equal(8, w.GenericCount);
        Assert.Equal(new[] { 9, 9, 10, 9, 10, 9, 9, 10 }, w.CountsAtPoints);
        Assert.Equal(w.CountsAtPoints, w.NullitiesAtPoints);                 // all eight semisimple
        // Q*_gap(4) = 1/x₀, x₀ the root of the degree-12 pinned polynomial from the (2,2) block, 1.342243 bisected
        Assert.InRange(Math.Abs(w.HandoverQ!.Value - 1.342243), 0, 5e-7);
        Assert.True(w.GapSidesHold);
        var deg12 = Poly(256, 0, 1280, 0, -2240, 0, -1696, 0, 68, 0, 132, 0, 9);
        var b22 = w.DiagonalBlocks.Single(b => b.PKet == 2);
        Assert.NotNull(ExceptionalCouplingWitness.ExactDivide(b22.Coefficients, deg12));
        Assert.Contains(Math.Round(w.ExceptionalSet[0], 9), b22.Roots.Select(r => Math.Round(r, 9)));
    }

    /// <summary>The window count is pinned only together with count = nullity (semisimple, so the window is a law;
    /// three of the star's four points and K₄'s point are multiple roots in γ yet semisimple in λ, the gate's G5
    /// nullities). The ring's γ/J = 2 is defective, 16 eigenvalues in the window against a nullity of 12: a −1 here,
    /// where the else branch asserts the nullity 12 exactly and count > nullity as the bound; the exact (8, 4) is read
    /// in the ring's own test below.</summary>
    [Theory]
    [InlineData("ring", new[] { 1.378129, 2.0, 2.309401 }, new[] { 9, -1, 9 })]
    [InlineData("star", new[] { 0.681250, 1.154701, 1.732051, 2.309401 }, new[] { 12, 10, 10, 9 })]
    [InlineData("complete", new[] { 2.309401 }, new[] { 11 })]
    public void N4_OtherGraphs_PointsAndCounts(string topology, double[] points, int[] counts)
    {
        var w = new ExceptionalCouplingWitness(4, topology);
        Assert.Equal(points.Length, w.ExceptionalSet.Count);
        for (int i = 0; i < points.Length; i++) Assert.Equal(points[i], w.ExceptionalSet[i], 6);
        Assert.Equal(8, w.GenericCount);
        for (int i = 0; i < counts.Length; i++)
            if (counts[i] >= 0)
            {
                Assert.Equal(counts[i], w.CountsAtPoints[i]);
                Assert.Equal(counts[i], w.NullitiesAtPoints[i]);
            }
            else
            {
                Assert.Equal(12, w.NullitiesAtPoints[i]);                    // the gate's nullity, read at γ/J = 2 exactly
                Assert.True(w.CountsAtPoints[i] > w.NullitiesAtPoints[i],
                    $"the ring's γ/J = 2 is defective: count {w.CountsAtPoints[i]} should exceed the nullity 12");
            }
    }

    [Fact]
    public void N3_RingIsK3_AndStarIsTheChain_SharingTheirPinnedTables()
    {
        var ring = new ExceptionalCouplingWitness(3, "ring");
        Assert.True(ring.HasPinnedTable);
        Assert.True(ring.PinnedCoversAllRoots);
        Assert.Equal(new[] { Math.Round(Math.Sqrt(3), 10) }, ring.ExceptionalSet.Select(r => Math.Round(r, 10)));
        var star = new ExceptionalCouplingWitness(3, "star");
        Assert.True(star.HasPinnedTable);
        Assert.True(star.PinnedCoversAllRoots);
        Assert.Equal(2, star.ExceptionalSet.Count);
        Assert.Throws<ArgumentException>(() => new ExceptionalCouplingWitness(2, "ring"));
    }

    [Fact]
    public void N4Ring_AtGammaOverJ2_TheRootIsDoubleInGamma_AndTheModesAreDefective_Exactly()
    {
        var w = new ExceptionalCouplingWitness(4, "ring");
        var b11 = w.DiagonalBlocks.Single(b => b.PKet == 1);
        Assert.Equal(1, b11.DistinctPositiveRoots);
        Assert.Equal(2, b11.PositiveRootsWithMultiplicity);
        Assert.Equal(2.0, b11.Roots[0], 10);
        // the gate's eigensolver read 16 eigenvalues against a nullity of 12 at this point (8 extra algebraic, 4
        // geometric over the generic 8); here the same numbers are exact over ℚ(i), no window
        var e = Assert.Single(w.ExactAtRationalPoints);
        Assert.Equal((new BigInteger(2), BigInteger.One), (e.Numerator, e.Denominator));
        Assert.Equal((8, 4), (e.Algebraic, e.Geometric));
    }

    [Fact]
    public void N3Chain_HasNoRationalPoint_AndN4Chain_Neither()
    {
        Assert.Empty(new ExceptionalCouplingWitness(3, "chain").ExactAtRationalPoints);
        Assert.Empty(new ExceptionalCouplingWitness(4, "chain").ExactAtRationalPoints);
    }

    // ---------------------------------------------------------------- the exact helpers themselves

    [Fact]
    public void ExactDivide_IsExactOrNull()
    {
        // (x² − 3)(x − 2) = x³ − 2x² − 3x + 6
        var prod = Poly(6, -3, -2, 1);
        Assert.Equal(Poly(-2, 1), ExceptionalCouplingWitness.ExactDivide(prod, Poly(-3, 0, 1)));
        Assert.Null(ExceptionalCouplingWitness.ExactDivide(prod, Poly(-2, 0, 1)));          // x² − 2 does not divide
        Assert.Null(ExceptionalCouplingWitness.ExactDivide(Poly(1, 1), Poly(1, 2)));         // 2x + 1 ∤ x + 1 over ℤ
        Assert.Equal(Poly(1, 1), ExceptionalCouplingWitness.ExactDivide(Poly(1, 3, 2), Poly(1, 2)));  // (2x+1)(x+1)
    }

    [Fact]
    public void PositiveRootCounting_AndIsolation_OnKnownPolynomials()
    {
        // (x − 2)²·(x + 1)·(x² + 1): one positive root, multiplicity 2
        var p = Poly(-2, 1);
        p = ExceptionalCouplingWitness.Multiply(p, Poly(-2, 1));
        p = ExceptionalCouplingWitness.Multiply(p, Poly(1, 1));
        p = ExceptionalCouplingWitness.Multiply(p, Poly(1, 0, 1));
        Assert.Equal(1, ExceptionalCouplingWitness.DistinctPositiveRoots(p));
        Assert.Equal(2, ExceptionalCouplingWitness.PositiveRootsWithMultiplicity(p));
        Assert.Equal(new[] { 2.0 }, ExceptionalCouplingWitness.IsolatePositiveRoots(p).Select(r => Math.Round(r, 10)));
        // x·(x² − 3)·(x² − 2): two positive roots, zero excluded
        var q = ExceptionalCouplingWitness.Multiply(Poly(0, 1), ExceptionalCouplingWitness.Multiply(Poly(-3, 0, 1), Poly(-2, 0, 1)));
        Assert.Equal(2, ExceptionalCouplingWitness.DistinctPositiveRoots(q));
        var roots = ExceptionalCouplingWitness.IsolatePositiveRoots(q);
        Assert.Equal(Math.Sqrt(2), roots[0], 10);
        Assert.Equal(Math.Sqrt(3), roots[1], 10);
        // no positive root at all
        Assert.Equal(0, ExceptionalCouplingWitness.DistinctPositiveRoots(Poly(1, 0, 1)));
        Assert.Empty(ExceptionalCouplingWitness.IsolatePositiveRoots(Poly(1, 0, 1)));
    }

    // ---------------------------------------------------------------- the claim

    [Fact]
    public void Claim_Is_Registered_With_Three_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ExceptionalCouplingSetClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<ExceptionalCouplingSetClaim>().Tier);
        var ancestors = registry.AncestorsOf<ExceptionalCouplingSetClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F50WeightOneDegeneracyPi2Inheritance), ancestors);
        Assert.Contains(typeof(JointPopcountSectors), ancestors);
    }

    [Fact]
    public void Claim_Statement_CarriesTheFences()
    {
        var s = ExceptionalCouplingSetClaim.Shared.Name;
        Assert.Contains("FINITE", s);
        Assert.Contains("real eigenvalues only", s);
        Assert.Contains("rung k = 1 only", s);
        Assert.Contains("not a theorem consequence", s);
        Assert.Contains("Pauli book", s);
    }
}
