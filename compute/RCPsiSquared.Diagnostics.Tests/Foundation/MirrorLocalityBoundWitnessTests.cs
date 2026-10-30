using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The mirror's locality bound under a temperature gradient (THE_CUBE_SITE_BY_SITE finding 1, the Z
/// corner). Every number here is exact: the integer generator, the rational steady state verified by L·ρ = 0,
/// the rank over ℚ. The N = 2 letters are the second route's (sympy over ℚ(i), block G of
/// simulations/cube_site_by_site_findings.py): ¼·II + 5/56·IZ + 27/224·ZI + 67/1568·ZZ + 1/112·(XY − YX).</summary>
public class MirrorLocalityBoundWitnessTests
{
    private static BigRational Q(long n, long d) => new(n, d);

    private static BigRational Letter(MirrorLocalityBoundWitness w, string label) =>
        w.Letters().Single(l => l.Label == label).Value;

    [Fact]
    public void TwoTemperatures_N2_TheSteadyStateIsExactAndOfFullRank()
    {
        var w = new MirrorLocalityBoundWitness(n: 2, down: new long[] { 3, 2 }, up: new long[] { 1, 1 });
        Assert.True(w.ZCornerExact);
        Assert.Equal(7, w.C);
        Assert.Equal(0, w.PrimesSkipped);
        Assert.Single(MirrorLocalityBoundWitness.NullspaceModP(w.TwiceGenerator, MirrorLocalityBoundWitness.Primes[0]));
        Assert.True(w.RankIsExact);
        Assert.Equal(4, w.SchmidtRank);
        Assert.Equal(4, w.MaxRank);
        // the letters, against the sympy route
        Assert.Equal(6, w.Letters().Count());
        Assert.Equal(Q(1, 4), Letter(w, "II"));
        Assert.Equal(Q(5, 56), Letter(w, "IZ"));
        Assert.Equal(Q(27, 224), Letter(w, "ZI"));
        Assert.Equal(Q(67, 1568), Letter(w, "ZZ"));
        Assert.Equal(Q(1, 112), Letter(w, "XY"));      // the current: the antisymmetric pair +XY, −YX
        Assert.Equal(Q(-1, 112), Letter(w, "YX"));
    }

    [Fact]
    public void TwoTemperatures_N2_Mirrored_SwapsTheSitesAndReversesTheCurrent()
    {
        var w = new MirrorLocalityBoundWitness(n: 2, down: new long[] { 2, 3 }, up: new long[] { 1, 1 });
        Assert.Equal(4, w.SchmidtRank);
        Assert.Equal(Q(27, 224), Letter(w, "IZ"));
        Assert.Equal(Q(5, 56), Letter(w, "ZI"));
        Assert.Equal(Q(-1, 112), Letter(w, "XY"));
        Assert.Equal(Q(1, 112), Letter(w, "YX"));
    }

    [Fact]
    public void OneTemperature_N2_IsAProduct_Rank1()
    {
        var w = new MirrorLocalityBoundWitness(n: 2, down: new long[] { 3, 3 }, up: new long[] { 1, 1 });
        Assert.True(w.ZCornerExact);
        Assert.True(w.RankIsExact);
        Assert.Equal(1, w.SchmidtRank);
        // the Gibbs product (I + κZ)/2 ⊗ (I + κZ)/2, κ = (3 − 1)/(3 + 1) = ½: c_IZ·c_ZI = c_II·c_ZZ and no current
        Assert.Equal(Letter(w, "II") * Letter(w, "ZZ"), Letter(w, "IZ") * Letter(w, "ZI"));
        Assert.DoesNotContain(w.Letters(), l => l.Label == "XY" || l.Label == "YX");
        Assert.Equal(Q(1, 8), Letter(w, "IZ"));
    }

    [Fact]
    public void NoHopping_TwoTemperatures_IsAProduct_Rank1()
    {
        var w = new MirrorLocalityBoundWitness(n: 2, down: new long[] { 3, 2 }, up: new long[] { 1, 1 }, j: 0);
        Assert.True(w.RankIsExact);
        Assert.Equal(1, w.SchmidtRank);
        Assert.DoesNotContain(w.Letters(), l => l.Label == "XY" || l.Label == "YX");
    }

    [Fact]
    public void N3_Gradient_WithZZ_FullRankAcrossBothCuts_AndOneTemperatureIsRank1()
    {
        foreach (int cut in new[] { 1, 2 })
        {
            var w = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 2, 1 }, up: new long[] { 1, 1, 1 }, j: 1, delta: 1, cut: cut);
            Assert.True(w.ZCornerExact);
            Assert.True(w.RankIsExact);
            Assert.Equal(4, w.SchmidtRank);
        }
        var one = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 3, 3 }, up: new long[] { 1, 1, 1 }, j: 1, delta: 1, cut: 1);
        Assert.True(one.RankIsExact);
        Assert.Equal(1, one.SchmidtRank);
    }

    [Fact]
    public void N3_ACutAcrossTheEqualRatioBondAlone_IsStillFullRank()
    {
        // (3, 3, 2): bond 0-1 joins equal ratios, bond 1-2 carries the gradient; the cut 1|2 crosses only the equal one
        var w = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 3, 2 }, up: new long[] { 1, 1, 1 }, cut: 1);
        Assert.True(w.RankIsExact);
        Assert.Equal(4, w.SchmidtRank);
        var other = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 3, 2 }, up: new long[] { 1, 1, 1 }, cut: 2);
        Assert.Equal(4, other.SchmidtRank);
    }

    [Fact]
    public void N3_TheZZTermIsPinnedByValue_XZYIsOddInDelta()
    {
        // at N = 2 the ZZ term commutes with the whole steady state (its letters are Δ-free); at N = 3 the letter XZY
        // appears only with Δ ≠ 0 and flips sign with it, while XYI is Δ-even (sympy over ℚ(i), the method of block G)
        var plus = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 2, 1 }, up: new long[] { 1, 1, 1 }, j: 1, delta: 1);
        var minus = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 2, 1 }, up: new long[] { 1, 1, 1 }, j: 1, delta: -1);
        var zero = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 2, 1 }, up: new long[] { 1, 1, 1 }, j: 1, delta: 0);
        Assert.Equal(Q(5385638735, 11804814373464), Letter(plus, "XZY"));
        Assert.Equal(Q(-5385638735, 11804814373464), Letter(minus, "XZY"));
        Assert.DoesNotContain(zero.Letters(), l => l.Label == "XZY");
        Assert.Equal(Q(46727008493, 7869876248976), Letter(plus, "XYI"));
        Assert.Equal(Q(46727008493, 7869876248976), Letter(minus, "XYI"));
        Assert.Equal(Q(43493, 6998448), Letter(zero, "XYI"));
    }

    [Fact]
    public void N4_Gradient_FullRankAtTheMiddleCut_16of16()
    {
        var edge = new MirrorLocalityBoundWitness(n: 4, down: new long[] { 6, 5, 4, 3 }, up: new long[] { 2, 2, 2, 2 }, cut: 1);
        Assert.True(edge.ZCornerExact);
        Assert.True(edge.RankIsExact);
        Assert.Equal(4, edge.SchmidtRank);
        var middle = new MirrorLocalityBoundWitness(n: 4, down: new long[] { 6, 5, 4, 3 }, up: new long[] { 2, 2, 2, 2 }, cut: 2);
        Assert.True(middle.RankIsExact);
        Assert.True(middle.PrimesUsed > 2 && middle.PrimesUsed <= MirrorLocalityBoundWitness.MaxPrimes);
        Assert.Equal(16, middle.SchmidtRank);
        Assert.Equal(16, middle.MaxRank);
    }

    [Fact]
    public void TheGeneratorPreservesTheTrace_TheIdentityRowVanishes()
    {
        // trace preservation: the I⊗N row of L is zero (no string feeds the identity coefficient)
        var w = new MirrorLocalityBoundWitness(n: 3, down: new long[] { 3, 2, 1 }, up: new long[] { 1, 1, 1 }, delta: 1);
        var M = w.TwiceGenerator;
        for (int s = 0; s < M.GetLength(1); s++) Assert.Equal(0, M[0, s]);
    }

    [Fact]
    public void Children_AreAllLive()
    {
        var w = new MirrorLocalityBoundWitness();
        var children = ((IInspectable)w).Children.ToList();
        Assert.Equal(5, children.Count);
        Assert.All(children, c => Assert.Equal(NodeProvenance.Live, c.Provenance));
        Assert.Contains("= 4 of 4 (exact over ℚ)", w.Summary);
        Assert.Equal(2, w.PrimesUsed);
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorLocalityBoundWitness(n: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorLocalityBoundWitness(n: 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorLocalityBoundWitness(n: 2, cut: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorLocalityBoundWitness(n: 2, down: new long[] { 3 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MirrorLocalityBoundWitness(n: 2, down: new long[] { 0, 2 }, up: new long[] { 0, 1 }));
    }

    [Fact]
    public void ThePrimesAreTestedNotRecalled()
    {
        Assert.Equal(MirrorLocalityBoundWitness.MaxPrimes, MirrorLocalityBoundWitness.Primes.Length);
        Assert.All(MirrorLocalityBoundWitness.Primes, p => Assert.True(p < (1L << 30) && p > (1L << 29)));
        Assert.Equal(MirrorLocalityBoundWitness.Primes.Distinct().Count(), MirrorLocalityBoundWitness.Primes.Length);
        foreach (long p in MirrorLocalityBoundWitness.Primes)
            for (long d = 2; d * d <= p; d++) Assert.NotEqual(0, p % d);
    }

    [Fact]
    public void RationalReconstruction_RecoversSmallFractions()
    {
        var m = (System.Numerics.BigInteger)MirrorLocalityBoundWitness.Primes[0] * MirrorLocalityBoundWitness.Primes[1];
        // x ≡ 5/56 (mod m), m the product of the two primes (composite, so the inverse is by Euclid)
        var x = 5 * ModInverse(56, m) % m;
        var q = MirrorLocalityBoundWitness.RationalReconstruct(x, m);
        Assert.Equal(Q(5, 56), q!.Value);
        var y = (m - 27 * ModInverse(224, m) % m) % m;          // −27/224
        Assert.Equal(Q(-27, 224), MirrorLocalityBoundWitness.RationalReconstruct(y, m)!.Value);
        // a residue that no fraction inside Wang's bound explains does not close (3/10^21 has none; checked, not
        // guaranteed by the size of d alone, since most residues DO have some small fraction behind them)
        var big = System.Numerics.BigInteger.Parse("1000000000000000000000");
        var z = 3 * ModInverse(big, m) % m;
        Assert.Null(MirrorLocalityBoundWitness.RationalReconstruct(z, m));
    }

    private static System.Numerics.BigInteger ModInverse(System.Numerics.BigInteger a, System.Numerics.BigInteger m)
    {
        System.Numerics.BigInteger g = m, x = 0, y = 1, r = a % m;
        if (r < 0) r += m;
        (System.Numerics.BigInteger r0, System.Numerics.BigInteger r1, System.Numerics.BigInteger s0, System.Numerics.BigInteger s1) = (g, r, x, y);
        while (!r1.IsZero)
        {
            var qq = r0 / r1;
            (r0, r1) = (r1, r0 - qq * r1);
            (s0, s1) = (s1, s0 - qq * s1);
        }
        Assert.Equal(System.Numerics.BigInteger.One, r0);
        return ((s0 % m) + m) % m;
    }

    // ---------------------------------------------------------------- the claim

    [Fact]
    public void Claim_Is_Registered_With_Its_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<MirrorLocalityBoundClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<MirrorLocalityBoundClaim>().Tier);
        var ancestors = registry.AncestorsOf<MirrorLocalityBoundClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PalindromeTwoEndCountClaim), ancestors);
        Assert.Contains(typeof(F84ThermalAmplitudeDampingPi2Inheritance), ancestors);
    }
}
