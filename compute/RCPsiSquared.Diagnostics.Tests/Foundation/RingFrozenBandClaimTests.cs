using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below gate for F170 (<see cref="RingFrozenBandClaim"/>, <see cref="RingFrozenBandWitness"/>).
/// Every count is an exact rank modulo primes, so every comparison is exact. Beside the theorem's numbers
/// sit reads that can fail: zeros off the locus and off the root, the ladder residual across and within the
/// two rings, the overfull middle rung at N = 4, the gcd excess at N = 9, and a deliberately wrong Ψ (its
/// Jordan-Wigner signs dropped), which shows which reads can see Ψ at all.</summary>
public sealed class RingFrozenBandClaimTests
{
    private static readonly Lazy<RingFrozenBandWitness> W3 = new(() => new RingFrozenBandWitness(3));
    private static readonly Lazy<RingFrozenBandWitness> W5 = new(() => new RingFrozenBandWitness(5));

    private static RingFrozenBandWitness Witness(int n) => n == 3 ? W3.Value : W5.Value;

    private static long[] Locus(int n)
    {
        long[] seed = { 29, -13, 7, 19 };
        var g = Enumerable.Repeat(256L, n).ToArray();
        for (int i = 0; i < n / 2; i++) { g[i] += seed[i % 4]; g[n - 1 - i] -= seed[i % 4]; }
        return g;
    }

    private static long[][] SignDropped(long[][] psi) =>
        psi.Select(row => row.Select(Math.Abs).ToArray()).ToArray();

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    public void Every_Rung_Carries_The_Corner_And_Nothing_New_Is_Seeded_Above_It(int n)
    {
        var w = Witness(n);
        foreach (var ring in new[] { w.SpinRing, w.ReversedSpinRing })
        {
            Assert.Equal(n - 1, ring.Count);
            foreach (var r in ring)
            {
                Assert.Equal(n / 2, r.Locus);
                Assert.Equal(n / 2, r.SecondCoupling);
                Assert.Equal(n - 1, r.Uniform);
                Assert.Equal(r.P == 1 ? n / 2 : 0, r.LocusLowestWeight);
                Assert.Equal(r.P == 1 ? n - 1 : 0, r.UniformLowestWeight);
                // Algebraic multiplicity: dim ker A² meets the kernel count, so there is no Jordan chain.
                Assert.Equal(n / 2, r.LocusAlgebraic);
                Assert.Equal(n - 1, r.UniformAlgebraic);
            }
        }
        Assert.True(w.Holds, w.Summary);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(6)]
    public void Ladder_Commutes_Across_The_Two_Rings_And_Not_Within_One(int n)
    {
        long[] g = Locus(n);
        for (int p = 2; p < n; p++)
        {
            Assert.Equal(0L, RingFrozenBandWitness.LadderResidual(n, p, +1, -1, 768, g));
            Assert.Equal(0L, RingFrozenBandWitness.LadderResidual(n, p, -1, +1, 768, g));
            Assert.True(RingFrozenBandWitness.LadderResidual(n, p, +1, +1, 768, g) > 0, $"N={n} p={p}");
        }
    }

    [Fact]
    public void Witness_Ladder_Fields_Read_Zero_Across_And_Positive_Within()
    {
        Assert.Equal(0L, W5.Value.LadderAcrossRings);
        Assert.True(W5.Value.LadderWithinOneRing > 0);
    }

    [Fact]
    public void A_Psi_Without_Jordan_Wigner_Signs_Fails_The_Ladder_And_The_Even_Control_But_Not_The_Odd_Zeros()
    {
        long[] g5 = Locus(5);
        for (int p = 2; p < 5; p++)
        {
            var wrong = SignDropped(RingFrozenBandWitness.Psi(5, p));
            // The odd-N lowest-weight zeros cannot see the error: a frozen space of dimension two meets the
            // kernel of almost any map only in zero.
            Assert.Equal(0, RingFrozenBandWitness.LowestWeightCount(5, p, +1, 768, g5, 256, wrong));
            // The ladder can: 2J per entry, J = 768 grid units.
            Assert.Equal(1536L, RingFrozenBandWitness.LadderResidual(5, p, +1, -1, 768, g5, wrong));
        }
        long[] g4 = Locus(4);
        Assert.Equal(2, RingFrozenBandWitness.LowestWeightCount(4, 2, +1, 768, g4, 256));
        Assert.Equal(1, RingFrozenBandWitness.LowestWeightCount(4, 2, +1, 768, g4, 256,
            SignDropped(RingFrozenBandWitness.Psi(4, 2))));
    }

    [Fact]
    public void Corner_Sweep_Meets_The_Floor_And_Both_Falsifiers_Read_Zero()
    {
        var w = W3.Value;
        Assert.Equal(8, w.CornerSweep.Count);
        foreach (var c in w.CornerSweep)
        {
            Assert.Equal(c.N / 2, c.Locus);
            Assert.Equal(c.N - 1, c.Uniform);
            Assert.Equal(0, c.OffLocus);
            Assert.Equal(0, c.WrongRoot);
        }
    }

    [Fact]
    public void Even_Control_Reads_An_Overfull_Middle_Rung_At_Four()
    {
        var w = W3.Value;
        Assert.Equal(new[] { 2, 4, 2 }, w.EvenControl.Total);
        Assert.Equal(new[] { 2, 2, 0 }, w.EvenControl.LowestWeight);
    }

    [Fact]
    public void Single_Pair_Profiles_At_Nine_Read_The_Gcd_Excess()
    {
        Assert.Equal(new[] { 4, 5, 4, 4 }, W3.Value.SinglePairCorners);
        Assert.Equal(new[] { 4, 5, 4, 4 },
            Enumerable.Range(0, 4).Select(l => RingFrozenBandWitness.SinglePairPrediction(9, l)));
    }

    [Fact]
    public void Resonance_Count_Is_Zero_Exactly_At_Odd_N()
    {
        for (int n = 3; n <= 41; n++)
            foreach (int s in new[] { 1, -1 })
            {
                int count = RingFrozenBandWitness.ResonantPairCount(n, s);
                if (n % 2 == 1) Assert.Equal(0, count);
                else Assert.True(count > 0, $"N={n} s={s}: no resonant pair at even N");
                Assert.Equal(n % 2 == 1, RingFrozenBandClaim.NoResonance(n));
            }
    }

    [Fact]
    public void Sector_Identity_Holds_For_Both_Spin_Rings_From_Three_To_Eight()
    {
        for (int n = 3; n <= 8; n++)
        {
            Assert.True(RingFrozenBandWitness.SectorIdentityHolds(n, 768, +1), $"N={n} wrap +1");
            Assert.True(RingFrozenBandWitness.SectorIdentityHolds(n, 768, -1), $"N={n} wrap −1");
        }
        Assert.True(W5.Value.SectorIdentity);
        Assert.True(W5.Value.PeriodicDiffersAtTwo);
    }

    [Fact]
    public void Cell_Floor_Exceeds_One_Above_The_Seed_At_Odd_N()
    {
        for (int n = 5; n <= 41; n += 2) Assert.True(RingFrozenBandClaim.CellFloorExceedsOne(n), $"N={n}");
        Assert.Equal((6, 5), RingFrozenBandClaim.CellFloor(5, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => RingFrozenBandClaim.CellFloor(6, 2));
        Assert.Equal(2, RingFrozenBandClaim.RungCount(5, uniform: false));
        Assert.Equal(4, RingFrozenBandClaim.RungCount(5, uniform: true));
        Assert.Throws<ArgumentOutOfRangeException>(() => RingFrozenBandClaim.RungCount(6, uniform: false));
    }

    [Fact]
    public void Witness_Refuses_Even_N_The_Zero_Mean_Stratum_And_Rates_The_Profile_Cannot_Keep_Positive()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new RingFrozenBandWitness(4));
        var zero = Assert.Throws<ArgumentOutOfRangeException>(() => new RingFrozenBandWitness(5, gammaBar: 0.0));
        Assert.Contains("zero-mean", zero.Message);
        var small = Assert.Throws<ArgumentOutOfRangeException>(() => new RingFrozenBandWitness(5, gammaBar: 0.03));
        Assert.Contains("positive", small.Message);
    }

    [Fact]
    public void Claim_Is_Registered_With_Its_Two_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<RingFrozenBandClaim>();
        Assert.Same(registry.Get<FrozenDivisorClaim>(), claim.FrozenDivisor);
        Assert.Same(registry.Get<SpectatorIntertwinerClaim>(), claim.Spectator);
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }
}
