using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

/// <summary>Unit tests for <see cref="F1PalindromeOrbitPairing"/>: the Π-orbit partition of
/// the joint-popcount sectors. Verifies the Π-image permutation rule, its order-4 property,
/// the Π² = X⊗N identity, Π-fixedness, the distinct-spectral-class count, and that the
/// <see cref="F1PalindromeOrbitPairing.PartitionByPiOrbit{TSector}"/> partition covers every
/// sector exactly once with the correct follower kind.</summary>
public sealed class F1PalindromeOrbitPairingTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var x = new F1PalindromeOrbitPairing(new SymmetryFamilyInventory());
        Assert.Equal(Tier.Tier1Derived, x.Tier);
    }

    [Fact]
    public void Constructor_NullInventory_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F1PalindromeOrbitPairing(null!));
    }

    [Fact]
    public void DisplayNameAndSummary_AreNonEmpty()
    {
        var x = new F1PalindromeOrbitPairing(new SymmetryFamilyInventory());
        Assert.False(string.IsNullOrWhiteSpace(x.DisplayName));
        Assert.False(string.IsNullOrWhiteSpace(x.Summary));
    }

    // ----------------------------------------------------------------------
    // PiImage: the whole-sector permutation rule (p_c, p_r) ↦ (N − p_r, p_c)
    // ----------------------------------------------------------------------

    [Theory]
    // Hand-computed: PiImage(N, p_c, p_r) = (N - p_r, p_c).
    [InlineData(3, 0, 0, 3, 0)]
    [InlineData(3, 1, 0, 3, 1)]
    [InlineData(3, 2, 3, 0, 2)]
    [InlineData(4, 1, 2, 2, 1)]
    [InlineData(4, 0, 0, 4, 0)]
    [InlineData(5, 2, 3, 2, 2)]
    [InlineData(6, 3, 3, 3, 3)]
    public void PiImage_MapsCorrectly(int N, int pCol, int pRow, int expectedPCol, int expectedPRow)
    {
        var (pcOut, prOut) = F1PalindromeOrbitPairing.PiImage(N, pCol, pRow);
        Assert.Equal(expectedPCol, pcOut);
        Assert.Equal(expectedPRow, prOut);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void PiImage_AppliedFourTimes_IsIdentity(int N)
    {
        for (int pc = 0; pc <= N; pc++)
        for (int pr = 0; pr <= N; pr++)
        {
            var label = (PCol: pc, PRow: pr);
            for (int k = 0; k < 4; k++)
                label = F1PalindromeOrbitPairing.PiImage(N, label.PCol, label.PRow);
            Assert.Equal((pc, pr), label);
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void PiImage_AppliedTwice_EqualsXNPairSector(int N)
    {
        // Π² = X⊗N: PiImage∘PiImage must equal XGlobalChargeConjugationPairing.PairSector.
        for (int pc = 0; pc <= N; pc++)
        for (int pr = 0; pr <= N; pr++)
        {
            var once = F1PalindromeOrbitPairing.PiImage(N, pc, pr);
            var twice = F1PalindromeOrbitPairing.PiImage(N, once.PCol, once.PRow);
            var xn = XGlobalChargeConjugationPairing.PairSector(N, pc, pr);
            Assert.Equal(xn, twice);
        }
    }

    // ----------------------------------------------------------------------
    // IsPiFixed: only the centre sector at even N
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(4, 2, 2, true)]
    [InlineData(6, 3, 3, true)]
    [InlineData(4, 1, 3, false)]
    [InlineData(4, 2, 1, false)]
    [InlineData(3, 1, 2, false)]
    public void IsPiFixed_DetectsCentreSectorOnly(int N, int pCol, int pRow, bool expected)
    {
        Assert.Equal(expected, F1PalindromeOrbitPairing.IsPiFixed(N, pCol, pRow));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void IsPiFixed_CountsExactlyOneAtEvenN_ZeroAtOddN(int N)
    {
        int fixedCount = 0;
        for (int pc = 0; pc <= N; pc++)
        for (int pr = 0; pr <= N; pr++)
            if (F1PalindromeOrbitPairing.IsPiFixed(N, pc, pr))
                fixedCount++;
        Assert.Equal(N % 2 == 0 ? 1 : 0, fixedCount);
    }

    // ----------------------------------------------------------------------
    // DistinctSpectralClasses: accounts for (N+1)² and is ≤ the X⊗N count
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(10)]
    public void DistinctSpectralClasses_AccountsForAllSectorsAndIsAtMostXNCount(int N)
    {
        int distinct = F1PalindromeOrbitPairing.DistinctSpectralClasses(N);
        int total = (N + 1) * (N + 1);
        int piFixed = N % 2 == 0 ? 1 : 0;

        // Π-orbits: piFixed singletons + (total - piFixed)/4 orbits of size 4.
        Assert.Equal(0, (total - piFixed) % 4);
        Assert.Equal(piFixed + (total - piFixed) / 4, distinct);

        // The Π-orbit groups 4 sectors where X⊗N groups 2, so it is never more classes.
        int xnDistinct = XGlobalChargeConjugationPairing.DistinctSpectralClasses(N);
        Assert.True(distinct <= xnDistinct,
            $"N={N}: Π-orbit classes {distinct} should be ≤ X⊗N classes {xnDistinct}.");
    }

    // ----------------------------------------------------------------------
    // PartitionByPiOrbit: cover, disjointness, follower-kind correctness
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void PartitionByPiOrbit_PrimariesAndFollowers_CoverAllSectorsDisjointly(int N)
    {
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (primaries, followerToPrimary) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow));

        // Primaries and followers are disjoint and together cover every sector index.
        var primarySet = primaries.ToHashSet();
        var followerSet = followerToPrimary.Keys.ToHashSet();
        Assert.Empty(primarySet.Intersect(followerSet));
        Assert.Equal(sectors.Count, primarySet.Count + followerSet.Count);
        var covered = primarySet.Union(followerSet).ToHashSet();
        for (int i = 0; i < sectors.Count; i++)
            Assert.Contains(i, covered);

        // The primary count equals the distinct-spectral-class count.
        Assert.Equal(F1PalindromeOrbitPairing.DistinctSpectralClasses(N), primaries.Count);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void PartitionByPiOrbit_EveryFollowerKind_IsCorrect(int N)
    {
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (_, followerToPrimary) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow));

        foreach (var (followerIdx, follower) in followerToPrimary)
        {
            var followerLabel = (sectors[followerIdx].PCol, sectors[followerIdx].PRow);
            var primaryLabel = (sectors[follower.PrimaryIndex].PCol, sectors[follower.PrimaryIndex].PRow);

            // Count how many PiImage applications take the primary to this follower.
            var walk = primaryLabel;
            int k = 0;
            for (int step = 1; step <= 4; step++)
            {
                walk = F1PalindromeOrbitPairing.PiImage(N, walk.Item1, walk.Item2);
                if (walk == followerLabel) { k = step; break; }
            }
            Assert.True(k is 1 or 2 or 3,
                $"N={N}: follower {followerLabel} is not reached from primary {primaryLabel} in 1..3 Π-steps.");

            var expectedKind = k == 2
                ? F1PalindromeOrbitPairing.F1FollowerKind.XnCopy
                : F1PalindromeOrbitPairing.F1FollowerKind.F1Reflect;
            Assert.Equal(expectedKind, follower.Kind);

            // Cross-check: the XnCopy follower is exactly the X⊗N partner of the primary.
            if (follower.Kind == F1PalindromeOrbitPairing.F1FollowerKind.XnCopy)
            {
                var xnPartner = XGlobalChargeConjugationPairing.PairSector(
                    N, primaryLabel.Item1, primaryLabel.Item2);
                Assert.Equal(xnPartner, followerLabel);
            }
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    public void PartitionByPiOrbit_OddN_EveryOrbitHasOnePrimaryAndThreeFollowers(int N)
    {
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (primaries, followerToPrimary) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow));

        // At odd N there is no Π-fixed sector; every primary owns exactly three followers,
        // and every follower kind appears: one XnCopy, two F1Reflect per orbit.
        var followersByPrimary = followerToPrimary.Values
            .GroupBy(f => f.PrimaryIndex)
            .ToDictionary(g => g.Key, g => g.ToList());
        foreach (int p in primaries)
        {
            Assert.True(followersByPrimary.ContainsKey(p), $"N={N}: primary {p} has no followers.");
            var fs = followersByPrimary[p];
            Assert.Equal(3, fs.Count);
            Assert.Equal(1, fs.Count(f => f.Kind == F1PalindromeOrbitPairing.F1FollowerKind.XnCopy));
            Assert.Equal(2, fs.Count(f => f.Kind == F1PalindromeOrbitPairing.F1FollowerKind.F1Reflect));
        }
    }

    [Fact]
    public void PartitionByPiOrbit_EvenN_CentreSectorIsPrimaryWithNoFollowers()
    {
        const int N = 4;
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (primaries, followerToPrimary) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow));

        int centreIdx = -1;
        for (int i = 0; i < sectors.Count; i++)
            if (sectors[i].PCol == N / 2 && sectors[i].PRow == N / 2)
                centreIdx = i;
        Assert.True(centreIdx >= 0);

        // The centre sector is a primary and is referenced by no follower.
        Assert.Contains(centreIdx, primaries);
        Assert.DoesNotContain(centreIdx, followerToPrimary.Keys);
        Assert.DoesNotContain(centreIdx, followerToPrimary.Values.Select(f => f.PrimaryIndex));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void PartitionByPiOrbit_WithSectorSize_SortsPrimariesDescendingBySize(int N)
    {
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (primaries, _) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow), s => s.Size);

        for (int i = 1; i < primaries.Count; i++)
            Assert.True(sectors[primaries[i - 1]].Size >= sectors[primaries[i]].Size,
                $"N={N}: primaries not sorted descending by size at index {i}.");
    }
}
