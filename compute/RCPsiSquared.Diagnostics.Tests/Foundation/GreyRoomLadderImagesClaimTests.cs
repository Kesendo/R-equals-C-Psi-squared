using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below gate for F171 (<see cref="GreyRoomLadderImagesClaim"/>,
/// <see cref="GreyRoomLadderImagesWitness"/>). Every identity is exact in integers; every image rank is an
/// exact lower bound met from above, every loss an exact rank, and every room an upper bound modulo primes,
/// so every comparison is exact in its stated direction. The literals are the page's (docs/THE_GREY_GLASS.md,
/// read there by the committed gate simulations/grey_glass_gate.py, its checks G7, G8, G10, G11, G12 and G13, and
/// by the sober twin GreyRoom; the odd ring's
/// empty (2, 1) floor by PROOF_RING_GAP_DOMINANCE; the odd ring's (1, 1) by F170's N − 1 at the uniform point),
/// asserted BESIDE the ladder theorem's closed form C(N, k), so a gate built only from the formula cannot pass
/// by construction. Beside the theorem's numbers sit reads that can fail: the losses on the rings (measured,
/// not proved, on the odd ring), the stagger dropped, the Heisenberg chain, a profile the reflection does not
/// keep, and the row k = 1 against the free-fermion floor witness on F86's own block builder. The reads that are
/// this witness's own and no page's, pinned so they cannot drift: 22 at a mirror-symmetric profile, the
/// Heisenberg 3 and 1 at two profiles, the skew chain's rooms all certified, the odd ring's (1, 1) at N = 3,
/// the stagger dropped lost on every route read.</summary>
public sealed class GreyRoomLadderImagesClaimTests
{
    private static readonly Dictionary<(int N, bool Ring), Lazy<GreyRoomLadderImagesWitness>> Cache = new();

    private static GreyRoomLadderImagesWitness W(int n, bool ring)
    {
        lock (Cache)
        {
            if (!Cache.TryGetValue((n, ring), out var w))
                Cache[(n, ring)] = w = new Lazy<GreyRoomLadderImagesWitness>(() => new GreyRoomLadderImagesWitness(n, ring));
            return w.Value;
        }
    }

    private static GreyRoomLadderImagesWitness.RoomRead Read(int n, bool ring, int p, int q, int k) =>
        Read(W(n, ring), p, q, k);

    private static GreyRoomLadderImagesWitness.RoomRead Read(GreyRoomLadderImagesWitness w, int p, int q, int k) =>
        w.Rooms.Single(r => r.P == p && r.Q == q && r.K == k);

    private static bool StrictlyInside(int n, GreyRoomLadderImagesWitness.RoomRead r) =>
        Math.Abs(r.P - r.Q) < r.K && r.K < Math.Min(r.P + r.Q, 2 * n - r.P - r.Q);

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Open_Chain_Every_Image_Lies_On_Its_Plane_Is_One_To_One_And_Is_Kept_Exactly(int n)
    {
        var w = W(n, false);
        Assert.True(w.Holds, string.Join("; ", w.Failures));
        Assert.True(w.ReadsAsThePage, string.Join("; ", w.OpenReads));
        Assert.True(w.Uniform);
        Assert.True(w.Connected);
        Assert.True(w.ChainReflectionKeepsCouplings);
        Assert.True(w.LightIsOneNumberPerPlane);
        Assert.True(w.LaddersCommute);
        Assert.Equal((n + 1) * (n + 1), w.Rooms.Select(r => (r.P, r.Q)).Distinct().Count());
        foreach (var r in w.Rooms)
        {
            Assert.True(r.ImageOnPlane, $"N={n} ({r.P},{r.Q}) k={r.K}: image leaves its plane");
            Assert.Equal(r.LadderDimension, r.ImageRankModP);
            Assert.True(r.KeptExactly, $"N={n} ({r.P},{r.Q}) k={r.K}: identity fails on the chain");
            Assert.False(r.LostByRank);
            Assert.Equal(r.FillsBlock, r.BlockDimension == r.LadderDimension);
            // the room meets the ladder bound everywhere but at the middle height k = N/2 strictly inside a window,
            // where the reflection builds two operators more (the page's G8, measured at N = 3 to 6; read here from
            // N = 2, where no window holds its middle height strictly inside)
            bool middle = n % 2 == 0 && r.K == n / 2 && StrictlyInside(n, r);
            Assert.Equal(r.LadderDimension + (middle ? 2 : 0), r.RoomUpperBound);
            Assert.Equal(!middle, r.Certified);
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Block_1_1_Of_The_Uniform_Open_Chain_Holds_Exactly_The_Two_Ladder_Images(int n)
    {
        var two = Read(n, false, 1, 1, 2);
        var zero = Read(n, false, 1, 1, 0);
        Assert.Equal(GreyRoomLadderImagesClaim.LadderDimension(n, 2), two.RoomUpperBound);
        Assert.True(two.Certified);
        Assert.Equal(1, zero.RoomUpperBound);
        Assert.True(zero.Certified);
        Assert.Equal(GreyRoomLadderImagesClaim.UniformOpenChainBlock11Room(n), two.RoomUpperBound + zero.RoomUpperBound);
    }

    [Fact]
    public void Eleven_At_N_Five() => Assert.Equal(11, GreyRoomLadderImagesClaim.UniformOpenChainBlock11Room(5));

    [Fact]
    public void The_Middle_Height_Reads_The_Reflections_Two_More_Beside_The_Kept_Image()
    {
        // the page's G8: 8 against 6 at N = 4, 22 against 20 at N = 6, a read (the products are proved grey, what
        // they add to the images is measured)
        var r4 = Read(4, false, 2, 2, 2);
        Assert.Equal(6, r4.LadderDimension);
        Assert.Equal(8, r4.RoomUpperBound);
        Assert.True(r4.KeptExactly);
        Assert.False(r4.Certified);
        foreach (var (p, q) in new[] { (3, 2), (2, 3), (4, 3), (3, 4) })
        {
            var r6 = Read(6, false, p, q, 3);
            Assert.Equal(20, r6.LadderDimension);
            Assert.Equal(22, r6.RoomUpperBound);
            Assert.True(r6.KeptExactly);
        }
        // the middle height at the bottom or the top of a window holds the image alone
        Assert.True(Read(6, false, 5, 2, 3).Certified);
        Assert.True(Read(6, false, 2, 1, 3).Certified);
        Assert.True(Read(4, false, 3, 1, 2).Certified);
    }

    [Fact]
    public void A_Profile_The_Reflection_Does_Not_Keep_Loses_The_Two_More_And_A_Mirror_Symmetric_One_Keeps_Them()
    {
        // the page's G8: at couplings the reflection does not keep the room falls to the ladder count; at
        // mirror-symmetric couplings the two more stay. Every image is still kept exactly: F171's chain
        // statement holds at any real couplings
        var skew = new GreyRoomLadderImagesWitness(4, bondCouplings: new[] { 0.75, 0.25, 1.25 });
        Assert.True(skew.Holds, string.Join("; ", skew.Failures));
        Assert.False(skew.Uniform);
        Assert.False(skew.ChainReflectionKeepsCouplings);
        Assert.All(skew.Rooms, r => Assert.True(r.KeptExactly));
        Assert.All(skew.Rooms, r => Assert.True(r.Certified, $"({r.P},{r.Q}) k={r.K} not certified at the skew profile"));
        Assert.Equal(6, Read(skew, 2, 2, 2).RoomUpperBound);
        // (1, 1) at this profile: certified exactly at this input, beyond the proof's uniform chain
        Assert.True(Read(skew, 1, 1, 2).Certified);
        Assert.True(Read(skew, 1, 1, 0).Certified);

        var mirror = new GreyRoomLadderImagesWitness(6, bondCouplings: new[] { 0.5, 1.25, 0.75, 1.25, 0.5 });
        Assert.True(mirror.Holds, string.Join("; ", mirror.Failures));
        Assert.True(mirror.ChainReflectionKeepsCouplings);
        Assert.Equal(22, Read(mirror, 3, 2, 3).RoomUpperBound);
        Assert.Equal(22, Read(mirror, 4, 3, 3).RoomUpperBound);
        Assert.True(Read(mirror, 2, 2, 2).Certified);
        // the Heisenberg (1, 1) at these profiles: on the mirror-symmetric chain of six F140's floor ⌊N/2⌋ = 3 is a
        // theorem (a frozen mode at exactly −4γ lies in the room) and the mod-p room meets it, certified; the chain
        // of four the reflection does not keep has no floor and reads 1 against 2, a read pinned so it cannot
        // drift; the ZZ term loses both ladders' images at either profile (both connected)
        Assert.Equal(3, mirror.HeisenbergBlock11Room);
        Assert.True(mirror.HeisenbergFloorHolds);
        Assert.True(mirror.HeisenbergBlock11Certified);
        Assert.Equal(1, skew.HeisenbergBlock11Room);
        Assert.False(skew.HeisenbergBlock11Certified);
        Assert.True(skew.ReadsAsThePage, string.Join("; ", skew.OpenReads));
        Assert.True(mirror.ReadsAsThePage, string.Join("; ", mirror.OpenReads));
        Assert.True(mirror.HeisenbergEtaImageLost);
        Assert.True(mirror.HeisenbergSpinImageLost);
        Assert.True(skew.HeisenbergEtaImageLost);
        Assert.True(skew.HeisenbergSpinImageLost);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Row_K_1_Is_The_Free_Fermion_Floor_On_F86s_Own_Block_Builder(int n)
    {
        // two block builders, one Krylov routine: the (n, n + 1) blocks cut from the Pauli chain here, F86's
        // CoherenceBlock there; the literal n guards the routine itself
        var floor = new FreeFermionFloorWitness(n);
        for (int lower = 0; lower < n; lower++)
        {
            var r = Read(n, false, lower, lower + 1, 1);
            Assert.Equal(n, r.RoomUpperBound);
            Assert.True(r.Certified);
            Assert.Equal(floor.Blocks[lower].FloorUpperBoundModP, r.RoomUpperBound);
            Assert.True(floor.Blocks[lower].FloorCertified);
        }
    }

    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    public void Even_Ring_Keeps_Even_Heights_Exactly_And_Loses_Odd_Heights_Off_The_Corners(int n)
    {
        var w = W(n, true);
        Assert.True(w.Holds, string.Join("; ", w.Failures));
        Assert.True(w.ReadsAsThePage, string.Join("; ", w.OpenReads));
        Assert.True(w.FluxRuleDock);
        Assert.True(w.WrapSignLawHolds);
        foreach (var r in w.Rooms)
        {
            // the rule, hand-spelled: even k kept, odd k only in a corner (the losses exact at N = 4, 6, 8)
            bool expect = r.K % 2 == 0 || r.FillsBlock;
            Assert.Equal(expect, r.RuleProvedKept);
            Assert.True(r.ImageOnPlane);
            Assert.Equal(r.LadderDimension, r.ImageRankModP);
            if (expect)
            {
                Assert.True(r.KeptExactly, $"N={n} ring ({r.P},{r.Q}) k={r.K}: identity fails");
                // each ladder step turns the flux, so the start's wrap sign is (−1)^(s+r) relative to the target's;
                // the starts (0, 0) and (N, 0) cannot hop, so either sign closes their identity
                if (!r.FillsBlock && r.K != 0 && r.K != n)
                    Assert.Equal((r.SpinSteps + r.EtaSteps) % 2 == 0 ? 1 : -1, r.SourceWrapSign);
                Assert.True(r.RoomUpperBound >= r.LadderDimension);
            }
            else
            {
                Assert.True(r.LostByRank, $"N={n} ring ({r.P},{r.Q}) k={r.K}: loss not shown by rank");
                Assert.False(r.KeptExactly);
            }
        }
    }

    [Fact]
    public void Even_Ring_Rooms_Read_The_Pages_Numbers()
    {
        // N = 4: two more than the ladder images with the mirror products of both bond axes (9), from the antipodal pairs (G10, G11)
        Assert.Equal(11, Read(4, true, 2, 2, 2).RoomUpperBound);
        Assert.Equal(6, Read(4, true, 2, 2, 2).LadderDimension);
        // (1, 1): the odd half with the U-even frozen modes beside it, 6 + 2 and 15 + 2 (the page's G7)
        Assert.Equal(8, Read(4, true, 1, 1, 2).RoomUpperBound);
        Assert.Equal(17, Read(6, true, 1, 1, 2).RoomUpperBound);
        // N = 6, (3, 2) at k = 3: the image of (3, 0) is lost, and the room is the three bond axes' mirror products alone (G10)
        var lost = Read(6, true, 3, 2, 3);
        Assert.True(lost.LostByRank);
        Assert.Equal(20, lost.LadderDimension);
        Assert.Equal(6, lost.RoomUpperBound);
        Assert.Equal(0, GreyRoomLadderImagesClaim.LadderBound(6, 3, 2, 3, ring: true, wrapBondNonzero: true));
        // the ring of six's rooms above the ladders in (2, 2) and (3, 3), the page's Open list: upper bounds, the
        // same Krylov elimination as the sober twin's on a Hamiltonian built another way
        Assert.Equal(19, Read(6, true, 2, 2, 2).RoomUpperBound);
        Assert.Equal(15, Read(6, true, 2, 2, 2).LadderDimension);
        Assert.Equal(18, Read(6, true, 3, 3, 2).RoomUpperBound);
        Assert.Equal(18, Read(6, true, 3, 3, 4).RoomUpperBound);
    }

    [Fact]
    public void Even_Ring_At_A_Profile_Loses_The_Odd_Heights_Off_The_Corners_Too()
    {
        // the proof's §6: off the corners the loss at odd k is a coupling-free map's, so it holds at every nonzero
        // wrap bond, exact at N = 4
        var w = new GreyRoomLadderImagesWitness(4, ring: true, bondCouplings: new[] { 0.75, 0.25, 1.25, 0.5 });
        Assert.True(w.Holds, string.Join("; ", w.Failures));
        Assert.True(w.WrapBondNonzero);
        Assert.True(Read(w, 2, 1, 1).LostByRank);
        Assert.True(Read(w, 3, 2, 1).LostByRank);
        Assert.True(Read(w, 2, 2, 2).KeptExactly);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    public void Odd_Ring_Keeps_An_Image_When_It_Comes_From_A_Corner_By_The_Ladder_Of_Its_Parity_And_Loses_The_Rest_Here(int n)
    {
        var w = W(n, true);
        Assert.True(w.Holds, string.Join("; ", w.Failures));
        Assert.True(w.ReadsAsThePage, string.Join("; ", w.OpenReads));   // the odd ring's losses are reads, as the page measured them
        Assert.True(w.FluxRuleDock);
        Assert.True(w.WrapSignLawHolds);
        foreach (var r in w.Rooms)
        {
            // the rule, hand-spelled: kept from a corner by the ladder of k's parity; that the rest is lost is
            // measured at N = 3, 5, 7 (the page's G13), not proved, and read exactly here at the uniform ring
            bool expect = r.FillsBlock ||
                          (r.K % 2 == 0 ? r.SpinSteps == 0 || r.SpinSteps == r.K : r.EtaSteps == 0 || r.EtaSteps == n - r.K);
            Assert.Equal(expect, r.RuleProvedKept);
            Assert.True(r.ImageOnPlane);
            Assert.Equal(r.LadderDimension, r.ImageRankModP);
            if (!expect)
            {
                Assert.True(r.LostByRank, $"N={n} ring ({r.P},{r.Q}) k={r.K}: loss not shown by rank");
                continue;
            }
            Assert.True(r.KeptExactly, $"N={n} ring ({r.P},{r.Q}) k={r.K}: identity fails");
            if (r.FillsBlock || r.K == 0 || r.K == n) continue;
            // the wrap-sign law: (−1)^(s+r) on a single-ladder route, −(−1)^(s+r) on a mixed route through a corner
            int steps = (r.SpinSteps + r.EtaSteps) % 2 == 0 ? 1 : -1;
            bool mixed = r.SpinSteps > 0 && r.EtaSteps > 0;
            Assert.Equal(mixed ? -steps : steps, r.SourceWrapSign);
            Assert.Equal(GreyRoomLadderImagesWitness.WrapSignLaw(n, r.SpinSteps, r.EtaSteps), r.SourceWrapSign);
        }
    }

    [Fact]
    public void Odd_Ring_Rooms_Read_The_Pages_Numbers()
    {
        // (1, 1) at the distance 2: the image of (2, 0), one spin step at even k, is lost, and the room is an upper
        // bound meeting F170's N − 1 frozen modes of the uniform point (the page reads N − 1 at N = 5, 7, 9; N = 3 is
        // this witness's read, 2 = F170's uniform count there)
        foreach (int n in new[] { 3, 5 })
        {
            var r = Read(n, true, 1, 1, 2);
            Assert.Equal(n - 1, r.RoomUpperBound);
            Assert.Equal(RingFrozenBandClaim.RungCount(n, uniform: true), r.RoomUpperBound);
            Assert.True(r.LostByRank);
            Assert.False(r.RuleProvedKept);
        }
        // at k = 1 nothing is kept off the four outer blocks: the ring's Δn = 1 floor at N = 5 is theirs alone
        // (PROOF_RING_GAP_DOMINANCE, a second route for that row)
        var middle = Read(5, true, 2, 1, 1);
        Assert.True(middle.LostByRank);
        Assert.Equal(0, middle.RoomUpperBound);
        Assert.True(Read(5, true, 2, 1, 3).KeptExactly);
        Assert.True(Read(5, true, 1, 0, 1).FillsBlock);
    }

    [Fact]
    public void A_Ring_With_Its_Wrap_Bond_At_Zero_Is_The_Open_Chain()
    {
        var w = new GreyRoomLadderImagesWitness(5, ring: true, bondCouplings: new[] { 0.75, 0.75, 0.75, 0.75, 0.0 });
        Assert.True(w.Holds, string.Join("; ", w.Failures));
        Assert.False(w.WrapBondNonzero);
        Assert.All(w.Rooms, r => Assert.True(r.RuleProvedKept && (r.KeptExactly || r.FillsBlock)));
        Assert.Equal(Read(5, false, 2, 1, 1).RoomUpperBound, Read(w, 2, 1, 1).RoomUpperBound);
        Assert.Equal(10, Read(w, 1, 1, 2).RoomUpperBound);
    }

    [Fact]
    public void Controls_The_Stagger_Dropped_And_The_Heisenberg_Chain_Lose_The_Images_By_Rank()
    {
        foreach (int n in new[] { 2, 3, 4, 5, 6 })
        {
            var w = W(n, false);
            // the unstaggered spin ladder on every route with a spin step that does not fill the corner (0, k): at
            // N = 2 the one route, (2, 0) into (1, 1), reads lost (the room there is one line and the unstaggered
            // vector another); every route read is lost, a read
            Assert.True(w.StaggerDroppedImagesRead > 0);
            Assert.Equal(w.StaggerDroppedImagesRead, w.StaggerDroppedImagesLost);
            // the uniform chain is palindromic: F140's floor ⌊N/2⌋ is a theorem for the room, the page reads exactly
            // ⌊N/2⌋ at N = 3 to 8, and the two meet: certified
            Assert.Equal(n / 2, w.HeisenbergBlock11Room);
            Assert.True(w.HeisenbergFloorHolds);
            Assert.True(w.HeisenbergBlock11Certified);
            // at N = 2 the ZZ energy is a function of the particle number, so both images survive on the
            // Heisenberg chain; from N = 3 on the ZZ term loses both on a chain with no zero bond, a read
            Assert.Equal(n >= 3, w.HeisenbergEtaImageLost);
            Assert.Equal(n >= 3, w.HeisenbergSpinImageLost);
        }
        // with zero bonds making equal dimers the one-particle ZZ diagonal is constant and the spin image survives on
        // the Heisenberg chain; nothing is asked there, and the theorems still hold
        var dimers = new GreyRoomLadderImagesWitness(4, bondCouplings: new[] { 0.75, 0.0, 0.75 });
        Assert.False(dimers.Connected);
        Assert.True(dimers.Holds, string.Join("; ", dimers.Failures));
        Assert.True(dimers.ReadsAsThePage, string.Join("; ", dimers.OpenReads));
        Assert.False(dimers.HeisenbergSpinImageLost);
        Assert.True(W(4, true).StaggerDroppedImagesRead > 0);
        Assert.Equal(W(4, true).StaggerDroppedImagesRead, W(4, true).StaggerDroppedImagesLost);
        Assert.True(W(5, true).StaggerDroppedImagesRead > 0);
        Assert.Equal(W(5, true).StaggerDroppedImagesRead, W(5, true).StaggerDroppedImagesLost);
    }

    [Fact]
    public void The_Two_Ladders_Commute_As_Integer_Matrices_To_N_Seven()
    {
        for (int n = 2; n <= 7; n++) Assert.True(GreyRoomLadderImagesWitness.LaddersCommuteAt(n), $"N={n}");
    }

    [Fact]
    public void Single_Blocks_Past_The_Sweep_Read_The_Page()
    {
        // (1, 1) at N = 7 and 8: C(N, 2) + 1, certified (the page's G7 to N = 12)
        var seven = GreyRoomLadderImagesWitness.ReadRoom(7, false, 1, 1, 2);
        Assert.Equal(21, seven.RoomUpperBound);
        Assert.True(seven.Certified);
        var eight = GreyRoomLadderImagesWitness.ReadRoom(8, false, 1, 1, 2);
        Assert.Equal(28, eight.RoomUpperBound);
        Assert.True(eight.Certified);
        Assert.True(GreyRoomLadderImagesWitness.ReadRoom(7, false, 1, 1, 0).Certified);
        Assert.True(GreyRoomLadderImagesWitness.ReadRoom(8, false, 1, 1, 0).Certified);
        // the row k = 1 in the middle block at N = 7, the free-fermion floor
        var floor = GreyRoomLadderImagesWitness.ReadRoom(7, false, 3, 4, 1);
        Assert.Equal(7, floor.RoomUpperBound);
        Assert.True(floor.Certified);
    }

    [Fact]
    public void Block_4_2_At_N_Eight_Reads_72_Against_70()
    {
        // the page's G8: the ladder image of (4, 0) with the mirror products fills the room of (4, 2) at the
        // distance 4, 72 against C(8, 4) = 70; read here as the upper bound beside the kept image
        var r = GreyRoomLadderImagesWitness.ReadRoom(8, false, 4, 2, 4);
        Assert.Equal(70, r.LadderDimension);
        Assert.True(r.KeptExactly);
        Assert.Equal(72, r.RoomUpperBound);
    }

    [Fact]
    public void Guards_Hold()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(2, ring: true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(7));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(4, gamma: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(4, j: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(4, j: 3.0));
        Assert.Throws<ArgumentException>(() => new GreyRoomLadderImagesWitness(4, bondCouplings: new[] { 1.0, 1.0 }));
        Assert.Throws<ArgumentException>(() => new GreyRoomLadderImagesWitness(4, bondCouplings: new[] { 0.0, 0.0, 0.0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoomLadderImagesWitness(4, bondCouplings: new[] { 1.0, 5.0, 1.0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GreyRoomLadderImagesWitness.ReadRoom(9, false, 1, 1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => GreyRoomLadderImagesWitness.ReadRoom(4, false, 1, 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => GreyRoomLadderImagesWitness.ReadRoom(4, false, 1, 1, 2, j: 3.0));
    }

    [Fact]
    public void The_Rule_And_The_Window_Are_The_Pages()
    {
        Assert.Equal(new[] { 0, 2 }, GreyRoomLadderImagesClaim.Window(6, 1, 1));
        Assert.Equal(new[] { 1, 3, 5 }, GreyRoomLadderImagesClaim.Window(6, 3, 2));
        Assert.Equal(new[] { 0, 2, 4 }, GreyRoomLadderImagesClaim.Window(6, 4, 4));
        Assert.Equal((1, 1), GreyRoomLadderImagesClaim.LadderSteps(6, 3, 2, 3));
        Assert.Equal((0, 2), GreyRoomLadderImagesClaim.LadderSteps(6, 4, 2, 2));
        Assert.Equal((2, 0), GreyRoomLadderImagesClaim.LadderSteps(6, 0, 2, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => GreyRoomLadderImagesClaim.LadderSteps(6, 3, 3, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => GreyRoomLadderImagesClaim.LadderSteps(6, 1, 1, 4));   // above the window's top
        // the chain keeps everything, a ring with its wrap bond at zero is the chain
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 3, 2, 3, ring: false, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 3, 2, 3, ring: true, wrapBondNonzero: false));
        // the even ring
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 2, 2, 2, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 3, 3, 4, ring: true, wrapBondNonzero: true));
        Assert.False(GreyRoomLadderImagesClaim.ImageProvedKept(6, 3, 2, 3, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 3, 0, 3, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(6, 6, 3, 3, ring: true, wrapBondNonzero: true));
        // the odd ring
        Assert.False(GreyRoomLadderImagesClaim.ImageProvedKept(5, 2, 1, 1, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(5, 2, 1, 3, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(5, 1, 0, 1, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(5, 0, 1, 1, ring: true, wrapBondNonzero: true));
        Assert.False(GreyRoomLadderImagesClaim.ImageProvedKept(5, 1, 1, 2, ring: true, wrapBondNonzero: true));
        Assert.True(GreyRoomLadderImagesClaim.ImageProvedKept(5, 1, 3, 2, ring: true, wrapBondNonzero: true));    // s = k = 2: from the corner (0, 2)
        Assert.False(GreyRoomLadderImagesClaim.ImageProvedKept(5, 1, 1, 1, ring: true, wrapBondNonzero: true));   // off the window
        // the transposes agree
        for (int p = 0; p <= 5; p++)
            for (int q = 0; q <= 5; q++)
                foreach (int k in GreyRoomLadderImagesClaim.Window(5, p, q))
                    Assert.Equal(GreyRoomLadderImagesClaim.ImageProvedKept(5, p, q, k, true, true),
                                 GreyRoomLadderImagesClaim.ImageProvedKept(5, q, p, k, true, true));
    }

    [Fact]
    public void Claim_Is_Registered_With_Its_Five_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GreyRoomLadderImagesClaim>();
        Assert.Same(registry.Get<SpectatorIntertwinerClaim>(), claim.Spectator);
        Assert.Same(registry.Get<PinnedBlockFloorClaim>(), claim.Pinned);
        Assert.Same(registry.Get<ChiralKClaim>(), claim.ChiralK);
        Assert.Same(registry.Get<FrozenDivisorClaim>(), claim.FrozenDivisor);
        Assert.Same(registry.Get<RingFrozenBandClaim>(), claim.RingBand);
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("F171", claim.Name);
        Assert.Contains("measured", claim.Name);
        Assert.Contains("PROOF_GREY_GLASS_LADDER_IMAGES", claim.Anchor);
    }
}
