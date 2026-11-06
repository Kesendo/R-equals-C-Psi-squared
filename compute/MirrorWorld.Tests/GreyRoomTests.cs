using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the grey room (adopted 2026-11-06 from docs/THE_GREY_GLASS.md and
// docs/proofs/PROOF_GREY_GLASS_LADDER_IMAGES.md, F171). Everything here is a rank by elimination over GF(p) at
// two primes: no eigensolver anywhere, the Seed genre. The literals are the page's numbers, read there by the
// committed gates simulations/grey_glass_gate.py (G7, G8, G10) and asserted BESIDE the ladder theorem's closed
// form C(N, k), so a gate built only from the formula cannot pass by construction. Where the page states no
// number the test asserts only what must hold between the two bounds.
public class GreyRoomTests
{
    static readonly World W = new();

    static GreyRoom Chain(int n, bool zz = false) => new(W, n, Enumerable.Repeat(1L, n - 1).ToArray(), ring: false, heisenberg: zz);
    static GreyRoom RingOf(int n, bool zz = false) => new(W, n, Enumerable.Repeat(1L, n).ToArray(), ring: true, heisenberg: zz);

    [Fact]
    public void The_Window_Runs_From_The_Difference_To_The_Smaller_Sum_In_Steps_Of_Two()
    {
        var g = Chain(6);
        Assert.Equal(new[] { 0, 2 }, g.Window(1, 1));
        Assert.Equal(new[] { 1, 3, 5 }, g.Window(3, 2));
        Assert.Equal(new[] { 2 }, g.Window(2, 0));
        Assert.Equal(new[] { 0, 2, 4 }, g.Window(4, 4));   // 2N - p - q = 4 caps the top
        Assert.Equal(4 * 6, Enumerable.Range(0, 7).Sum(p => Enumerable.Range(0, 7).Count(q => g.IsCorner(p, q))));
    }

    // the (1, 1) block of the uniform open XY chain: the two ladder images and nothing else, C(N, 2) + 1,
    // proved at every N (the page's G7 counts it at N = 3 to 12).
    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    public void Uniform_Open_Chain_Block_1_1_Holds_Exactly_The_Two_Ladder_Images(int n)
    {
        var g = Chain(n);
        Assert.Equal(GreyRoom.Choose(n, 2), g.Room(1, 1, 2));
        Assert.Equal(1, g.Room(1, 1, 0));
        Assert.True(g.Certified(1, 1, 2));
    }

    [Fact]
    public void Eleven_At_N_Five() => Assert.Equal(11, Chain(5).Room(1, 1, 2) + Chain(5).Room(1, 1, 0));

    // on every other block of the chain the room meets the ladder bound at N = 3 to 6, everywhere but at the
    // height N/2 strictly inside the window (G8): the two blocks where it does not are pinned below.
    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)]
    public void Uniform_Open_Chain_Rooms_Meet_The_Ladder_Bound_Off_The_Middle_Height(int n)
    {
        var g = Chain(n);
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
                foreach (int k in g.Window(p, q))
                {
                    int room = g.Room(p, q, k), ladder = g.LadderBound(p, q, k);
                    Assert.True(room >= ladder, $"({p},{q}) k={k}: room {room} below the proved bound {ladder}");
                    bool middle = n % 2 == 0 && k == n / 2 && Math.Abs(p - q) < k && k < Math.Min(p + q, 2 * n - p - q);
                    if (!middle) Assert.Equal(ladder, room);
                }
    }

    [Fact]
    public void The_Reflection_Builds_Two_More_At_The_Middle_Height()
    {
        var g4 = Chain(4);
        Assert.Equal(8, g4.Room(2, 2, 2));          // 8 against the ladders' 6 (G8)
        Assert.Equal(6, g4.LadderBound(2, 2, 2));
        Assert.Equal(4, g4.Products(2, 2, 2));        // a product and its transpose land together in the middle block; at N = 4 the transposes fall into the ladder span, so 8 = 6 + 2
        var g6 = Chain(6);
        Assert.Equal(22, g6.Room(3, 2, 3));         // 22 against 20, and its transpose
        Assert.Equal(22, g6.Room(2, 3, 3));
        Assert.Equal(22, g6.Room(4, 3, 3));         // (4, 3) is the fold of (3, 2) up to the transpose
        Assert.Equal(20, g6.LadderBound(3, 2, 3));
        Assert.Equal(2, g6.Products(3, 2, 3));
        Assert.Equal(2, g6.Products(3, 4, 3));
        Assert.Equal(0, g6.Products(3, 2, 1));
        Assert.Equal(0, g6.Products(3, 3, 2));     // (3, 3): neither side is 2j, and the height is not N/2
        Assert.Equal(0, g6.Products(3, 3, 3));     // off the window of (3, 3), which is {0, 2, 4, 6}
    }

    [Fact]
    public void Couplings_The_Reflection_Does_Not_Keep_Fall_To_The_Ladder_Count()
    {
        var g = new GreyRoom(W, 4, new long[] { 1, 2, 3 });
        Assert.False(g.ReflectionKeepsCouplings());
        Assert.Equal(0, g.Products(2, 2, 2));
        Assert.Equal(6, g.Room(2, 2, 2));
        Assert.True(g.Certified(2, 2, 2));
    }

    // the even ring: the ladder images at even k are kept at any real couplings, at odd k only in a corner; the
    // room at k = N/2 in (N/2, 2) is 11 at N = 4 and the bond axes' 6 at N = 6 on the uniform ring (G10).
    [Fact]
    public void Even_Ring_Keeps_Even_Heights_And_Odd_Heights_Only_In_The_Corners()
    {
        var g6 = RingOf(6);
        Assert.True(g6.Kept(2, 2, 2));
        Assert.True(g6.Kept(3, 3, 4));
        Assert.False(g6.Kept(3, 2, 3));
        Assert.True(g6.Kept(3, 0, 3));
        Assert.True(g6.Kept(6, 3, 3));
        Assert.Equal(6, g6.Room(3, 2, 3));          // the mirror products alone, below C(6, 3) = 20
        Assert.Equal(0, g6.LadderBound(3, 2, 3));
        Assert.True(g6.Room(2, 2, 2) >= g6.LadderBound(2, 2, 2));
        Assert.Equal(15, g6.LadderBound(2, 2, 2));
        var g4 = RingOf(4);
        Assert.Equal(11, g4.Room(2, 2, 2));         // two more than the ladder images with the mirror products of both bond axes (9), from the antipodal pairs
    }

    // the odd ring: an image is kept exactly from a corner by the ladder its parity allows; at k = 1 nothing is
    // kept off the four outer blocks, which is the ring's Delta n = 1 floor (PROOF_RING_GAP_DOMINANCE).
    [Theory]
    [InlineData(5)] [InlineData(7)]
    public void Odd_Ring_Keeps_An_Image_Only_From_A_Corner(int n)
    {
        var g = RingOf(n);
        Assert.False(g.Kept(2, 1, 1));              // r = 1, neither 0 nor N - 1
        Assert.True(g.Kept(2, 1, 3));               // r = 0 at odd k
        Assert.True(g.Kept(1, 0, 1));
        Assert.True(g.Kept(0, 1, 1));
        Assert.Equal(0, g.Room(2, 1, 1));           // the floor sits in the four outer blocks alone
        Assert.True(g.Room(2, 1, 3) >= GreyRoom.Choose(n, 3));
        Assert.False(g.Kept(1, 1, 2));              // S- of (2, 0), one step at even k: not from a corner
    }

    // the block (1, 1) elsewhere (the page's paragraph, G7): the odd ring holds its N - 1 frozen modes, the even
    // ring the odd half with the U-even frozen modes beside it, 2, 2 and 4 of them at N = 4, 6, 8; the open
    // Heisenberg chain holds floor(N/2), the frozen modes of the uniform point.
    [Theory]
    [InlineData(5, 4)] [InlineData(7, 6)]
    public void Odd_Ring_Block_1_1_Is_Its_Frozen_Modes(int n, int room) =>
        Assert.Equal(room, RingOf(n).Room(1, 1, 2));

    [Theory]
    [InlineData(4, 8)] [InlineData(6, 17)]
    public void Even_Ring_Block_1_1_Is_The_Odd_Half_With_The_U_Even_Frozen_Modes(int n, int room) =>
        Assert.Equal(room, RingOf(n).Room(1, 1, 2));

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    public void Open_Heisenberg_Chain_Block_1_1_Holds_Floor_N_Over_2(int n)
    {
        var g = Chain(n, zz: true);
        Assert.Equal(n / 2, g.Room(1, 1, 2));
        Assert.Equal(1, g.Room(1, 1, 0));
        Assert.Equal(0, g.LadderBound(1, 1, 2));    // an XY theorem: the ZZ term loses the spin ladder
    }

    [Theory]
    [InlineData(4)] [InlineData(5)] [InlineData(6)]
    public void Heisenberg_Ring_Block_1_1_Is_The_XY_Rings(int n) =>
        Assert.Equal(RingOf(n).Room(1, 1, 2), RingOf(n, zz: true).Room(1, 1, 2));

    // the reading is scale-free and sign-free in the coupling (one coupling, the XY book)
    [Fact]
    public void The_Room_Is_Scale_And_Sign_Free()
    {
        foreach (long j in new[] { 1L, 3L, -1L, 7L })
        {
            var g = new GreyRoom(W, 5, Enumerable.Repeat(j, 4).ToArray());
            Assert.Equal(10, g.Room(1, 1, 2));
            Assert.Equal(10, g.Room(2, 1, 3));
        }
    }

    // a ring whose wrap bond is zero is the open chain: F171 keeps every image there, and the rooms are the chain's
    [Fact]
    public void A_Ring_With_Its_Wrap_Bond_At_Zero_Is_The_Open_Chain()
    {
        var r6 = new GreyRoom(W, 6, new long[] { 1, 1, 1, 1, 1, 0 }, ring: true);
        Assert.True(r6.Kept(3, 2, 3));
        Assert.Equal(22, r6.Room(3, 2, 3));
        Assert.Equal(Chain(6).Room(3, 2, 3), r6.Room(3, 2, 3));
        var r5 = new GreyRoom(W, 5, new long[] { 1, 1, 1, 1, 0 }, ring: true);
        Assert.True(r5.Kept(2, 1, 1));
        Assert.Equal(5, r5.Room(2, 1, 1));
    }

    [Fact]
    public void A_Zero_Bond_Grows_The_Room_And_Is_Not_Fenced_Away()
    {
        // two open halves: the commutant of a disconnected hopping is larger, and the count simply reports it
        var g = new GreyRoom(W, 4, new long[] { 1, 0, 1 });
        Assert.True(g.Room(1, 1, 2) >= GreyRoom.Choose(4, 2));
        Assert.True(g.Room(1, 1, 0) > 1);
    }

    // a reading of this object alone, on no gate of the page (one route, the mod-p upper bound): the uniform ring
    // of six holds more than the ladders in its other coloured blocks too; pinned so it cannot drift silently,
    // and the page's Open list carries it as a reading
    [Fact]
    public void Ring_Of_Six_Reads_More_Than_The_Ladders_In_The_Blocks_2_2_And_3_3()
    {
        var g = RingOf(6);
        Assert.Equal(19, g.Room(2, 2, 2));
        Assert.Equal(15, g.LadderBound(2, 2, 2));
        Assert.Equal(18, g.Room(3, 3, 2));
        Assert.Equal(18, g.Room(3, 3, 4));
        Assert.Equal(15, g.LadderBound(3, 3, 4));
    }

    [Fact]
    public void Guards_Hold()
    {
        Assert.Throws<ArgumentException>(() => new GreyRoom(W, 1, Array.Empty<long>()));
        Assert.Throws<ArgumentException>(() => new GreyRoom(W, 4, new long[] { 1, 1 }));
        Assert.Throws<ArgumentException>(() => new GreyRoom(W, 4, new long[] { 1, 1, 1 }, ring: true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoom(W, 3, new long[] { 1, 1L << 41 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Chain(4).Room(5, 0, 1));
        Assert.Equal(0, Chain(4).Room(1, 1, 1));   // a height off the window holds no cell
        Assert.Equal(new[] { "room", "images", "kept", "products" }, Chain(4).Own);
        Assert.Throws<ArgumentException>(() => new GreyRoom(W, 2, new long[] { 1, 1 }, ring: true));
        Assert.Throws<ArgumentOutOfRangeException>(() => new GreyRoom(W, 9, Enumerable.Repeat(1L, 8).ToArray()));
        Assert.Equal(new[] { "x", "y", "z" }, Chain(4).Inherited);
    }
}
