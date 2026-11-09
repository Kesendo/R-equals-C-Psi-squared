using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the carrier (built 2026-11-08). The expected sheets and turns are the census's,
// simulations/results/far_element_sheets_census.txt (part 1, the carriers of PROOF_ONE_DEPHASED_SEAT),
// the sheets read there on framework's far_element_shape over Q(i) and the turns by the census's own turns() on that
// element, a route sharing no code with this object: sheets A 4, B 2, C 1, D 3, E 0, E's corner 0, the N = 4 family 2;
// turns A X+ Y- Z-, B Z-, C Z-, D, E and E's corner none, the family X+. C's twelve strings are
// blind_seat_sheet_gate.py's census row.
public class CarrierTests
{
    static readonly World Frame = new();

    static string One(int n, int s, char c) { var a = new string('I', n).ToCharArray(); a[s] = c; return new string(a); }
    static string Two(int n, int i, int j, char c) { var a = new string('I', n).ToCharArray(); a[i] = c; a[j] = c; return new string(a); }

    static EndCount Row(int n, (int, int, long)[] bonds, (int, long, long)[] fields)
    {
        var h = new List<(string, long)>();
        foreach (var (i, j, w) in bonds)
            if (w != 0) foreach (char c in "XYZ") h.Add((Two(n, i, j, c), w));
        foreach (var (s, x, y) in fields)
        {
            if (x != 0) h.Add((One(n, s, 'X'), x));
            if (y != 0) h.Add((One(n, s, 'Y'), y));
        }
        return new EndCount(Frame, n, h, new[] { One(n, 0, 'Z') });
    }

    static (int, int, long)[] Tri(long du, long dv, long uv) => new[] { (0, 1, du), (0, 2, dv), (1, 2, uv) };

    static EndCount CensusRow(string row) => row switch
    {
        "A" => Row(3, Tri(1, 2, 3), new[] { (0, 1L, 0L), (1, 2L, 0L), (2, 3L, 0L) }),
        "B" => Row(3, Tri(2, 2, 3), new[] { (0, 1L, 1L), (1, 3L, 4L), (2, 4L, 3L) }),
        "C" => Row(3, Tri(3, 6, -2), new[] { (1, 3L, 4L), (2, 5L, 0L) }),
        "D" => Row(3, Tri(1, 0, 3), new[] { (0, 1L, 0L), (1, -1L, -2L), (2, 0L, 2L) }),
        "E" => Row(3, Tri(1, 2, 6), new[] { (1, 8L, 0L), (2, -5L, 3L) }),
        "E's corner" => Row(3, Tri(2, -2, 0), new[] { (1, 1L, 4L), (2, 4L, 1L) }),
        "N = 4 family" => Row(4, new[] { (1, 2, 2L), (0, 1, 1L), (0, 2, -1L), (1, 3, -1L), (2, 3, 1L) },
            new[] { (1, -1L, 0L), (2, 1L, 0L), (3, 0L, 3L) }),
        _ => throw new ArgumentException(row),
    };

    static string TurnText(IReadOnlyDictionary<char, int?> t)
    {
        var homogeneous = "XYZ".Where(p => t[p] is not null).Select(p => $"{p}{(t[p] == 1 ? '+' : '-')}").ToList();
        return homogeneous.Count == 0 ? "none" : string.Join(" ", homogeneous);
    }

    [Theory]
    [InlineData("A", 4, "X+ Y- Z-")]
    [InlineData("B", 2, "Z-")]
    [InlineData("C", 1, "Z-")]
    [InlineData("D", 3, "none")]
    [InlineData("E", 0, "none")]
    [InlineData("E's corner", 0, "none")]
    [InlineData("N = 4 family", 2, "X+")]
    public void The_Census_Carriers_Have_Their_Sheets_And_Turns(string row, int sheets, string turns)
    {
        var c = Carrier.Of(CensusRow(row));
        Assert.NotNull(c);
        Assert.Equal((sheets, turns), (c!.Sheets, TurnText(c.Turns)));
    }

    // A constant term of H moves nothing in L, so it moves no reading. On this row (W = YZZ, a colouring) ZZZ and XII
    // reverse every term of H, and a constant counted among the terms would drop both from the subgroup and count a
    // sheet more: the sheets are 3 with 5·III and without it, as framework's far_element_shape reads both rows.
    [Fact]
    public void A_Constant_Term_Moves_No_Reading()
    {
        var h = new List<(string, long)> { ("YIZ", -1), ("YXX", 4), ("YII", -3), ("YYY", -5), ("ZYI", -3) };
        var plain = Carrier.Of(new EndCount(Frame, 3, h, new[] { "ZII" }))!;
        var shifted = Carrier.Of(new EndCount(Frame, 3, h.Append(("III", 5L)).ToList(), new[] { "ZII" }))!;
        Assert.Equal((3, 3), (plain.Sheets, shifted.Sheets));
    }

    // Docked, not inherited: the points read off the masks are LetterCube's coordinates, string for string, on every
    // census carrier, and the letter counts (k_Z the letters X or Y, k_X the letters Y or Z, k_Y the letters X or Z).
    [Fact]
    public void The_Points_Are_The_Cube_Coordinates()
    {
        int Count(string s, string letters) => s.Count(letters.Contains);
        foreach (var row in new[] { "A", "B", "C", "D", "E", "E's corner", "N = 4 family" })
        {
            var c = Carrier.Of(CensusRow(row))!;
            Assert.Equal(c.Element.Select(t => LetterCube.Coordinates(t.Letters)), c.Points);
            Assert.Equal(c.Element.Select(t => new LetterCube.Point(Count(t.Letters, "XY"), Count(t.Letters, "YZ"),
                                                                    Count(t.Letters, "XZ"))), c.Points);
        }
    }

    // A turn is a reading of the points: the sign of the turn by Z^N is (-1)^{k_Z}, so a Z-turn that reverses W means
    // every height is odd, and C's strings sit at the heights 1 and 3.
    [Fact]
    public void C_Has_Twelve_Strings_Of_Odd_Height()
    {
        var c = Carrier.Of(CensusRow("C"))!;
        Assert.Equal(12, c.Element.Count);
        Assert.Equal(new[] { 1, 3 }, c.Heights);
        Assert.All(c.Points, p => Assert.Equal(1, p.KZ % 2));
    }

    // No carrier where the row does not pair: C's couplings detuned (near 1, far 0), the u, v singlet at J_du = J_dv
    // with parallel fields of one length and a field on d off them (near 2, far 1); and none where the ends are larger
    // than one: the triangle with no fields, coloured by X^3 and by Y^3 alike, both commuting with every bond, and a
    // row with one colouring, XXX under X fields of one strength, whose ends are two (the counts are read before the
    // colouring is).
    [Fact]
    public void A_Row_Without_One_Carrier_Has_None()
    {
        var detuned = Row(3, Tri(3, 6, -3), new[] { (1, 3L, 4L), (2, 5L, 0L) });
        Assert.Equal((1, 0), detuned.UpperCounts());
        Assert.Null(Carrier.Of(detuned));
        var singlet = Row(3, Tri(2, 2, 3), new[] { (0, 1L, 2L), (1, 3L, 1L), (2, 3L, 1L) });
        Assert.Equal((2, 1), singlet.UpperCounts());
        Assert.Null(Carrier.Of(singlet));
        var bare = Row(3, Tri(1, 2, 3), Array.Empty<(int, long, long)>());
        var (near, far) = bare.UpperCounts();
        Assert.True(near == far && far > 1, $"{near} {far}");
        Assert.Null(Carrier.Of(bare));
        var oneColour = Row(3, Tri(2, 2, 3), new[] { (0, 1L, 0L), (1, 1L, 0L), (2, 1L, 0L) });
        Assert.Equal((2, 2), oneColour.UpperCounts());
        Assert.Equal(new[] { "XXX" }, oneColour.Colourings().Select(s => s.ToString(3)));
        Assert.Null(Carrier.Of(oneColour));
    }

    // Past the lift (one jump, N = 6: a lit span of 2048 strings) the row's one colouring is its far element: the
    // Heisenberg chain with distinct X fields, coloured by X^6, one string at the height 6 that every turn keeps.
    [Fact]
    public void Past_The_Lift_A_Single_Colouring_Is_The_Carrier()
    {
        var chain = Row(6, Enumerable.Range(0, 5).Select(i => (i, i + 1, 1L)).ToArray(),
            Enumerable.Range(0, 6).Select(i => (i, (long)(i + 1), 0L)).ToArray());
        Assert.Equal((1, 1), chain.UpperCounts());
        Assert.Null(chain.FarElement());
        var c = Carrier.Of(chain)!;
        Assert.Equal(new[] { ("XXXXXX", 1L) }, c.Element);
        Assert.Equal("X+ Y+ Z+", TurnText(c.Turns));
        Assert.Equal(new[] { 6 }, c.Heights);
    }

    // The guard: the group is found by enumerating every string, so a row past MaxSites is refused, not read.
    [Fact]
    public void Past_MaxSites_The_Carrier_Refuses()
    {
        const int n = Carrier.MaxSites + 1;
        var h = Enumerable.Range(0, n - 1).SelectMany(i => "XYZ".Select(c => (Two(n, i, i + 1, c), 1L))).ToList();
        var count = new EndCount(Frame, n, h, Enumerable.Range(0, n).Select(s => One(n, s, 'Z')).ToArray());
        Assert.Throws<ArgumentOutOfRangeException>(() => Carrier.Of(count));
    }

    // the two buckets: the carrier owns its readings and hangs on the end count, whose readings arrive as inheritance;
    // the jumps it reads are the count's.
    [Fact]
    public void Ontology_The_Carrier_Owns_Its_Readings_And_Hangs_On_The_End_Count()
    {
        var count = CensusRow("C");
        Assert.Equal(new[] { "ZII" }, count.Jumps);
        var c = Carrier.Of(count)!;
        Assert.Equal(new[] { "points", "turns", "sheets" }, c.Own);
        Assert.IsType<EndCount>(c.Parent);
        Assert.Contains("near", c.Inherited);
        Assert.Contains("far", c.Inherited);
        Assert.Contains("x", c.Inherited);
    }
}
