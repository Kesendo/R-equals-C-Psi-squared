using System.Numerics;

namespace MirrorWorld;

// The carrier (built 2026-11-08 beside docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md, the paragraphs "The
// carriers' locality, and their sheets" and "The sheet and blindness, at N = 3 and N = 4", and its census
// simulations/far_element_sheets_census.py; the name is the proof's word, narrowed as the fence below says). Where a
// row pairs through ONE far element, near = far = 1 in F158's count, that element W carries the palindrome: every
// mirror onto -L - 2*sigma sends the identity onto W's line. EndCount finds W and certifies it; this object holds W as
// a thing and reads what can be read off it, two things on the letter cube and one beside it, none of which needs an
// eigensolver:
//
//   points  every string's place (k_Z, k_X, k_Y) on the cube, the popcounts of its x mask, its z mask and x XOR z,
//           read off PauliString's masks here and docked in the tests against LetterCube's coordinates and the
//           letter counts, the way EndCount docks what it meets (the cube is not this object's parent: the element is
//           the end count's);
//   turns   for P = X, Y, Z whether W is homogeneous under the turn by P^N (conjugation by P on every site, the sign
//           (-1)^{k_P} on a string, F63's three turns), and the sign it then puts on every string of W: +1 where the
//           turn keeps W and -1 where it reverses W; the proof reads two such turns as sheets, C's one sheet (Z^N
//           reverses C's W, F61's parity) and one of the N = 4 family's two (X^N keeps it);
//   sheets  the Pauli strings g with g W g = +-W form a group, and its GF(2) dimension beyond the subgroup that the
//           jumps and the strings keeping every term of H or reversing every one (g H g = +-H) generate is the number
//           of W's sheets, the carriers' paragraph's count. A sheet can show as a turn (C's sheet Z_uZ_v is Z^N
//           modulo the jump, the parity of the height), but the count is site-resolved GF(2) algebra beside the cube:
//           the N = 4 family's other sheet, IIIY, is one letter on one site and no turn, and D's group, three sheets
//           beyond its subgroup, holds no turn at all. The far end being one-dimensional, every string of that
//           subgroup keeps or reverses W, so the subgroup lies inside the group; the object checks this theorem, a
//           failure being a bug. A constant term of H moves nothing in L and is left out of the subgroup's test, as
//           framework's far_element_shape leaves it out.
//
// What stays outside: whether blindness brings a sheet is the main repo's question, read on finite grids by the
// census; the rule of "The sheet and blindness" for far elements of odd height lives in the main repo's battery
// (TrianglePalindromeLocusWitness, node "the sheet and blindness"); the locality half of the carriers' paragraph
// stays outside too, its operator Schmidt ranks in Diagnostics' OperatorSchmidtRank, the triangle witness and
// framework's far_element_shape, and which site permutations make W a product in the census. This object reads, it
// does not decide; it neither searches for a mirror (the Sextant's) nor reads site-local mirrors (the Theodolite's).
//
// Words, fenced at the door. A CARRIER here is the far element of one row whose far end is one-dimensional, narrower
// than the proof's word, which calls every invertible far element a carrier (it calls XXX one at every coupling of the
// triangle without fields, which XXX and YYY both colour, so its far end is larger than one and it has none here).
// The word has many other owners, in this world and in the main repo, typed and untyped, and none of them is meant;
// it is in check_spent_words. A SHEET is a homogeneity of W, not F127's sheet lattice, a Riemann sheet of the EP
// documents or a sheet of the quantum-classical fold. A TURN is F63's, conjugation by P on every site, the map that
// LetterCube's table names the turn by P^N, MirrorGroup holds as Ad_{Z^N} and Ad_{X^N}, EndCount reads among its
// cube rotations as the pi rotations about the letter axes, and the Lattice runs for X^N as its double turn; the
// world uses the word for other maps as well (LetterCube's half-turns permute the letters, GammaFold's turns act on
// gamma), and none of them is meant. The HEIGHT is k_Z, as the proof says it; LetterCube's Height is the rate a
// dissipator charges a string, 2·gamma·k_Z under Z dephasing at rate gamma.
public sealed class Carrier : GameObject
{
    /// <summary>The largest N read: the group is found by enumerating all 4^N strings.</summary>
    public const int MaxSites = 10;

    public int N { get; }
    public IReadOnlyList<(string Letters, long Coefficient)> Element { get; }
    readonly IReadOnlyList<PauliString> hamiltonian;   // H's terms without a constant
    readonly IReadOnlyList<PauliString> jumps;
    readonly Lazy<(int Group, int Subgroup)> dims;

    // parent = the end count: the element is what its far end holds, inherited, and the setting with it.
    Carrier(EndCount count, IReadOnlyList<(string Letters, long Coefficient)> element) : base(count)
    {
        N = count.N;
        Element = Array.AsReadOnly(element.ToArray());
        hamiltonian = count.Terms.Select(t => PauliString.Parse(t.Letters)).Where(s => !s.IsIdentity).ToList();
        jumps = count.Jumps.Select(PauliString.Parse).ToList();
        dims = new Lazy<(int, int)>(ReadDims);
    }

    /// <summary>The carrier of a row, where the row pairs through one far element: the count's upper bounds are
    /// near = far = 1 and either its lifted far vector certifies as an exact invertible element or, where it does not
    /// lift (a lit span past EndCount.MaxLiftColumns, N = 6 to 8 under one jump, where EndCount's span cap ends the
    /// count, or a vector past the reconstruction), the row's one colouring is that element. Null otherwise: a row
    /// that does not pair has no carrier, a larger far end leaves the identity's image open, and a far element that
    /// neither lifts nor is a single colouring is not read. Throws past MaxSites, and EndCount's span cap throws first
    /// where the jumps are too few (every one-jump row from N = 9).</summary>
    public static Carrier? Of(EndCount count)
    {
        if (count.N > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(count), count.N, $"the carrier reads up to {MaxSites} sites.");
        if (count.UpperCounts() != (1, 1)) return null;
        if (count.FarElement() is { Check.Certifies: true } f) return new Carrier(count, f.Element);
        var colours = count.Colourings();
        return colours.Count == 1 ? new Carrier(count, new[] { (colours[0].ToString(count.N), 1L) }) : null;
    }

    public IReadOnlyList<LetterCube.Point> Points => Element.Select(t => PauliString.Parse(t.Letters))
        .Select(s => new LetterCube.Point(BitOperations.PopCount(s.X), BitOperations.PopCount(s.Z),
                                          BitOperations.PopCount(s.X ^ s.Z))).ToList();

    /// <summary>The heights k_Z the strings of W take, sorted and distinct.</summary>
    public IReadOnlyList<int> Heights => Points.Select(p => p.KZ).Distinct().OrderBy(k => k).ToList();

    /// <summary>For P in X, Y, Z: the sign (-1)^{k_P} every string of W shares, +1 where the turn by P^N keeps W and
    /// -1 where it reverses W, or null where W is not homogeneous under it.</summary>
    public IReadOnlyDictionary<char, int?> Turns
    {
        get
        {
            var pts = Points;
            int? Turn(Func<LetterCube.Point, int> k)
            {
                var signs = pts.Select(p => k(p) % 2 == 0 ? 1 : -1).Distinct().ToList();
                return signs.Count == 1 ? signs[0] : null;
            }
            return new Dictionary<char, int?> { ['X'] = Turn(p => p.KX), ['Y'] = Turn(p => p.KY), ['Z'] = Turn(p => p.KZ) };
        }
    }

    /// <summary>The number of W's sheets: dim of { g : g W g = +-W } minus dim of the subgroup the jumps and the
    /// strings with g H g = +-H generate, both over GF(2).</summary>
    public int Sheets => dims.Value.Group - dims.Value.Subgroup;

    static bool Homogeneous(PauliString g, IEnumerable<PauliString> set) =>
        set.Select(s => PauliString.Commute(g, s)).Distinct().Count() <= 1;

    (int, int) ReadDims()
    {
        var w = Element.Select(t => PauliString.Parse(t.Letters)).ToList();
        var group = new List<PauliString>();
        var subgroup = new List<PauliString>(jumps);
        for (ulong x = 0; x < 1UL << N; x++)
        for (ulong z = 0; z < 1UL << N; z++)
        {
            var g = new PauliString(x, z);
            if (Homogeneous(g, w)) group.Add(g);
            if (Homogeneous(g, hamiltonian)) subgroup.Add(g);
        }
        if (!subgroup.All(g => Homogeneous(g, w)))
            throw new InvalidOperationException(
                "a string of the subgroup (a jump, or a string keeping or reversing H) does not keep or reverse W: " +
                "the far end cannot be one-dimensional, so this is a bug.");
        return (Dim(group), Dim(subgroup));
    }

    // GF(2) dimension of the span of strings, each a vector (x mask, z mask) of 2N bits.
    int Dim(IEnumerable<PauliString> strings)
    {
        var basis = new List<ulong>();
        foreach (var g in strings)
        {
            ulong v = g.X | (g.Z << N);
            foreach (var b in basis) v = Math.Min(v, v ^ b);
            if (v != 0) { basis.Add(v); basis.Sort((a, c) => c.CompareTo(a)); }
        }
        return basis.Count;
    }

    // left: what the carrier itself produces.
    public override IReadOnlyList<string> Own => new[] { "points", "turns", "sheets" };
}
