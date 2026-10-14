using System.Numerics;

namespace MirrorWorld;

// The letter cube (adopted 2026-10-14 from docs/THE_ONE_SQUARE.md §7 and §8, the cube of the three one-site squares;
// gates simulations/diagonal_square_gate.py S8 and S9, the Python primitive simulations/framework/letter_cube.py).
//
// Each dephasing letter P in {Z, X, Y} has its own one-site square, and a letter sits on P's anti-diagonal exactly
// when it anticommutes with P. Give each letter one coordinate per square: I = (0,0,0), X = (1,0,1), Y = (1,1,0),
// Z = (0,1,1) in the order (k_Z, k_X, k_Y). On N sites a Pauli string sits at
//
//     (k_Z, k_X, k_Y) = (n_X + n_Y, n_Y + n_Z, n_X + n_Z),
//
// k_P the number of sites whose letter anticommutes with P. Read off the masks of PauliString: k_Z = |x|,
// k_X = |z|, k_Y = |x XOR z|. The strings fill the even lattice points of the tetrahedron with vertices I^N, X^N,
// Y^N, Z^N, C(N+3, 3) of them; w = (k_Z + k_X + k_Y)/2 is the number of non-identity letters. The three
// dephasing diagonals Q_P = sum_l P_l . P_l are the three heights: Q_P(sigma) = (N - 2 k_P) sigma. So dephasing
// along every letter, rate gamma_{l,P} on site l, charges a string the height
//
//     2 sum_l sum_P gamma_{l,P} [the letter on site l anticommutes with P],
//
// which with the same rates on every site is 2(gamma_Z k_Z + gamma_X k_X + gamma_Y k_Y): the dissipator is a
// height and moves nothing. The four one-site rates pair under a half-turn about P's axis (I with P, the other
// two with each other) exactly when gamma_P = 0. R, the right multiplication by X^N (I <-> X, Y <-> Z), carries
// (k_Z, k_X, k_Y) to (N - k_Z, k_X, N - k_Y): the populations' edge k_Z = 0 onto the far edge k_Z = N. The
// centre (N/2, N/2, N/2) holds a string only when 4 | N. The parity cube of the hardness bloc is this cube read
// mod 2 in its first two coordinates, with the third read as half a sum: y_par = (k_Z + k_X - k_Y)/2 mod 2,
// which is no coordinate's parity.
//
// The Hamiltonian moves strings (§8, gate S9). A placed Pauli term T gives [T, sigma] = 0 when the two commute and
// 2 i^k Target when they anticommute, T sigma = i^k Target with Target the phase-free product, one string; a two-site term changes w by exactly one, a bond P (x) P steps
// +-2 along one of the two axes other than P, and a one-site field keeps w and steps on a face diagonal of the plane
// k_P = const. The Pauli dissipator never moves a string, so every move a string makes on the cube is the
// Hamiltonian's.
//
// A jump that is not a Pauli string moves strings too (§9, the amplitude-damping paragraph, gates M8, M9). Any
// one-site jump J at rate r acts on the four letters through its Pauli-transfer matrix T[a, b] = Tr(P_a D(P_b))/2,
// D(rho) = r (J rho J^dagger - {J^dagger J, rho}/2): its diagonal is a height (minus each letter's rate), its
// off-diagonal entries are arrows, a move from letter b to letter a with a coefficient. sigma^- = |0><1| at rate r is
// the lights r/4 along X and along Y (heights r/2, r/2, r on X, Y, Z) plus one arrow I -> Z with coefficient +r;
// sigma^+ carries the arrow -r; at equal loss and gain the arrows cancel and the lights remain. Computed exactly in
// Gaussian integers: 4T is integral for an integral J and rate. The move table of §9 asked of a given system is not
// here.
//
// Deliberately outside: spectra, eigenvalue counts on a plane, centroid trajectories (paths, the main repo's).
//
// Parent: the World. The three axes are "anticommutes with the frame's z, x, y"; the frame is what this object
// inherits.
//
// Words, fenced at the door: "cube" here is the integer letter cube; the parity cube is Hardness.Cube and its
// three bits, MirrorGroup's "cube" the characters, EndCount.CubeRotations the 24 frame rotations. "Height" is
// the rate a Pauli dissipator charges a string, not a level of H. "Edge" is an edge of the tetrahedron, not
// Crack's band edge, Lattice's bridge between worlds or a bond of a graph. "Half-turn" permutes the four letters
// (I with P, the other two with each other), not EndCount.HalfTurnAxis's pi rotations, which sign letters. "Centre"
// is the cube's point (N/2, N/2, N/2), not BlindSeat's seat, Divisor's room or Atmosphere's block centres.
public sealed class LetterCube : GameObject
{
    public LetterCube(World world) : base(world) { }

    // left: what the cube itself produces.
    public override IReadOnlyList<string> Own => new[] { "coordinates", "height", "step", "arrow" };

    public readonly record struct Point(int KZ, int KX, int KY)
    {
        public int Weight => (KZ + KX + KY) / 2;
        public override string ToString() => $"({KZ}, {KX}, {KY})";
    }

    // (k_Z, k_X, k_Y): how many letters anticommute with Z, with X, with Y.
    public static Point Coordinates(PauliString s)
        => new(BitOperations.PopCount(s.X), BitOperations.PopCount(s.Z), BitOperations.PopCount(s.X ^ s.Z));

    public static Point Coordinates(string letters) => Coordinates(PauliString.Parse(letters));

    // w: the number of non-identity letters.
    public static int Weight(PauliString s) => BitOperations.PopCount(s.X | s.Z);

    // every string on n sites, in mask order
    public static IEnumerable<PauliString> Strings(int n)
    {
        if (n < 1 || n > 12)
            throw new ArgumentOutOfRangeException(nameof(n), n, "1 to 12 sites (4^n strings are enumerated)");
        ulong size = 1UL << n;
        for (ulong x = 0; x < size; x++)
            for (ulong z = 0; z < size; z++)
                yield return new PauliString(x, z);
    }

    // the lattice points the strings reach, with how many strings sit at each
    public static IReadOnlyDictionary<Point, long> Points(int n)
    {
        var count = new Dictionary<Point, long>();
        foreach (var s in Strings(n))
        {
            var p = Coordinates(s);
            count[p] = count.TryGetValue(p, out var c) ? c + 1 : 1;
        }
        return count;
    }

    // the even lattice points of the tetrahedron scaled by n: n_X, n_Y, n_Z >= 0 with n_X + n_Y + n_Z <= n
    public static IReadOnlyList<Point> TetrahedronPoints(int n)
    {
        var pts = new List<Point>();
        for (int nx = 0; nx <= n; nx++)
            for (int ny = 0; nx + ny <= n; ny++)
                for (int nz = 0; nx + ny + nz <= n; nz++)
                    pts.Add(new Point(nx + ny, ny + nz, nx + nz));
        return pts;
    }

    // C(n+3, 3), the number of even lattice points of the scaled tetrahedron
    public static long TetrahedronCount(int n) => (long)(n + 1) * (n + 2) * (n + 3) / 6;

    // the strings at the centre (n/2, n/2, n/2): n_X = n_Y = n_Z = n/4, so n!/((n/4)!)^4 when 4 | n, else none
    public static BigInteger CentreCount(int n)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "n >= 1");
        if (n % 4 != 0) return BigInteger.Zero;
        BigInteger f(int m) { BigInteger r = 1; for (int i = 2; i <= m; i++) r *= i; return r; }
        var q = f(n / 4);
        return f(n) / (q * q * q * q);
    }

    // The height under Pauli dephasing with integer rates per site and letter: rates[l] = (gamma_Z, gamma_X,
    // gamma_Y) on site l. A string decays at this rate: L(sigma) = -Height(sigma) sigma.
    public static long Height(PauliString s, (long Z, long X, long Y)[] rates)
    {
        if (rates.Length < 1 || rates.Length > PauliString.MaxSites)
            throw new ArgumentOutOfRangeException(nameof(rates), rates.Length, $"one rate triple per site, 1 to {PauliString.MaxSites} sites");
        ulong outside = ~((1UL << rates.Length) - 1);
        if (((s.X | s.Z) & outside) != 0)
            throw new ArgumentException("the string has letters beyond the sites the rates cover", nameof(s));
        long h = 0;
        for (int l = 0; l < rates.Length; l++)
        {
            bool x = ((s.X >> l) & 1) == 1, z = ((s.Z >> l) & 1) == 1;
            if (x) h += rates[l].Z;            // X or Y anticommutes with Z
            if (z) h += rates[l].X;            // Y or Z anticommutes with X
            if (x ^ z) h += rates[l].Y;        // X or Z anticommutes with Y
        }
        return 2 * h;
    }

    // the same rates on every site: 2(gamma_Z k_Z + gamma_X k_X + gamma_Y k_Y)
    public static long Height(Point p, long gammaZ, long gammaX, long gammaY)
        => 2 * (gammaZ * p.KZ + gammaX * p.KX + gammaY * p.KY);

    // The half-turn about P's axis on one site: flip the two coordinates other than k_P, k -> 1 - k, and read the
    // letter at the image point. It sends I to P and the other two letters to each other.
    public static char HalfTurn(char letter, char axis)
    {
        var p = Coordinates(letter.ToString());
        var q = axis switch
        {
            'Z' => new Point(p.KZ, 1 - p.KX, 1 - p.KY),
            'X' => new Point(1 - p.KZ, p.KX, 1 - p.KY),
            'Y' => new Point(1 - p.KZ, 1 - p.KX, p.KY),
            _ => throw new ArgumentException($"'{axis}' is not an axis letter", nameof(axis)),
        };
        return "IXYZ".Single(c => Coordinates(c.ToString()) == q);
    }

    // With one-site rates (gamma_Z, gamma_X, gamma_Y) the half-turn about P's axis is a pairing of the four letters'
    // costs (every letter and its image sum to the same) exactly when gamma_P = 0. Returns the axes whose half-turn
    // pairs.
    public static IReadOnlyList<char> PairingHalfTurns(long gammaZ, long gammaX, long gammaY)
    {
        long cost(char letter) => Height(Coordinates(letter.ToString()), gammaZ, gammaX, gammaY);
        var axes = new List<char>();
        foreach (char p in "ZXY")
        {
            var sums = "IXYZ".Select(c => cost(c) + cost(HalfTurn(c, p))).Distinct().Count();
            if (sums == 1) axes.Add(p);
        }
        return axes;
    }

    // ---- the steps (§8) ----

    // The move a placed term makes on a string: null when they commute (the commutator is zero), else the string
    // the commutator is proportional to and the step (dk_Z, dk_X, dk_Y).
    public static (PauliString Target, Point Step)? Step(PauliString term, PauliString s)
    {
        if (PauliString.Commute(term, s)) return null;
        var (target, _) = PauliString.Multiply(term, s);
        var a = Coordinates(s);
        var b = Coordinates(target);
        return (target, new Point(b.KZ - a.KZ, b.KX - a.KX, b.KY - a.KY));
    }

    // [term, sigma] = 2 i^PhasePower Target when they anticommute (term sigma = i^PhasePower Target, sigma term =
    // -term sigma), null when they commute. Enough to rebuild every entry of ad_H in the string basis.
    public static (PauliString Target, int PhasePower)? CommutatorCoefficient(PauliString term, PauliString s)
    {
        if (PauliString.Commute(term, s)) return null;
        return PauliString.Multiply(term, s);
    }

    // every step the placed terms make on any string of n sites
    public static IReadOnlySet<Point> Moves(IEnumerable<PauliString> terms, int n)
    {
        var ts = terms.ToList();
        var steps = new HashSet<Point>();
        foreach (var s in Strings(n))
            foreach (var t in ts)
                if (Step(t, s) is { } r) steps.Add(r.Step);
        return steps;
    }

    // the pairs of distinct lattice points the placed terms join (a move that lands on the point it started from,
    // possible for terms on three or more sites, gives no edge), each pair ordered
    public static IReadOnlySet<(Point, Point)> Edges(IEnumerable<PauliString> terms, int n)
    {
        var ts = terms.ToList();
        var edges = new HashSet<(Point, Point)>();
        foreach (var s in Strings(n))
            foreach (var t in ts)
                if (Step(t, s) is { } r)
                {
                    var a = Coordinates(s);
                    var b = Coordinates(r.Target);
                    if (a == b) continue;
                    edges.Add(Order(a, b));
                }
        return edges;
    }

    static (Point, Point) Order(Point a, Point b)
        => (a.KZ, a.KX, a.KY).CompareTo((b.KZ, b.KX, b.KY)) <= 0 ? (a, b) : (b, a);

    // the placed terms of a template on every window of an open chain of n sites: a k-letter template sits on sites
    // i..i+k-1, a one-letter template on every site
    public static IReadOnlyList<PauliString> ChainTerms(IEnumerable<string> templates, int n)
    {
        var placed = new List<PauliString>();
        foreach (var t in templates)
            for (int i = 0; i + t.Length <= n; i++)
                placed.Add(PauliString.Parse(new string('I', i) + t + new string('I', n - t.Length - i)));
        return placed;
    }

    // ---- the arrows (§9, gates M8, M9) ----

    // A Gaussian integer, the exact arithmetic of the transfer matrices.
    public readonly record struct Gauss(long Re, long Im)
    {
        public static readonly Gauss Zero = new(0, 0);
        public static readonly Gauss One = new(1, 0);
        public static readonly Gauss I = new(0, 1);
        public static Gauss operator +(Gauss a, Gauss b) => new(a.Re + b.Re, a.Im + b.Im);
        public static Gauss operator -(Gauss a, Gauss b) => new(a.Re - b.Re, a.Im - b.Im);
        public static Gauss operator -(Gauss a) => new(-a.Re, -a.Im);
        public static Gauss operator *(Gauss a, Gauss b) => new(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);
        public static Gauss operator *(long k, Gauss a) => new(k * a.Re, k * a.Im);
        public Gauss Conj => new(Re, -Im);
        public bool IsZero => Re == 0 && Im == 0;
        public override string ToString() => Im == 0 ? $"{Re}" : Re == 0 ? $"{Im}i" : $"{Re}{(Im < 0 ? "-" : "+")}{Math.Abs(Im)}i";
    }

    static readonly string Letters = "IXYZ";

    // the four letters as 2x2 Gaussian matrices
    public static Gauss[,] LetterMatrix(char letter) => letter switch
    {
        'I' => new[,] { { Gauss.One, Gauss.Zero }, { Gauss.Zero, Gauss.One } },
        'X' => new[,] { { Gauss.Zero, Gauss.One }, { Gauss.One, Gauss.Zero } },
        'Y' => new[,] { { Gauss.Zero, -Gauss.I }, { Gauss.I, Gauss.Zero } },
        'Z' => new[,] { { Gauss.One, Gauss.Zero }, { Gauss.Zero, -Gauss.One } },
        _ => throw new ArgumentException($"'{letter}' is not a Pauli letter", nameof(letter)),
    };

    // sigma^- = |0><1| and sigma^+ = |1><0|, Z|0> = |0>
    public static Gauss[,] Lowering => new[,] { { Gauss.Zero, Gauss.One }, { Gauss.Zero, Gauss.Zero } };
    public static Gauss[,] Raising => new[,] { { Gauss.Zero, Gauss.Zero }, { Gauss.One, Gauss.Zero } };

    static Gauss[,] Mul2(Gauss[,] a, Gauss[,] b)
    {
        var c = new Gauss[2, 2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                c[i, j] = a[i, 0] * b[0, j] + a[i, 1] * b[1, j];
        return c;
    }

    static Gauss[,] Dagger2(Gauss[,] a) => new[,] { { a[0, 0].Conj, a[1, 0].Conj }, { a[0, 1].Conj, a[1, 1].Conj } };

    // 4 T[a, b] = 2 Tr(P_a D(P_b)) with D(rho) = r (J rho J^dagger - {J^dagger J, rho}/2), letters in the order
    // I, X, Y, Z; integral for an integral J and rate.
    public static Gauss[,] TransferTimesFour(Gauss[,] jump, long rate)
    {
        var jd = Dagger2(jump);
        var jdj = Mul2(jd, jump);
        var t = new Gauss[4, 4];
        for (int b = 0; b < 4; b++)
        {
            var pb = LetterMatrix(Letters[b]);
            var sandwich = Mul2(jump, Mul2(pb, jd));
            var left = Mul2(jdj, pb);
            var right = Mul2(pb, jdj);
            var d2 = new Gauss[2, 2];                         // 2 D(P_b) / rate, integral
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    d2[i, j] = 2 * sandwich[i, j] - left[i, j] - right[i, j];
            for (int a = 0; a < 4; a++)
            {
                var pa = LetterMatrix(Letters[a]);
                var prod = Mul2(pa, d2);
                t[a, b] = rate * (prod[0, 0] + prod[1, 1]);   // Tr(P_a 2D(P_b)) = 4 T[a, b]
            }
        }
        return t;
    }

    // several jumps on one site add
    public static Gauss[,] TransferTimesFour(IEnumerable<(Gauss[,] Jump, long Rate)> jumps)
    {
        var t = new Gauss[4, 4];
        foreach (var (j, r) in jumps)
        {
            var tj = TransferTimesFour(j, r);
            for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) t[a, b] = t[a, b] + tj[a, b];
        }
        return t;
    }

    // The height of each letter, times four: minus the transfer's diagonal (I, X, Y, Z).
    public static Gauss[] HeightsTimesFour(Gauss[,] transferTimesFour)
        => Enumerable.Range(0, 4).Select(a => -transferTimesFour[a, a]).ToArray();

    // The arrows: every off-diagonal entry, a move from letter From to letter To with its coefficient times four.
    public static IReadOnlyList<(char From, char To, Gauss CoefficientTimesFour)> Arrows(Gauss[,] transferTimesFour)
    {
        var arrows = new List<(char, char, Gauss)>();
        for (int b = 0; b < 4; b++)
            for (int a = 0; a < 4; a++)
                if (a != b && !transferTimesFour[a, b].IsZero)
                    arrows.Add((Letters[b], Letters[a], transferTimesFour[a, b]));
        return arrows;
    }

    // R, the right multiplication by X^N, letter by letter (I <-> X, Y <-> Z), up to phase
    public static PauliString RightByXN(PauliString s, int n)
    {
        if (n < 1 || n > PauliString.MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"1 to {PauliString.MaxSites} sites");
        ulong all = (1UL << n) - 1;
        if (((s.X | s.Z) & ~all) != 0)
            throw new ArgumentException($"the string has letters beyond site {n - 1}", nameof(s));
        return new PauliString(s.X ^ all, s.Z);
    }

    // R's action on the cube read as a map of points: (k_Z, k_X, k_Y) -> (N - k_Z, k_X, N - k_Y)
    public static Point RAsPoint(Point p, int n)
    {
        if (n < 1 || p.KZ > n || p.KX > n || p.KY > n)
            throw new ArgumentOutOfRangeException(nameof(n), n, "a point of the cube on n >= 1 sites");
        return new(n - p.KZ, p.KX, n - p.KY);
    }

    // the parity cube (bit_a, bit_b, y_par) read off a point: k_Z mod 2, k_X mod 2, half a sum mod 2
    public static (int A, int B, int Y) Parities(Point p) => (p.KZ % 2, p.KX % 2, ((p.KZ + p.KX - p.KY) / 2) % 2);
}
