namespace MirrorWorld;

using static MirrorWorld.LetterCube;

// The theodolite (built 2026-10-20 beside experiments/THE_CUBE_SITE_BY_SITE.md and its Python twin
// simulations/framework/diagnostics/cube_letters.py; the name is Tom's). Why the name: a surveyor sets a theodolite up
// at one station after another along a traverse and reads two angles at each. This object reads, site by site, a
// mirror of the generator that is a product over the sites: at every site a turn about Z and a thermal weight (an
// imaginary turn), and from site to site the shape of the traverse, straight (one map on every site), a fixed
// deflection (the map turned by the same angle from station to station, the way a curve is staked out), or repeating.
//
// What it reads. A sandwich E(rho) = A·rho^(T)·B with A and B products over the sites and every site factor monomial:
// diagonal, diag(1, x), or anti-diagonal, [[0, x], [1, 0]]. On the letters such a factor is a pure turn: it multiplies
// sigma^+ by a power of its turn x or swaps sigma^+ and sigma^-, and keeps or flips Z. So every condition the mirror
// must meet is binomial, a product of turns equal to a constant read off the rates and the couplings (Binomials):
//   per site      E_l·L_l° = T(L_l°)·E_l, with T(M) = -M for the Pi-reflection -L - c, -M^dagger for the
//                 reflection -L^dagger - c and M for a copy (E·L = L·E, a symmetry), from the site's heights,
//                 arrows and field;
//   per support   A_S·h_S^(T)·A_S^-1 = s·h_S and B_S^-1·h_S^(T)·B_S = s·h_S, the left and the right factor apart,
//                 s = -1 for (Pi, plain), (reflection, transposed) and (copy, transposed) and +1 otherwise: the
//                 Hamiltonian part of the sector split, exactly, since a bond's terms are traceless on every site and
//                 no constant enters. A copy's bonds read like a reflection's (-L_H^dagger = L_H for a Hermitian H);
//                 only its sites differ.
// The plain all-diagonal copies are left out: they form a group holding the identity, which the solver would return,
// so a copy's silence means no copy of the other kinds, not no symmetry (the turns about Z and the moves by Z^N, the
// turn rho -> Z^N·rho·Z^N and the shift rho -> rho·Z^N, are in the group left out). The first copy read so is U =
// Sigma·rho^T·Sigma on the XY chain, Sigma the sign on every second site, the involution
// docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md splits the block (1, 1) by.
// A reading is the reduced binomial system with the kinds of the factors; it never holds a number. On a chain built
// from one repeated cell the same equations with a step per site, x_l = x_0·r^l and y_l = y_0·s^l, read every length
// at once: two positions of a window fix the step, so a traverse holds on every open chain. Uniform: one map on every
// site. Turned: the map turned about Z by an angle from site to site, z with |z| = 1, pinned by z^k = Num/Den with
// |Num| = |Den|. Periodic: kinds and turns repeating with a period p. Stepped: the turns step, but no turn of the map
// carries one site onto the next.
//
// Relabelled (ReadRelabelled): the sandwich after a permutation P_pi of the sites, E·P_pi, for every pi that carries
// the setting's supports onto themselves (on a chain the reversal, on a ring its rotations and reflections). P_pi
// changes no string's letter counts and so no height; it moves each site's table to its image, and a reflection turns
// the letters of a bond's string around (XY on (l, l+1) becomes YX), so the same binomial reading, asked between the
// relabelled setting and the setting itself, decides it. It reads, under loss alone, a bond with ZZ and s, a, w all
// nonzero at rates the reversal keeps: the reversal times a uniform turn, then S_XY
// (experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md).
//
// Sectors (ReadSectors): an unwatched site where H uses exactly one letter (so H conserves it) is a control: turned onto Z, its bit is a number H never changes, so H = sum over the controls' bits sigma of
// |sigma><sigma| ⊗ H_sigma, and a bond to a control is a field on its neighbour whose sign the bit sets
// (docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md Theorem 4). The sandwich may then depend on the bits: A =
// sum_sigma |pi_A(sigma)><sigma| ⊗ A_sigma and B likewise, every A_sigma, B_tau a product of monomial factors on the
// other sites with their own turns x_l.sigma, y_l.tau, each site's kind the same in every sector, the controls kept
// or flipped. No dissipator touches a control,
// so L is a sum of blocks (sigma, tau), each a generator on the other sites with H_sigma on the ket and H_tau on the
// bra, and every block is read site by site as above: a site's table with one field on the ket and another on the bra,
// a support's left equation with the ket's Hamiltonian and its right one with the bra's. The blocks' constants (the
// terms on the controls alone) must agree with the relation: e_sigma - e_tau the same on the target block for a
// reflection or a copy, its negative for the Pi-reflection. It reads -3XY with fields -2X under Z dephasing at (0, 3),
// where no product reads: in sector s site 1 sees the field (-2, -3s), perpendicular to its watching, and its turn
// follows the sector, the conditioned turn exp(-i·theta·X_0·Z_1), a frame in the sense of
// experiments/THE_PALINDROME_AS_A_COLOURING.md (it commutes with every jump), read without its angle. The frames of
// the other sites are asked blind, 3^k of them for k sites that are not controls.
//
// Parent: the letter cube, whose letters, heights, arrows and one-site moves it reads with; the sextant's sibling,
// which asks the cube's moves on the whole generator, uniform on every site. Its own are the letters (the binomial
// equations), the reading and every N.
//
// Deliberately outside: the numbers (the Python twin takes the roots), site maps that are no such sandwich (the Python
// twin's cube search and machine), bits a sandwich depends on other than a conserved unwatched letter's and kinds that
// change with a control's bit (ReadSectors),
// frames other than one Clifford axis per site carried onto Z (ReadRotated, ReadSectors), weights
// along X or Y, relabellings and turned frames in EveryN, a relabelling and a turned frame together, and a silence's
// meaning: when the letters are silent, nothing of this form may exist or the element may not be such a sandwich, and
// the theodolite does not decide which.
//
// Words, fenced at the door: a "turn" here is a site factor's entry x, a turn about Z (a real angle) or a thermal weight
// (an imaginary one), not the letter cube's turn by X^N and not GammaFold's per-site turn gamma -> -gamma; a "kind" is
// diagonal or anti-diagonal; a "reading" is the theodolite's answer for one setting, not the sextant's Reading of a
// move; a "traverse" is its answer for every N; a "setting" is the system read, as in the sextant.
public sealed class Theodolite : GameObject
{
    public Theodolite(LetterCube cube) : base(cube) { }

    // left: what the theodolite itself produces.
    public override IReadOnlyList<string> Own => new[] { "letters", "reading", "every N" };

    public enum Relation { Pi, Reflection, Copy }

    // the branches asked: null both (plain first), or one of them
    static bool[] Branches(bool? transposed) => transposed is { } t ? new[] { t } : new[] { false, true };

    // the plain all-diagonal copies, the group that holds the identity (left out, see the head comment)
    static bool Gauge(Relation rel, bool transposed, IEnumerable<(Kind A, Kind B)> kinds)
        => rel == Relation.Copy && !transposed && kinds.All(k => k == (Kind.Diagonal, Kind.Diagonal));
    public enum Kind { Diagonal, AntiDiagonal }
    public enum Shape { Uniform, Turned, Periodic, Stepped }

    static readonly (Kind, Kind)[] KindPairs =
    {
        (Kind.Diagonal, Kind.Diagonal), (Kind.Diagonal, Kind.AntiDiagonal),
        (Kind.AntiDiagonal, Kind.Diagonal), (Kind.AntiDiagonal, Kind.AntiDiagonal),
    };

    // A setting: n sites, H as Pauli strings with integer coefficients (a one-letter string is a field), one-site jumps
    // with integer rates. Rational rates and couplings are scaled to integers first: every equation the theodolite reads
    // is homogeneous in one site's table or one support's terms, so a common factor changes none of them.
    public sealed class Setting
    {
        public int N { get; }
        public IReadOnlyList<(PauliString Term, long Coefficient)> H { get; }
        public IReadOnlyList<(int Site, Gauss[,] Jump, long Rate)> Jumps { get; }

        public Setting(int n, IEnumerable<(PauliString, long)> h, IEnumerable<(int, Gauss[,], long)> jumps)
        {
            if (n < 1 || n > 12) throw new ArgumentOutOfRangeException(nameof(n), n, "1 to 12 sites");
            N = n;
            H = h.Where(t => !t.Item1.IsIdentity && t.Item2 != 0).ToList();
            Jumps = jumps.ToList();
            ulong inside = (1UL << n) - 1;
            foreach (var (term, _) in H)
                if (((term.X | term.Z) & ~inside) != 0) throw new ArgumentException($"a term of H reaches past site {n - 1}", nameof(h));
            foreach (var (site, _, _) in Jumps)
                if (site < 0 || site >= n) throw new ArgumentOutOfRangeException(nameof(jumps), site, $"sites run from 0 to {n - 1}");
        }

        // H from strings, site 0 first ("XIY" puts X on site 0 and Y on site 2)
        public static Setting Of(int n, IEnumerable<(string Letters, long Coefficient)> h, IEnumerable<(int, Gauss[,], long)> jumps)
            => new(n, h.Select(t => (PauliString.Parse(t.Letters.PadRight(n, 'I')), t.Coefficient)), jumps);

        internal int[] Support(PauliString t) => Enumerable.Range(0, N).Where(l => t.Letter(l) != 'I').ToArray();

        // The same system seen in a turned frame, times four: site l's axis axes[l] carried onto Z by U_l, so the
        // generator is 4·Ad_U·L·Ad_U^-1. U_X = (X + Z)/√2 (X <-> Z, Y -> -Y), U_Y = (1 - iX)/√2 (Y -> Z, Z -> -Y),
        // U_Z = 1. A jump J on a turned site becomes M·J·M^dagger, M = √2·U integral, so its rate stays; every other
        // rate and every coefficient of H takes the factor four.
        public Setting Rotated(string axes)
        {
            if (axes.Length != N || axes.Any(a => a is not ('X' or 'Y' or 'Z')))
                throw new ArgumentException("one axis X, Y or Z per site", nameof(axes));
            var h = H.Select(t =>
            {
                var letters = new char[N];
                long sign = 1;
                for (int l = 0; l < N; l++)
                {
                    var (letter, sg) = Turn(axes[l], t.Term.Letter(l));
                    letters[l] = letter;
                    sign *= sg;
                }
                return (PauliString.Parse(new string(letters)), 4 * sign * t.Coefficient);
            });
            var jumps = Jumps.Select(j => axes[j.Site] == 'Z'
                ? (j.Site, j.Jump, 4 * j.Rate)
                : (j.Site, Sandwich(TurnMatrix(axes[j.Site]), j.Jump), j.Rate));
            return new Setting(N, h, jumps);
        }

        // the letter P becomes sign·P' under U·P·U^dagger
        static (char, long) Turn(char axis, char p) => (axis, p) switch
        {
            ('X', 'X') => ('Z', 1), ('X', 'Z') => ('X', 1), ('X', 'Y') => ('Y', -1),
            ('Y', 'Y') => ('Z', 1), ('Y', 'Z') => ('Y', -1), ('Y', 'X') => ('X', 1),
            _ => (p, 1),
        };

        static Gauss[,] TurnMatrix(char axis) => axis == 'X'
            ? new[,] { { Gauss.One, Gauss.One }, { Gauss.One, -Gauss.One } }                    // X + Z
            : new[,] { { Gauss.One, new Gauss(0, -1) }, { new Gauss(0, -1), Gauss.One } };     // 1 - iX

        static Gauss[,] Sandwich(Gauss[,] m, Gauss[,] j)
        {
            var o = new Gauss[2, 2];
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                {
                    var v = Gauss.Zero;
                    for (int k = 0; k < 2; k++)
                        for (int q = 0; q < 2; q++)
                            v = v + m[a, k] * j[k, q] * m[b, q].Conj;                            // M·J·M^dagger
                    o[a, b] = v;
                }
            return o;
        }

        // the same system with site l carried to site pi[l]: P_pi·L·P_pi^-1 is this setting's generator
        public Setting Relabelled(IReadOnlyList<int> pi)
        {
            if (pi.Count != N || pi.Distinct().Count() != N || pi.Any(v => v < 0 || v >= N))
                throw new ArgumentException("a relabelling is a permutation of the sites", nameof(pi));
            var h = H.Select(t =>
            {
                var letters = new char[N];
                for (int l = 0; l < N; l++) letters[pi[l]] = t.Term.Letter(l);
                return (PauliString.Parse(new string(letters)), t.Coefficient);
            });
            return new Setting(N, h, Jumps.Select(j => (pi[j.Site], j.Jump, j.Rate)));
        }
    }

    // One site of a chain repeated: its jumps, its fields (letter, coefficient) and its windows (letters, offsets from
    // the window's first site, coefficient), e.g. ("XZX", {0, 1, 2}, 3).
    public sealed class Cell
    {
        public IReadOnlyList<(Gauss[,] Jump, long Rate)> Jumps { get; }
        public IReadOnlyList<(char Letter, long Coefficient)> Fields { get; }
        public IReadOnlyList<(string Letters, int[] Offsets, long Coefficient)> Windows { get; }

        public Cell(IEnumerable<(Gauss[,], long)> jumps, IEnumerable<(char, long)>? fields = null,
            IEnumerable<(string, int[], long)>? windows = null)
        {
            Jumps = jumps.ToList();
            Fields = (fields ?? Enumerable.Empty<(char, long)>()).ToList();
            Windows = (windows ?? Enumerable.Empty<(string, int[], long)>()).ToList();
        }

        // identity letters dropped, offsets shifted to start at 0, one-letter windows moved to the fields
        public Cell Normalized()
        {
            var fields = new List<(char, long)>(Fields);
            var windows = new List<(string, int[], long)>();
            foreach (var (letters, off, c) in Windows)
            {
                var kept = letters.Select((ch, i) => (ch, o: off[i])).Where(t => t.ch != 'I').ToList();
                if (kept.Count == 0) continue;
                if (kept.Count == 1) { fields.Add((kept[0].ch, c)); continue; }
                int m = kept.Min(t => t.o);
                windows.Add((new string(kept.Select(t => t.ch).ToArray()), kept.Select(t => t.o - m).ToArray(), c));
            }
            return new Cell(Jumps, fields, windows);
        }

        public int Span => Windows.Count == 0 ? 1 : Windows.Max(w => w.Offsets.Max()) + 1;

        public Setting Chain(int length)
        {
            var h = new List<(PauliString, long)>();
            for (int l = 0; l < length; l++)
                foreach (var (letter, c) in Fields)
                    h.Add((OnSites(length, new[] { l }, new[] { letter }), c));
            foreach (var (letters, off, c) in Windows)
                for (int start = 0; start + off.Max() < length; start++)
                    h.Add((OnSites(length, off.Select(o => start + o).ToArray(), letters.ToCharArray()), c));
            var jumps = Enumerable.Range(0, length).SelectMany(l => Jumps.Select(j => (l, j.Jump, j.Rate)));
            return new Setting(length, h, jumps);
        }

        static PauliString OnSites(int n, int[] sites, char[] letters)
        {
            var s = Enumerable.Repeat('I', n).ToArray();
            for (int i = 0; i < sites.Length; i++) s[sites[i]] = letters[i];
            return PauliString.Parse(new string(s));
        }
    }

    // ---------------------------------------------------------------- the letters: one site, one support

    sealed class SilentException : Exception
    {
        public SilentException(string where) : base(where) { }
    }

    static Gauss[,] Id2 => LetterMatrix('I');

    static Gauss[,] Mul2(Gauss[,] a, Gauss[,] b)
    {
        var c = new Gauss[2, 2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                c[i, j] = a[i, 0] * b[0, j] + a[i, 1] * b[1, j];
        return c;
    }

    static Gauss[,] Dagger(Gauss[,] a)
    {
        int d = a.GetLength(0);
        var t = new Gauss[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                t[i, j] = a[j, i].Conj;
        return t;
    }

    static Gauss[,] Transpose(Gauss[,] a)
    {
        int d = a.GetLength(0);
        var t = new Gauss[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                t[i, j] = a[j, i];
        return t;
    }

    // rho -> a·rho·b on one site, in cells (|0><0|, |0><1|, |1><0|, |1><1|, row-major): kron(a, b^T)
    static Gauss[,] Lr(Gauss[,] a, Gauss[,] b)
    {
        var m = new Gauss[4, 4];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++)
                        m[2 * i + j, 2 * k + l] = a[i, k] * b[l, j];
        return m;
    }

    // one site's generator times two (integral): the field -i[h, .] and every jump rate·D[J]; with braFields the bra
    // carries its own field, -i(h·rho - rho·h')
    static Gauss[,] SiteTimesTwo(IEnumerable<(char Letter, long Coefficient)> fields, IEnumerable<(Gauss[,] Jump, long Rate)> jumps,
        IEnumerable<(char Letter, long Coefficient)>? braFields = null)
    {
        var t = new Gauss[4, 4];
        Gauss[,] Field(IEnumerable<(char Letter, long Coefficient)> fs)
        {
            var h = new Gauss[2, 2];
            foreach (var (letter, c) in fs)
            {
                var p = LetterMatrix(letter);
                for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) h[i, j] = h[i, j] + c * p[i, j];
            }
            return h;
        }
        var hl = Lr(Field(fields), Id2);
        var hr = Lr(Id2, Field(braFields ?? fields));
        var minusTwoI = new Gauss(0, -2);
        for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) t[r, c] = t[r, c] + minusTwoI * (hl[r, c] - hr[r, c]);
        foreach (var (jump, rate) in jumps)
        {
            var jd = Dagger(jump);
            var jdj = Mul2(jd, jump);
            var sandwich = Lr(jump, jd);
            var left = Lr(jdj, Id2);
            var right = Lr(Id2, jdj);
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    t[r, c] = t[r, c] + rate * (2 * sandwich[r, c] - left[r, c] - right[r, c]);
        }
        return t;
    }

    // the traceless part times eight: 8·(L - tr L/4) = 4·(2L) - tr(2L)·1, integral
    static Gauss[,] TracelessTimesEight(Gauss[,] twice)
    {
        var tr = twice[0, 0] + twice[1, 1] + twice[2, 2] + twice[3, 3];
        var m = new Gauss[4, 4];
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                m[r, c] = 4 * twice[r, c] - (r == c ? tr : Gauss.Zero);
        return m;
    }

    static Gauss[,] Target(Gauss[,] m, Relation rel)
    {
        if (rel == Relation.Copy) return m;
        var t = rel == Relation.Pi ? m : Dagger(m);
        int d = t.GetLength(0);
        var o = new Gauss[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) o[i, j] = -t[i, j];
        return o;
    }

    // a factor's nonzero entries (row, column, exponent of its turn): diag(1, x) or [[0, x], [1, 0]]
    static (int Row, int Col, int Exp)[] FactorEntries(Kind k) => k == Kind.Diagonal
        ? new[] { (0, 0, 0), (1, 1, 1) }
        : new[] { (0, 1, 1), (1, 0, 0) };

    // rho -> A·rho^(T)·B on one site, in cells: per entry whether it is nonzero and its exponents of (x, y)
    static (bool On, int Ex, int Ey)[,] SiteMap(Kind ka, Kind kb, bool transposed)
    {
        var e = new (bool, int, int)[4, 4];
        foreach (var (i, k, ex) in FactorEntries(ka))
            foreach (var (l, j, ey) in FactorEntries(kb))
                e[2 * i + j, 2 * k + l] = (true, ex, ey);
        if (!transposed) return e;
        int[] perm = { 0, 2, 1, 3 };
        var t = new (bool, int, int)[4, 4];
        for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) t[r, c] = e[r, perm[c]];
        return t;
    }

    // The site equation E·L° = T(L°)·E for every piece, as binomials in (x, y); null when an entry is a single
    // monomial (no nonzero solution); throws past two monomials (silent).
    static List<Binomial>? SiteEquations(IReadOnlyList<(Gauss[,] From, Gauss[,] To)> pieces, Relation rel, Kind ka, Kind kb,
        bool transposed, string x, string y, string origin)
    {
        var e = SiteMap(ka, kb, transposed);
        var outs = new List<Binomial>();
        foreach (var (l8, to) in pieces)
        {
            var tg = Target(to, rel);
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                {
                    var terms = new Dictionary<(int, int), Gauss>();
                    for (int k = 0; k < 4; k++)
                    {
                        if (e[r, k].On && !l8[k, c].IsZero) AddTerm(terms, (e[r, k].Ex, e[r, k].Ey), l8[k, c]);
                        if (e[k, c].On && !tg[r, k].IsZero) AddTerm(terms, (e[k, c].Ex, e[k, c].Ey), -tg[r, k]);
                    }
                    var nz = terms.Where(kv => !kv.Value.IsZero).ToList();
                    if (nz.Count == 0) continue;
                    if (nz.Count == 1) return null;
                    if (nz.Count > 2) throw new SilentException($"{origin}: {nz.Count} monomials");
                    var ((ex1, ey1), c1) = (nz[0].Key, nz[0].Value);
                    var ((ex2, ey2), c2) = (nz[1].Key, nz[1].Value);
                    var exps = new Dictionary<string, int>();
                    if (ex1 != ex2) exps[x] = ex1 - ex2;
                    if (ey1 != ey2) exps[y] = ey1 - ey2;
                    outs.Add(new Binomial(exps, GaussBig.Of(-c2), GaussBig.Of(c1), origin));
                }
        }
        return outs;
    }

    static void AddTerm(Dictionary<(int, int), Gauss> terms, (int, int) m, Gauss v)
        => terms[m] = terms.TryGetValue(m, out var old) ? old + v : v;

    // the Hamiltonian of the terms whose support is exactly these sites, on them in this order (first site the most
    // significant bit)
    static Gauss[,] HamiltonianOn(Setting s, int[] support)
    {
        int k = support.Length, d = 1 << k;
        var h = new Gauss[d, d];
        foreach (var (term, c) in s.H)
        {
            if (!s.Support(term).SequenceEqual(support)) continue;
            var m = new Gauss[1, 1] { { Gauss.One } };
            foreach (int site in support) m = Kron(m, LetterMatrix(term.Letter(site)));
            for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) h[i, j] = h[i, j] + c * m[i, j];
        }
        return h;
    }

    static Gauss[,] Kron(Gauss[,] a, Gauss[,] b)
    {
        int ra = a.GetLength(0), rb = b.GetLength(0);
        var m = new Gauss[ra * rb, ra * rb];
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < ra; j++)
                for (int k = 0; k < rb; k++)
                    for (int l = 0; l < rb; l++)
                        m[i * rb + k, j * rb + l] = a[i, j] * b[k, l];
        return m;
    }

    // The support equation, left (A_S·h'^(T)·A_S^-1 = s·h) or right (B_S^-1·h'^(T)·B_S = s·h), as binomials in the
    // support's turns, h' the source's terms and h the target's (the same setting unless the sites are relabelled).
    // A monomial factor sends basis state i to pi(i) with the monomial m_i, so every entry of the
    // conjugated h is one monomial times one entry of h: (A h A^-1)_ij = m_i h_pi(i)pi(j) / m_j and
    // (B^-1 h B)_ij = h_pi(i)pi(j) m_pi(j) / m_pi(i). Null when an entry has no nonzero solution.
    static List<Binomial>? SupportEquations(Gauss[,] hFrom, Gauss[,] h, Relation rel, bool transposed, IReadOnlyList<Kind> kinds,
        IReadOnlyList<string> turns, bool left, string origin)
    {
        int k = kinds.Count, d = 1 << k;
        int sign = (rel == Relation.Pi && !transposed) || (rel != Relation.Pi && transposed) ? -1 : 1;
        var ht = transposed ? Transpose(hFrom) : hFrom;
        int Bit(int i, int p) => (i >> (k - 1 - p)) & 1;
        int Image(int i)
        {
            int o = i;
            for (int p = 0; p < k; p++)
                if (kinds[p] == Kind.AntiDiagonal) o ^= 1 << (k - 1 - p);
            return o;
        }
        int[] Exponent(int i)
        {
            var e = new int[k];
            for (int p = 0; p < k; p++)
                e[p] = (kinds[p] == Kind.Diagonal) == (Bit(i, p) == 1) ? 1 : 0;
            return e;
        }
        var outs = new List<Binomial>();
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                var conj = ht[Image(i), Image(j)];
                var rhs = sign * h[i, j];
                if (conj.IsZero && rhs.IsZero) continue;
                if (conj.IsZero || rhs.IsZero) return null;
                var top = left ? Exponent(i) : Exponent(Image(j));
                var bottom = left ? Exponent(j) : Exponent(Image(i));
                var exps = new Dictionary<string, int>();
                for (int p = 0; p < k; p++)
                    if (top[p] != bottom[p]) exps[turns[p]] = top[p] - bottom[p];
                if (exps.Count == 0)
                {
                    if (conj == rhs) continue;
                    return null;
                }
                outs.Add(new Binomial(exps, GaussBig.Of(rhs), GaussBig.Of(conj), origin));
            }
        return outs;
    }

    static IReadOnlyList<Gauss[,]> SitePieces(Setting s, int l, bool rateFree)
    {
        var fields = s.H.Where(t => s.Support(t.Term).SequenceEqual(new[] { l }))
            .Select(t => (t.Term.Letter(l), t.Coefficient)).ToList();
        var jumps = s.Jumps.Where(j => j.Site == l).Select(j => (j.Jump, j.Rate)).ToList();
        return Pieces(fields, jumps, rateFree);
    }

    // the whole site table, or with rateFree every channel at rate 1 and the fields apart (the element must hold at
    // every rate of every channel)
    static IReadOnlyList<Gauss[,]> Pieces(IReadOnlyList<(char, long)> fields, IReadOnlyList<(Gauss[,] Jump, long Rate)> jumps, bool rateFree)
    {
        if (!rateFree) return new[] { TracelessTimesEight(SiteTimesTwo(fields, jumps)) };
        var pieces = jumps.Select(j => TracelessTimesEight(SiteTimesTwo(Array.Empty<(char, long)>(), new[] { (j.Jump, 1L) }))).ToList();
        if (fields.Count > 0) pieces.Add(TracelessTimesEight(SiteTimesTwo(fields, Array.Empty<(Gauss[,], long)>())));
        return pieces;
    }

    static List<int[]> Supports(Setting s)
        => s.H.Select(t => s.Support(t.Term)).Where(sup => sup.Length >= 2)
            .GroupBy(sup => string.Join(",", sup)).Select(g => g.First()).OrderBy(sup => string.Join(",", sup), StringComparer.Ordinal)
            .ToList();

    // ---------------------------------------------------------------- reading a setting

    // Relabelling: null for the site-local sandwich; otherwise pi, site l carried to site pi[l] before the sandwich.
    public sealed record Reading(bool Transposed, IReadOnlyList<(Kind A, Kind B)> Kinds, Binomials.Reduced System,
        IReadOnlyList<Binomial> Equations, IReadOnlyList<int>? Relabelling = null, string? Axes = null,
        IReadOnlyList<int>? Controls = null)
    {
        public IEnumerable<string> Lines()
        {
            yield return $"rho -> {(Transposed ? "A·rho^T·B" : "A·rho·B")} with " +
                         (Controls is { } cs ? $"A = sum over the bits sigma of the controls [{string.Join(", ", cs)}] of their kept or flipped bits ⊗ A_sigma, B likewise" : "A = ⊗A_l, B = ⊗B_l") +
                         (Relabelling is { } p ? $", after the sites are relabelled l -> [{string.Join(", ", p)}]" : "") +
                         (Axes is { } a ? $", in the frame where site l's axis {a} (letter l of the word) is Z" : "");
            for (int l = 0; l < Kinds.Count; l++)
                yield return Controls is { } c && c.Contains(l)
                    ? $"site {l}: a control, its bit {(Kinds[l].A == Kind.Diagonal ? "kept" : "flipped")} by A, {(Kinds[l].B == Kind.Diagonal ? "kept" : "flipped")} by B"
                    : Controls is not null
                        ? $"site {l}: A_{l} {Name(Kinds[l].A, $"x{l}.sigma")}, B_{l} {Name(Kinds[l].B, $"y{l}.tau")}"
                        : $"site {l}: A_{l} {Name(Kinds[l].A, $"x{l}")}, B_{l} {Name(Kinds[l].B, $"y{l}")}";
            foreach (var line in System.Lines()) yield return $"  {line}";
        }
    }

    static string Name(Kind k, string turn) => k == Kind.Diagonal ? $"diag(1, {turn})" : $"[[0, {turn}], [1, 0]]";

    static string X(int l) => $"x{l}";
    static string Y(int l) => $"y{l}";

    // A monomial sandwich carrying L onto T(L) - c, read from the letters, or null when the letters are silent.
    // branch: null asks the plain branch and then the transposed one, false or true only that one.
    public static Reading? Read(Setting s, Relation rel, bool? branch = null) => Between(s, s, null, rel, branch, same: true);

    // A sandwich in a turned frame: per site one axis (X, Y or Z) is carried onto Z by a one-site Clifford U_l (the
    // cube's h_zx for X, the quarter turn about X for Y, nothing for Z), and the letters are read there, so the element
    // is E·Ad_U with Ad_U = ⊗U_l·(.)·U_l^dagger, a similarity. X dephasing read on its own axis becomes Z dephasing and
    // a field in the YZ plane a field in the XY plane, the theodolite's own ground. Axes is the word of the axes, site
    // 0 first; the all-Z word is Read itself and is not asked. A site's axis is asked only where that site's table can
    // be carried onto itself in that frame, and the first word that reads is returned; null when none does. A site
    // with no field and no jump fits every axis, so a silence costs 3^k readings for k such sites.
    public static Reading? ReadRotated(Setting s, Relation rel, bool? branch = null)
    {
        int n = s.N;
        var fit = new List<char>[n];
        for (int l = 0; l < n; l++) fit[l] = new List<char>();
        foreach (char a in "ZXY")
        {
            var turned = s.Rotated(new string(a, n));
            for (int l = 0; l < n; l++)
            {
                var pieces = SitePieces(turned, l, rateFree: false);
                if (SiteFits(pieces, pieces, rel, branch)) fit[l].Add(a);
            }
        }
        IEnumerable<string> Words(int l)
        {
            if (l == n) { yield return ""; yield break; }
            foreach (char a in fit[l])
                foreach (var rest in Words(l + 1)) yield return a + rest;
        }
        foreach (var word in Words(0))
        {
            if (word.All(a => a == 'Z')) continue;
            if (Read(s.Rotated(word), rel, branch) is { } r) return r with { Axes = word };
        }
        return null;
    }

    // A sandwich conditioned on the bits of the controls (the head comment's Sectors), asked in the plain frame of the
    // other sites first and then in every turned one; null when the setting has no control or nothing reads.
    public static Reading? ReadSectors(Setting s, Relation rel, bool? branch = null)
    {
        int n = s.N;
        var controls = new List<int>();
        var letter = new char[n];
        for (int u = 0; u < n; u++)
        {
            if (s.Jumps.Any(j => j.Site == u && j.Rate != 0)) continue;
            var used = s.H.Select(t => t.Term.Letter(u)).Where(c => c != 'I').Distinct().ToList();
            if (used.Count != 1) continue;
            controls.Add(u);
            letter[u] = used[0];
        }
        if (controls.Count == 0 || controls.Count == n) return null;
        var rest = Enumerable.Range(0, n).Where(l => !controls.Contains(l)).ToArray();
        for (int code = 0; code < (int)Math.Pow(3, rest.Length); code++)
        {
            var word = new char[n];
            foreach (int u in controls) word[u] = letter[u];
            for (int i = 0, c = code; i < rest.Length; i++, c /= 3) word[rest[i]] = "ZXY"[c % 3];
            var w = new string(word);
            bool plain = w.All(a => a == 'Z');
            if (SectorRead(plain ? s : s.Rotated(w), controls, rest, rel, branch) is { } r) return plain ? r : r with { Axes = w };
        }
        return null;
    }

    static Reading? SectorRead(Setting t, List<int> controls, int[] rest, Relation rel, bool? branch)
    {
        int n = t.N, k = controls.Count, m = rest.Length, nb = 1 << k;
        // the sector settings on the other sites, and their constants
        var sectors = new Setting[nb];
        var constant = new long[nb];
        for (int sg = 0; sg < nb; sg++)
        {
            var terms = new Dictionary<string, long>();
            foreach (var (term, c) in t.H)
            {
                long sign = 1;
                for (int i = 0; i < k; i++)
                {
                    char ch = term.Letter(controls[i]);
                    if (ch is not ('I' or 'Z')) throw new InvalidOperationException($"a control's letter is not on Z: {term.ToString(n)}");
                    if (ch == 'Z' && ((sg >> i) & 1) == 1) sign = -sign;
                }
                var letters = new string(rest.Select(term.Letter).ToArray());
                if (letters.All(ch => ch == 'I')) { constant[sg] += sign * c; continue; }
                terms[letters] = (terms.TryGetValue(letters, out long o) ? o : 0) + sign * c;
            }
            var jumps = t.Jumps.Where(j => rest.Contains(j.Site)).Select(j => (Array.IndexOf(rest, j.Site), j.Jump, j.Rate));
            sectors[sg] = new Setting(m, terms.Select(kv => (PauliString.Parse(kv.Key), kv.Value)), jumps);
        }
        string Bits(int sg) => new string(Enumerable.Range(0, k).Select(i => (char)('0' + ((sg >> i) & 1))).ToArray());
        string Xs(int i, int sg) => $"x{rest[i]}.{Bits(sg)}";
        string Ys(int i, int sg) => $"y{rest[i]}.{Bits(sg)}";
        List<(char, long)> Fields(Setting st, int i) => st.H.Where(h => st.Support(h.Term).SequenceEqual(new[] { i }))
            .Select(h => (h.Term.Letter(i), h.Coefficient)).ToList();
        var supports = sectors.SelectMany(Supports).GroupBy(sup => string.Join(",", sup)).Select(g => g.First())
            .OrderBy(sup => string.Join(",", sup), StringComparer.Ordinal).ToList();
        var doneAt = Enumerable.Range(0, m).Select(i => supports.Where(sup => sup.Max() == i).ToList()).ToList();
        var hams = new Dictionary<(int, string), Gauss[,]>();
        Gauss[,] Ham(int sg, int[] sup)
        {
            var key = (sg, string.Join(",", sup));
            if (!hams.TryGetValue(key, out var h)) hams[key] = h = HamiltonianOn(sectors[sg], sup);
            return h;
        }
        foreach (bool transposed in Branches(branch))
            for (int code = 0; code < 1 << (2 * k); code++)
            {
                int flipA = 0, flipB = 0;
                for (int i = 0; i < k; i++)
                {
                    if (((code >> (2 * i)) & 1) == 1) flipA |= 1 << i;
                    if (((code >> (2 * i + 1)) & 1) == 1) flipB |= 1 << i;
                }
                // every block (ket sg, bra tg): the sectors A and B sit on, and the target block
                var blocks = new List<(int Sg, int Tg, int A, int B, int Ts, int Tt)>();
                bool fits = true;
                for (int sg = 0; sg < nb && fits; sg++)
                    for (int tg = 0; tg < nb && fits; tg++)
                    {
                        int a = transposed ? tg : sg, b = transposed ? sg : tg;
                        int ts = a ^ flipA, tt = b ^ flipB;
                        long src = constant[sg] - constant[tg], tgt = constant[ts] - constant[tt];
                        fits = rel == Relation.Pi ? src == -tgt : src == tgt;
                        blocks.Add((sg, tg, a, b, ts, tt));
                    }
                if (!fits) continue;
                var allowed = new List<List<((Kind A, Kind B) K, List<Binomial> Eqs)>>();
                for (int i = 0; i < m; i++)
                {
                    var site = t.Jumps.Where(j => j.Site == rest[i]).Select(j => (j.Jump, j.Rate)).ToList();
                    var opts = new List<((Kind, Kind), List<Binomial>)>();
                    foreach (var (ka, kb) in KindPairs)
                    {
                        var eqs = new List<Binomial>();
                        bool ok = true;
                        foreach (var bl in blocks)
                        {
                            var from = TracelessTimesEight(SiteTimesTwo(Fields(sectors[bl.Sg], i), site, Fields(sectors[bl.Tg], i)));
                            var to = TracelessTimesEight(SiteTimesTwo(Fields(sectors[bl.Ts], i), site, Fields(sectors[bl.Tt], i)));
                            try
                            {
                                var b = SiteEquations(new[] { (from, to) }, rel, ka, kb, transposed, Xs(i, bl.A), Ys(i, bl.B), $"site {rest[i]}");
                                if (b is null) { ok = false; break; }
                                eqs.AddRange(b);
                            }
                            catch (SilentException) { ok = false; break; }
                        }
                        if (ok) opts.Add(((ka, kb), eqs));
                    }
                    allowed.Add(opts);
                }
                if (allowed.Any(o => o.Count == 0)) continue;
                var choice = new int[m];

                List<Binomial>? Equations(int upto)
                {
                    var eqs = new List<Binomial>();
                    for (int i = 0; i <= upto; i++)
                    {
                        eqs.AddRange(allowed[i][choice[i]].Eqs);
                        foreach (var sup in doneAt[i])
                        {
                            string where = $"support {{{string.Join(", ", sup.Select(q => rest[q]))}}}";
                            foreach (var (a, ts) in blocks.Select(bl => (bl.A, bl.Ts)).Distinct())
                            {
                                var e = SupportEquations(Ham(a, sup), Ham(ts, sup), rel, transposed, sup.Select(q => allowed[q][choice[q]].K.A).ToList(),
                                    sup.Select(q => Xs(q, a)).ToList(), left: true, $"{where}, left");
                                if (e is null) return null;
                                eqs.AddRange(e);
                            }
                            foreach (var (b, tt) in blocks.Select(bl => (bl.B, bl.Tt)).Distinct())
                            {
                                var e = SupportEquations(Ham(b, sup), Ham(tt, sup), rel, transposed, sup.Select(q => allowed[q][choice[q]].K.B).ToList(),
                                    sup.Select(q => Ys(q, b)).ToList(), left: false, $"{where}, right");
                                if (e is null) return null;
                                eqs.AddRange(e);
                            }
                        }
                    }
                    return eqs;
                }

                List<string> Unknowns(int upto) => Enumerable.Range(0, upto + 1)
                    .SelectMany(i => Enumerable.Range(0, nb).SelectMany(sg => new[] { Xs(i, sg), Ys(i, sg) })).ToList();

                bool Rec(int i)
                {
                    if (i == m)                                           // the plain copies holding the identity left out
                        return !(rel == Relation.Copy && !transposed && flipA == 0 && flipB == 0
                                 && choice.Select((ci, q) => allowed[q][ci].K).All(kk => kk == (Kind.Diagonal, Kind.Diagonal)));
                    for (int ci = 0; ci < allowed[i].Count; ci++)
                    {
                        choice[i] = ci;
                        var eqs = Equations(i);
                        if (eqs is null || Binomials.Solve(eqs, Unknowns(i)) is null) continue;
                        if (Rec(i + 1)) return true;
                    }
                    return false;
                }

                if (!Rec(0)) continue;
                var all = Equations(m - 1)!;
                var reduced = Binomials.Solve(all, Unknowns(m - 1))!;
                var kinds = new (Kind A, Kind B)[n];
                for (int i = 0; i < k; i++)
                    kinds[controls[i]] = (((flipA >> i) & 1) == 1 ? Kind.AntiDiagonal : Kind.Diagonal, ((flipB >> i) & 1) == 1 ? Kind.AntiDiagonal : Kind.Diagonal);
                for (int i = 0; i < m; i++) kinds[rest[i]] = allowed[i][choice[i]].K;
                return new Reading(transposed, kinds, reduced, all, Controls: controls);
            }
        return null;
    }

    // A sandwich between two settings on the same sites, E·L(from)·E^-1 = T(L(to)) - c, c fixed by the traces (for a
    // copy, c = (Tr L(to) - Tr L(from))/4^N, the spectra the same up to that shift). The per-site constants cancel in
    // the sum, so the binomial reading of each site table and each support decides it as for one setting; the plain
    // diagonal copies count here, since between two settings they are not the identity's group. The first use: F137's
    // exchange of loss and gain, X on every site carrying L(H, loss) onto L(X^N·H·X^N, gain).
    public static Reading? ReadBetween(Setting from, Setting to, Relation rel, bool? branch = null)
    {
        if (from.N != to.N) throw new ArgumentException("the two settings sit on the same number of sites", nameof(to));
        return Between(from, to, null, rel, branch, same: false);
    }

    // A sandwich after a relabelling of the sites, E·P_pi carrying L onto T(L) - c, P_pi moving site l to pi[l]: the
    // sandwich then carries L of the relabelled setting onto T(L) of the setting itself. The relabellings asked are
    // the automorphisms of the setting's supports other than the identity (on a chain the reversal, on a ring its
    // rotations and reflections), and the first that reads is returned; null when none does. A site l may go to v only
    // when site l's table can be carried onto site v's alone (a necessary condition, so nothing is lost), which keeps
    // the enumeration small where the supports alone do not (fields only, a complete graph).
    public static Reading? ReadRelabelled(Setting s, Relation rel, bool? branch = null)
    {
        int n = s.N;
        var pieces = Enumerable.Range(0, n).Select(l => SitePieces(s, l, rateFree: false)).ToList();
        var fits = new bool[n, n];
        for (int l = 0; l < n; l++)
            for (int v = 0; v < n; v++)
                fits[l, v] = SiteFits(pieces[l], pieces[v], rel, branch);
        foreach (var pi in Automorphisms(s, (l, v) => fits[l, v]))
            if (Between(s.Relabelled(pi), s, pi, rel, branch, same: false) is { } r) return r;
        return null;
    }

    // some kinds and branch for which the site equation carrying one site's table onto another's has a solution
    static bool SiteFits(IReadOnlyList<Gauss[,]> from, IReadOnlyList<Gauss[,]> to, Relation rel, bool? branch)
    {
        var pairs = from.Zip(to).ToList();
        foreach (bool transposed in Branches(branch))
            foreach (var (ka, kb) in KindPairs)
            {
                try
                {
                    var eqs = SiteEquations(pairs, rel, ka, kb, transposed, "x", "y", "site");
                    if (eqs is not null && Binomials.Solve(eqs, new[] { "x", "y" }) is not null) return true;
                }
                catch (SilentException) { }
            }
        return false;
    }

    // the permutations pi of the sites, not the identity, that carry the set of supports onto itself, enumerated
    // lazily; allowed(l, v) may forbid sending site l to v
    public static IEnumerable<int[]> Automorphisms(Setting s, Func<int, int, bool>? allowed = null)
    {
        int n = s.N;
        var supports = Supports(s);
        var keys = new HashSet<string>(supports.Select(sup => string.Join(",", sup)));
        var pi = new int[n];
        var used = new bool[n];
        IEnumerable<int[]> Rec(int l)
        {
            if (l == n)
            {
                if (Enumerable.Range(0, n).Any(q => pi[q] != q)) yield return (int[])pi.Clone();
                yield break;
            }
            for (int v = 0; v < n; v++)
            {
                if (used[v] || (allowed is not null && !allowed(l, v))) continue;
                pi[l] = v;
                used[v] = true;
                bool ok = supports.Where(sup => sup.Max() == l)
                    .All(sup => keys.Contains(string.Join(",", sup.Select(q => pi[q]).OrderBy(q => q))));
                if (ok)
                    foreach (var p in Rec(l + 1)) yield return p;
                used[v] = false;
            }
        }
        return Rec(0);
    }

    static Reading? Between(Setting from, Setting s, int[]? relabelling, Relation rel, bool? branch, bool same)
    {
        var branches = Branches(branch);
        int n = s.N;
        var pieces = Enumerable.Range(0, n)
            .Select(l => (IReadOnlyList<(Gauss[,], Gauss[,])>)SitePieces(from, l, rateFree: false).Zip(SitePieces(s, l, rateFree: false)).ToList())
            .ToList();
        var supports = Supports(s).Concat(Supports(from)).GroupBy(sup => string.Join(",", sup)).Select(g => g.First())
            .OrderBy(sup => string.Join(",", sup), StringComparer.Ordinal).ToList();
        var hams = supports.ToDictionary(sup => string.Join(",", sup), sup => HamiltonianOn(s, sup));
        var hamsFrom = supports.ToDictionary(sup => string.Join(",", sup), sup => HamiltonianOn(from, sup));
        var doneAt = Enumerable.Range(0, n).Select(l => supports.Where(sup => sup.Max() == l).ToList()).ToList();
        foreach (bool transposed in branches)
        {
            var allowed = new List<List<((Kind A, Kind B) K, List<Binomial> Eqs)>>();
            for (int l = 0; l < n; l++)
            {
                var opts = new List<((Kind, Kind), List<Binomial>)>();
                foreach (var (ka, kb) in KindPairs)
                {
                    try
                    {
                        var b = SiteEquations(pieces[l], rel, ka, kb, transposed, X(l), Y(l), $"site {l}");
                        if (b is not null) opts.Add(((ka, kb), b));
                    }
                    catch (SilentException) { }
                }
                allowed.Add(opts);
            }
            if (allowed.Any(o => o.Count == 0)) continue;
            var choice = new int[n];

            List<Binomial>? Equations(int upto)
            {
                var eqs = new List<Binomial>();
                for (int l = 0; l <= upto; l++)
                {
                    eqs.AddRange(allowed[l][choice[l]].Eqs);
                    foreach (var sup in doneAt[l])
                    {
                        var h = hams[string.Join(",", sup)];
                        var hf = hamsFrom[string.Join(",", sup)];
                        var a = SupportEquations(hf, h, rel, transposed, sup.Select(q => allowed[q][choice[q]].K.A).ToList(),
                            sup.Select(X).ToList(), left: true, $"support {{{string.Join(", ", sup)}}}, left");
                        var b = SupportEquations(hf, h, rel, transposed, sup.Select(q => allowed[q][choice[q]].K.B).ToList(),
                            sup.Select(Y).ToList(), left: false, $"support {{{string.Join(", ", sup)}}}, right");
                        if (a is null || b is null) return null;
                        eqs.AddRange(a);
                        eqs.AddRange(b);
                    }
                }
                return eqs;
            }

            List<string> Unknowns(int upto) => Enumerable.Range(0, upto + 1).SelectMany(l => new[] { X(l), Y(l) }).ToList();

            bool Rec(int l)
            {
                if (l == n) return !same || !Gauge(rel, transposed, Enumerable.Range(0, n).Select(q => allowed[q][choice[q]].K));
                for (int ci = 0; ci < allowed[l].Count; ci++)
                {
                    choice[l] = ci;
                    var eqs = Equations(l);
                    if (eqs is null || Binomials.Solve(eqs, Unknowns(l)) is null) continue;
                    if (Rec(l + 1)) return true;
                }
                return false;
            }

            if (Rec(0))
            {
                var eqs = Equations(n - 1)!;
                var reduced = Binomials.Solve(eqs, Unknowns(n - 1))!;
                return new Reading(transposed, Enumerable.Range(0, n).Select(l => allowed[l][choice[l]].K).ToList(), reduced, eqs, relabelling);
            }
        }
        return null;
    }

    // ---------------------------------------------------------------- every N

    public sealed record Traverse(Shape Shape, int Period, bool Transposed, IReadOnlyList<(Kind A, Kind B)> Kinds,
        Binomials.Reduced System, (int K, GaussBig Num, GaussBig Den)? Turn)
    {
        public IEnumerable<string> Lines()
        {
            string form = Transposed ? "A·rho^T·B" : "A·rho·B";
            yield return Shape switch
            {
                Shape.Uniform => $"uniform: rho -> {form}, one map on every site",
                Shape.Turned => $"turned: rho -> {form}, the map turned about Z from site to site by z with " +
                                (Turn is { } t ? $"z^{t.K} = {Binomial.Ratio(t.Num, t.Den)}" : "any z"),
                Shape.Periodic => $"periodic, period {Period}: rho -> {form}",
                _ => $"stepped: rho -> {form}, the turns step from site to site without a turn of the map",
            };
            for (int c = 0; c < Kinds.Count; c++)
                yield return $"  class {c}: A {Name(Kinds[c].A, "x")}, B {Name(Kinds[c].B, "y")}";
            foreach (var line in System.Lines()) yield return $"  {line}";
        }
    }

    static List<Binomial>? ChainEquations(Cell cell, Relation rel, bool transposed, IReadOnlyList<(Kind A, Kind B)> kinds,
        int length, IReadOnlyList<Gauss[,]> pieces)
    {
        var s = cell.Chain(length);
        var eqs = new List<Binomial>();
        for (int l = 0; l < length; l++)
        {
            var b = SiteEquations(pieces.Select(q => (q, q)).ToList(), rel, kinds[l].A, kinds[l].B, transposed, X(l), Y(l), $"site {l}");
            if (b is null) return null;
            eqs.AddRange(b);
        }
        foreach (var sup in Supports(s))
        {
            var h = HamiltonianOn(s, sup);
            var a = SupportEquations(h, h, rel, transposed, sup.Select(q => kinds[q].A).ToList(), sup.Select(X).ToList(), true,
                $"support {{{string.Join(", ", sup)}}}, left");
            var b = SupportEquations(h, h, rel, transposed, sup.Select(q => kinds[q].B).ToList(), sup.Select(Y).ToList(), false,
                $"support {{{string.Join(", ", sup)}}}, right");
            if (a is null || b is null) return null;
            eqs.AddRange(a);
            eqs.AddRange(b);
        }
        return eqs;
    }

    // rewrite the binomials in new unknowns: each old unknown is a monomial in the new ones
    static List<Binomial> Substitute(IEnumerable<Binomial> eqs, Func<string, IReadOnlyDictionary<string, int>> rule)
    {
        var outs = new List<Binomial>();
        foreach (var b in eqs)
        {
            var e = new Dictionary<string, int>();
            foreach (var (u, k) in b.Exponents)
                foreach (var (v, m) in rule(u))
                    e[v] = (e.TryGetValue(v, out int old) ? old : 0) + k * m;
            outs.Add(b with { Exponents = e.Where(kv => kv.Value != 0).ToDictionary(kv => kv.Key, kv => kv.Value) });
        }
        return outs;
    }

    static Binomial Pin(string u, int k, string origin) => new(new Dictionary<string, int> { [u] = k }, GaussBig.One, GaussBig.One, origin);

    static Binomial Relate(string a, int ka, string b, int kb, string origin)
        => new(new Dictionary<string, int> { [a] = ka, [b] = kb }, GaussBig.One, GaussBig.One, origin);

    // The relations a turn about Z by z (|z| = 1) imposes on the steps: Ad_u, u = diag(1, z), leaves a diagonal left
    // factor of a plain sandwich and an anti-diagonal one of a transposed sandwich unchanged, and multiplies the others'
    // turns by z^-2 (plain, either side) or z^2 (transposed, left) and z^-2 (transposed, right).
    static List<Binomial> TurnRelations(bool transposed, (Kind A, Kind B) k)
    {
        var rel = new List<Binomial>();
        if (!transposed)
        {
            rel.Add(k.A == Kind.Diagonal ? Pin("r", 1, "turn") : Relate("r", 1, "z", 2, "turn"));
            rel.Add(k.B == Kind.Diagonal ? Pin("s", 1, "turn") : Relate("s", 1, "z", 2, "turn"));
        }
        else
        {
            rel.Add(k.A == Kind.Diagonal ? Relate("r", 1, "z", -2, "turn") : Pin("r", 1, "turn"));
            rel.Add(k.B == Kind.Diagonal ? Relate("s", 1, "z", 2, "turn") : Pin("s", 1, "turn"));
        }
        return rel;
    }

    // A chain built from one repeated cell, read from the letters for every length (open chains), or null when the
    // letters are silent. rateFree: the element must hold at every rate of every channel. branch: as in Read.
    public static Traverse? EveryN(Cell cell, Relation rel, bool rateFree = false, int maxPeriod = 4, bool? branch = null)
    {
        var branches = Branches(branch);
        cell = cell.Normalized();
        int span = cell.Span;
        var pieces = Pieces(cell.Fields, cell.Jumps, rateFree);
        var allowed = new Dictionary<bool, List<(Kind, Kind)>>();
        foreach (bool transposed in branches)
        {
            allowed[transposed] = new List<(Kind, Kind)>();
            foreach (var kp in KindPairs)
            {
                try
                {
                    if (SiteEquations(pieces.Select(q => (q, q)).ToList(), rel, kp.Item1, kp.Item2, transposed, "u", "v", "site") is not null)
                        allowed[transposed].Add(kp);
                }
                catch (SilentException) { }
            }
        }
        int length = span + 1;
        var geo = new List<(bool Tr, (Kind A, Kind B) K, List<Binomial> Eqs)>();
        IReadOnlyDictionary<string, int> Geometric(string u)
        {
            int l = int.Parse(u[1..]);
            var (root, step) = u[0] == 'x' ? ("X0", "r") : ("Y0", "s");
            return l == 0 ? new Dictionary<string, int> { [root] = 1 } : new Dictionary<string, int> { [root] = 1, [step] = l };
        }
        foreach (bool transposed in branches)
            foreach (var kp in allowed[transposed])
            {
                if (Gauge(rel, transposed, new[] { kp })) continue;
                try
                {
                    var eqs = ChainEquations(cell, rel, transposed, Enumerable.Repeat(kp, length).ToList(), length, pieces);
                    if (eqs is not null) geo.Add((transposed, kp, Substitute(eqs, Geometric)));
                }
                catch (SilentException) { }
            }
        var steps = new[] { "X0", "Y0", "r", "s" };
        IReadOnlyDictionary<string, int> NoStep(string u) => u is "r" or "s"
            ? new Dictionary<string, int>()
            : new Dictionary<string, int> { [u] = 1 };
        foreach (var (tr, kp, g) in geo)
        {
            // uniform: the steps are 1, so the system is read in the turns alone
            var red = Binomials.Solve(Substitute(g, NoStep), new[] { "X0", "Y0" });
            if (red is not null) return new Traverse(Shape.Uniform, 1, tr, new[] { kp }, red, null);
        }
        Traverse? stepped = null;
        foreach (var (tr, kp, g) in geo)
        {
            var red = Binomials.Solve(g, steps);
            if (red is null) continue;
            var turned = Binomials.Solve(g.Concat(TurnRelations(tr, kp)), new[] { "X0", "Y0", "r", "s", "z" });
            if (turned is not null)
            {
                var alone = turned.Alone("z");
                if (alone is null || alone.Value.Num.Norm == alone.Value.Den.Norm)
                    return new Traverse(Shape.Turned, 1, tr, new[] { kp }, turned, alone);
            }
            stepped ??= new Traverse(Shape.Stepped, 1, tr, new[] { kp }, red, null);
        }
        for (int p = 2; p <= maxPeriod; p++)
        {
            int len = p + span - 1;
            var classes = Enumerable.Range(0, p).SelectMany(c => new[] { $"P{c}", $"Q{c}" }).ToList();
            IReadOnlyDictionary<string, int> Periodic(string u)
            {
                int l = int.Parse(u[1..]);
                return new Dictionary<string, int> { [$"{(u[0] == 'x' ? 'P' : 'Q')}{l % p}"] = 1 };
            }
            foreach (bool transposed in branches)
                foreach (var pattern in Patterns(allowed[transposed], p))
                {
                    if (Gauge(rel, transposed, pattern)) continue;
                    try
                    {
                        var eqs = ChainEquations(cell, rel, transposed, Enumerable.Range(0, len).Select(l => pattern[l % p]).ToList(), len, pieces);
                        if (eqs is null) continue;
                        var red = Binomials.Solve(Substitute(eqs, Periodic), classes);
                        if (red is not null) return new Traverse(Shape.Periodic, p, transposed, pattern, red, null);
                    }
                    catch (SilentException) { }
                }
        }
        return stepped;
    }

    static IEnumerable<IReadOnlyList<(Kind, Kind)>> Patterns(IReadOnlyList<(Kind, Kind)> options, int p)
    {
        if (p == 0)
        {
            yield return Array.Empty<(Kind, Kind)>();
            yield break;
        }
        foreach (var head in options)
            foreach (var tail in Patterns(options, p - 1))
                yield return new[] { head }.Concat(tail).ToList();
    }
}
