namespace MirrorWorld;

using static MirrorWorld.LetterCube;

// The sextant (built 2026-10-14 beside docs/THE_ONE_SQUARE.md §9; the name is Tom's). Why the name: a sextant brings two
// images into coincidence through two mirrors and reads off the angle at which they meet. This object sets the image
// S L S^-1 of a move beside L and beside -L^dagger - c and -L - c, reports whether they coincide, and reads off c, the
// offset at which they do. It reads a setting against the letter cube: given H in Pauli strings with integer coefficients and one-site jumps with integer rates (1 to 6
// sites), it holds the generator exactly in the string basis, its Hamiltonian part and its dissipator apart, and asks
// every move of the cube the relation form of §9's two questions on the full generator: a copy (S L S^-1 = L), the
// reflection of §9 (S L S^-1 = -L^dagger - c, F158's far end), the reflection of F1's Pi (S L S^-1 = -L - c), or
// neither, c fixed by the trace; the readings on L_H and on the dissipator annotate the verdict. Its search asks the
// whole one-site closure, uniform on every site. It reproduces §9 under Z dephasing (R and the shift by Y^N reflect,
// the transpose followed by R is Pi_Z), under loss (a copy keeps the arrow, S_XY reflects) and in a thermal bath,
// where the search finds S_XY times the copies as the only reflections by itself.
//
// Parent: the letter cube. Everything the sextant reads with (coordinates, heights, steps, arrows, the moves) is the
// cube's and inherited; its own are the verdict and the search.
//
// Deliberately outside: spectra (the verdicts are exact identities of the generator, no eigensolver), site
// permutations (EndCount.SymmetryElement holds the half-turn composed with one), non-uniform moves beyond what a
// caller builds by hand.
//
// Words, fenced at the door: "reflection" here is the generator relation, not the chain's site reversal or Seed's
// reflection seat; "copy" is a move commuting with L; "setting" is the system read, not a configuration of MirrorWorld's
// running engine; "verdict" and Sextant.Reading are a move's reading against a setting, not EndCount's palindrome
// verdict and its Reading.
public sealed class Sextant : GameObject
{
    public Sextant(LetterCube cube) : base(cube) { }

    // left: what the sextant itself produces.
    public override IReadOnlyList<string> Own => new[] { "verdict", "search" };


    // A setting: n sites, H as Pauli strings with integer coefficients, and one-site jumps (site, matrix, integer
    // rate). Its generator is held four times over in the string basis, the Hamiltonian part and the dissipator
    // apart, entries Gaussian integers.
    public sealed class Setting
    {
        public int N { get; }
        public IReadOnlyList<(PauliString Term, long Coefficient)> H { get; }
        public IReadOnlyList<(int Site, Gauss[,] Jump, long Rate)> Jumps { get; }
        public IReadOnlyDictionary<(PauliString Row, PauliString Col), Gauss> HamiltonianTimesFour { get; }
        public IReadOnlyDictionary<(PauliString Row, PauliString Col), Gauss> DissipatorTimesFour { get; }

        public Setting(int n, IEnumerable<(PauliString, long)> h, IEnumerable<(int, Gauss[,], long)> jumps)
        {
            if (n < 1 || n > 6) throw new ArgumentOutOfRangeException(nameof(n), n, "1 to 6 sites (4^n strings)");
            N = n;
            H = h.ToList();
            Jumps = jumps.ToList();
            var lh = new Dictionary<(PauliString, PauliString), Gauss>();
            foreach (var s in LetterCube.Strings(n))
                foreach (var (t, c) in H)
                    if (LetterCube.CommutatorCoefficient(t, s) is { } r)
                    {
                        // L_H(s) = -i c [T, s] = -i c 2 i^k Target; times four: -8 c i^(k+1) Target
                        Add(lh, (r.Target, s), (-8 * c) * IPower(r.PhasePower + 1));
                    }
            HamiltonianTimesFour = Clean(lh);
            var ld = new Dictionary<(PauliString, PauliString), Gauss>();
            foreach (var group in Jumps.GroupBy(j => j.Site))
            {
                int l = group.Key;
                if (l < 0 || l >= n) throw new ArgumentOutOfRangeException(nameof(jumps), l, "a jump's site lies on the system");
                var t4 = LetterCube.TransferTimesFour(group.Select(j => (j.Jump, j.Rate)));
                foreach (var s in LetterCube.Strings(n))
                {
                    int b = LetterCube.LetterIndex(s, l);
                    for (int a = 0; a < 4; a++)
                        if (!t4[a, b].IsZero) Add(ld, (LetterCube.WithLetter(s, l, a), s), t4[a, b]);
                }
            }
            DissipatorTimesFour = Clean(ld);
        }

        // H placed on every window of an open chain
        public static Setting Chain(int n, IEnumerable<(string Template, long Coefficient)> bonds, IEnumerable<(int, Gauss[,], long)> jumps)
        {
            var h = new List<(PauliString, long)>();
            foreach (var (t, c) in bonds)
                foreach (var placed in LetterCube.ChainTerms(new[] { t }, n)) h.Add((placed, c));
            return new Setting(n, h, jumps);
        }

        // the same jump at the same rate on every site
        public static IEnumerable<(int, Gauss[,], long)> OnEverySite(int n, Gauss[,] jump, long rate)
            => Enumerable.Range(0, n).Select(l => (l, jump, rate));
    }

    static Gauss IPower(int k) => (((k % 4) + 4) % 4) switch { 0 => Gauss.One, 1 => Gauss.I, 2 => -Gauss.One, _ => -Gauss.I };

    static void Add(Dictionary<(PauliString, PauliString), Gauss> m, (PauliString, PauliString) key, Gauss v)
        => m[key] = m.TryGetValue(key, out var old) ? old + v : v;

    static Dictionary<(PauliString, PauliString), Gauss> Clean(Dictionary<(PauliString, PauliString), Gauss> m)
        => m.Where(kv => !kv.Value.IsZero).ToDictionary(kv => kv.Key, kv => kv.Value);

    public enum OnHamiltonian { Keeps, Flips, Breaks }
    public enum OnDissipator { Keeps, Reflects, Breaks }

    // The verdict on the full L: a copy (S L S^-1 = L), the reflection of §9 (S L S^-1 = -L^dagger - c, F158's far
    // end), the reflection of F1's Pi (S L S^-1 = -L - c, the shape Pi_Z = R D has), or neither. Both reflections pair
    // the spectrum about -c/2. The verdicts are exclusive by precedence, copy before reflection before Pi-reflection:
    // with no jumps every copy is also a reflection with c = 0, and with a Hermitian L (Pauli jumps, no H) every
    // reflection is also a Pi-reflection, and each is reported once, as the first that holds.
    public enum Verdict { Copy, Reflection, PiReflection, Neither }

    // The two questions of §9, asked of a move. The verdict is read on the full L; the readings on L_H (keeps it, flips
    // it to -L_H, or neither) and on the dissipator (keeps it, reflects it to -D^dagger - c, or neither) annotate it, and
    // agree with it for traceless jumps: a jump with a trace hides part of its dissipator in a commutator (1 + iZ at
    // rate 1 with a field Z is pure Z dephasing), so the parts can break where the whole is kept. c is fixed by the
    // trace, c = -2 Tr(L)/4^n, and the commutator part is traceless.
    public sealed record Reading(string Move, OnHamiltonian OnH, OnDissipator OnD, Verdict Verdict, long MinusTwoTraceTimesFour, long Dim)
    {
        public string C => Verdict is Verdict.Reflection or Verdict.PiReflection || OnD == OnDissipator.Reflects
            ? Fraction(MinusTwoTraceTimesFour, 4 * Dim) : "";
    }

    static string Fraction(long num, long den)
    {
        if (num == 0) return "0";
        long g = Gcd(Math.Abs(num), den);
        return den / g == 1 ? $"{num / g}" : $"{num / g}/{den / g}";
    }

    static long Gcd(long a, long b) { while (b != 0) (a, b) = (b, a % b); return a; }

    static Dictionary<(PauliString, PauliString), Gauss> Conjugate(IReadOnlyDictionary<(PauliString Row, PauliString Col), Gauss> x, LetterCube.Move m)
    {
        var y = new Dictionary<(PauliString, PauliString), Gauss>();
        foreach (var ((row, col), v) in x)
        {
            var (ri, rp) = m.Apply(row);
            var (ci, cp) = m.Apply(col);
            Add(y, (ri, ci), rp * v * cp.Conj);           // the phases are units: the inverse is the conjugate
        }
        return Clean(y);
    }

    static bool SameMatrix(IReadOnlyDictionary<(PauliString, PauliString), Gauss> a, IReadOnlyDictionary<(PauliString, PauliString), Gauss> b)
        => a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out var v) && v == kv.Value);

    public static Reading Ask(Setting sys, LetterCube.Move m)
    {
        if (m.Sites.Length != sys.N) throw new ArgumentException("the move acts on the system's sites", nameof(m));
        long dim = 1L << (2 * sys.N);
        var h = sys.HamiltonianTimesFour;
        var hm = Conjugate(h, m);
        var onH = SameMatrix(hm, h) ? OnHamiltonian.Keeps
                : SameMatrix(hm, h.ToDictionary(kv => kv.Key, kv => -kv.Value)) ? OnHamiltonian.Flips
                : OnHamiltonian.Breaks;
        var d = sys.DissipatorTimesFour;
        var dm = Conjugate(d, m);
        long traceTimesFour = 0;
        foreach (var ((r, c), v) in d) if (r == c) traceTimesFour += v.Re;
        // c = -2 Tr(D)/dim: compare dim (S D S^-1) with -dim D^dagger + 2 Tr(D) 1, everything times four
        var lhs = dm.ToDictionary(kv => kv.Key, kv => dim * kv.Value);
        var rhs = new Dictionary<(PauliString, PauliString), Gauss>();
        foreach (var ((r, c), v) in d) Add(rhs, (c, r), (-dim) * v.Conj);
        foreach (var s in LetterCube.Strings(sys.N)) Add(rhs, (s, s), new Gauss(2 * traceTimesFour, 0));
        var onD = SameMatrix(dm, d) ? OnDissipator.Keeps
                : SameMatrix(lhs, Clean(rhs)) ? OnDissipator.Reflects
                : OnDissipator.Breaks;
        // the verdict on the full L = L_H + D (Tr L_H = 0, so Tr L = Tr D)
        var l = new Dictionary<(PauliString, PauliString), Gauss>();
        foreach (var kv in h) Add(l, kv.Key, kv.Value);
        foreach (var kv in d) Add(l, kv.Key, kv.Value);
        var full = Clean(l);
        var lm = Conjugate(full, m);
        var lmScaled = lm.ToDictionary(kv => kv.Key, kv => dim * kv.Value);
        var dagger = new Dictionary<(PauliString, PauliString), Gauss>();
        var minus = new Dictionary<(PauliString, PauliString), Gauss>();
        foreach (var ((r, c), v) in full)
        {
            Add(dagger, (c, r), (-dim) * v.Conj);
            Add(minus, (r, c), (-dim) * v);
        }
        foreach (var s in LetterCube.Strings(sys.N))
        {
            Add(dagger, (s, s), new Gauss(2 * traceTimesFour, 0));
            Add(minus, (s, s), new Gauss(2 * traceTimesFour, 0));
        }
        var verdict = SameMatrix(lm, full) ? Verdict.Copy
                    : SameMatrix(lmScaled, Clean(dagger)) ? Verdict.Reflection
                    : SameMatrix(lmScaled, Clean(minus)) ? Verdict.PiReflection
                    : Verdict.Neither;
        return new Reading(m.Name, onH, onD, verdict, -2 * traceTimesFour, dim);
    }

    // The sextant's search: every uniform one-site map of the closure, asked both questions; every map whose verdict
    // is not Neither (copies, reflections and Pi-reflections), each once modulo a global phase.
    public static IReadOnlyList<(LetterCube.SiteMove Map, Reading Reading)> Search(Setting sys)
    {
        var found = new List<(LetterCube.SiteMove, Reading)>();
        foreach (var m in LetterCube.Closure())
        {
            var r = Ask(sys, LetterCube.Move.Uniform(m, sys.N, m.Describe()));
            if (r.Verdict != Verdict.Neither) found.Add((m, r));
        }
        return found;
    }

}
