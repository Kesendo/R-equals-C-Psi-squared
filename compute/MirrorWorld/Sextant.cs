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
// Its decision has two exact answers and two that are not. WITNESSED: an element of the closure reflects (a yes, exact
// at this N), or there are no jumps, where -L^dagger = L and the identity reflects with c = 0, or every jump is a Pauli
// letter at a positive rate (F158's hypothesis) and EndCount's verdict is an exact palindrome; there the element can
// sit outside the closure, a site-dependent string (EndCount's colouring) or a product with a turned letter such as
// X⊗Z⊗(Y − 2X), whose angle the couplings set. OBSTRUCTED: an exact broken EndCount verdict there (F158's necessity),
// or an odd centred trace Tr((L + c/2)^m) not zero modulo a prime, so no element of any kind, inside the search or
// outside it, reflects (a no, exact). SYMMETRIC-MOD-P: every odd trace up to the dimension vanishes modulo the two fixed
// primes and no element of the closure reflects; evidence that the spectrum is symmetric, not a proof (the integers
// could be multiples of both primes), and even a symmetric spectrum is a pairing, not yet an element (Jordan blocks).
// DM bonds under loss read so, where F137 proves the pairing by its own route. OPEN: the traces were read below the
// dimension and none obstructs. A failed search alone is never a no.
//
// The traces are the moment route the repo already holds: Hardness.OddPowerSums (F117's trace face, odd power sums
// of the recentred A + gamma Q under Z dephasing), F120's moment tower (docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md)
// and the third-moment row of experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md. Its fence holds here too: the traces
// decide the pairing of the complex spectrum, not whether the rates (real parts) alone pair. The reflection logic uses
// that L is real in the Pauli basis (it preserves Hermiticity), so its spectrum is closed under conjugation.
//
// Parent: the letter cube. Everything the sextant reads with (coordinates, heights, steps, arrows, the moves) is the
// cube's and inherited; its own are the verdict, the search and the decision. Under F158's hypothesis the decision asks
// EndCount (a sibling under the frame, not a parent), so the zone left SYMMETRIC-MOD-P holds mainly jumps outside it,
// such as loss with DM bonds, and Pauli-letter settings whose end count stays a rank reading.
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
    // F158's hypothesis on the jumps: every one a Pauli letter at a positive rate. There the end count decides the
    // palindrome about -sigma, and its exact readings are the sextant's exact answers; a rank reading is not taken.
    static Decision? EndCountReading(Setting s)
    {
        var letters = new List<string>();
        foreach (var (site, jump, rate) in s.Jumps)
        {
            char? letter = "XYZ".Cast<char?>().FirstOrDefault(c => SameMatrix(jump, LetterMatrix(c!.Value)));
            if (letter is null || rate <= 0) return null;
            letters.Add(new string(Enumerable.Range(0, s.N).Select(l => l == site ? letter.Value : 'I').ToArray()));
        }
        var count = new EndCount(new World(), s.N, s.H.Select(h => (h.Term.ToString(s.N), h.Coefficient)).ToList(), letters);
        var reading = count.Verdict();
        if (!EndCount.IsExact(reading)) return null;
        long c = 2 * s.Jumps.Sum(j => j.Rate);
        return EndCount.IsPalindrome(reading)
            ? new Decision(Palindrome.Witnessed, $"F158, the end count: {reading} (c = {c}; no element of the closure reflects, the element sits outside it)")
            : new Decision(Palindrome.Obstructed, $"F158, the end count: {reading}");
    }

    static bool SameMatrix(Gauss[,] a, Gauss[,] b)
    {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    public Sextant(LetterCube cube) : base(cube) { }

    // left: what the sextant itself produces.
    public override IReadOnlyList<string> Own => new[] { "verdict", "search", "decision" };


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
            ulong inside = (1UL << n) - 1;
            foreach (var (term, _) in H)
                if (((term.X | term.Z) & ~inside) != 0) throw new ArgumentException($"a term of H reaches past site {n - 1}", nameof(h));
            foreach (var (site, _, _) in Jumps)
                if (site < 0 || site >= n) throw new ArgumentOutOfRangeException(nameof(jumps), site, $"sites run from 0 to {n - 1}");
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

    // ---- the no: the odd centred traces ----
    //
    // A witness proves a palindrome; the search failing proves nothing. The no comes from an invariant every reflecting
    // element would force: if any similarity carries L to -L^dagger - c or to -L - c, the spectrum is symmetric about
    // -c/2, so every odd power of M = L + c/2 has trace zero. One odd trace that is not zero rules out every element,
    // inside the search space or outside it. Computed exactly modulo the two primes of ModP (both 1 mod 4, so i is a
    // number there), on 4 dim M, which has Gaussian-integer entries: a trace that is not zero modulo a prime is not zero.
    // When every odd trace up to the dimension vanishes modulo both primes, the characteristic polynomial is even or odd
    // (Newton's identities), so the spectrum is symmetric about -c/2, read modulo two primes.

    // Tr((4 dim M)^m) modulo p for every odd m up to maxPower, the first that is not zero
    public static (int Power, long Prime)? OddTraceObstruction(Setting s, int maxPower)
    {
        var strings = LetterCube.Strings(s.N).ToList();
        int dim = strings.Count;
        var index = new Dictionary<PauliString, int>();
        for (int i = 0; i < dim; i++) index[strings[i]] = i;
        long trace4 = 0;
        var full = new Dictionary<(PauliString, PauliString), Gauss>();
        foreach (var kv in s.HamiltonianTimesFour) Add(full, kv.Key, kv.Value);
        foreach (var kv in s.DissipatorTimesFour) Add(full, kv.Key, kv.Value);
        foreach (var ((r, c), v) in full) if (r == c) trace4 += v.Re;
        foreach (long p in ModP.Primes)
        {
            long root = ModP.SqrtMinusOne(p);
            long G(Gauss g) => ModP.Mod(ModP.Mod(g.Re, p) + ModP.MulMod(ModP.Mod(g.Im, p), root, p), p);
            // the columns of 4 dim M = dim (4L) - (4 Tr L) 1, sparse
            var cols = new List<(int Row, long Value)>[dim];
            for (int j = 0; j < dim; j++) cols[j] = new List<(int, long)>();
            var diag = new long[dim];
            foreach (var ((r, c), v) in full)
            {
                long x = ModP.MulMod(ModP.Mod(dim, p), G(v), p);
                if (r == c) diag[index[r]] = ModP.AddMod(diag[index[r]], x, p);
                else cols[index[c]].Add((index[r], x));
            }
            long shift = ModP.Mod(-trace4, p);
            for (int i = 0; i < dim; i++) cols[i].Add((i, ModP.AddMod(diag[i], shift, p)));
            // P = M^m, dense, advanced by P <- P M column by column
            var P = new long[dim, dim];
            for (int j = 0; j < dim; j++) foreach (var (r, v) in cols[j]) P[r, j] = ModP.AddMod(P[r, j], v, p);
            for (int m = 1; m <= maxPower; m++)
            {
                if (m > 1)
                {
                    var next = new long[dim, dim];
                    for (int j = 0; j < dim; j++)
                        foreach (var (k, v) in cols[j])
                            for (int i = 0; i < dim; i++)
                                if (P[i, k] != 0) next[i, j] = ModP.AddMod(next[i, j], ModP.MulMod(P[i, k], v, p), p);
                    P = next;
                }
                if (m % 2 == 0) continue;
                long tr = 0;
                for (int i = 0; i < dim; i++) tr = ModP.AddMod(tr, P[i, i], p);
                if (tr != 0) return (m, p);
            }
        }
        return null;
    }

    public enum Palindrome { Witnessed, Obstructed, SymmetricModP, Open }

    // The sextant's answer to "does the spectrum pair about -c/2?" (the class header has the levels): WITNESSED
    // (exact), OBSTRUCTED (exact), SYMMETRIC-MOD-P (evidence: the odd traces vanish modulo two primes, no element of
    // the closure reflects), OPEN (read below the dimension). By default the traces are read up to the dimension when
    // it is at most 256 (N <= 4) and up to 64 beyond, which leaves the answer OPEN there unless a trace obstructs.
    public sealed record Decision(Palindrome Kind, string Detail);

    public static Decision Decide(Setting s, int? maxPower = null)
    {
        if (s.DissipatorTimesFour.Count == 0)
            return new Decision(Palindrome.Witnessed, "no jumps: -L^dagger = L, the identity reflects with c = 0");
        var witness = Search(s).FirstOrDefault(f => f.Reading.Verdict is Verdict.Reflection or Verdict.PiReflection);
        if (witness.Map is not null)
            return new Decision(Palindrome.Witnessed, $"{witness.Reading.Verdict}: {witness.Map.Describe()} (c = {witness.Reading.C})");
        if (EndCountReading(s) is { } ended)
            return ended;
        int dim = 1 << (2 * s.N);
        int upTo = maxPower ?? (dim <= 256 ? dim : 64);
        var obstruction = OddTraceObstruction(s, upTo);
        if (obstruction is { } o)
            return new Decision(Palindrome.Obstructed, $"Tr(M^{o.Power}) is not zero modulo {o.Prime}");
        return upTo >= dim
            ? new Decision(Palindrome.SymmetricModP, $"every odd trace up to m = {dim} vanishes modulo two primes (evidence, not proof); no element of the closure reflects")
            : new Decision(Palindrome.Open, $"no odd trace up to m = {upTo} obstructs, below the dimension {dim}; no element of the closure reflects");
    }
}
