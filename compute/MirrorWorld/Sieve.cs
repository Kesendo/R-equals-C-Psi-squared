namespace MirrorWorld;

using static MirrorWorld.LetterCube;

// The sieve (built 2026-11-10; Tom's word for what it does is narrowing down). A whole family of cells is passed
// through the theodolite's every-N reading first: what it reads holds on every open chain (and on the rings its closed
// traverse names), a proof for every N in milliseconds, and nothing more is asked of it. Only the rest goes to the end
// count (F158, EndCount), at fixed N, exactly, where the cost grows with 2^N. Per cell and geometry (chain, ring) the
// N read are split into the covered ones (a route of the reading holds there) and the uncovered ones, where the reading
// is silent and the end count is asked; the verdict:
//   Read      the reading covers some N read, and every uncovered N read is broken (on the chain the open traverse
//             covers every N, so nothing is asked);
//   Broken    no N read is covered, and every N read is broken;
//   Residue   every uncovered N read is palindromic: an element the every-N reading does not see (no product of its
//             site factors, or a relabelling, a turned frame or a longer period, which it leaves out);
//   Mixed     the uncovered N read are palindromic at some N and broken at others.
// The yes side is exact for every N; the other three are exact at the N read and say nothing past them. A silence is
// not a no: the residue is where it isn't, and an element understood there is a reader the sieve can gain.
//
// The first census (ClaudeTasks/2026-11-10-everyn-census, private): every cell of one to three of the nine two-letter
// bond terms, signs ±, under Z dephasing, N up to 6. On the chain the reading and the end count split it with nothing
// between. On the ring the residue is XZ ± ZY and YZ ± ZX, with or without ZZ: broken on every chain read and paired on
// every ring read, and for +XZ+ZY the far element is no product (XXX + YYY at three sites). CZ on every ring bond,
// diagonal and so commuting with the dephasing, carries the ring Hamiltonian of XZ ± ZY (± ZZ) onto that of
// ZX ± YZ (± ZZ), so up to the sign of H the two families are one.
//
// Parent: the theodolite, whose every-N reading is the sieve's first pass; the end count is asked, not inherited.
// Only Pauli-letter jumps (X, Y or Z, any rate), since the end count reads Pauli-string jumps; the rates do not enter
// its kernels.
public sealed class Sieve : GameObject
{
    public Sieve(Theodolite theodolite) : base(theodolite) { }

    // left: what the sieve itself produces.
    public override IReadOnlyList<string> Own => new[] { "sieve", "residue" };

    public enum Verdict { Read, Broken, Residue, Mixed }

    // One cell through the sieve on one geometry: the every-N reading (null when silent), and per N read a letter:
    // E covered by the reading (no end count asked), P/B palindromic/broken by the end count exactly, p/b by its ranks.
    public sealed record Pass(Theodolite.Traverse? Traverse, int From, string Letters, Verdict Verdict);

    public sealed record Entry(string Name, Theodolite.Cell Cell, Pass Chain, Pass Ring);

    static readonly Theodolite.Relation[] Pairing = { Theodolite.Relation.Pi, Theodolite.Relation.Reflection };

    public static Entry Run(string name, Theodolite.Cell cell, int maxN)
    {
        var norm = cell.Normalized();
        JumpLetters(norm.Jumps.Select(j => j.Jump));                                   // Pauli letters only, before any work
        Theodolite.Traverse? open = null, closed = null;
        foreach (var rel in Pairing)
        {
            open ??= Theodolite.EveryN(norm, rel);
            closed ??= Theodolite.EveryN(norm, rel, closed: true);
        }
        var chain = Through(open, Math.Max(3, norm.Span), maxN, n => open is not null, norm.Chain);
        var ring = Through(closed, norm.SmallestRing, maxN, n => closed is not null && closed.Closes(n), norm.Ring);
        return new Entry(name, cell, chain, ring);
    }

    static Pass Through(Theodolite.Traverse? t, int from, int maxN, Func<int, bool> covered, Func<int, Theodolite.Setting> build)
    {
        var letters = new System.Text.StringBuilder();
        for (int n = from; n <= maxN; n++)
        {
            if (covered(n)) { letters.Append('E'); continue; }
            var v = Count(build(n)).Verdict();
            letters.Append(EndCount.IsPalindrome(v) ? (EndCount.IsExact(v) ? 'P' : 'p') : (EndCount.IsExact(v) ? 'B' : 'b'));
        }
        string s = letters.ToString();
        bool anyCovered = s.Contains('E');
        bool anyPal = s.Any(c => c is 'P' or 'p'), anyBroken = s.Any(c => c is 'B' or 'b');
        Verdict verdict = anyPal && anyBroken ? Verdict.Mixed
            : anyPal ? Verdict.Residue
            : anyCovered ? Verdict.Read
            : Verdict.Broken;
        return new Pass(t, from, s, verdict);
    }

    // the end count of a setting, its jumps read off as Pauli strings (each must be X, Y or Z on one site)
    public static EndCount Count(Theodolite.Setting s)
    {
        int n = s.N;
        var h = s.H.Select(t => (t.Term.ToString(n), t.Coefficient)).ToList();
        var jumps = s.Jumps.Select(j => (j.Site, Letter: JumpLetters(new[] { j.Jump })[0]))
            .Select(j => new string('I', j.Site) + j.Letter + new string('I', n - j.Site - 1)).Distinct().ToList();
        return new EndCount(new World(), n, h, jumps);
    }

    static string JumpLetters(IEnumerable<Gauss[,]> jumps)
    {
        var letters = new List<char>();
        foreach (var j in jumps)
        {
            char? hit = "XYZ".Cast<char?>().FirstOrDefault(p => Same(j, LetterMatrix(p!.Value)));
            if (hit is null) throw new ArgumentException("the sieve reads Pauli-letter jumps only (X, Y or Z)");
            letters.Add(hit.Value);
        }
        return new string(letters.ToArray());
    }

    static bool Same(Gauss[,] a, Gauss[,] b)
    {
        for (int i = 0; i < 2; i++)
            for (int k = 0; k < 2; k++)
                if (a[i, k] != b[i, k]) return false;
        return true;
    }

    static readonly string[] BondLetters = { "XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ" };

    // Every cell of one to maxTerms of the nine two-letter bond terms on (0, 1), the first coefficient +1 and the
    // others ±1, named like "+XZ-ZY", with the jumps jumpLetters at rate 1 on every site.
    public static IEnumerable<(string Name, Theodolite.Cell Cell)> BondFamily(int maxTerms, string jumpLetters)
    {
        var jumps = jumpLetters.Select(j => (LetterMatrix(j), 1L)).ToList();
        for (int size = 1; size <= maxTerms; size++)
            foreach (var combo in Combos(BondLetters.Length, size))
                for (int signs = 0; signs < 1 << (size - 1); signs++)
                {
                    var terms = combo.Select((t, i) => (BondLetters[t], i == 0 || ((signs >> (i - 1)) & 1) == 0 ? 1L : -1L)).ToList();
                    yield return (string.Concat(terms.Select(t => (t.Item2 > 0 ? "+" : "-") + t.Item1)),
                        new Theodolite.Cell(jumps, null, terms.Select(t => (t.Item1, new[] { 0, 1 }, t.Item2))));
                }
    }

    static IEnumerable<int[]> Combos(int n, int k)
    {
        if (k == 0) { yield return Array.Empty<int>(); yield break; }
        for (int first = 0; first <= n - k; first++)
            foreach (var rest in Combos(n - first - 1, k - 1))
                yield return new[] { first }.Concat(rest.Select(r => r + first + 1)).ToArray();
    }
}
