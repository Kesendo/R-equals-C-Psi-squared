using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for Theorem 8 of the complement connection (<c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>,
/// typed <see cref="PalindromeComplementConnectionClaim.PlanarHoppingPalindrome"/>): hopping bonds J(XX + YY) + D(XY − YX) on any
/// graph, a diagonal part unchanged under flipping every bit, one-site fields a·X + b·Y + c·Z, every site dephased along Z. The rule,
/// read on the block K of H on the states with at most one excitation (every cycle of K a real product of entries, no Z field), is held
/// against Theorem 2 walked on the bitstrings (<see cref="ComplementConnectionGraph"/>, every entry a Gaussian integer), which never
/// forms K: random rows on random graphs, rows built from a gauge (pairing as built, one entry turned in sign on every second), rows with repeated,
/// reversed, cancelling and zero bonds and repeated field sites, the uniform twisted rings, pure DM on named graphs. Three
/// controls must disagree: the rule with the fields left out of K, with the D parts left out, and with the fields' entries
/// a + ib in place of a − ib.</summary>
public sealed class PlanarHoppingWitness : IInspectable
{
    public sealed record FamilyReading(string Name, int Rows, int Pairing, int Mismatches, string? FirstMismatch)
    {
        public bool Agrees => Mismatches == 0;
    }

    public IReadOnlyList<FamilyReading> Families { get; }
    public int FieldsIgnoredMismatches { get; }
    public int TwistsIgnoredMismatches { get; }
    public int ConjugatedFieldsMismatches { get; }
    public bool AllAgree => Families.All(f => f.Agrees);
    public bool ControlsBreak => FieldsIgnoredMismatches > 0 && TwistsIgnoredMismatches > 0 && ConjugatedFieldsMismatches > 0;

    /// <summary>Theorem 2 on n sites for the row: palindrome iff the near and far counts agree.</summary>
    public static bool PairsByTheGraph(int n, IReadOnlyList<(int I, int J, long Hop, long Dm)> bonds,
        IReadOnlyList<(int Site, long X, long Y, long Z)> fields, IReadOnlyList<(int I, int J, long Zz)> diagonal)
    {
        var terms = new List<(PauliMask, BigInteger)>();
        string On(params (int Site, char Letter)[] letters)
        {
            var s = Enumerable.Repeat('I', n).ToArray();
            foreach (var (site, letter) in letters) s[site] = letter;
            return new string(s);
        }
        foreach (var (i, j, hop, dm) in bonds)
        {
            if (hop != 0) { terms.Add((PauliMask.Parse(On((i, 'X'), (j, 'X'))), hop)); terms.Add((PauliMask.Parse(On((i, 'Y'), (j, 'Y'))), hop)); }
            if (dm != 0) { terms.Add((PauliMask.Parse(On((i, 'X'), (j, 'Y'))), dm)); terms.Add((PauliMask.Parse(On((i, 'Y'), (j, 'X'))), -(BigInteger)dm)); }
        }
        foreach (var (site, x, y, z) in fields)
        {
            if (x != 0) terms.Add((PauliMask.Parse(On((site, 'X'))), x));
            if (y != 0) terms.Add((PauliMask.Parse(On((site, 'Y'))), y));
            if (z != 0) terms.Add((PauliMask.Parse(On((site, 'Z'))), z));
        }
        foreach (var (i, j, zz) in diagonal) terms.Add((PauliMask.Parse(On((i, 'Z'), (j, 'Z'))), zz));
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, terms));
        return near == far;
    }

    public PlanarHoppingWitness(int rows = 1500, int seed = 20261111)
    {
        var families = new List<FamilyReading>();
        int fieldsIgnored = 0, twistsIgnored = 0, conjugated = 0;
        FamilyReading Read(string name, IEnumerable<(string Label, int N, (int, int, long, long)[] Bonds, (int, long, long, long)[] Fields, (int, int, long)[] Diag)> rowsIn)
        {
            int count = 0, pairing = 0, bad = 0;
            string? first = null;
            foreach (var (label, n, b, f, d) in rowsIn)
            {
                bool rule = PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, b, f);
                bool graph = PairsByTheGraph(n, b, f, d);
                count++;
                if (graph) pairing++;
                if (rule != graph) { bad++; first ??= label; }
                if (graph != PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, b, f, PalindromeComplementConnectionClaim.PlanarHoppingReading.FieldsIgnored)) fieldsIgnored++;
                if (graph != PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, b, f, PalindromeComplementConnectionClaim.PlanarHoppingReading.TwistsIgnored)) twistsIgnored++;
                if (graph != PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, b, f, PalindromeComplementConnectionClaim.PlanarHoppingReading.FieldsConjugated)) conjugated++;
            }
            return new FamilyReading(name, count, pairing, bad, first);
        }

        // random rows: random graphs on 2 to 7 sites, each bond hopping, DM or both, a ZZ weight from {−2 … 2} on a third
        // of the bonds, fields of random direction in the plane on about half the sites, a Z part on one field in twenty
        var rng = new Random(seed);
        IEnumerable<(string, int, (int, int, long, long)[], (int, long, long, long)[], (int, int, long)[])> RandomRows()
        {
            long[] hops = { 1, -1, 2 }, dms = { 1, -1, 2, 3 }, xs = { 0, 1, -1, 2 }, ys = { 0, 0, 1, -1, 2 };
            for (int r = 0; r < rows; r++)
            {
                int n = rng.Next(2, 8);
                var pairs = (from i in Enumerable.Range(0, n) from j in Enumerable.Range(i + 1, n - i - 1) select (i, j)).ToList();
                for (int k = pairs.Count - 1; k > 0; k--) { int m = rng.Next(k + 1); (pairs[k], pairs[m]) = (pairs[m], pairs[k]); }   // Fisher-Yates
                var bonds = new List<(int, int, long, long)>();
                foreach (var (i, j) in pairs.Take(rng.Next(1, Math.Min(pairs.Count, n + 2) + 1)))
                {
                    int kind = rng.Next(3);
                    bonds.Add((i, j, kind == 1 ? 0 : hops[rng.Next(hops.Length)], kind == 0 ? 0 : dms[rng.Next(dms.Length)]));
                }
                var fields = Enumerable.Range(0, n).Where(_ => rng.Next(2) == 0)
                    .Select(l => (l, xs[rng.Next(xs.Length)], ys[rng.Next(ys.Length)], rng.Next(20) == 0 ? 1L : 0L)).ToArray();
                var diag = bonds.Where(_ => rng.Next(3) == 0).Select(b => (b.Item1, b.Item2, (long)rng.Next(-2, 3))).ToArray();
                yield return ($"random row {r} (N = {n})", n, bonds.ToArray(), fields, diag);
            }
        }
        families.Add(Read("random graphs on 2 to 7 sites, hopping and DM bonds, ZZ, fields (a Z part on one in twenty)", RandomRows()));

        // rows built from a gauge: Gaussian site phases u_l, every bond 2(J − iD) = 2r·ū_i·u_j and every field a − ib = r'·u_l, so
        // that K is real after the gauge and the row pairs; on every second one D or one b is turned in sign, which breaks the row
        // where that entry lies on a cycle of K and is neither real nor purely imaginary (a zero entry, or a turn that only
        // negates the entry, changes no cycle product's reality)
        IEnumerable<(string, int, (int, int, long, long)[], (int, long, long, long)[], (int, int, long)[])> GaugeRows()
        {
            var phases = new (long Re, long Im)[] { (1, 0), (0, 1), (1, 1), (2, 1), (1, 2), (3, 4), (1, -1), (2, -1) };
            for (int r = 0; r < 300; r++)
            {
                int n = rng.Next(2, 8);
                var u = Enumerable.Range(0, n).Select(_ => phases[rng.Next(phases.Length)]).ToArray();
                var pairs = (from i in Enumerable.Range(0, n) from j in Enumerable.Range(i + 1, n - i - 1) select (i, j)).ToList();
                var chosen = pairs.Where(_ => rng.Next(2) == 0).ToList();
                if (chosen.Count == 0) chosen.Add(pairs[0]);
                var bonds = new List<(int, int, long, long)>();
                foreach (var (i, j) in chosen)
                {
                    long s = rng.Next(2) == 0 ? 1 : -1;
                    long re = u[i].Re * u[j].Re + u[i].Im * u[j].Im, im = u[i].Re * u[j].Im - u[i].Im * u[j].Re;   // ū_i·u_j
                    bonds.Add((i, j, s * re, -s * im));
                }
                var fields = Enumerable.Range(0, n).Where(_ => rng.Next(2) == 0).Select(l => { long s = rng.Next(2) == 0 ? 1 : -1; return (l, s * u[l].Re, -s * u[l].Im, 0L); }).ToList();
                if (r % 2 == 1)
                {
                    if (rng.Next(2) == 0 || fields.Count == 0) { int k = rng.Next(bonds.Count); var b = bonds[k]; bonds[k] = (b.Item1, b.Item2, b.Item3, -b.Item4); }
                    else { int k = rng.Next(fields.Count); var f = fields[k]; fields[k] = (f.Item1, f.Item2, -f.Item3, 0L); }
                }
                yield return ($"gauge row {r} (N = {n})", n, bonds.ToArray(), fields.ToArray(), Array.Empty<(int, int, long)>());
            }
        }
        families.Add(Read("rows built from a gauge on 2 to 7 sites, one entry turned on every second", GaugeRows()));

        // repeated, reversed, cancelling and zero bonds, repeated field sites (Z parts that cancel among them)
        IEnumerable<(string, int, (int, int, long, long)[], (int, long, long, long)[], (int, int, long)[])> EdgeRows()
        {
            for (int r = 0; r < 200; r++)
            {
                int n = rng.Next(2, 6);
                var bonds = new List<(int, int, long, long)>();
                for (int k = rng.Next(1, 2 * n); k > 0; k--)
                {
                    int i = rng.Next(n), j = rng.Next(n - 1); if (j >= i) j++;
                    long hop = rng.Next(-2, 3), dm = rng.Next(-2, 3);
                    bonds.Add((i, j, hop, dm));
                    if (rng.Next(4) == 0) bonds.Add((i, j, -hop, -dm));                 // cancels the bond just placed
                    if (rng.Next(4) == 0) bonds.Add((j, i, hop, dm));                    // the same pair read backwards
                }
                var fields = new List<(int, long, long, long)>();
                for (int k = rng.Next(0, n + 2); k > 0; k--)
                {
                    int l = rng.Next(n);
                    long z = rng.Next(5) == 0 ? rng.Next(-1, 2) : 0;
                    fields.Add((l, rng.Next(-1, 3), rng.Next(-1, 2), z));
                    if (z != 0 && rng.Next(2) == 0) fields.Add((l, 0, 0, -z));          // a Z part cancelled on the same site
                }
                yield return ($"edge row {r} (N = {n})", n, bonds.ToArray(), fields.ToArray(), Array.Empty<(int, int, long)>());
            }
        }
        families.Add(Read("repeated, reversed, cancelling and zero bonds, repeated field sites, on 2 to 5 sites", EdgeRows()));

        // the uniform twisted rings: t = 2(J − iD) on every bond, with ZZ; pairs exactly where N·arg t is a multiple of π
        IEnumerable<(string, int, (int, int, long, long)[], (int, long, long, long)[], (int, int, long)[])> Rings()
        {
            foreach (var (j, d) in new (long, long)[] { (1, 0), (0, 1), (1, 1), (1, -1), (3, 4), (2, 1) })
                for (int n = 3; n <= 10; n++)
                {
                    var bonds = Enumerable.Range(0, n).Select(k => (k, (k + 1) % n, j, d)).ToArray();
                    var diag = Enumerable.Range(0, n).Select(k => (k, (k + 1) % n, 1L)).ToArray();
                    yield return ($"ring N = {n}, J = {j}, D = {d}", n, bonds, Array.Empty<(int, long, long, long)>(), diag);
                }
        }
        families.Add(Read("the uniform twisted rings, J + iD in six directions, N = 3 to 10", Rings()));

        // pure DM on named graphs: bipartite or not
        IEnumerable<(string, int, (int, int, long, long)[], (int, long, long, long)[], (int, int, long)[])> Named()
        {
            var graphs = new (string Name, int N, (int, int)[] Edges)[]
            {
                ("triangle", 3, new[] { (0, 1), (1, 2), (0, 2) }),
                ("K4", 4, new[] { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) }),
                ("square", 4, new[] { (0, 1), (1, 2), (2, 3), (0, 3) }),
                ("pentagon", 5, new[] { (0, 1), (1, 2), (2, 3), (3, 4), (0, 4) }),
                ("the three-cube", 8, new[] { (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7) }),
                ("K_{2,3}", 5, new[] { (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4) }),
                ("the bowtie", 5, new[] { (0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4) }),
            };
            foreach (var (name, n, edges) in graphs)
                yield return ($"pure DM on {name}", n, edges.Select(e => (e.Item1, e.Item2, 0L, 1L)).ToArray(), Array.Empty<(int, long, long, long)>(), Array.Empty<(int, int, long)>());
        }
        families.Add(Read("pure DM on named graphs (pairs exactly on the bipartite ones)", Named()));

        Families = families;
        FieldsIgnoredMismatches = fieldsIgnored;
        TwistsIgnoredMismatches = twistsIgnored;
        ConjugatedFieldsMismatches = conjugated;
    }

    public string DisplayName => "PlanarHoppingWitness (Theorem 8 of the complement connection)";

    public string Summary =>
        $"Theorem 8 on {Families.Count} families: {Families.Sum(f => f.Rows)} rows against Theorem 2 on the bitstrings, " +
        $"{Families.Sum(f => f.Mismatches)} disagreeing; the controls (the fields left out of K, the D parts left out, the fields' " +
        $"entries conjugated) disagree {FieldsIgnoredMismatches}, {TwistsIgnoredMismatches}, {ConjugatedFieldsMismatches} times";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var f in Families)
                yield return new InspectableNode(
                    displayName: f.Name,
                    summary: $"{f.Rows} rows, {f.Pairing} pairing, " +
                             (f.Agrees ? "every one agreeing with Theorem 2" : $"{f.Mismatches} DISAGREEING, first {f.FirstMismatch}"),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the controls, which must disagree",
                summary: $"the rule with the fields left out of K: {FieldsIgnoredMismatches} rows; with the D parts left out: " +
                         $"{TwistsIgnoredMismatches}; with the fields' entries a + ib in place of a − ib: {ConjugatedFieldsMismatches}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
