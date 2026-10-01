using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the complement connection (claim <see cref="PalindromeComplementConnectionClaim"/>,
/// proof <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>): with one single-site Pauli jump on
/// every site, turn each jump to Z and the far end of F158 becomes the flat sections of
///
/// <code>
///     d_y / d_x = H_xy / H_x̄ȳ        (x̄ the bitwise complement)
/// </code>
///
/// on H's hopping graph over bitstrings. The near count is the number of components of the hopping graph
/// Γ_H, the far count the number of GOOD components of its union Γ with the complement image (equal
/// moduli on every edge, H_xx = H_x̄x̄, trivial holonomy), and the palindrome holds exactly when every
/// component of Γ is good (Theorem 2).
///
/// <para><b>A different road from the gate.</b> The sober base's <c>EndCount.ComplementConnection</c> and
/// <c>EndCount.SectorConnection</c> are the gate; this witness shares none of their code. It turns the
/// letters with one proper rotation per site on <see cref="PauliMask"/> and hands the strings to the Core's
/// <see cref="ComplementConnectionGraph"/> (shared with <c>PalindromeSoftCertifier.DecideAtN</c>), which reads
/// each matrix entry off the mask form P(x, z) = i^|x∧z|·X^x·Z^z (so ⟨y ⊕ x|P|y⟩ = i^|x∧z|·(−1)^|z∧y|), holds the section as a
/// quotient of two <see cref="GaussianInteger"/>s and checks every edge by cross-multiplication, with no
/// division and no float. Beside it, it reads <see cref="PalindromeStringSpanWitness"/> on the same row,
/// which never forms a basis state: two constructions meeting.</para>
///
/// <para><b>The theorems it reads.</b> The Hamiltonian is the one of <c>twoend</c> and
/// <c>twoendstrings</c>: bonds of weight 10 on a connected topology (Heisenberg by default, or XX + YY,
/// or ZZ), letter fields of weight 3. With every site dephased under Heisenberg bonds every row lies in
/// Theorem 1's class (one common axis, turned to Z by a global rotation) or in Theorem 3's (mixed axes),
/// and the witness prints the rule's prediction beside the graph's verdict: under one axis a, the
/// palindrome holds exactly when no field lies along a and the fields do not use both other letters;
/// under mixed axes, exactly when two axes occur and every field lies along the third letter, the near
/// end one-dimensional. XX + YY bonds under Z on every site lie in Theorem 1's class as well.</para>
///
/// <para><b>Undephased sites.</b> Where every undephased site keeps a letter (every term touching it uses
/// that one letter there: at ZZ bonds with its field absent or along Z), the witness turns the kept letter
/// to Z as well, splits H into the sectors of their bits and reads Theorem 4: near = Σ_{σ,τ} hom(H_τ, H_σ),
/// far = Σ_{σ,τ} hom(H_τ, H̄_σ), each hom the good components of the pair graph, walked the same way as the
/// flat sections, with the cross-sector terms (σ ≠ τ) printed apart, since they are what the undephased
/// sites add. Where an undephased site keeps no letter (a Heisenberg or XX + YY bond touches it, or a ZZ
/// bond and an X or Y field) the blocks over the dephased bits are matrices and no graph reads the row; the string route decides it,
/// and under Heisenberg bonds with exactly one undephased site Theorem 5's rule is printed beside it: the
/// palindrome holds exactly when some letter is no dephased site's axis and every field, the undephased
/// site's included, lies along it.</para>
///
/// <para>Args: <c>--N</c> (2..<see cref="MaxN"/>, default 3), <c>--deph</c> (X/Y/Z or '.' per site, at
/// least one jump, default all Z), <c>--field</c> (X/Y/Z or '.' per site, default none),
/// <c>--topology</c> (chain|ring|complete, default chain), <c>--bonds</c> (XYZ|XY|ZZ, default XYZ).</para></summary>
public sealed class ComplementConnectionWitness : IInspectable
{
    /// <summary>A cost guard: the graph has 2^N vertices.</summary>
    public const int MaxN = 12;

    private const long Bond = 10, Field = 3;

    private readonly int _n;
    private readonly string _deph, _field, _topology, _bonds;
    private readonly (int A, int B)[] _edges;
    private readonly char?[] _kept;

    public ComplementConnectionWitness(int n, string? deph = null, string? field = null, string? topology = null,
                                       string? bonds = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"--root complement builds a graph on 2^N bitstrings and is guarded at N in 2..{MaxN}; got {n}. " +
                "The theorems carry no N, only this witness does.");
        _n = n;
        _deph = Letters(deph, n, 'Z', nameof(deph));
        if (_deph.All(c => c == '.'))
            throw new ArgumentException("--deph names no jump; F158 needs at least one.", nameof(deph));
        _field = Letters(field, n, '.', nameof(field));
        _bonds = PalindromeStringSpanWitness.ParseBonds(bonds);
        _topology = (topology ?? "chain").ToLowerInvariant();
        _edges = _topology switch
        {
            "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray(),
            "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray(),
            "complete" => (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
            _ => throw new ArgumentException($"--topology takes chain, ring or complete; got \"{topology}\"."),
        };
        _kept = Enumerable.Range(0, n).Select(KeptLetter).ToArray();
    }

    private static string Letters(string? spec, int n, char fallback, string arg)
    {
        if (string.IsNullOrWhiteSpace(spec)) return new string(fallback, n);
        string s = spec.Trim().ToUpperInvariant();
        if (s.Length != n)
            throw new ArgumentException($"--{arg} takes exactly one letter per site (N = {n}); got \"{s}\".");
        foreach (char c in s)
            if (!(c is 'X' or 'Y' or 'Z' or '.' or 'I'))
                throw new ArgumentException($"--{arg} letter '{c}' is not one of X Y Z or '.'.");
        return new string(s.Select(c => c == 'I' ? '.' : c).ToArray());
    }

    /// <summary>The letter an undephased site keeps: the one letter every term touching it uses there (Z
    /// when none touches it), or null when two letters meet at it. Dephased sites keep nothing here.</summary>
    private char? KeptLetter(int site)
    {
        if (_deph[site] != '.') return null;
        var letters = new HashSet<char>();
        if (_edges.Any(e => e.A == site || e.B == site))
            foreach (char p in _bonds) letters.Add(p);
        if (_field[site] != '.') letters.Add(_field[site]);
        return letters.Count switch { 0 => 'Z', 1 => letters.First(), _ => null };
    }

    /// <summary>The number of sites that carry no jump.</summary>
    public int Undephased => _deph.Count(c => c == '.');

    /// <summary>True when a graph reads the row: every site dephased (Theorem 2) or every undephased site
    /// keeping a letter (Theorem 4). Otherwise the string route decides it.</summary>
    public bool GraphReads => Enumerable.Range(0, _n).All(s => _deph[s] != '.' || _kept[s] is not null);

    public string Name => "complement";

    public string DisplayName =>
        $"the complement connection live (N = {_n}, {_topology}, bonds {PalindromeStringSpanWitness.BondName(_bonds)}, " +
        $"dephasing {_deph}, fields {_field})";

    public string Summary
    {
        get
        {
            var r = Read();
            string counts = r.Near is { } near && r.Far is { } far
                ? string.Format(CultureInfo.InvariantCulture, "near {0}, far {1} on the graph", near, far)
                : "no graph reads this row, the string route decides";
            return string.Format(CultureInfo.InvariantCulture, "{0}: {1}; {2}{3}",
                counts, r.Palindrome ? "palindrome" : "broken", r.TheoremClass,
                r.Predicted is { } p ? $" predicts {(p ? "palindrome" : "broken")}" : "");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>What one inspect recomputes. Near and Far are the graph's counts (null where no graph reads
    /// the row); with every site dephased Near is the number of hopping components, Far the number of good
    /// components of the union graph, which has UnionComponents components. Sectors is the number of
    /// sectors of the kept letters (1 with every site dephased), CrossNear and CrossFar the part of the two
    /// counts from pairs of distinct sectors. Palindrome is the graph's verdict where it reads the row and
    /// the string route's otherwise; Exact says whether that verdict is exact. Predicted is the closed rule
    /// of Theorem 1, 3 or 5 where one applies.</summary>
    public sealed record Reading(
        int? Near, int? Far, int? UnionComponents, int Sectors, int CrossNear, int CrossFar,
        bool Palindrome, bool Exact, string TheoremClass, bool? Predicted);

    private Reading? _reading;
    private PalindromeStringSpanWitness.Reading? _strings;

    private PalindromeStringSpanWitness.Reading Strings() =>
        _strings ??= new PalindromeStringSpanWitness(_n, _deph, _field, _topology, bonds: _bonds).Read();

    public Reading Read()
    {
        if (_reading is not null) return _reading;
        var (cls, predicted) = Prediction();
        if (!GraphReads)
        {
            var s = Strings();
            return _reading = new Reading(null, null, null, 0, 0, 0,
                PalindromeStringSpanWitness.IsPalindrome(s.Verdict), PalindromeStringSpanWitness.IsExact(s.Verdict),
                cls, predicted);
        }

        var h = ComplementConnectionGraph.Columns(_n, TurnedTerms().Select(t => (t.S, new BigInteger(t.C))));
        int dephased = 0, undephased = 0;
        for (int s = 0; s < _n; s++)
            if (_deph[s] != '.') dephased |= 1 << s; else undephased |= 1 << s;
        foreach (var (column, row) in h.Select((row, column) => (column, row)))
            foreach (int target in row.Keys)
                if (((target ^ column) & undephased) != 0)
                    throw new InvalidOperationException(
                        $"H moves an undephased bit ({column} to {target}) although every undephased site keeps a letter");

        var stringsOfS = SubMasks(dephased).ToList();
        var sectors = SubMasks(undephased).ToList();
        int near = 0, far = 0, crossNear = 0, crossFar = 0, union = 0;
        foreach (int sigma in sectors)
            foreach (int tau in sectors)
            {
                // hom(H_τ, H_σ): B = H_σ; and hom(H_τ, H̄_σ): B = H_σ read on complemented S bits
                Func<int, int, GaussianInteger> a = (x, y) => ComplementConnectionGraph.Entry(h, x | tau, y | tau);
                Func<int, int, GaussianInteger> b = (x, y) => ComplementConnectionGraph.Entry(h, x | sigma, y | sigma);
                Func<int, int, GaussianInteger> bBar = (x, y) =>
                    ComplementConnectionGraph.Entry(h, (x ^ dephased) | sigma, (y ^ dephased) | sigma);
                // the candidate neighbours of x: the rows of its columns in H_τ, H_σ and H̄_σ, on the bits of S
                Func<int, IEnumerable<int>> nearCandidates = x =>
                    h[x | tau].Keys.Concat(h[x | sigma].Keys).Select(y => y & dephased);
                Func<int, IEnumerable<int>> farCandidates = x =>
                    h[x | tau].Keys.Select(y => y & dephased)
                        .Concat(h[(x ^ dephased) | sigma].Keys.Select(y => (y ^ dephased) & dephased));
                var (_, nearGood) = ComplementConnectionGraph.Hom(stringsOfS, a, b, nearCandidates);
                var (farComponents, farGood) = ComplementConnectionGraph.Hom(stringsOfS, a, bBar, farCandidates);
                near += nearGood;
                far += farGood;
                if (sigma != tau) { crossNear += nearGood; crossFar += farGood; }
                else if (sectors.Count == 1) union = farComponents;
            }
        return _reading = new Reading(near, far, sectors.Count == 1 ? union : null, sectors.Count, crossNear, crossFar,
            near == far, true, cls, predicted);
    }

    private static IEnumerable<int> SubMasks(int mask)
    {
        for (int s = mask; ; s = (s - 1) & mask)
        {
            yield return s;
            if (s == 0) yield break;
        }
    }

    // ---- the turn to Z, and the entries ----

    /// <summary>The proper rotation of one site's letters that sends its jump letter to Z, as
    /// (image letter, sign) for X, Y, Z: identity for Z; for X the Hadamard (X ↔ Z, Y → −Y); for Y the
    /// cycle X → Y → Z → X. Both are rotations (determinant +1), so they are conjugations by one-site
    /// unitaries and keep every kernel dimension.</summary>
    public static (char Letter, int Sign) Turn(char jump, char letter) => (jump, letter) switch
    {
        ('Z', _) => (letter, 1),
        ('X', 'X') => ('Z', 1),
        ('X', 'Z') => ('X', 1),
        ('X', 'Y') => ('Y', -1),
        ('Y', 'X') => ('Y', 1),
        ('Y', 'Y') => ('Z', 1),
        ('Y', 'Z') => ('X', 1),
        _ => throw new ArgumentException($"no turn for jump {jump}, letter {letter}"),
    };

    /// <summary>The letter turned to Z on each site: its jump where it is dephased, the kept letter where
    /// it is not.</summary>
    private char TurnLetter(int site) => _deph[site] != '.' ? _deph[site] : _kept[site] ?? 'Z';

    private IEnumerable<(PauliMask S, long C)> TurnedTerms()
    {
        var terms = new List<(char[] L, long C)>();
        foreach (var (a, b) in _edges)
            foreach (char p in _bonds)
            {
                var l = Enumerable.Repeat('I', _n).ToArray();
                l[a] = l[b] = p;
                terms.Add((l, Bond));
            }
        for (int s = 0; s < _n; s++)
            if (_field[s] != '.')
            {
                var l = Enumerable.Repeat('I', _n).ToArray();
                l[s] = _field[s];
                terms.Add((l, Field));
            }
        foreach (var (l, c) in terms)
        {
            int sign = 1;
            var turned = new PauliLetter[_n];
            for (int s = 0; s < _n; s++)
            {
                if (l[s] == 'I') continue;
                var (t, g) = Turn(TurnLetter(s), l[s]);
                sign *= g;
                turned[s] = PauliLetterExtensions.FromSymbol(t);
            }
            yield return (PauliMask.FromLetters(turned), sign * c);
        }
    }

    // ---- the theorems' rules ----

    private (string Class, bool? Predicted) Prediction()
    {
        var axes = _deph.Where(c => c != '.').Distinct().ToList();
        var fields = _field.Where(c => c != '.').Distinct().ToList();
        string bonds = PalindromeStringSpanWitness.BondName(_bonds);
        if (Undephased == 0)
        {
            if (_bonds == "XYZ" && axes.Count > 1)
                return ($"Theorem 3 ({axes.Count} axes)",
                    PalindromeComplementConnectionClaim.MixedAxesPalindrome(axes, fields));
            if (axes.Count == 1 && (_bonds == "XYZ" || axes[0] == 'Z' && _bonds == "XY"))
            {
                char a = axes[0];
                var others = "XYZ".Where(c => c != a).ToList();
                return ($"Theorem 1 (one common axis {a}, bonds {bonds})",
                    !fields.Contains(a) && !(fields.Contains(others[0]) && fields.Contains(others[1])));
            }
            return ($"Theorem 2 alone (bonds {bonds}, no closed rule)", null);
        }
        if (GraphReads)
            return ($"Theorem 4 ({Undephased} undephased site{(Undephased == 1 ? "" : "s")} keeping " +
                    $"{string.Join(", ", Enumerable.Range(0, _n).Where(s => _deph[s] == '.').Select(s => $"{_kept[s]} at {s}"))})",
                null);
        if (_bonds == "XYZ" && Undephased == 1)
            return ("Theorem 5 (one undephased site, Heisenberg bonds)",
                PalindromeComplementConnectionClaim.OneUndephasedSitePalindrome(axes, fields));
        return ($"outside the theorems ({Undephased} undephased sites, bonds {bonds}, a site keeping no letter)", null);
    }

    // ---- the children ----

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            if (!GraphReads)
                yield return new InspectableNode("why no graph reads this row",
                    summary: "An undephased site here keeps no letter: the terms touching it use two letters or more " +
                             "there (" + (_bonds == "Z" ? "its field leaves the bond's Z" : "the bond itself carries two or more") +
                             "), so H does not split into sectors of its bit, and its blocks over the dephased bits are " +
                             "matrices that do not commute. The string route decides the row." +
                             (r.Predicted is null
                                 ? " No theorem gives a rule here (the matrix-valued case, open)."
                                 : " Theorem 5 reaches it by other means, and its rule is printed beside the reading."));
            else if (r.Sectors == 1)
            {
                yield return new InspectableNode("the hopping graph, and its union with the complement image",
                    summary: string.Format(CultureInfo.InvariantCulture,
                        "Each jump turned to Z by a proper rotation of its site's letters; H read on the 2^N = {0} " +
                        "bitstrings exactly over the Gaussian integers. The hopping graph has {1} components (the " +
                        "near count, dim ker L), its union with the complement image {2}.",
                        1 << _n, r.Near, r.UnionComponents));
                yield return new InspectableNode("the flat sections (Theorem 2)",
                    summary: string.Format(CultureInfo.InvariantCulture,
                        "{0} of the {1} union components are good: equal moduli on every edge, H_xx = H_x̄x̄, and a " +
                        "section d_y/d_x = H_xy/H_x̄ȳ that closes on every edge, checked by cross-multiplication with " +
                        "no division. That is the far count, dim ker(L + 2σ), and the verdict: {2}.",
                        r.Far, r.UnionComponents, r.Palindrome ? "every component good, the palindrome holds" : "broken"));
            }
            else
                yield return new InspectableNode("the sectors of the kept letters (Theorem 4)",
                    summary: string.Format(CultureInfo.InvariantCulture,
                        "Each jump and each kept letter turned to Z; H moves no undephased bit (checked on every " +
                        "entry), so it splits into {0} sectors H_σ on the 2^{1} bitstrings of the dephased sites. " +
                        "Over the {2} ordered pairs of sectors the near count is Σ hom(H_τ, H_σ) = {3} and the far " +
                        "count Σ hom(H_τ, H̄_σ) = {4}, each hom the good components of the pair graph, walked like the " +
                        "flat sections. The pairs of distinct sectors give {5} near and {6} far: what the undephased " +
                        "sites add. The verdict: {7}.",
                        r.Sectors, _n - Undephased, r.Sectors * r.Sectors, r.Near, r.Far, r.CrossNear, r.CrossFar,
                        r.Palindrome ? "the two sums agree, the palindrome holds" : "the two sums differ, broken"));

            if (r.Predicted is { } predicted)
                yield return new InspectableNode("the theorem's rule, beside the reading",
                    summary: string.Format(CultureInfo.InvariantCulture,
                        "{0}: the rule predicts {1}, the {2} reads {3}{4}. {5}{6}",
                        r.TheoremClass, predicted ? "a palindrome" : "a break", GraphReads ? "graph" : "string route",
                        r.Palindrome ? "a palindrome" : "a break", r.Exact ? "" : " (a rank reading, not exact)",
                        predicted == r.Palindrome ? "They meet." : "THEY DO NOT MEET, READ IT.",
                        Undephased == 0 && _bonds == "XYZ" && _deph.Distinct().Count() > 1
                            ? r.Near == 1
                                ? " Under mixed axes Lemma A puts the near count at 1, and it is 1."
                                : $" Under mixed axes Lemma A puts the near count at 1; HERE IT IS {r.Near}, READ IT."
                            : ""));

            yield return new InspectableNode("the string route, beside this one",
                summary: MeetTheStringRoute(r));

            yield return new InspectableNode("a break beside the canonical row",
                summary: ComparisonRow());

            yield return new InspectableNode("what this witness does NOT decide",
                summary: "Undephased sites that keep no letter are read by the string route, not by a graph; with " +
                         "two or more of them under Heisenberg bonds no theorem gives a rule (F138 (b)'s SWAP rows pair " +
                         "without a colouring). Nor jumps that are not single-site letters, nor Hamiltonians other " +
                         "than this family's bonds (weight 10) and letter fields (weight 3): Theorems 2 and 4 hold for " +
                         "any real combination of Pauli strings, but this witness builds only that family. The guard " +
                         "on N is the cost of 2^N vertices and carries no physics.");
        }
    }

    /// <summary>The string route's near span has 2^(N + undephased) strings and is refused past
    /// 2^<see cref="PalindromeStringSpanWitness.MaxSpanBits"/>.</summary>
    private bool StringRouteFits => _n + Undephased <= PalindromeStringSpanWitness.MaxSpanBits;

    private string MeetTheStringRoute(Reading r)
    {
        if (!StringRouteFits)
            return string.Format(CultureInfo.InvariantCulture,
                "Not run: the string route's span has 2^{0} strings, past its guard of 2^{1}; the graph reads this " +
                "row alone.", _n + Undephased, PalindromeStringSpanWitness.MaxSpanBits);
        var s = Strings();
        if (r.Near is not { } near || r.Far is not { } far)
            return string.Format(CultureInfo.InvariantCulture,
                "twoendstrings counts the two ends as commutators on spans of Pauli strings and reads near {0}, far {1}, " +
                "{2}, {3}; this row's verdict is that one.{4}",
                s.NearUpper, s.FarUpper, s.Verdict,
                PalindromeStringSpanWitness.IsExact(s.Verdict) ? "an exact reading" : "a rank reading",
                s.Colourings.Count > 0 ? $" Its colourings: {string.Join(", ", s.Colourings)}." : " No colouring.");
        bool meet = s.NearUpper == near && s.FarUpper == far
                    && PalindromeStringSpanWitness.IsPalindrome(s.Verdict) == r.Palindrome;
        return string.Format(CultureInfo.InvariantCulture,
            "twoendstrings counts the two ends as commutators on spans of Pauli strings and reads near {0}, far {1}, " +
            "{2}; the graph reads {3} and {4}. {5} The two share no construction: one ranks commutators on strings " +
            "and never forms a basis state, the other walks basis states and never forms a commutator.",
            s.NearUpper, s.FarUpper, s.Verdict, near, far,
            meet ? "They meet." : "THEY DO NOT MEET, READ IT.");
    }

    private string ComparisonRow()
    {
        var canonical = new ComplementConnectionWitness(_n, new string('Z', _n), null, _topology).Read();
        string mixed = _n >= 3 ? "XYZ" + new string('Z', _n - 3) : "XY";
        var threeAxes = new ComplementConnectionWitness(_n, mixed, null, _topology).Read();
        return string.Format(CultureInfo.InvariantCulture,
            "Heisenberg bonds, Z on every site, no field: hopping {0}, good {1} (the popcount shells, N + 1 = {2}). " +
            "Beside it, dephasing {3} with no field: hopping {4}, good {5}, {6}. {7}",
            canonical.Near, canonical.Far, PalindromeTwoEndCountClaim.CanonicalChainCount(_n),
            mixed, threeAxes.Near, threeAxes.Far, threeAxes.Palindrome ? "a palindrome" : "broken",
            (_n >= 3, threeAxes.Palindrome, threeAxes.Near) switch
            {
                (true, false, 1) => "Three axes leave no letter to colour with, and the far end is empty (Theorem 3).",
                (false, true, 1) => "Two axes and no field: the third letter Z colours both sites (Theorem 3).",
                _ => "THIS IS NOT WHAT THEOREM 3 SAYS, READ IT.",
            });
    }
}
