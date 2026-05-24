using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Shared k-body Pauli-pair enumeration: unordered (term1, term2) pairs at
/// given k_body such that both terms are non-trivial (not all-I), share the same
/// Klein index, and share the same Y-parity. Letters drawn from {I, X, Y, Z}.
///
/// <para>Single source of truth replacing previously duplicated copies in
/// <c>F104KBodyTrichotomyVerificationTests</c>, <c>F105KBodyTrichotomyVerificationTestsN5K3</c>,
/// and <c>F87Z2CubedEnumerationN5K3Tool</c>. F106 (k=4) and future F107+ regimes
/// call this same helper with different k values.</para>
///
/// <para>At k=3: 294 pairs. At k=4: 4248 pairs (per F105+ exploration map).
/// Generic k uses base-4 indexing over 4^k letter sequences; outer loop is
/// 4^k × 4^k = 16^k iterations before filters (k=3: 4096, k=4: 65536; both tiny).</para></summary>
public static class Z2HomogeneousKBodyEnumeration
{
    private static readonly PauliLetter[] AllLetters =
        { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    public static List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> Enumerate(int k)
    {
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k), $"k must be positive; got {k}");

        // Materialize all 4^k k-letter sequences via base-4 indexing.
        int total = 1;
        for (int i = 0; i < k; i++) total *= 4;

        var allTerms = new List<PauliLetter[]>(total);
        for (int idx = 0; idx < total; idx++)
        {
            var seq = new PauliLetter[k];
            int x = idx;
            for (int j = 0; j < k; j++)
            {
                seq[j] = AllLetters[x & 3];
                x >>= 2;
            }
            allTerms.Add(seq);
        }

        var seen = new HashSet<string>();
        var items = new List<(PauliTerm, PauliTerm, (int, int), int)>(total);

        foreach (var t1Letters in allTerms)
        {
            if (t1Letters.All(l => l == PauliLetter.I)) continue;
            var term1 = new PauliTerm(t1Letters, Complex.One);

            foreach (var t2Letters in allTerms)
            {
                if (t2Letters.All(l => l == PauliLetter.I)) continue;
                var term2 = new PauliTerm(t2Letters, Complex.One);

                if (term1.KleinIndex != term2.KleinIndex) continue;
                if (term1.YParity != term2.YParity) continue;

                string k1 = LettersKey(t1Letters);
                string k2 = LettersKey(t2Letters);
                string key = string.Compare(k1, k2, StringComparison.Ordinal) <= 0
                    ? k1 + "|" + k2
                    : k2 + "|" + k1;
                if (!seen.Add(key)) continue;

                items.Add((term1, term2, term1.KleinIndex, term1.YParity));
            }
        }

        return items;
    }

    private static string LettersKey(IReadOnlyList<PauliLetter> letters) =>
        string.Concat(letters.Select(LetterChar));

    private static char LetterChar(PauliLetter l) => l switch
    {
        PauliLetter.I => 'I',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        PauliLetter.Z => 'Z',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };
}
