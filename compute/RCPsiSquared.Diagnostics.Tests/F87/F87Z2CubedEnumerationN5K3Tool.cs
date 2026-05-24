using System.Numerics;
using System.Text.Json;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F105 one-shot enumeration tool: re-classify the 294 Z₂³-homogeneous +
/// Y-par-homogeneous k=3 Pauli pairs at N=5 (same enumeration as F103, larger chain)
/// via F104's <see cref="PauliPairTrichotomy.Classify"/> k-body overload, and export
/// a 4×3×2×3 count grid to <c>simulations/results/f87_z2cubed_split_n5_k3_counts.json</c>
/// for human inspection and to feed the F105 Claim's frozen counts.
///
/// <para>Runtime: ~3h batch at N=5 (per-call ~12s, 882 classifications). Tagged
/// SLOW_F105_BATCH so it skips by default in CI; re-run manually via
/// <c>--filter "Category=SLOW_F105_BATCH"</c>.</para></summary>
[Trait("Category", "SLOW_F105_BATCH")]
public class F87Z2CubedEnumerationN5K3Tool
{
    private static readonly PauliLetter[] AllLetters =
        { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    [Fact]
    public void EnumerateClassifyAndExportToJson()
    {
        var chain = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05);

        var items = EnumerateZ2HomogeneousK3();
        Assert.Equal(294, items.Count);

        var counts = ClassifyAndGroup(chain, items);

        var json = SerializeToJson(counts);

        string outDir = Path.Combine(
            FindRepoRoot(), "simulations", "results");
        Directory.CreateDirectory(outDir);
        string outPath = Path.Combine(outDir, "f87_z2cubed_split_n5_k3_counts.json");
        File.WriteAllText(outPath, json);

        // Sanity: total classifications across grid = 882 (294 pairs × 3 dephase letters).
        int total = counts.Values.Sum();
        Assert.Equal(882, total);
    }

    private static List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> EnumerateZ2HomogeneousK3()
    {
        var seen = new HashSet<string>();
        var items = new List<(PauliTerm, PauliTerm, (int, int), int)>();
        foreach (var a in AllLetters)
            foreach (var b in AllLetters)
                foreach (var c in AllLetters)
                {
                    var t1Letters = new[] { a, b, c };
                    if (t1Letters.All(l => l == PauliLetter.I)) continue;
                    var term1 = new PauliTerm(t1Letters, Complex.One);

                    foreach (var d in AllLetters)
                        foreach (var e in AllLetters)
                            foreach (var f in AllLetters)
                            {
                                var t2Letters = new[] { d, e, f };
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

    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> ClassifyAndGroup(
        ChainSystem chain,
        List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> items)
    {
        // Parallelize over (pair × dephase letter). PauliPairTrichotomy.Classify is pure
        // (constructs fresh H/L/M per call); ChainSystem and PauliTerm are immutable records.
        // 882 independent classifications saturate all available cores via PLINQ.
        var classifications = items
            .SelectMany(item => DephaseLetters.Select(d => (item, dephase: d)))
            .AsParallel()
            .Select(x =>
            {
                var templates = new[] { x.item.Term1, x.item.Term2 };
                var cls = PauliPairTrichotomy.Classify(chain, templates, dephaseLetter: x.dephase);
                return (x.item.Klein, Dephase: DephaseChar(x.dephase), x.item.YPar, Cls: cls);
            })
            .ToList();

        return classifications
            .GroupBy(c => (c.Klein, c.Dephase, c.YPar, c.Cls))
            .ToDictionary(g => g.Key, g => g.Count());
    }

    private static string SerializeToJson(Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> counts)
    {
        var grid = new List<object>();
        foreach (var kvp in counts.OrderBy(k => k.Key.Klein.A).ThenBy(k => k.Key.Klein.B)
                                  .ThenBy(k => k.Key.Dephase).ThenBy(k => k.Key.YPar).ThenBy(k => k.Key.Cls.ToString()))
        {
            grid.Add(new
            {
                klein_a = kvp.Key.Klein.A,
                klein_b = kvp.Key.Klein.B,
                dephase = kvp.Key.Dephase.ToString(),
                y_par = kvp.Key.YPar,
                trichotomy = kvp.Key.Cls.ToString().ToLowerInvariant(),
                count = kvp.Value,
            });
        }
        return JsonSerializer.Serialize(new
        {
            n = 5,
            k = 3,
            total_pairs = 294,
            total_classifications = 882,
            grid,
        }, new JsonSerializerOptions { WriteIndented = true });
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

    private static char DephaseChar(PauliLetter l) => l switch
    {
        PauliLetter.Z => 'Z',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    private static string FindRepoRoot()
    {
        var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
        while (dir != null && !File.Exists(Path.Combine(dir.FullName, "README.md")) &&
                                !Directory.Exists(Path.Combine(dir.FullName, ".git")))
            dir = dir.Parent;
        if (dir == null)
            throw new InvalidOperationException("Could not find repo root from " + Directory.GetCurrentDirectory());
        return dir.FullName;
    }
}
