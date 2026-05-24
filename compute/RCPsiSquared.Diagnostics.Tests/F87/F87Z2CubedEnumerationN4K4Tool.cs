using System.Text.Json;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F106 one-shot enumeration tool: classify the 4248 Z₂³-homogeneous +
/// Y-par-homogeneous k=4 Pauli pairs at N=4 (new enumeration vs F103/F105's 294
/// at k=3) via F104's <see cref="PauliPairTrichotomy.Classify"/> k-body overload,
/// export a 4×3×2×3 count grid to
/// <c>simulations/results/f87_z2cubed_split_n4_k4_counts.json</c> for human
/// inspection and to feed the F106 Claim's frozen counts.
///
/// <para>Runtime: ~2-3min PLINQ at N=4 (per-call ~0.2s sequential, 12744
/// classifications; PLINQ on 24 cores saturates). Tagged SLOW_F106_BATCH so it
/// skips by default in CI; re-run manually via
/// <c>--filter "Category=SLOW_F106_BATCH"</c>.</para></summary>
[Trait("Category", "SLOW_F106_BATCH")]
public class F87Z2CubedEnumerationN4K4Tool
{
    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    [Fact]
    public void EnumerateClassifyAndExportToJson()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);

        var items = Z2HomogeneousKBodyEnumeration.Enumerate(4);
        Assert.Equal(4248, items.Count);

        var counts = ClassifyAndGroup(chain, items);

        var json = SerializeToJson(counts);

        string outDir = Path.Combine(FindRepoRoot(), "simulations", "results");
        Directory.CreateDirectory(outDir);
        string outPath = Path.Combine(outDir, "f87_z2cubed_split_n4_k4_counts.json");
        File.WriteAllText(outPath, json);

        // Sanity: total classifications across grid = 12744 (4248 pairs × 3 dephase letters).
        int total = counts.Values.Sum();
        Assert.Equal(12744, total);
    }

    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> ClassifyAndGroup(
        ChainSystem chain,
        List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> items)
    {
        // Parallelize over (pair × dephase letter). PauliPairTrichotomy.Classify is pure;
        // ChainSystem and PauliTerm are immutable records. 12744 classifications saturate
        // 24 cores via PLINQ.
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
            n = 4,
            k = 4,
            total_pairs = 4248,
            total_classifications = 12744,
            grid,
        }, new JsonSerializerOptions { WriteIndented = true });
    }

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
