using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Theorem 8 of the complement connection: the rule on H's block over the states with at most one excitation agrees
/// with Theorem 2 walked on the bitstrings on every row the witness reads, and the three controls (the fields left out of K,
/// the D parts left out, the fields' entries conjugated) disagree, so the battery can fail.</summary>
public class PlanarHoppingWitnessTests
{
    [Fact]
    public void EveryFamilyAgreesWithTheoremTwo_AndTheControlsDisagree()
    {
        var w = new PlanarHoppingWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.ControlsBreak, $"controls {w.FieldsIgnoredMismatches}, {w.TwistsIgnoredMismatches}, {w.ConjugatedFieldsMismatches}");
        Assert.Equal(new[] { (1500, 610), (300, 225), (200, 109), (48, 16), (7, 3) }, w.Families.Select(f => (f.Rows, f.Pairing)).ToArray());   // measured
        Assert.Equal((354, 525, 137), (w.FieldsIgnoredMismatches, w.TwistsIgnoredMismatches, w.ConjugatedFieldsMismatches));   // measured, quoted on the page
    }

    // The corollaries by hand, the rule read far past the bitstrings' reach.
    [Fact]
    public void TheRuleReadsThePagesCorollaries()
    {
        static (int, int, long, long)[] Ring(int n, long j, long d) => Enumerable.Range(0, n).Select(k => (k, (k + 1) % n, j, d)).ToArray();
        var none = System.Array.Empty<(int, long, long, long)>();
        for (int n = 3; n <= 40; n++)
        {
            Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, Ring(n, 1, 0), none));                  // Theorem 1
            Assert.Equal(n % 2 == 0, PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, Ring(n, 0, 1), none));    // pure DM: bipartite
            Assert.Equal(n % 4 == 0, PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, Ring(n, 1, 1), none));    // arg t = −π/4
            Assert.False(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(n, Ring(n, 3, 4), none));                // arg t/π irrational
        }
        // a chain carries no cycle: XX + YY bonds with DM parts (D = 3 here) pair without fields
        var chain = Enumerable.Range(0, 5).Select(k => (k, k + 1, 2L, 3L)).ToArray();
        Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(6, chain, none));
        // a DM bond with X fields on both sites breaks (P = 1·(−2i)·1, not real); X on one and Y on the other pairs (P = 1·(−2i)·i)
        var dm = new[] { (0, 1, 0L, 1L) };
        Assert.False(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(2, dm, new[] { (0, 1L, 0L, 0L), (1, 1L, 0L, 0L) }));
        Assert.False(PlanarHoppingWitness.PairsByTheGraph(2, dm, new[] { (0, 1L, 0L, 0L), (1, 1L, 0L, 0L) }, System.Array.Empty<(int, int, long)>()));
        Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(2, dm, new[] { (0, 1L, 0L, 0L), (1, 0L, 1L, 0L) }));
        Assert.True(PlanarHoppingWitness.PairsByTheGraph(2, dm, new[] { (0, 1L, 0L, 0L), (1, 0L, 1L, 0L) }, System.Array.Empty<(int, int, long)>()));
        // one field X + Y alone on a site of the triangle pairs: outside Theorem 1's letters, inside Theorem 8
        var triangle = new[] { (0, 1, 1L, 0L), (1, 2, 1L, 0L), (0, 2, 1L, 0L) };
        Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(3, triangle, new[] { (0, 1L, 1L, 0L) }));
        Assert.True(PlanarHoppingWitness.PairsByTheGraph(3, triangle, new[] { (0, 1L, 1L, 0L) }, System.Array.Empty<(int, int, long)>()));
        // a Z field breaks; two Z parts on one site that cancel do not
        Assert.False(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(3, triangle, new[] { (1, 0L, 0L, 1L) }));
        Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(3, triangle, new[] { (1, 0L, 0L, 1L), (1, 0L, 0L, -1L) }));
        Assert.False(PlanarHoppingWitness.PairsByTheGraph(3, triangle, new[] { (1, 0L, 0L, 1L) }, System.Array.Empty<(int, int, long)>()));
        Assert.True(PlanarHoppingWitness.PairsByTheGraph(3, triangle, new[] { (1, 0L, 0L, 1L), (1, 0L, 0L, -1L) }, System.Array.Empty<(int, int, long)>()));
        // coefficients past long's half range, the entries formed in BigInteger: P = 2((2^62 + 1)² + 1), real, pairs
        Assert.True(PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(2, new[] { (0, 1, (1L << 62) + 1, 1L) },
            new[] { (0, (1L << 62) + 1, -1L, 0L), (1, 1L, 0L, 0L) }));
        Assert.True(PlanarHoppingWitness.PairsByTheGraph(2, new[] { (0, 1, (1L << 62) + 1, 1L) },
            new[] { (0, (1L << 62) + 1, -1L, 0L), (1, 1L, 0L, 0L) }, System.Array.Empty<(int, int, long)>()));
        Assert.Throws<System.ArgumentException>(() => PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(2, new[] { (0, 5, 1L, 0L) }, new[] { (0, 0L, 0L, 1L) }));
        Assert.Throws<System.ArgumentException>(() => PalindromeComplementConnectionClaim.PlanarHoppingPalindrome(2, new[] { (0, 1, 1L, 0L) }, new[] { (7, 0L, 0L, 1L) }));
    }
}
