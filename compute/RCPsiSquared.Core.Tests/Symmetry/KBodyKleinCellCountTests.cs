using RCPsiSquared.Core.Pauli;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The k-body Klein cells counted by enumeration against the letter-character closed form that
/// <c>Pi2OpenQuestions</c> states: each letter carries (bit_a, bit_b), I (0,0), X (1,0), Y (1,1), Z (0,1); with
/// identities every cell holds 4^(k−1) strings (the all-identity string counted in (0,0)), without identity letters the (0,0) cell holds (3^k + 3(−1)^k)/4 and each
/// other cell (3^k − (−1)^k)/4. The cell of a string is read through
/// <see cref="PauliLetter"/>'s own BitA/BitB.</summary>
public class KBodyKleinCellCountTests
{
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void KleinCellCounts_FollowTheLetterCharacters(int k)
    {
        var withI = new long[2, 2];
        var full = new long[2, 2];
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        for (long code = 0; code < 1L << (2 * k); code++)
        {
            int a = 0, b = 0; bool hasI = false;
            for (int i = 0; i < k; i++)
            {
                var l = letters[(int)((code >> (2 * i)) & 3)];
                a ^= l.BitA(); b ^= l.BitB();
                if (l == PauliLetter.I) hasI = true;
            }
            withI[a, b]++;
            if (!hasI) full[a, b]++;
        }
        long pow3 = (long)Math.Pow(3, k), sign = k % 2 == 0 ? 1 : -1;
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
            {
                Assert.Equal(1L << (2 * (k - 1)), withI[a, b]);
                long expected = (a, b) == (0, 0) ? (pow3 + 3 * sign) / 4 : (pow3 - sign) / 4;
                Assert.Equal(expected, full[a, b]);
            }
    }

    [Fact]
    public void TwoBody_FullSupport_IsTheTableOfNineBilinears()
    {
        // XX, YY, ZZ in (0,0); two strings in each of the other three cells
        var counts = new Dictionary<(int, int), int>();
        foreach (var p in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
            foreach (var q in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
            {
                var cell = (p.BitA() ^ q.BitA(), p.BitB() ^ q.BitB());
                counts[cell] = counts.GetValueOrDefault(cell) + 1;
            }
        Assert.Equal(3, counts[(0, 0)]);
        Assert.Equal(2, counts[(0, 1)]);
        Assert.Equal(2, counts[(1, 0)]);
        Assert.Equal(2, counts[(1, 1)]);
    }
}
