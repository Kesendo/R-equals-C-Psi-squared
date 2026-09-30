using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

/// <summary>PauliMask judged against <see cref="PauliString.Build"/>, the dense route that never sees
/// a mask. The entries are 0, ±1, ±i, so the dense products are exact in double and compared
/// with ==, no tolerance.</summary>
public class PauliMaskTests
{
    static readonly PauliLetter[] Letters = { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    static IEnumerable<PauliLetter[]> All(int n)
    {
        for (int k = 0; k < 1 << (2 * n); k++)
            yield return Enumerable.Range(0, n).Select(l => Letters[(k >> (2 * l)) & 3]).ToArray();
    }

    static Complex IPow(int k) => k switch { 0 => Complex.One, 1 => Complex.ImaginaryOne, 2 => -Complex.One, _ => -Complex.ImaginaryOne };

    static bool Equal(MathNet.Numerics.LinearAlgebra.Matrix<Complex> a, MathNet.Numerics.LinearAlgebra.Matrix<Complex> b)
    {
        for (int i = 0; i < a.RowCount; i++)
            for (int j = 0; j < a.ColumnCount; j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    [Fact]
    public void Product_And_Phase_Match_The_Dense_Matrices_On_Every_Pair_At_Two_Sites()
    {
        var all = All(2).ToList();
        foreach (var a in all)
            foreach (var b in all)
            {
                var (c, k) = PauliMask.Multiply(PauliMask.FromLetters(a), PauliMask.FromLetters(b));
                var dense = PauliString.Build(a) * PauliString.Build(b);
                Assert.True(Equal(dense, PauliString.Build(c.ToLetters(2)) * IPow(k)),
                    $"{PauliMask.FromLetters(a).ToString(2)}·{PauliMask.FromLetters(b).ToString(2)}");
            }
    }

    [Fact]
    public void Commutation_Matches_The_Dense_Commutator_And_Fixes_The_Phase_Parity_At_Three_Sites()
    {
        var all = All(3).ToList();
        var rng = new Random(20260930);
        for (int t = 0; t < 300; t++)
        {
            var a = all[rng.Next(all.Count)];
            var b = all[rng.Next(all.Count)];
            var (ma, mb) = (PauliString.Build(a), PauliString.Build(b));
            bool commute = Equal(ma * mb, mb * ma);
            Assert.Equal(commute, PauliMask.Commute(PauliMask.FromLetters(a), PauliMask.FromLetters(b)));
            Assert.Equal(commute ? 0 : 1, PauliMask.Multiply(PauliMask.FromLetters(a), PauliMask.FromLetters(b)).PhasePower % 2);
        }
    }

    [Fact]
    public void A_Long_Mixed_Product_Is_The_Product_Of_Its_One_Site_Products()
    {
        // 40 sites, past what a dense check reaches: each site's letter and phase from the 2x2
        // matrices, and the total phase their product.
        var rng = new Random(40);
        for (int t = 0; t < 30; t++)
        {
            var a = Enumerable.Range(0, 40).Select(_ => Letters[rng.Next(4)]).ToArray();
            var b = Enumerable.Range(0, 40).Select(_ => Letters[rng.Next(4)]).ToArray();
            Complex phase = 1;
            var letters = new PauliLetter[40];
            for (int l = 0; l < 40; l++)
            {
                var m = PauliString.Build(new[] { a[l] }) * PauliString.Build(new[] { b[l] });
                var hit = Letters.SelectMany(c => Enumerable.Range(0, 4).Select(k => (c, k)))
                    .Single(ck => Equal(m, PauliString.Build(new[] { ck.c }) * IPow(ck.k)));
                letters[l] = hit.c;
                phase *= IPow(hit.k);
            }
            var (prod, kk) = PauliMask.Multiply(PauliMask.FromLetters(a), PauliMask.FromLetters(b));
            Assert.Equal(letters, prod.ToLetters(40));
            Assert.Equal(phase, IPow(kk));
        }
    }

    [Fact]
    public void Parse_And_Print_Are_Inverse_And_Use_PauliLetters_Own_Bits()
    {
        foreach (var s in All(3)) Assert.Equal(s, PauliMask.FromLetters(s).ToLetters(3));
        var y = PauliMask.Parse("IYZX");
        Assert.Equal(0b1010UL, y.X);   // X parts at sites 1 (Y) and 3 (X)
        Assert.Equal(0b0110UL, y.Z);   // Z parts at sites 1 (Y) and 2 (Z)
        Assert.Equal("IYZX", y.ToString(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => y.ToString(3));
        Assert.Throws<ArgumentException>(() => PauliMask.Parse("XQ"));
    }
}
