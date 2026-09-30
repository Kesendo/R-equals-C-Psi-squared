using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for PauliString, the exact string algebra EndCount is written in. Every gate judges
// the masks against dense 2^N matrices built letter by letter from the four 2x2 Pauli matrices, a route
// that never sees a mask, a popcount or the phase formula. The entries are 0, +-1, +-i, so the dense
// products are exact in double and compared with ==, no tolerance.
public class PauliStringTests
{
    static readonly Complex I1 = Complex.ImaginaryOne;

    static Complex[,] Letter(char c) => c switch
    {
        'I' => new Complex[,] { { 1, 0 }, { 0, 1 } },
        'X' => new Complex[,] { { 0, 1 }, { 1, 0 } },
        'Y' => new Complex[,] { { 0, -I1 }, { I1, 0 } },
        'Z' => new Complex[,] { { 1, 0 }, { 0, -1 } },
        _ => throw new ArgumentException(c.ToString()),
    };

    // site 0 is the most significant tensor factor
    internal static Complex[,] Dense(string letters)
    {
        Complex[,] m = { { 1 } };
        foreach (char c in letters) m = Kron(m, Letter(c));
        return m;
    }

    static Complex[,] Kron(Complex[,] a, Complex[,] b)
    {
        int ra = a.GetLength(0), rb = b.GetLength(0);
        var o = new Complex[ra * rb, ra * rb];
        for (int i = 0; i < ra; i++) for (int j = 0; j < ra; j++)
            for (int k = 0; k < rb; k++) for (int l = 0; l < rb; l++)
                o[i * rb + k, j * rb + l] = a[i, j] * b[k, l];
        return o;
    }

    internal static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        int n = a.GetLength(0);
        var o = new Complex[n, n];
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++)
        {
            if (a[i, k] == Complex.Zero) continue;
            for (int j = 0; j < n; j++) o[i, j] += a[i, k] * b[k, j];
        }
        return o;
    }

    static bool Equal(Complex[,] a, Complex[,] b)
    {
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (a[i, j] != b[i, j]) return false;
        return true;
    }

    static Complex[,] Scale(Complex[,] a, Complex s)
    {
        int n = a.GetLength(0);
        var o = new Complex[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) o[i, j] = a[i, j] * s;
        return o;
    }

    static IEnumerable<string> AllStrings(int n)
    {
        const string L = "IXYZ";
        for (int k = 0; k < (1 << (2 * n)); k++)
        {
            var c = new char[n];
            for (int l = 0; l < n; l++) c[l] = L[(k >> (2 * l)) & 3];
            yield return new string(c);
        }
    }

    static Complex IPow(int k) => k switch { 0 => 1, 1 => I1, 2 => -1, _ => -I1 };

    [Fact]
    public void TheProductAndItsPhase_MatchTheDenseMatrices_OnEveryPairAtTwoSites()
    {
        var all = AllStrings(2).ToList();
        foreach (var a in all)
            foreach (var b in all)
            {
                var (prod, k) = PauliString.Multiply(PauliString.Parse(a), PauliString.Parse(b));
                Assert.True(Equal(Mul(Dense(a), Dense(b)), Scale(Dense(prod.ToString(2)), IPow(k))), $"{a}·{b}");
            }
    }

    [Fact]
    public void CommutationMatchesTheDenseCommutator_AndFixesThePhaseParity_AtThreeSites()
    {
        var all = AllStrings(3).ToList();
        var rng = new Random(20260930);
        for (int t = 0; t < 400; t++)
        {
            string a = all[rng.Next(all.Count)], b = all[rng.Next(all.Count)];
            var (pa, pb) = (PauliString.Parse(a), PauliString.Parse(b));
            bool denseCommute = Equal(Mul(Dense(a), Dense(b)), Mul(Dense(b), Dense(a)));
            Assert.Equal(denseCommute, PauliString.Commute(pa, pb));
            var (_, k) = PauliString.Multiply(pa, pb);
            Assert.Equal(denseCommute ? 0 : 1, k % 2);      // commuting: real phase; anticommuting: +-i
        }
    }

    [Fact]
    public void EveryStringSquaresToTheIdentity_WithPhasePlusOne()
    {
        foreach (var s in AllStrings(3))
        {
            var (prod, k) = PauliString.Multiply(PauliString.Parse(s), PauliString.Parse(s));
            Assert.True(prod.IsIdentity);
            Assert.Equal(0, k);
        }
    }

    [Fact]
    public void ParseAndPrint_AreInverse_AndRefuseForeignLetters()
    {
        foreach (var s in AllStrings(3)) Assert.Equal(s, PauliString.Parse(s).ToString(3));
        Assert.Throws<ArgumentException>(() => PauliString.Parse("XQZ"));
        Assert.Throws<ArgumentException>(() => PauliString.Parse(""));
        Assert.Throws<ArgumentOutOfRangeException>(() => PauliString.Parse("IIX").ToString(2));
    }

    // Mixed letters at 32 sites against the product taken site by site with the 2x2 matrices: the
    // phase is the product of the 32 one-site phases, a check a wrong popcount formula cannot pass by
    // landing on 0 mod 4.
    [Fact]
    public void A_Mixed_ThirtyTwoSite_Product_Is_The_Product_Of_Its_Sites()
    {
        const string L = "IXYZ";
        var rng = new Random(32);
        for (int t = 0; t < 50; t++)
        {
            var a = new string(Enumerable.Range(0, 32).Select(_ => L[rng.Next(4)]).ToArray());
            var b = new string(Enumerable.Range(0, 32).Select(_ => L[rng.Next(4)]).ToArray());
            Complex phase = 1;
            var letters = new char[32];
            for (int l = 0; l < 32; l++)
            {
                var m = Mul(Letter(a[l]), Letter(b[l]));
                var found = L.SelectMany(c => new[] { 0, 1, 2, 3 }.Select(k => (c, k)))
                    .Single(ck => Equal(m, Scale(Letter(ck.c), IPow(ck.k))));
                letters[l] = found.c;
                phase *= IPow(found.k);
            }
            var (prod, kk) = PauliString.Multiply(PauliString.Parse(a), PauliString.Parse(b));
            Assert.Equal(new string(letters), prod.ToString(32));
            Assert.Equal(phase, IPow(kk));
        }
    }

    [Fact]
    public void ThirtyTwoSites_AreTheBound()
    {
        var s = PauliString.Parse(new string('Y', 32));
        var (prod, k) = PauliString.Multiply(s, PauliString.Parse(new string('X', 32)));
        Assert.Equal(new string('Z', 32), prod.ToString(32));
        Assert.Equal(0, k);                                   // (-i)^32 = 1: YX = -iZ at every site
        Assert.Throws<ArgumentException>(() => PauliString.Parse(new string('X', 33)));
        Assert.Throws<ArgumentOutOfRangeException>(() => s.ToString(33));
        Assert.Throws<ArgumentOutOfRangeException>(() => s.Letter(32));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PauliString(0, 1UL << 32));
    }
}
