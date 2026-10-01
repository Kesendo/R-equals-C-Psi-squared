using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The complement connection carried into the certifier (<see cref="PalindromeSoftCertifier.DecideAtN"/>,
/// Theorem 2 of <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>). Judged against the spectral
/// authority <see cref="PauliPairTrichotomy"/>, which builds L and pairs its eigenvalues, and against the
/// certifier's own N-free soft certificates, which must hold at every N.</summary>
public class PalindromeDecideAtNTests
{
    private const PauliLetter I = PauliLetter.I, X = PauliLetter.X, Y = PauliLetter.Y, Z = PauliLetter.Z;
    private static readonly PauliLetter[] Letters = { I, X, Y, Z };
    private static PauliTerm T(params PauliLetter[] ls) => new(ls, Complex.One);

    private static IEnumerable<PauliTerm> TwoBody() =>
        from a in Letters from b in Letters
        where a != I && b != I
        select T(a, b);

    // Every pair of two-body templates (and every single one) at N = 4: the graph's verdict is the
    // authority's (palindromic = truly or soft, broken = hard) on every set, both verdicts occurring.
    [Fact]
    public void The_Graph_Meets_The_Spectral_Authority_On_Every_Two_Body_Pair_At_N4()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var templates = TwoBody().ToList();
        int sets = 0, broken = 0, palindromic = 0;
        for (int i = 0; i < templates.Count; i++)
            for (int j = i; j < templates.Count; j++)
            {
                var terms = i == j ? new[] { templates[i] } : new[] { templates[i], templates[j] };
                var (pal, _, _) = PalindromeSoftCertifier.ComplementConnectionAtN(terms, 4);
                var authority = PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: Z);
                Assert.True(pal == (authority != TrichotomyClass.Hard),
                    $"{string.Join("+", terms.Select(t => t.Label))}: graph {pal}, authority {authority}");
                if (pal) palindromic++; else broken++;
                sets++;
            }
        Assert.Equal(9 * 10 / 2, sets);
        Assert.True(broken > 5 && palindromic > 5, $"broken {broken}, palindromic {palindromic}");
    }

    // Three-body pairs at N = 4, the k-body path of the authority: every third pair of the diagonal-cell
    // Mixed templates the hard sweep uses.
    [Fact]
    public void The_Graph_Meets_The_Spectral_Authority_On_Three_Body_Pairs_At_N4()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var cell = PalindromeHardSweepTests.DiagonalCellMixedTerms(3).ToList();
        int sets = 0, broken = 0;
        for (int i = 0; i < cell.Count; i++)
            for (int j = i + 1; j < cell.Count; j += 3)
            {
                var terms = new[] { cell[i], cell[j] };
                var (pal, _, _) = PalindromeSoftCertifier.ComplementConnectionAtN(terms, 4);
                var authority = PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: Z);
                Assert.True(pal == (authority != TrichotomyClass.Hard),
                    $"{cell[i].Label}+{cell[j].Label}: graph {pal}, authority {authority}");
                if (!pal) broken++;
                sets++;
            }
        Assert.True(sets > 20 && broken > 0, $"sets {sets}, broken {broken}");
    }

    // Decide's N-free verdicts stand, and where it is Undetermined the graph decides at this N.
    [Fact]
    public void DecideAtN_Keeps_Decide_And_Fills_Its_Undetermined_Rows()
    {
        var frustrated = new[] { T(X, X, X), T(X, X, Y), T(Y, X, X) };
        Assert.Equal(PalindromeSoftCertifier.Decision.Undetermined, PalindromeSoftCertifier.Decide(frustrated, 4).Verdict);
        var d = PalindromeSoftCertifier.DecideAtN(frustrated, 4);
        var authority = PauliPairTrichotomy.Classify(new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05), frustrated, dephaseLetter: Z);
        Assert.Equal(authority == TrichotomyClass.Hard ? PalindromeSoftCertifier.Decision.Hard : PalindromeSoftCertifier.Decision.Soft,
            d.Verdict);
        Assert.True(d.SoftStrategy == PalindromeSoftCertifier.SoftStrategy.ComplementConnection
                    || d.HardStrategy == PalindromeSoftCertifier.HardStrategy.ComplementConnection, d.Reason);

        var xy = new[] { T(X, X), T(Y, Y) };
        Assert.Equal(PalindromeSoftCertifier.Decide(xy, 4), PalindromeSoftCertifier.DecideAtN(xy, 4));
        var hardPair = new[] { T(X, X, Z), T(X, Z, X) };                    // Decide's N-free Hard passes through
        Assert.Equal(PalindromeSoftCertifier.HardStrategy.DiagonalCellValuation,
            PalindromeSoftCertifier.DecideAtN(hardPair, 4).HardStrategy);
    }

    // An N-free soft certificate claims the palindrome at every N; the graph must agree wherever one is given.
    [Fact]
    public void Every_N_Free_Soft_Certificate_Holds_On_The_Graph()
    {
        int certified = 0;
        var templates = TwoBody().ToList();
        foreach (int n in new[] { 3, 5 })
            for (int i = 0; i < templates.Count; i++)
                for (int j = i + 1; j < templates.Count; j++)
                {
                    var terms = new[] { templates[i], templates[j] };
                    if (!PalindromeSoftCertifier.Certify(terms, n).Certified) continue;
                    certified++;
                    Assert.True(PalindromeSoftCertifier.ComplementConnectionAtN(terms, n).Palindrome,
                        $"{templates[i].Label}+{templates[j].Label} at N = {n}: certified soft, broken on the graph");
                }
        Assert.True(certified > 10, $"certified {certified}");
    }

    // Past the Liouvillian: the XY chain at N = 14 under Z has the popcount shells at both ends, N + 1 each.
    // A template slides over every window, so a one-site Z template is a uniform Z field: its diagonal
    // h·(N − 2|x|) equals its complement's only on the half-filled shell, the one good component left at N = 6.
    [Fact]
    public void The_XY_Chain_At_N14_Has_The_Popcount_Shells_At_Both_Ends()
    {
        Assert.Equal((true, 15, 15), PalindromeSoftCertifier.ComplementConnectionAtN(new[] { T(X, X), T(Y, Y) }, 14));
        var withField = new[] { T(X, X), T(Y, Y), new PauliTerm(new[] { Z }, new Complex(0.3, 0)) };
        Assert.Equal((false, 7, 1), PalindromeSoftCertifier.ComplementConnectionAtN(withField, 6));
    }

    // The coefficients are scaled to integers by one power of two, exactly: a decimal like 0.1 keeps every bit.
    [Fact]
    public void Dyadic_Scaling_Is_Exact()
    {
        Assert.Equal(new BigInteger[] { 2, 1, -12, 0 }, ComplementConnectionGraph.DyadicIntegers(new[] { 0.5, 0.25, -3.0, 0.0 }));
        var odd = ComplementConnectionGraph.DyadicIntegers(new[] { -0.1, double.Epsilon, 1.0 });   // a negative shift, a subnormal
        Assert.Equal(BigInteger.One << 1074, odd[2]);
        Assert.Equal(BigInteger.One, odd[1]);
        Assert.Equal(-(new BigInteger(3602879701896397) << (1074 - 55)), odd[0]);
        var tenth = ComplementConnectionGraph.DyadicIntegers(new[] { 0.1, 1.0 });
        Assert.Equal(BigInteger.One << 55, tenth[1]);                   // 0.1 = m·2^−55 with m odd, the smallest scale
        Assert.False(tenth[0].IsEven);
        Assert.True((double)tenth[0] / (double)tenth[1] == 0.1);        // both exactly representable, the quotient exact
        Assert.Throws<ArgumentException>(() => PalindromeSoftCertifier.ComplementConnectionAtN(
            new[] { new PauliTerm(new[] { X, Y }, new Complex(0, 1)) }, 3));
    }
}
