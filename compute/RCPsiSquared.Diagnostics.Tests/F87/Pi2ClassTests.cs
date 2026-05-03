using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class Pi2ClassTests
{
    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);

    [Fact]
    public void Empty_IsTruly()
    {
        Assert.Equal(Pi2Class.Truly, Array.Empty<PauliPairBondTerm>().Pi2ClassOf());
    }

    [Fact]
    public void XXplusYY_IsTruly()
    {
        var terms = new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) };
        Assert.Equal(Pi2Class.Truly, terms.Pi2ClassOf());
    }

    [Fact]
    public void Heisenberg_XX_YY_ZZ_IsTruly()
    {
        var terms = new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.Y, PauliLetter.Y),
            Term(PauliLetter.Z, PauliLetter.Z),
        };
        Assert.Equal(Pi2Class.Truly, terms.Pi2ClassOf());
    }

    [Fact]
    public void XYplusYX_IsPi2OddPure()
    {
        var terms = new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) };
        Assert.Equal(Pi2Class.Pi2OddPure, terms.Pi2ClassOf());
    }

    [Fact]
    public void XYAlone_IsPi2OddPure()
    {
        var terms = new[] { Term(PauliLetter.X, PauliLetter.Y) };
        Assert.Equal(Pi2Class.Pi2OddPure, terms.Pi2ClassOf());
    }

    [Fact]
    public void YZplusZY_IsPi2EvenNonTruly()
    {
        var terms = new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) };
        Assert.Equal(Pi2Class.Pi2EvenNonTruly, terms.Pi2ClassOf());
    }

    [Fact]
    public void TrulyMixedWithEvenNonTruly_IsPi2EvenNonTruly()
    {
        // XX (truly) + YZ (Π²-even non-truly) — both Π²-even, but YZ breaks truly criterion.
        var terms = new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Z) };
        Assert.Equal(Pi2Class.Pi2EvenNonTruly, terms.Pi2ClassOf());
    }

    [Fact]
    public void XXplusXY_IsMixed()
    {
        // XX (Π²-even, truly) + XY (Π²-odd) → Mixed.
        var terms = new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) };
        Assert.Equal(Pi2Class.Mixed, terms.Pi2ClassOf());
    }

    [Fact]
    public void XYplusYZ_IsMixed()
    {
        // XY (Π²-odd) + YZ (Π²-even non-truly) → Mixed (both classes coexist).
        var terms = new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.Z) };
        Assert.Equal(Pi2Class.Mixed, terms.Pi2ClassOf());
    }

    [Theory]
    [InlineData(Pi2Class.Truly, TrichotomyClass.Truly)]
    [InlineData(Pi2Class.Pi2OddPure, TrichotomyClass.Soft)]
    [InlineData(Pi2Class.Pi2EvenNonTruly, TrichotomyClass.Soft)]
    [InlineData(Pi2Class.Mixed, TrichotomyClass.Hard)]
    public void PredictedTrichotomy_MatchesEmpiricalCanonicalCases(Pi2Class pi2, TrichotomyClass expected)
    {
        Assert.Equal(expected, pi2.PredictedTrichotomy());
    }

    [Theory]
    [MemberData(nameof(CanonicalCrossWalkCases))]
    public void Algebraic_AndSpectral_Agree_OnAllCanonicalWitnesses(
        string label, PauliPairBondTerm[] terms, Pi2Class expectedPi2, TrichotomyClass expectedSpectral)
    {
        // Both lenses must produce consistent classifications on the F87 canonical witnesses.
        // The algebraic test (no L build) gives Pi2Class; the spectral test (full L + Spec)
        // gives TrichotomyClass; the cross-walk maps the first to the second.
        var actualPi2 = ((IReadOnlyList<PauliPairBondTerm>)terms).Pi2ClassOf();
        Assert.Equal(expectedPi2, actualPi2);

        var actualSpectral = PauliPairTrichotomy.Classify(MakeChain(3), terms);
        Assert.Equal(expectedSpectral, actualSpectral);

        Assert.Equal(actualSpectral, actualPi2.PredictedTrichotomy());
    }

    public static IEnumerable<object[]> CanonicalCrossWalkCases() => new[]
    {
        new object[]
        {
            "XX+YY",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) },
            Pi2Class.Truly,
            TrichotomyClass.Truly,
        },
        new object[]
        {
            "XX+YY+ZZ (Heisenberg)",
            new[]
            {
                Term(PauliLetter.X, PauliLetter.X),
                Term(PauliLetter.Y, PauliLetter.Y),
                Term(PauliLetter.Z, PauliLetter.Z),
            },
            Pi2Class.Truly,
            TrichotomyClass.Truly,
        },
        new object[]
        {
            "XY+YX (bond-flip soft)",
            new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) },
            Pi2Class.Pi2OddPure,
            TrichotomyClass.Soft,
        },
        new object[]
        {
            "YZ+ZY (EQ-030 soft)",
            new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) },
            Pi2Class.Pi2EvenNonTruly,
            TrichotomyClass.Soft,
        },
        new object[]
        {
            "XX+XY (mixed hard)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) },
            Pi2Class.Mixed,
            TrichotomyClass.Hard,
        },
    };
}
