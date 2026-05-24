using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using Xunit;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliKBodyChainExtensionsTests
{
    [Fact]
    public void ChainKBody_SingleK3TermAtN3_ReturnsExpectedMatrix()
    {
        // N=3, single XYZ template, k=3 fills the entire chain with exactly one
        // sliding position (l=0). Result must equal PauliString.Build([X, Y, Z]).
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One),
        };
        var actual = templates.ChainKBody(N: 3);
        var expected = PauliString.Build(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z });
        Assert.True((actual - expected).FrobeniusNorm() < 1e-12,
            $"single-position k=3 chain build does not match X⊗Y⊗Z; norm error = {(actual - expected).FrobeniusNorm()}");
    }

    [Fact]
    public void ChainKBody_SingleK3TermAtN4_SlidingWindowGivesTwoPositions()
    {
        // N=4, single XYZ template, k=3 has 2 sliding positions: l=0 → XYZI; l=1 → IXYZ.
        // Result must equal sum of those two tensor products.
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One),
        };
        var actual = templates.ChainKBody(N: 4);
        var xyzi = PauliString.Build(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z, PauliLetter.I });
        var ixyz = PauliString.Build(new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z });
        var expected = xyzi + ixyz;
        Assert.True((actual - expected).FrobeniusNorm() < 1e-12,
            $"sliding-window k=3 chain build does not match XYZI + IXYZ; norm error = {(actual - expected).FrobeniusNorm()}");
    }

    [Fact]
    public void ChainKBody_MultipleTemplates_SumsContributions()
    {
        // N=4, two distinct templates: XYZ and ZYX, each sliding over 2 positions.
        // Result must equal sum of each template's individual sliding window.
        var t1 = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One);
        var t2 = new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Y, PauliLetter.X }, Complex.One);

        var combined = new[] { t1, t2 }.ChainKBody(N: 4);
        var first = new[] { t1 }.ChainKBody(N: 4);
        var second = new[] { t2 }.ChainKBody(N: 4);

        Assert.True((combined - (first + second)).FrobeniusNorm() < 1e-12,
            "ChainKBody on multiple templates must equal the sum of individual template results");
    }

    [Fact]
    public void ChainKBody_ThrowsOnKGreaterThanN()
    {
        // k=4 template on N=3 chain cannot fit; must throw ArgumentException.
        var templates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z, PauliLetter.X }, Complex.One),
        };
        Assert.Throws<ArgumentException>(() => templates.ChainKBody(N: 3));
    }
}
