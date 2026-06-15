using System;
using System.Linq;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F117 Pascal-Gram positivity witness: the live SOS recompute must reproduce the exact
/// integer coefficients the Python verifier (f87_pascal_gram_positivity.py BLOCK 3) pins on all
/// five branch representatives (d = 1, 3, 5), each strictly positive (a sum of squares), and the
/// mod-4 selection rule must single out the firing class for deg ≤ 3.</summary>
public class PascalGramPositivityWitnessTests
{
    [Fact]
    public void SosCoefficient_MatchesExactCrt_AndIsStrictlyPositive_AllFiveCases()
    {
        foreach (var c in PascalGramPositivityWitness.Cases)
        {
            double got = PascalGramPositivityWitness.ComputeCoefficient(c);
            Assert.True(got > 0, $"{c.Name}: SOS coefficient {got} must be > 0 (sum of squares)");
            Assert.True(Math.Abs(got - c.Expected) <= 1e-6 * c.Expected,
                $"{c.Name}: live SOS {got} != exact CRT {c.Expected}");
        }
    }

    [Theory]
    [InlineData("IXXZ+XIXZ", new[] { 1 })]
    [InlineData("K3", new[] { 3 })]
    [InlineData("FLUX", new[] { 3 })]
    [InlineData("MULTIZ", new[] { 3 })]
    [InlineData("IIXY+ZXZY", new[] { 1, 5 })]
    public void SelectionRule_AllowedClasses_ContainTheFiredClass(string name, int[] expectedAllowed)
    {
        var c = PascalGramPositivityWitness.Cases.Single(x => x.Name == name);
        int[] allowed = PascalGramPositivityWitness.AllowedClasses(c.MStar, c.Ell);
        Assert.Equal(expectedAllowed, allowed);
        Assert.Contains(c.D, allowed);
    }

    [Fact]
    public void Witness_SummaryReportsAllFiveReproducing_AndEveryCoefficientPositive()
    {
        var w = new PascalGramPositivityWitness();
        Assert.Contains($"{PascalGramPositivityWitness.Cases.Count}/{PascalGramPositivityWitness.Cases.Count}", w.Summary);
        for (int i = 0; i < PascalGramPositivityWitness.Cases.Count; i++)
            Assert.True(w.Coefficient(i) > 0, $"case {i} coefficient must be positive");
    }
}
