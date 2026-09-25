using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CpsiEnvelopeTheoremClaimTests
{
    [Fact]
    public void Claim_TierIsOpenQuestion()
    {
        // The statement is the surviving main-peak question; the refutations are computing members below.
        Assert.Equal(Tier.OpenQuestion, CpsiEnvelopeTheoremClaim.Shared.Tier);
    }

    [Fact]
    public void LocalField_FirstTwoMaximaRise_OnTheGammaZeroClosedForm()
    {
        // Golden-section search on the closed form, one hump on each side of the kink at t = π/2.
        var (first, second) = CpsiEnvelopeTheoremClaim.LocalFieldFirstTwoMaxima();

        // The refutation: a later local maximum above the earlier one.
        Assert.True(second.Cpsi > first.Cpsi,
            $"expected a rising pair, got {first.Cpsi:R} (t = {first.T:R}) then {second.Cpsi:R} (t = {second.T:R})");

        // The six-decimal readings the proof page carries come from exact expm propagation in
        // simulations/envelope_n2_rises.py, a different route; they are the exact roundings of these values.
        Assert.Equal(0.734893, Math.Round(first.Cpsi, 6));
        Assert.Equal(0.991268, Math.Round(second.Cpsi, 6));
        Assert.Equal(0.9545, Math.Round(first.T, 4));
        Assert.Equal(2.3280, Math.Round(second.T, 4));
    }

    [Fact]
    public void LocalField_TheTwoMaximaAreTheFirstTwoAndSuccessive_OnAnIndependentGrid()
    {
        // A uniform grid on [0, π] sees exactly two interior local maxima and one local minimum between
        // them, the kink at π/2 where |sin 2t| touches zero. So nothing precedes the first maximum and
        // nothing sits between the two: they are the first two, and successive.
        const int n = 100_000;
        double h = Math.PI / n;
        var values = Enumerable.Range(0, n + 1)
            .Select(k => CpsiEnvelopeTheoremClaim.LocalFieldCpsiAtGammaZero(k * h)).ToArray();
        var maxima = new List<double>();
        var minima = new List<double>();
        for (int k = 1; k < n; k++)
        {
            if (values[k] > values[k - 1] && values[k] >= values[k + 1]) maxima.Add(k * h);
            if (values[k] < values[k - 1] && values[k] <= values[k + 1]) minima.Add(k * h);
        }
        Assert.Equal(2, maxima.Count);
        Assert.Single(minima);
        Assert.True(Math.Abs(minima[0] - Math.PI / 2.0) <= h, $"kink found at {minima[0]:R}");

        // The golden-section locations sit inside the grid step that brackets each grid maximum.
        var (first, second) = CpsiEnvelopeTheoremClaim.LocalFieldFirstTwoMaxima();
        Assert.True(Math.Abs(maxima[0] - first.T) < h, $"{maxima[0]:R} vs {first.T:R}");
        Assert.True(Math.Abs(maxima[1] - second.T) < h, $"{maxima[1]:R} vs {second.T:R}");
    }

    [Fact]
    public void NumberConserving_MicroMaximumBetweenSplitCorners_RisesToTheNextMainPeak()
    {
        // cos(π/8)|00⟩ + sin(π/8)|01⟩ under XXX at γ = 0.5, in closed form. The (0,1)-coherence corner stays
        // at π/4, the block-coherence corner moves to π/ω, ω = 2√(4 − γ²); between them sits a maximum.
        var m = CpsiEnvelopeTheoremClaim.MicroMaximumExample();
        Assert.True(m.SecondCornerT > m.FirstCornerT, "damping splits the corners");

        // A strict interior maximum of the interval between the corners is a local maximum of the curve.
        Assert.True(m.Micro.T > m.FirstCornerT && m.Micro.T < m.SecondCornerT, $"micro at {m.Micro.T:R}");
        Assert.True(m.Micro.Cpsi > m.FirstCornerCpsi, $"{m.Micro.Cpsi:R} vs {m.FirstCornerCpsi:R}");
        Assert.True(m.Micro.Cpsi > m.SecondCornerCpsi, $"{m.Micro.Cpsi:R} vs {m.SecondCornerCpsi:R}");

        // The refutation of the literal statement: a later local maximum above an earlier one.
        Assert.True(m.NextMainPeak.Cpsi > m.Micro.Cpsi);

        // An independent grid on (π/4, π/2): exactly these two maxima, and one minimum between them at the
        // moved corner, so they are successive.
        const int n = 20_000;
        double a = Math.PI / 4.0, h = (Math.PI / 2.0 - a) / n;
        var values = Enumerable.Range(0, n + 1)
            .Select(k => CpsiEnvelopeTheoremClaim.MicroMaximumExampleCpsi(a + k * h)).ToArray();
        var maxima = new List<double>();
        var minima = new List<double>();
        for (int k = 1; k < n; k++)
        {
            if (values[k] > values[k - 1] && values[k] >= values[k + 1]) maxima.Add(a + k * h);
            if (values[k] < values[k - 1] && values[k] <= values[k + 1]) minima.Add(a + k * h);
        }
        Assert.Equal(2, maxima.Count);
        Assert.Single(minima);
        Assert.True(Math.Abs(minima[0] - m.SecondCornerT) <= h, $"minimum at {minima[0]:R}");
        Assert.True(Math.Abs(maxima[0] - m.Micro.T) < h && Math.Abs(maxima[1] - m.NextMainPeak.T) < h);

        // simulations/envelope_n2_rises.py evaluates the same closed form in Python and checks it against the
        // matrix exponential of the 16×16 Liouvillian at these four times (to 1·10⁻¹⁶); the readings it
        // prints are the exact roundings of these values, so this pins the C# port of the closed form.
        Assert.Equal(0.087079030, Math.Round(m.Micro.Cpsi, 9));
        Assert.Equal(0.789207, Math.Round(m.Micro.T, 6));
        Assert.Equal(0.0986733, Math.Round(m.NextMainPeak.Cpsi, 7));
        Assert.Equal(1.00012, Math.Round(m.NextMainPeak.T, 5));
    }

    [Fact]
    public void PointwiseRise_NeedsTheLocalHamiltonian_RatesAreExact()
    {
        // ρ₀ = ((I + X/2 + Z/2)/2) ⊗ |0⟩⟨0| with first-site Z dephasing (Z₁ρZ₁ − ρ)/4. P' and L₁' are exact
        // dyadic rationals; CΨ' is one correctly rounded division by 3 of the exact dyadic 1/2 (or −1/4),
        // the same double as 1.0/6.0 (or −1.0/12.0). So all three compare exactly.
        var withH = CpsiEnvelopeTheoremClaim.PointwiseRiseRatesAtZero(withHamiltonian: true);
        Assert.True(withH.Purity == -1.0 / 8.0, $"P'(0) = {withH.Purity:R}");
        Assert.True(withH.L1 == 3.0 / 4.0, $"L1'(0) = {withH.L1:R}");
        Assert.True(withH.Cpsi == 1.0 / 6.0, $"CPsi'(0) = {withH.Cpsi:R}");

        // The mutation: the same state under the dephasing alone falls. H = Y⊗I is load-bearing.
        var dephasingOnly = CpsiEnvelopeTheoremClaim.PointwiseRiseRatesAtZero(withHamiltonian: false);
        Assert.True(dephasingOnly.Purity == -1.0 / 8.0, $"P'(0) = {dephasingOnly.Purity:R}");
        Assert.True(dephasingOnly.L1 == -1.0 / 4.0, $"L1'(0) = {dephasingOnly.L1:R}");
        Assert.True(dephasingOnly.Cpsi == -1.0 / 12.0, $"CPsi'(0) = {dephasingOnly.Cpsi:R}");
    }

    [Fact]
    public void Claim_PublicDirectParentProperties_AreExactlyF25AndQuarter()
    {
        PropertyInfo[] directParents = typeof(CpsiEnvelopeTheoremClaim)
            .GetProperties(BindingFlags.Instance | BindingFlags.Public | BindingFlags.DeclaredOnly)
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType))
            .ToArray();

        Assert.Equal(2, directParents.Length);

        var actual = directParents
            .Select(property => (Name: property.Name, PropertyType: property.PropertyType))
            .OrderBy(parent => parent.Name, StringComparer.Ordinal)
            .ToArray();
        var expected = new[]
        {
            (Name: "F25", PropertyType: typeof(F25CPsiBellPlusPi2Inheritance)),
            (Name: "Quarter", PropertyType: typeof(QuarterAsBilinearMaxvalClaim)),
        };

        Assert.Equal(expected, actual);
    }

    [Fact]
    public void Build_SharesOneQuarterInstance_AcrossF25AndTheDirectEdge()
    {
        var c = CpsiEnvelopeTheoremClaim.Shared;
        // Build() threads ONE algebraic QuarterAsBilinearMaxvalClaim into both F25 and the direct edge.
        Assert.Same(c.Quarter, c.F25.Quarter);
    }
}
