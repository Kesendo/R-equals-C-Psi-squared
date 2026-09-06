using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The centre-line object inherits its x/y/z frame from the system and keeps the two
/// different spectral orbits separate: the linear F1 map fixes only λ = −σ, while F1 followed by
/// conjugation fixes the whole line Re λ = −σ.</summary>
public class SelfMirrorObjectTests
{
    private static SelfMirrorObject Build(
        int n,
        double gamma = 0.1,
        double coupling = 1.0,
        HamiltonianType hamiltonianType = HamiltonianType.XY)
    {
        var h = new ChainSystem(n, coupling, gamma, hamiltonianType, TopologyKind.Chain).BuildHamiltonian();
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", gamma)).ToList();
        return new SelfMirrorObject(new MirrorSystem(n, h, channels));
    }

    [Fact]
    public void Object_InheritsFrame_AndSitsAtMinusSigma()
    {
        var obj = Build(4, gamma: 0.1);
        Assert.Equal(4, obj.N);              // N inherited from the system, not owned
        Assert.Equal(0.4, obj.Sigma, 10);    // σ = Σγ = Nγ inherited
        Assert.Equal(-0.4, obj.Center, 10);  // the object sits at Re λ = −σ
    }

    /// <summary>The count the object reports, rebuilt from below rather than from its own
    /// formula. On the certified branch H = 0 and every rate is γ, so the coherence |i⟩⟨j| decays
    /// at −2γ·popcount(i ⊕ j) and σ = Nγ; sitting on the centre line Re λ = −σ therefore means
    /// popcount(i ⊕ j) = N/2 exactly. Enumerating the 4^N pairs and counting that condition is a
    /// different route to the same integer, and it is the route the tests gate against. At σ = 0
    /// the line is Re λ = 0 and every pair qualifies, so the count is the whole space.</summary>
    private static int EnumerateCentreLinePairs(int n, bool sigmaIsZero)
    {
        int dim = 1 << n;
        int count = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                if (sigmaIsZero) { count++; continue; }
                if ((n & 1) == 0 && System.Numerics.BitOperations.PopCount((uint)(i ^ j)) == n / 2)
                    count++;
            }
        return count;
    }

    [Theory]
    [InlineData(2, true)]
    [InlineData(4, true)]
    [InlineData(3, false)]
    [InlineData(5, false)]
    public void CompositeFixedLine_IsTheEnumeratedCentreLineCount(int n, bool populated)
    {
        int count = Build(n, coupling: 0.0).CompositeFixedLineCount;
        Assert.Equal(EnumerateCentreLinePairs(n, sigmaIsZero: false), count);
        Assert.Equal(populated, count > 0);   // odd N: popcount cannot reach the half-integer N/2
    }

    [Theory]
    [InlineData(2, 8)]      // C(2,1) = 2 masks of popcount 1, times 2^2 = 4 kets
    [InlineData(4, 96)]     // C(4,2) = 6 masks of popcount 2, times 2^4 = 16 kets
    [InlineData(6, 1280)]   // C(6,3) = 20 masks of popcount 3, times 2^6 = 64 kets
    public void CompositeFixedLine_MatchesTheClosedForm(int n, int expected)
    {
        // 2^N · C(N, N/2), the closed form the object carries, written out. Neither 4^N (16, 256,
        // 4096) nor 2^N (4, 16, 64) equals these, so the numbers separate the count from the two
        // trivial answers a wrong sector condition would give.
        Assert.Equal(expected, EnumerateCentreLinePairs(n, sigmaIsZero: false));
        if (n <= 4)
            Assert.Equal(expected, Build(n, coupling: 0.0).CompositeFixedLineCount);
        Assert.NotEqual(1 << (2 * n), expected);
        Assert.NotEqual(1 << n, expected);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void CompositeFixedLine_AtZeroSigma_IsTheWholeSpace(int n)
    {
        // σ = 0 puts the centre line at Re λ = 0, where the whole (dephasing-free) spectrum sits:
        // the count is 4^N, at odd N too, so the odd-N zero above is a statement about the rate
        // sector and not about odd N as such.
        var obj = Build(n, gamma: 0.0, coupling: 0.0);
        Assert.True(obj.IsFixedSetResolved);
        Assert.Equal(EnumerateCentreLinePairs(n, sigmaIsZero: true), obj.CompositeFixedLineCount);
        Assert.Equal(1 << (2 * n), obj.CompositeFixedLineCount);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void CertifiedBranch_TwoCountsCoincide_BecauseItsSpectrumIsReal(int n)
    {
        // The certified branch is H = 0 with a uniform rate. There L is the dephasing generator
        // alone, whose spectrum is real, so every centre-line mode already has Im λ = 0 and the
        // composite fixed LINE collapses onto the linear F1 fixed POINT. The equality of the two
        // counts is therefore a consequence, and this pins the consequence together with its
        // reason: if the spectrum here ever carried a nonzero frequency, the second assertion
        // would fail before the first one could mislead.
        var h = Matrix<Complex>.Build.Dense(1 << n, 1 << n, Complex.Zero);
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", 0.1)).ToList();
        var system = new MirrorSystem(n, h, channels);
        var obj = new SelfMirrorObject(system);

        Assert.True(obj.IsFixedSetResolved);

        double worstFrequency = system.Spectrum.Modes.Max(m => System.Math.Abs(m.OscillationFrequency));
        Assert.True(worstFrequency == 0.0,
            $"the certified branch must carry a real spectrum; largest |Im λ| = {worstFrequency:E3}");

        Assert.Equal(obj.CompositeFixedLineCount, obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void ANonzeroHamiltonianDoesCarryFrequencies_SoTheTwoCountsCouldDiffer()
    {
        // The control for the test above: on the branch the object refuses, Im λ ≠ 0 genuinely
        // occurs, so the coincidence just pinned is a property of the certified branch and not a
        // property of the two definitions.
        var system = new MirrorSystem(
            2,
            new ChainSystem(2, 1.0, 0.1, HamiltonianType.XY, TopologyKind.Chain).BuildHamiltonian(),
            Enumerable.Range(0, 2).Select(l => new ChannelRate($"q{l}", 0.1)).ToList());

        double worstFrequency = system.Spectrum.Modes.Max(m => System.Math.Abs(m.OscillationFrequency));
        Assert.True(worstFrequency > 1e-6,
            $"expected a genuinely complex spectrum off the certified branch, got |Im λ|max = {worstFrequency:E3}");
        Assert.False(new SelfMirrorObject(system).IsFixedSetResolved);
    }

    [Fact]
    public void N2_NonzeroHamiltonianDoesNotInventEitherExactMultiplicity()
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(
            new Complex[] { 0, 1e-9, 3e-9, 7e-9 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Theory]
    [InlineData(1e3)]
    [InlineData(1e4)]
    public void N2_Heisenberg_IntermediateCouplingsDoNotReportAFalseChangedMultiplicity(double coupling)
    {
        var obj = Build(2, gamma: 0.1, coupling, HamiltonianType.Heisenberg);

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N2_GenuinelyTinyNonzeroFrequenciesAreNotRoundedIntoLinearFixedPoints()
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(
            new Complex[] { 0, 1e-15, 3e-15, 7e-15 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N1_GenuinelyNearCentrePairIsNotRoundedOntoTheCompositeFixedLine()
    {
        const double delta = 5e-8;
        double hScale = 0.5 * Math.Sqrt(1.0 - delta * delta);
        var h = Matrix<Complex>.Build.DenseOfArray(new[,]
        {
            { Complex.Zero, new Complex(hScale, 0.0) },
            { new Complex(hScale, 0.0), Complex.Zero },
        });
        var obj = new SelfMirrorObject(new MirrorSystem(
            1, h, new[] { new ChannelRate("q0", 1.0) }));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N2_LargeHamiltonianScale_DoesNotWidenTheDecayAxisIntoTheWholeSpectrum()
    {
        var obj = Build(2, gamma: 0.1, coupling: 1e8, hamiltonianType: HamiltonianType.Heisenberg);

        Assert.False(obj.IsFixedSetResolved);
        Assert.Contains("UNRESOLVED", obj.Summary);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N2_LargeUnrelatedFrequency_DoesNotTurnResolvedSmallFrequenciesIntoZero()
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { 0, 1, 1e8, 1e8 + 2 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N2_DecayWindow_DoesNotTurnExactlyResolvedTinyFrequenciesIntoZero()
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(
            new Complex[] { 0, 1e-9, 3e-9, 7e-9 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N2_UnrelatedHugeBranch_DoesNotEraseAnExactTinyCenterLineFrequency()
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(
            new Complex[] { 0, 1e-9, 1e8, 1e8 + 2 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.False(obj.IsFixedSetResolved);
        Assert.Contains("UNRESOLVED", obj.Summary);
        Assert.Throws<InvalidOperationException>(() => obj.LinearF1FixedPointCount);
    }

    [Theory]
    [InlineData(3, HamiltonianType.XY)]
    [InlineData(3, HamiltonianType.Heisenberg)]
    public void FixedCounts_AreInvariantUnderCommonEnergyRescaling(
        int n,
        HamiltonianType hamiltonianType)
    {
        var reference = Build(n, gamma: 0.1, coupling: 1.0, hamiltonianType);
        var rescaled = Build(n, gamma: 1e-11, coupling: 1e-10, hamiltonianType);

        Assert.False(reference.IsFixedSetResolved);
        Assert.False(rescaled.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => reference.CompositeFixedLineCount);
        Assert.Throws<InvalidOperationException>(() => rescaled.CompositeFixedLineCount);
    }

    [Fact]
    public void LiveStrings_NameCompositeMap_AndFenceItFromLinearF1()
    {
        var obj = Build(2);
        string rendered = string.Join("\n", new[] { obj.DisplayName, obj.Summary }
            .Concat(obj.Children.Select(c => $"{c.DisplayName}\n{c.Summary}")));

        Assert.Contains("λ ↦ −2σ − conj(λ)", rendered);
        Assert.Contains("linear F1", rendered);
        Assert.Contains("λ = −σ", rendered);
        Assert.DoesNotContain("each its own mirror under λ ↦ −2σ − λ", rendered);
    }

    [Fact]
    public void LiveStrings_RenderSmallNonzeroScaleInScientificNotation()
    {
        var obj = Build(2, gamma: 1e-11, coupling: 1e-10);

        Assert.Contains("σ = 2E-11", obj.Summary);
        Assert.Contains("−σ = -2E-11", obj.DisplayName);
        Assert.DoesNotContain("σ = 0", obj.Summary);
    }

    [Fact]
    public void WeakDissipationAgainstUnitHamiltonian_IsSurfacedAsUnresolved()
    {
        var obj = Build(2, gamma: 1e-11, coupling: 1.0, hamiltonianType: HamiltonianType.Heisenberg);

        Assert.False(obj.IsFixedSetResolved);
        Assert.Contains("UNRESOLVED", obj.Summary);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
    }

    [Fact]
    public void ClosedSystemWithNumericallyAmbiguousZeroFrequenciesIsUnresolved()
    {
        var obj = Build(2, gamma: 0.0, coupling: 1.0, hamiltonianType: HamiltonianType.Heisenberg);

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
    }
}
