using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The centre-line object inherits its x/y/z frame from the system and keeps the two
/// different spectral orbits separate: the linear F1 map fixes only λ = −σ, while F1 followed by
/// conjugation fixes the whole line Re λ = −σ. Both counts are exact algebraic multiplicities.</summary>
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

    /// <summary>The H = 0 count rebuilt from below rather than from the object's route. With H = 0
    /// and every rate γ, the coherence |i⟩⟨j| decays at −2γ·popcount(i ⊕ j) and σ = Nγ; sitting on
    /// the centre line Re λ = −σ therefore means popcount(i ⊕ j) = N/2 exactly. Enumerating the 4^N
    /// pairs and counting that condition is a different route to the same integer. At σ = 0 the line
    /// is Re λ = 0 and every pair qualifies, so the count is the whole space.</summary>
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
        // 2^N · C(N, N/2), written out. Neither 4^N (16, 256, 4096) nor 2^N (4, 16, 64) equals these,
        // so the numbers separate the count from the two trivial answers a wrong sector condition
        // would give.
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
    public void ZeroHamiltonian_TwoCountsCoincide_BecauseItsSpectrumIsReal(int n)
    {
        // H = 0 with a uniform rate: L is the dephasing generator alone, its spectrum is real, so
        // every centre-line mode already has Im λ = 0 and the composite fixed LINE collapses onto
        // the linear F1 fixed POINT. The equality of the two counts is a consequence, and this pins
        // it together with its reason: if the spectrum here ever carried a nonzero frequency, the
        // second assertion would fail before the first one could mislead.
        var h = Matrix<Complex>.Build.Dense(1 << n, 1 << n, Complex.Zero);
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", 0.1)).ToList();
        var system = new MirrorSystem(n, h, channels);
        var obj = new SelfMirrorObject(system);

        Assert.True(obj.IsFixedSetResolved);

        double worstFrequency = system.Spectrum.Modes.Max(m => System.Math.Abs(m.OscillationFrequency));
        Assert.True(worstFrequency == 0.0,
            $"H = 0 must carry a real spectrum; largest |Im λ| = {worstFrequency:E3}");

        Assert.Equal(obj.CompositeFixedLineCount, obj.LinearF1FixedPointCount);
    }

    /// <summary>The exact counts at N = 2 under γ = 0.1, as the factored characteristic polynomial in
    /// μ = λ + σ gives them (sympy over Q, at γ = 1/10; the double nearest 0.1, which is what the object
    /// reads, gives the same counts): Heisenberg μ⁴(25μ² + 24)(μ² + 1)²(5μ − 1)³(5μ + 1)³ at J = 1, so 4 + 2 + 4 = 10 on the line and 4 at the point, and the same 10/4 at J = 10³, 10⁴, 10⁸
    /// (25μ² + 24999999 and μ² + 10⁶ at J = 10³); XY (25μ² + 99)(μ² + 1)⁴(5μ − 1)³(5μ + 1)³, 10 on the
    /// line and none at the point. A floating window tied to the decay scale read the Heisenberg
    /// point count as 4, 3, 0 across J = 1, 10³, 10⁴; these do not move.</summary>
    [Theory]
    [InlineData(HamiltonianType.Heisenberg, 1.0, 10, 4)]
    [InlineData(HamiltonianType.Heisenberg, 1e3, 10, 4)]
    [InlineData(HamiltonianType.Heisenberg, 1e4, 10, 4)]
    [InlineData(HamiltonianType.Heisenberg, 1e8, 10, 4)]
    [InlineData(HamiltonianType.XY, 1.0, 10, 0)]
    [InlineData(HamiltonianType.XY, 1e3, 10, 0)]
    public void N2_Chain_ExactCounts(HamiltonianType hamiltonianType, double coupling, int line, int point)
    {
        var obj = Build(2, gamma: 0.1, coupling, hamiltonianType);

        Assert.True(obj.IsFixedSetResolved);
        Assert.Equal(line, obj.CompositeFixedLineCount);
        Assert.Equal(point, obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void ANonzeroHamiltonianCarriesFrequencies_AndTheTwoCountsPart()
    {
        // The control for the H = 0 coincidence: the XY chain at N = 2 has Im λ ≠ 0, and none of its
        // ten centre-line modes sits at the point λ = −σ, so the coincidence above is a property of
        // H = 0 and not of the two definitions.
        var system = new MirrorSystem(
            2,
            new ChainSystem(2, 1.0, 0.1, HamiltonianType.XY, TopologyKind.Chain).BuildHamiltonian(),
            Enumerable.Range(0, 2).Select(l => new ChannelRate($"q{l}", 0.1)).ToList());

        double worstFrequency = system.Spectrum.Modes.Max(m => System.Math.Abs(m.OscillationFrequency));
        Assert.True(worstFrequency > 1e-6,
            $"expected a genuinely complex spectrum, got |Im λ|max = {worstFrequency:E3}");
        var obj = new SelfMirrorObject(system);
        Assert.Equal(10, obj.CompositeFixedLineCount);
        Assert.Equal(0, obj.LinearF1FixedPointCount);
    }

    /// <summary>A diagonal H makes L diagonal in the |i⟩⟨j| basis: λ = −i(E_i − E_j) − 2γ·popcount(i ⊕ j).
    /// At N = 2 and γ = 0.1 the line Re λ = −0.2 holds the eight popcount-1 pairs, and none of them is
    /// at the point, because every popcount-1 energy difference below is nonzero, however small
    /// (10⁻¹⁵) and however far beside a large one (10⁸).</summary>
    [Theory]
    [InlineData(0.0, 1e-9, 3e-9, 7e-9)]
    [InlineData(0.0, 1e-15, 3e-15, 7e-15)]
    [InlineData(0.0, 1.0, 1e8, 1e8 + 2)]
    [InlineData(0.0, 1e-9, 1e8, 1e8 + 2)]
    public void N2_DiagonalHamiltonian_EightOnTheLine_NoneAtThePoint(double e0, double e1, double e2, double e3)
    {
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { e0, e1, e2, e3 });
        var channels = new[] { new ChannelRate("q0", 0.1), new ChannelRate("q1", 0.1) };
        var obj = new SelfMirrorObject(new MirrorSystem(2, h, channels));

        Assert.True(obj.IsFixedSetResolved);
        Assert.Equal(8, obj.CompositeFixedLineCount);
        Assert.Equal(0, obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void N1_NearCentrePair_IsNotOnTheLine()
    {
        // H = h·X at γ = 1: the Y-Z pair has λ = −γ ± √(γ² − 4h²) = −1 ± δ with δ ≈ 5·10⁻⁸, real and
        // off the centre Re λ = −1 on both sides; the other two modes sit at 0 and −2. Nothing is on
        // the line and nothing at the point, and the exact route says so without a window.
        const double delta = 5e-8;
        double hScale = 0.5 * Math.Sqrt(1.0 - delta * delta);
        var h = Matrix<Complex>.Build.DenseOfArray(new[,]
        {
            { Complex.Zero, new Complex(hScale, 0.0) },
            { new Complex(hScale, 0.0), Complex.Zero },
        });
        var obj = new SelfMirrorObject(new MirrorSystem(
            1, h, new[] { new ChannelRate("q0", 1.0) }));

        Assert.True(obj.IsFixedSetResolved);
        Assert.Equal(0, obj.CompositeFixedLineCount);
        Assert.Equal(0, obj.LinearF1FixedPointCount);
    }

    [Theory]
    [InlineData(2, HamiltonianType.XY, 10, 0)]
    [InlineData(2, HamiltonianType.Heisenberg, 10, 4)]
    [InlineData(3, HamiltonianType.XY, 0, 0)]
    [InlineData(3, HamiltonianType.Heisenberg, 0, 0)]
    public void FixedCounts_AreInvariantUnderCommonEnergyRescaling(
        int n,
        HamiltonianType hamiltonianType,
        int line,
        int point)
    {
        // (γ, J) = (0.1, 1) against (10⁻¹¹, 10⁻¹⁰): the doubles are not exactly proportional, and the
        // counts agree at both scales.
        var reference = Build(n, gamma: 0.1, coupling: 1.0, hamiltonianType);
        var rescaled = Build(n, gamma: 1e-11, coupling: 1e-10, hamiltonianType);

        Assert.Equal(line, reference.CompositeFixedLineCount);
        Assert.Equal(point, reference.LinearF1FixedPointCount);
        Assert.Equal(line, rescaled.CompositeFixedLineCount);
        Assert.Equal(point, rescaled.LinearF1FixedPointCount);
    }

    /// <summary>The site convention: site l is Kronecker factor l of H and carries γ_l. With
    /// H = Z_2 + Z_3 (a field on the last two sites) and γ = (1, 1, 2, 4), σ = 8 and the line needs
    /// disagreement sets D with Σ_D γ = 4: D = {0, 1, 2} or {3}, 16 pairs each. Every such pair flips
    /// site 2 or site 3 alone, so its energy difference is ±2 and none is at the point: 32/0. The
    /// reversed profile (4, 2, 1, 1) puts the line on {0} and {1, 2, 3}; flipping site 0 never changes
    /// Z_2 + Z_3 (16 at the point), flipping sites 2 and 3 together leaves it unchanged on the 8 kets
    /// where it is 0: 32/24 (both values predicted in sympy first). A map from sites to bits that ran the wrong way would
    /// swap the two rows.</summary>
    [Theory]
    [InlineData(1.0, 1.0, 2.0, 4.0, 32, 0)]
    [InlineData(4.0, 2.0, 1.0, 1.0, 32, 24)]
    public void SiteOrder_FieldOnTheLastTwoSites_NonUniformRates(
        double g0, double g1, double g2, double g3, int line, int point)
    {
        const int n = 4;
        var energies = new Complex[1 << n];
        for (int a = 0; a < 1 << n; a++)
        {
            int z2 = 1 - 2 * ((a >> 1) & 1);   // site 2 is bit N−1−2 = 1
            int z3 = 1 - 2 * (a & 1);          // site 3 is bit 0
            energies[a] = z2 + z3;
        }
        var h = Matrix<Complex>.Build.DiagonalOfDiagonalArray(energies);
        var gammas = new[] { g0, g1, g2, g3 };
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", gammas[l])).ToList();
        var obj = new SelfMirrorObject(new MirrorSystem(n, h, channels));

        Assert.Equal(line, obj.CompositeFixedLineCount);
        Assert.Equal(point, obj.LinearF1FixedPointCount);
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
    public void WeakDissipationAgainstUnitHamiltonian_KeepsTheHeisenbergCounts()
    {
        // γ = 10⁻¹¹ against J = 1: eleven decades between the two scales, and the same 10/4.
        var obj = Build(2, gamma: 1e-11, coupling: 1.0, hamiltonianType: HamiltonianType.Heisenberg);

        Assert.Equal(10, obj.CompositeFixedLineCount);
        Assert.Equal(4, obj.LinearF1FixedPointCount);
    }

    /// <summary>γ = 0: the line is Re λ = 0 and L = −i[H, ·] is anti-Hermitian, so all 4^N modes are on
    /// it; the point λ = 0 holds Σ_k m_k², the squared degeneracies of H (Heisenberg N = 2: triplet
    /// and singlet, 3² + 1² = 10; XY N = 2: energies −1, 0, 0, +1, so 1 + 4 + 1 = 6).</summary>
    [Theory]
    [InlineData(HamiltonianType.Heisenberg, 16, 10)]
    [InlineData(HamiltonianType.XY, 16, 6)]
    public void ClosedSystem_AllOnTheLine_PointIsTheSquaredDegeneracies(
        HamiltonianType hamiltonianType, int line, int point)
    {
        var obj = Build(2, gamma: 0.0, coupling: 1.0, hamiltonianType: hamiltonianType);

        Assert.Equal(line, obj.CompositeFixedLineCount);
        Assert.Equal(point, obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void AboveTheCostBound_TheCountThrowsInsteadOfGuessing()
    {
        // N = 5 XY: the joint-popcount sector (2, 3) has C(5,2)·C(5,3) = 100 > 36 entries.
        var obj = Build(5);

        Assert.False(obj.IsFixedSetResolved);
        Assert.Throws<InvalidOperationException>(() => obj.CompositeFixedLineCount);
        Assert.Contains("not computed", obj.Summary);
    }

    /// <summary>The counting step alone, on polynomials whose roots are known by construction.
    /// μ⁴(25μ² + 24)(μ² + 1)² has 4 + 2 + 4 = 10 roots on the imaginary axis and 4 at zero; the
    /// factor (5μ − 1)³(5μ + 1)³ that completes the N = 2 Heisenberg charpoly adds none. The pair
    /// (μ − δ)(μ + δ) with δ = 1/2^25, and the complex pair μ² − 2δμ + (δ² + 1), lie off the axis by δ
    /// and must count 0.</summary>
    [Fact]
    public void ImaginaryAxisCount_OnPolynomialsWithKnownRoots()
    {
        static GaussianInteger[] Mul(GaussianInteger[] a, GaussianInteger[] b)
        {
            var r = Enumerable.Repeat(GaussianInteger.Zero, a.Length + b.Length - 1).ToArray();
            for (int i = 0; i < a.Length; i++)
                for (int j = 0; j < b.Length; j++) r[i + j] += a[i] * b[j];
            return r;
        }
        static GaussianInteger[] Pow(GaussianInteger[] a, int k)
        {
            var r = new[] { GaussianInteger.One };
            for (int i = 0; i < k; i++) r = Mul(r, a);
            return r;
        }
        GaussianInteger[] mu = { 0, 1 };
        var onAxis = Mul(Mul(Pow(mu, 4), new GaussianInteger[] { 24, 0, 25 }), Pow(new GaussianInteger[] { 1, 0, 1 }, 2));
        var heis = Mul(Mul(onAxis, Pow(new GaussianInteger[] { -1, 5 }, 3)), Pow(new GaussianInteger[] { 1, 5 }, 3));
        Assert.Equal(10, CentreLineExactCount.ImaginaryAxisRootCount(heis));
        Assert.Equal(4, CentreLineExactCount.OrderAtZero(heis));

        var d = BigInteger.One << 25;   // δ = 1/d; scaled by d: (dμ − 1)(dμ + 1) and d²μ² − 2dμ + (1 + d²)
        var offReal = Mul(new GaussianInteger[] { -1, d }, new GaussianInteger[] { 1, d });
        var offComplex = new GaussianInteger[] { BigInteger.One + d * d, -2 * d, d * d };
        Assert.Equal(0, CentreLineExactCount.ImaginaryAxisRootCount(offReal));
        Assert.Equal(0, CentreLineExactCount.ImaginaryAxisRootCount(offComplex));
        Assert.Equal(2, CentreLineExactCount.ImaginaryAxisRootCount(Mul(offComplex, new GaussianInteger[] { 1, 0, 1 })));
    }
}
