using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Polarity;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F155's live witness: the closed form meeting a directly built no-jump generator.
///
/// <para>This is the integration anchor for the whole carrier. The Core tests pin the closed
/// form's arithmetic against hand-derived literals; here the SECOND route runs, building
/// G_H = −i(Hρ − ρH†) as a 4^N × 4^N superoperator and reading the asymmetry off the palindrome
/// residual. If the theorem were wrong, or the Kronecker order or the σ convention mismatched,
/// this is where it would show.</para></summary>
public class PhysicalGeneratorPolarityBreakWitnessTests
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly ITestOutputHelper _output;

    public PhysicalGeneratorPolarityBreakWitnessTests(ITestOutputHelper output) => _output = output;

    // ============================================================
    // The two routes meet
    // ============================================================

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void EveryConfiguration_HasTheTwoRoutesAgreeing(int n)
    {
        var witness = new PhysicalGeneratorPolarityBreakWitness(n);

        foreach (var r in witness.Readings)
        {
            _output.WriteLine(
                $"{r.Case.Name} [Pi_{r.Case.DephaseLetter}] {r.Verdict}: measured "
              + $"{r.Measured.ToString("G17", Inv)}, closed form {r.ClosedForm.ToString("G17", Inv)}, "
              + $"residual {r.Residual.ToString("G6", Inv)}, ratio {r.NoiseRatio.ToString("G6", Inv)}, "
              + $"bit-exact {r.IsBitExact}");

            Assert.True(r.IsBitExact || r.NoiseRatio <= PhysicalGeneratorPolarityBreakWitness.AgreementSlack,
                $"'{r.Case.Name}' has the two routes disagreeing beyond the float noise: "
              + $"residual {r.Residual.ToString("G17", Inv)} is "
              + $"{r.NoiseRatio.ToString("G6", Inv)}× the floor eps·‖M_anti‖².");
        }
    }

    [Fact]
    public void AnOverallScaleCannotMoveTheRatio_SoASweepInItProvesNothing()
    {
        // This test exists to record WHY the obvious scaling sweep is not evidence, because an
        // earlier version of this file used exactly that sweep as its tolerance law. Scaling BOTH
        // halves by one factor s makes the asymmetry, the closed form and ‖M_anti‖² all
        // homogeneous of degree 2, so the ratio is homogeneous of degree ZERO and nothing but ulp
        // quantisation can move it. A gate whose quantity is provably invariant under the only
        // parameter it varies cannot fail, however comfortable its numbers look.
        var relativeErrors = new List<double>();
        foreach (double scale in new[] { 1e-3, 3e-2, 1.0, 30.0, 1e3 })
        {
            var reading = ReadScaled(scale, scale);
            _output.WriteLine(
                $"s={scale.ToString("G6", Inv)}: ratio {reading.NoiseRatio.ToString("G6", Inv)}, "
              + $"relative error {RelativeError(reading).ToString("G6", Inv)}");
            relativeErrors.Add(RelativeError(reading));
        }

        // What IS invariant, and what the sweep can therefore honestly assert: the relative error
        // stays at the ulp level throughout. The ratio is not asserted, because it is noise.
        Assert.True(relativeErrors.Max() < 1e-13,
            $"the worst relative error under a pure overall scale is "
          + $"{relativeErrors.Max().ToString("G6", Inv)}; degree-2 homogeneity of both routes "
          + "means anything above the ulp level is a real defect.");
    }

    [Fact]
    public void PullingTheTwoHalvesApart_DestroysTheMeasurement_AndTheWitnessSaysSo()
    {
        // The sweep that CAN fail, and the reason IsInformative exists. Scaling A by s and B by
        // 1/s holds the true asymmetry fixed (it is bilinear) while the noise floor grows like s²
        // (it is quadratic in the larger half). So the measured route loses digit after digit on
        // an unchanged answer. At s = 1e8 it returns exactly 0.0 against a closed form of 17.5168,
        // and the ‖M_anti‖² agreement ratio stays far inside any sane bound (measured 0.178 on
        // this fixture).
        // The agreement ratio is not wrong there; it is answering a different question. What must
        // hold is that the witness declines to call such a reading a confirmation.
        double referenceClosedForm = ReadScaled(1.0, 1.0).ClosedForm;
        bool sawAnUninformativeReading = false;

        foreach (double s in new[] { 1e0, 1e5, 1e7, 1e8, 1e9 })
        {
            var reading = ReadScaled(s, 1.0 / s);

            _output.WriteLine(
                $"s={s.ToString("G3", Inv)}: measured {reading.Measured.ToString("G12", Inv)}, "
              + $"closed {reading.ClosedForm.ToString("G12", Inv)}, "
              + $"ratio {reading.NoiseRatio.ToString("G4", Inv)}, "
              + $"signal/noise {reading.SignalToNoise.ToString("G4", Inv)}, "
              + $"informative {reading.IsInformative}, "
              + $"relative error {RelativeError(reading).ToString("G4", Inv)}");

            // The closed form is the invariant here: it must not move at all across the sweep.
            Assert.Equal(referenceClosedForm, reading.ClosedForm, 10);

            if (!reading.IsInformative)
            {
                sawAnUninformativeReading = true;
                continue;
            }

            // Where the witness DOES claim the measurement had teeth, the two numbers it reports
            // must jointly PREDICT how many digits survive, and they do: the relative error is
            // the agreement ratio divided by the signal-to-noise, because both are measured in
            // the same units of eps·‖M_anti‖². Asserting a flat 1e-10 here was wrong and this
            // sweep caught it: at s = 1e5 the reading is genuinely informative (signal/noise
            // 1.8e5) and genuinely carries a relative error of 8.7e-7, which is exactly
            // 0.154 / 1.8e5. The bound below says there is no error source beyond the modelled
            // one, which is a statement that can fail.
            double predicted =
                PhysicalGeneratorPolarityBreakWitness.AgreementSlack / reading.SignalToNoise;
            Assert.True(RelativeError(reading) <= predicted,
                $"at s={s.ToString("G3", Inv)} the measured value is off by "
              + $"{RelativeError(reading).ToString("G4", Inv)} relative, above the "
              + $"{predicted.ToString("G4", Inv)} the reported signal-to-noise "
              + $"{reading.SignalToNoise.ToString("G4", Inv)} accounts for. An error the error "
              + "model does not predict is the finding.");
        }

        Assert.True(sawAnUninformativeReading,
            "the sweep is supposed to reach a regime where the measurement carries no "
          + "information; if every reading is informative, the probe has stopped probing.");
    }

    [Fact]
    public void TheAsymmetryAndItsFloorAreBlindToSigma_WhileTheRetiredDenominatorIsNot()
    {
        // The argument for the denominator, made measurable. M = Π·L·Π⁻¹ + L + 2σ·I, and the
        // 2σ·I term lies in Ad_Π's +1 eigenspace, which M_anti annihilates. So the asymmetry and
        // ‖M_anti‖² cannot see σ at all, while ‖M‖² grows with it without bound. A denominator
        // that responds to a parameter the numerator is provably blind to is not a noise model,
        // and that is why ‖M‖² is not used here. The figures quoted in the witness doc comment
        // come from this test.
        var a = new List<PauliTerm>();
        var b = new List<PauliTerm>();
        double[] aCoefficients = { 0.31, -0.72, 0.14, 0.58, -0.23, 0.91, -0.45, 0.07,
                                   0.66, -0.19, 0.83, -0.37, 0.52, 0.28, -0.61, 0.44 };
        double[] bCoefficients = { -0.27, 0.63, -0.88, 0.11, 0.74, -0.36, 0.49, -0.05,
                                   0.22, 0.95, -0.53, 0.68, -0.41, 0.17, 0.79, -0.34 };
        for (long k = 0; k < 16; k++)
        {
            var letters = PauliIndex.FromFlat(k, 2);
            a.Add(new PauliTerm(letters, new Complex(aCoefficients[k], 0.0)));
            b.Add(new PauliTerm(letters, new Complex(bCoefficients[k], 0.0)));
        }

        var generator = NoJumpGenerator.Build(
            new PauliHamiltonian(2, a).ToMatrix(), new PauliHamiltonian(2, b).ToMatrix());

        double? firstAsymmetry = null;
        double? firstMAnti = null;
        var mNorms = new List<double>();

        foreach (double sigma in new[] { 0.0, 0.05, 1.0, 17.3 })
        {
            var c = PolarityCoordinates.Decompose(generator, 2, sigma);
            _output.WriteLine(
                $"σ={sigma.ToString("G4", Inv)}: asymmetry {c.Asymmetry.ToString("G17", Inv)}, "
              + $"‖M‖² {c.MNormSquared.ToString("G6", Inv)}, "
              + $"‖M_anti‖² {c.MAntiNormSquared.ToString("G6", Inv)}");

            // Blind means BIT-IDENTICAL. There is an exact route here, so there is nothing to
            // tolerate: the term the σ shift adds is annihilated, not merely made small.
            firstAsymmetry ??= c.Asymmetry;
            firstMAnti ??= c.MAntiNormSquared;
            Assert.Equal(firstAsymmetry.Value, c.Asymmetry);
            Assert.Equal(firstMAnti.Value, c.MAntiNormSquared);
            mNorms.Add(c.MNormSquared);
        }

        // And the retired denominator does move, by well over an order of magnitude, on a
        // measurement that has not changed by a single bit.
        Assert.True(mNorms.Max() / mNorms.Min() > 10.0,
            $"‖M‖² ran only from {mNorms.Min().ToString("G6", Inv)} to "
          + $"{mNorms.Max().ToString("G6", Inv)}; if it no longer responds to σ, the argument for "
          + "preferring ‖M_anti‖² needs restating rather than assuming.");
    }

    [Fact]
    public void TheSuperoperatorDecomposeOverload_GuardsItsArguments()
    {
        // The route F155 needs into the polarity layer. Its guards are exercised nowhere else.
        var generator = NoJumpGenerator.Build(
            new PauliHamiltonian(2, new[] { PauliTerm.SingleSite(2, 0, PauliLetter.Z, Complex.One) })
                .ToMatrix(),
            Matrix<Complex>.Build.Dense(4, 4));

        Assert.Throws<ArgumentNullException>(
            () => PolarityCoordinates.Decompose(null!, 2, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => PolarityCoordinates.Decompose(generator, 0, 0.0));
        // A 16×16 superoperator is N=2; asking for N=3 wants 64×64.
        Assert.Throws<ArgumentException>(
            () => PolarityCoordinates.Decompose(generator, 3, 0.0));
    }

    private static double RelativeError(PhysicalGeneratorPolarityBreakReading r) =>
        r.ClosedForm == 0.0 ? Math.Abs(r.Residual) : Math.Abs(r.Residual / r.ClosedForm);

    /// <summary>One fixed N=2 configuration carrying content on every one of the 16 Pauli strings,
    /// with A scaled by <paramref name="aScale"/> and B by <paramref name="bScale"/>. The
    /// coefficients are written out rather than drawn from a generator, so the sweep is the same
    /// arithmetic on every run.</summary>
    private static PhysicalGeneratorPolarityBreakReading ReadScaled(double aScale, double bScale)
    {
        double[] aCoefficients =
        {
            0.31, -0.72, 0.14, 0.58, -0.23, 0.91, -0.45, 0.07,
            0.66, -0.19, 0.83, -0.37, 0.52, 0.28, -0.61, 0.44,
        };
        double[] bCoefficients =
        {
            -0.27, 0.63, -0.88, 0.11, 0.74, -0.36, 0.49, -0.05,
            0.22, 0.95, -0.53, 0.68, -0.41, 0.17, 0.79, -0.34,
        };

        var a = new List<PauliTerm>();
        var b = new List<PauliTerm>();
        for (long k = 0; k < 16; k++)
        {
            var letters = PauliIndex.FromFlat(k, 2);
            a.Add(new PauliTerm(letters, new Complex(aCoefficients[k] * aScale, 0.0)));
            b.Add(new PauliTerm(letters, new Complex(bCoefficients[k] * bScale, 0.0)));
        }

        return PhysicalGeneratorPolarityBreakWitness.Read(
            new PhysicalGeneratorPolarityBreakCase(
                $"A×{aScale}, B×{bScale}", 2, a, b, PauliLetter.Z, "scaling probe"));
    }

    // ============================================================
    // Individual configurations, each with its own verdict
    // ============================================================

    [Fact]
    public void PositiveControl_ReadsExactlyMinusFourToTheNPlusOne()
    {
        // Without a configuration whose value is large and known, a witness reporting zeros
        // everywhere would look like a confirmation instead of a broken instrument.
        var readings = new PhysicalGeneratorPolarityBreakWitness(3).Readings;
        var control = readings.Single(r => r.Case.Name.Contains("positive control"));

        Assert.Equal(-256.0, control.ClosedForm);
        Assert.Equal("BROKEN", control.Verdict);
        Assert.True(control.IsBitExact || control.NoiseRatio <= 8.0);
    }

    [Fact]
    public void HermitianH_IsBalanced_AndTheClosedFormIsExactlyZero()
    {
        var readings = new PhysicalGeneratorPolarityBreakWitness(2).Readings;
        var hermitian = readings.Single(r => r.Case.Name.Contains("Hermitian H"));

        Assert.Equal(0.0, hermitian.ClosedForm);
        Assert.Equal("BALANCED", hermitian.Verdict);
    }

    [Fact]
    public void DisjointOddSupports_AreBalanced_ThoughBothHalvesCarryOddContent()
    {
        // The case that shows the sorting is one-directional: bit_b-odd content in A is
        // necessary and NOT sufficient.
        var readings = new PhysicalGeneratorPolarityBreakWitness(2).Readings;
        var disjoint = readings.Single(r => r.Case.Name.Contains("disjoint odd supports"));

        Assert.Equal(0.0, disjoint.ClosedForm);
        Assert.Equal("BALANCED", disjoint.Verdict);
        Assert.False(disjoint.IsStructurallySilent,
            "both halves carry bit_b-odd strings, so this is a genuine balance and not a silence");
    }

    [Fact]
    public void ThePTDimer_IsBalancedAtItsExceptionalPoint_AndTheDetunedOneIsNot()
    {
        var readings = new PhysicalGeneratorPolarityBreakWitness(2).Readings;
        var atEp = readings.Single(r => r.Case.Name.Contains("exceptional point"));
        var detuned = readings.Single(r => r.Case.Name.Contains("detuning"));

        Assert.Equal(0.0, atEp.ClosedForm);
        Assert.Equal("BALANCED", atEp.Verdict);

        Assert.NotEqual(0.0, detuned.ClosedForm);
        Assert.Equal("BROKEN", detuned.Verdict);

        // The PT label decided nothing: same gain-loss profile, opposite verdicts, and the only
        // difference is one detuning term in A.
        _output.WriteLine(
            $"at the EP {atEp.ClosedForm.ToString("G17", Inv)}, detuned "
          + $"{detuned.ClosedForm.ToString("G17", Inv)}");
    }

    [Fact]
    public void TheF113ConfigurationUnderPiX_IsSilentRatherThanBalanced()
    {
        var readings = new PhysicalGeneratorPolarityBreakWitness(2).Readings;
        var silent = readings.Single(r => r.Case.DephaseLetter == PauliLetter.X);

        Assert.True(silent.IsStructurallySilent,
            "every string is a Z, and bit_a(Z) = 0, so there is no Π_X²-odd content at all");
        Assert.Equal("SILENT", silent.Verdict);

        // The exact claim of proof §(g): the two halves are individually zero, not cancelling.
        Assert.Equal(0.0, silent.MAntiNormSquared);
    }

    [Fact]
    public void TheStructuralSilenceTest_SkipsZeroCoefficientTerms()
    {
        // A term with a zero coefficient is not in the operator at all. Reading the letters
        // without the coefficient would flip an exact verdict on a term that contributes
        // nothing, which is the trap IsStructurallyDegenerate records having walked into.
        var a = new[] { PauliTerm.SingleSite(2, 0, PauliLetter.X, Complex.One) };
        var bWithGhost = new[]
        {
            PauliTerm.SingleSite(2, 0, PauliLetter.X, Complex.One),
            PauliTerm.SingleSite(2, 1, PauliLetter.Z, Complex.Zero),
        };

        Assert.True(PhysicalGeneratorPolarityBreakWitness.CarriesNoOddContent(
            new PhysicalGeneratorPolarityBreakCase(
                "ghost term", 2, a, bWithGhost, PauliLetter.Z, "probe")));
    }

    // ============================================================
    // The witness contract
    // ============================================================

    [Fact]
    public void Children_AreAllLive_TheWitnessRecomputesEveryRead()
    {
        var witness = new PhysicalGeneratorPolarityBreakWitness(2);
        foreach (var child in witness.Children)
            Assert.Equal(NodeProvenance.Live, child.Provenance);
    }

    [Fact]
    public void Summary_ReportsTheCountsAndTheWorstRatio()
    {
        string summary = new PhysicalGeneratorPolarityBreakWitness(2).Summary;

        Assert.Contains("N=2", summary);
        Assert.Contains("bit-exactly", summary);
        Assert.Contains("eps·‖M_anti‖²", summary);
        Assert.Contains("could have refuted", summary);
    }

    [Theory]
    [InlineData(1, 4)]
    [InlineData(2, 8)]
    [InlineData(3, 5)]
    [InlineData(4, 6)]
    public void StandardSet_HasTheCountItsProseClaims(int n, int expected)
    {
        // The count is N-dependent and nothing pinned it, which is how five separate documents
        // came to say "six configurations" while the set returns eight at N=2 and five at N=3.
        // Four always run; the disjoint-supports pair needs two sites; the PT dimer and its
        // detuned partner are N=2 objects; the cancelling zero needs equally many gain and loss
        // ends, hence even N.
        Assert.Equal(expected, PhysicalGeneratorPolarityBreakWitness.StandardSet(n).Count);
    }

    [Fact]
    public void StandardSet_NamesAreUnique_SoAReadingCanBeIdentifiedByName()
    {
        // The per-configuration tests below select by name; a duplicate would make Single throw
        // and, worse, would let two different physics hide behind one label.
        for (int n = 1; n <= 4; n++)
        {
            var names = PhysicalGeneratorPolarityBreakWitness.StandardSet(n)
                .Select(c => $"{c.Name}|{c.DephaseLetter}").ToList();
            Assert.Equal(names.Count, names.Distinct().Count());
        }
    }

    [Fact]
    public void Constructor_RejectsNBelowOne() =>
        Assert.Throws<ArgumentOutOfRangeException>(
            () => new PhysicalGeneratorPolarityBreakWitness(0));

    [Fact]
    public void Constructor_RejectsNAboveTheGuard() =>
        Assert.Throws<ArgumentOutOfRangeException>(
            () => new PhysicalGeneratorPolarityBreakWitness(
                PhysicalGeneratorPolarityBreakWitness.MaxN + 1));

    [Fact]
    public void Read_RejectsNull() =>
        Assert.Throws<ArgumentNullException>(() => PhysicalGeneratorPolarityBreakWitness.Read(null!));

}
