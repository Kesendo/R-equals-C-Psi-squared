using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One named configuration for <see cref="PhysicalGeneratorPolarityBreakWitness"/>:
/// the two Hermitian halves of H = A + iB as Pauli terms, the palindrome convention to read them
/// in, and what the configuration is FOR.</summary>
/// <param name="Name">Short label, used as the child node's display name.</param>
/// <param name="N">Qubit count.</param>
/// <param name="ATerms">Pauli terms of A, the Hermitian part of H.</param>
/// <param name="BTerms">Pauli terms of B; an empty list is a Hermitian H.</param>
/// <param name="DephaseLetter">The palindrome convention; Z is the theorem of proof §(a).</param>
/// <param name="What">One sentence on what this configuration decides, shown in the summary.</param>
public sealed record PhysicalGeneratorPolarityBreakCase(
    string Name,
    int N,
    IReadOnlyList<PauliTerm> ATerms,
    IReadOnlyList<PauliTerm> BTerms,
    PauliLetter DephaseLetter,
    string What);

/// <summary>What one configuration reads when the two routes are run against each other.</summary>
/// <param name="Case">The configuration this reading came from.</param>
/// <param name="Measured">The asymmetry read off the directly built generator.</param>
/// <param name="ClosedForm">F155's closed form, evaluated from the Pauli coefficients alone.</param>
/// <param name="Residual">Measured minus closed form.</param>
/// <param name="NoiseFloor">eps·‖M_anti‖², the absolute float noise of the measured route.</param>
/// <param name="MAntiNormSquared">‖M_anti‖², the whole of what the asymmetry is a difference of.</param>
/// <param name="IsStructurallySilent">The EXACT letter-level verdict: every string in A and in B
/// is Π²-even on this axis, so the generator carries no Π²-odd content at all.</param>
public sealed record PhysicalGeneratorPolarityBreakReading(
    PhysicalGeneratorPolarityBreakCase Case,
    double Measured,
    double ClosedForm,
    double Residual,
    double NoiseFloor,
    double MAntiNormSquared,
    bool IsStructurallySilent)
{
    /// <summary>|residual| / (eps·‖M_anti‖²): how far the two routes are apart, in units of the
    /// measured route's own float noise. A value near 1 means they agree as closely as the
    /// arithmetic permits.
    ///
    /// <para><b>The denominator is ‖M_anti‖² and not ‖M‖², and that is not a preference.</b> The
    /// asymmetry is a difference of the two halves of M_anti, so ‖M_anti‖² is what it is a
    /// difference OF; ‖M‖² is strictly larger, being ‖M_zero‖² + ‖M_anti‖². The repo retired the
    /// larger denominator once already on exactly this ground
    /// (<c>docs/CAUGHT_ERRORS.md</c>, 2026-08-06). Here it is measurably wrong and not merely
    /// looser: the asymmetry is provably blind to σ, since the 2σ·I term lies in Ad_Π's +1
    /// eigenspace which M_anti annihilates, so ‖M_anti‖² is σ-independent while ‖M‖² grows with σ
    /// without bound. A denominator that responds to a parameter the numerator cannot see is not a
    /// noise model. Gated by
    /// <c>TheAsymmetryAndItsFloorAreBlindToSigma_WhileTheRetiredDenominatorIsNot</c>, which carries
    /// the measured figures.</para>
    ///
    /// <para><b>Passing this is NOT the same as the measurement being accurate</b>, which is why
    /// <see cref="IsInformative"/> exists beside it. This ratio bounds the ABSOLUTE error, and
    /// where the two halves nearly cancel the absolute noise can exceed the whole signal. It is
    /// <see cref="double.NaN"/> where the floor is exactly 0, i.e. where M_anti vanishes and
    /// there is nothing to compare at all.</para></summary>
    public double NoiseRatio => NoiseFloor == 0.0 ? double.NaN : Math.Abs(Residual) / NoiseFloor;

    /// <summary>Signal over noise: |closed form| / (eps·‖M_anti‖²). This is the question
    /// <see cref="NoiseRatio"/> cannot answer, namely how far the predicted value stands above the
    /// float noise the measurement carries.
    ///
    /// <para>The asymmetry is BILINEAR in (A, B) while ‖M_anti‖² is quadratic in whichever half is
    /// larger, so driving the two halves apart, A by s and B by 1/s, holds the true value fixed
    /// while the noise floor grows like s². Measured at N=2 on the fixed sixteen-string fixture of
    /// <c>PullingTheTwoHalvesApart_DestroysTheMeasurement_AndTheWitnessSaysSo</c>: at s = 1 the two
    /// routes agree to 1 ulp, and at s = 10⁸ the measured route returns exactly 0.0 against a
    /// closed form of 17.5168, a relative error of 1, while the agreement ratio sits at a
    /// comfortable 0.178. Nothing is wrong with that ratio; it is answering a different question,
    /// and a gate that reported "passed" in that regime would be confirming nothing.</para>
    ///
    /// <para>It is 0 for a BALANCED configuration, where the predicted value IS zero, and that is
    /// not a defect: see <see cref="IsInformative"/> for why such a reading is refutable
    /// anyway.</para></summary>
    public double SignalToNoise =>
        NoiseFloor == 0.0 ? double.NaN : Math.Abs(ClosedForm) / NoiseFloor;

    /// <summary>Whether the measured route could have REFUTED the closed form, which is the only
    /// thing that makes an agreement worth reporting. Three cases, and they are three because the
    /// question is genuinely different in each.
    ///
    /// <list type="bullet">
    /// <item>No polarity content at all (‖M_anti‖² exactly 0, the SILENT case): nothing could have
    /// come out either way. NOT informative.</item>
    /// <item>The closed form predicts exactly ZERO on a non-empty M_anti: the measurement could
    /// have returned a nonzero asymmetry and did not, so it had teeth even though
    /// <see cref="SignalToNoise"/> is 0. Informative, and the question that decides it is
    /// <see cref="NoiseRatio"/> alone.</item>
    /// <item>The closed form predicts a NONZERO value: the measurement can separate it from zero
    /// only if it stands above its own noise, so this needs <see cref="SignalToNoise"/> above 1.
    /// Informative exactly then.</item>
    /// </list>
    ///
    /// <para>An earlier version wrote the middle case as <c>|| IsBitExact</c>. That reached the
    /// same verdict on today's configurations by the wrong route: it let the AGREEMENT itself count
    /// as evidence of the ability to disagree, which is circular, and it would have called a
    /// nonzero prediction informative merely because the subtraction happened to land on 0.0.
    /// Nothing in the threshold is tuned: 1 is where a value stops exceeding its own noise.</para></summary>
    public bool IsInformative =>
        MAntiNormSquared != 0.0
        && (ClosedForm == 0.0 || SignalToNoise > 1.0);

    /// <summary>True when the two routes agree BIT-EXACTLY. Reported, never assumed: which zeros
    /// are exact is a property of the configuration, not of the theorem.</summary>
    public bool IsBitExact => Residual == 0.0;

    /// <summary>SILENT where there is no polarity content to compare at all (the third answer
    /// beside the other two, and a "BALANCED" reading there would confirm nothing), BALANCED
    /// where the closed form is exactly zero on a non-empty M_anti, BROKEN otherwise.</summary>
    public string Verdict =>
        IsStructurallySilent || MAntiNormSquared == 0.0 ? "SILENT"
        : ClosedForm == 0.0 ? "BALANCED"
        : "BROKEN";
}

/// <summary>The live lab for F155 (<see cref="PhysicalGeneratorPolarityBreak"/>): the polarity
/// break of the physical (no-jump) generator ρ ↦ −i(Hρ − ρH†) for H = A + iB.
///
/// <para><b>What is derived and what is recomputed here.</b> The closed form
/// 4^(N+1)·Σ_{bit_b odd}(−1)^#Z·a_σ·b_σ is DERIVED, anchored in
/// <c>docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md</c>. What this witness does at
/// inspect time is run the OTHER route: build G_H as a 4^N × 4^N superoperator through
/// <see cref="NoJumpGenerator.Build"/>, push it through the palindrome residual and the ±i
/// refinement (<see cref="PolarityCoordinates.Decompose(MathNet.Numerics.LinearAlgebra.Matrix{Complex}, int, double, PauliLetter)"/>),
/// and read the asymmetry off the matrices. Two independent computations meeting; neither is a
/// lookup of the other.</para>
///
/// <para><b>The tolerance is a law, not a number, and it takes TWO numbers to state.</b> One route
/// is exact arithmetic on Pauli coefficients and the other is a chain of dense matrix products, so
/// there is no exact route to the comparison and a threshold that merely passes would be
/// worthless. Each reading therefore reports |residual| / (eps·‖M_anti‖²), the agreement in units
/// of the measured route's own noise, AND |closed form| / (eps·‖M_anti‖²), the signal in the same
/// units. The second is what decides whether the first means anything: the asymmetry is bilinear
/// in (A, B) while the noise floor is quadratic in whichever half is larger, so pulling the halves
/// apart can leave a comfortable agreement ratio on a measurement that has lost every digit.
/// Where the signal is at or under the noise the reading is marked NOT INFORMATIVE
/// (<see cref="PhysicalGeneratorPolarityBreakReading.IsInformative"/>) instead of passing quietly.
/// Whether a given zero is bit-exact is REPORTED
/// (<see cref="PhysicalGeneratorPolarityBreakReading.IsBitExact"/>), never assumed.</para>
///
/// <para><b>What the denominator is, and one thing it is not.</b> ‖M_anti‖² is what the asymmetry
/// is a difference OF. ‖M‖² would be strictly larger and is the form this repo retired on that
/// exact ground (<c>docs/CAUGHT_ERRORS.md</c>, 2026-08-06); here it is also σ-sensitive on a
/// numerator that is provably σ-blind, so it would loosen the floor by a factor the physics cannot
/// see.</para>
///
/// <para><b>σ is pinned to 0 and that is a decision.</b> The palindrome residual is
/// M = Π·L·Π⁻¹ + L + 2σ·I, and G_H carries no dephasing dissipator, so its σ is zero. The
/// chain-facing polarity paths take σ from a chain instead, which is why this witness goes
/// through the superoperator overload.</para>
///
/// <para>Reached by <c>dotnet run --project compute/RCPsiSquared.Cli -- inspect --root f155</c>;
/// pass <c>--N</c> to move the configurations that scale.</para>
/// </summary>
public sealed class PhysicalGeneratorPolarityBreakWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest N this witness will build. G_H is 4^N × 4^N dense and the residual
    /// costs several products of that size, so N=5 is already 1024×1024 per factor and N=6 would
    /// be 4096×4096: affordable, but nothing here needs it, and an inspect root must not be able
    /// to detonate on a stray <c>--N</c>.</summary>
    public const int MaxN = 5;

    /// <summary>How many multiples of the float noise floor the two routes may stand apart before
    /// the headline summary stops counting them as agreeing.
    ///
    /// <para><b>This is slack, not the law</b>, and it is worth being blunt about that. It is
    /// inherited from <c>EPS_RATIO_BOUND</c> in
    /// <c>simulations/f155_polarity_break_bilinear.py</c>, where it plays the same role, and
    /// neither side derives it: a residual accumulated over several dense 4^N × 4^N products can
    /// exceed a single-rounding floor by a small factor, and 8 is a generous ceiling on that factor
    /// rather than a computed bound. What carries the actual weight is the LAW, that the relative
    /// error is predicted by <see cref="PhysicalGeneratorPolarityBreakReading.NoiseRatio"/> divided
    /// by <see cref="PhysicalGeneratorPolarityBreakReading.SignalToNoise"/> and by nothing else,
    /// which <c>PullingTheTwoHalvesApart_DestroysTheMeasurement_AndTheWitnessSaysSo</c> tests
    /// across nine decades and which a ceiling alone could never establish.</para></summary>
    public const double AgreementSlack = 8.0;

    private readonly Lazy<IReadOnlyList<PhysicalGeneratorPolarityBreakReading>> _readings;

    public int N { get; }

    public PhysicalGeneratorPolarityBreakWitness(int n = 2)
    {
        if (n < 1)
            throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (n > MaxN)
            throw new ArgumentOutOfRangeException(
                nameof(n), $"N must be ≤ {MaxN}; got {n}. The generator is 4^N × 4^N dense.");
        N = n;
        _readings = new Lazy<IReadOnlyList<PhysicalGeneratorPolarityBreakReading>>(
            () => StandardSet(N).Select(Read).ToArray());
    }

    /// <summary>The readings, computed on first access and cached for the lifetime of this
    /// instance. A fresh instance recomputes from scratch; nothing here is banked.</summary>
    public IReadOnlyList<PhysicalGeneratorPolarityBreakReading> Readings => _readings.Value;

    /// <summary>The configurations, each one already a gate of
    /// <c>simulations/f155_polarity_break_bilinear.py</c> or its sibling script.
    ///
    /// <para><b>The count depends on N</b>, and the authority is
    /// <c>StandardSet_HasTheCountItsProseClaims</c> rather than any sentence here: four at N=1,
    /// eight at N=2, five at odd N above 1, six at even N above 2. The disjoint-supports case needs
    /// two sites, the PT dimer and its detuned partner are N=2 objects, and the cancelling zero
    /// needs an even number of sites to have equally many gain and loss ends.</para>
    ///
    /// <para>They do not all discriminate the same amount, and it is worth being exact about
    /// that. At N=2, three reach a NONZERO value (the F113 reduction, the positive control, the detuned
    /// dimer). Three reach zero through an EMPTY sum, which on the closed-form side is one branch
    /// three times (B = 0, disjoint supports, the dimer at its EP) though they differ on the
    /// measured side, where ‖M_anti‖² is 20, 272 and 16 at N=2. One reaches zero through a
    /// CANCELLATION, which is the only zero here that the CLOSED FORM reaches by subtraction rather
    /// than by an empty sum; the measured route could have returned nonzero on any of the four.
    /// One is silent.</para></summary>
    public static IReadOnlyList<PhysicalGeneratorPolarityBreakCase> StandardSet(int n)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");

        var cases = new List<PhysicalGeneratorPolarityBreakCase>();

        // (1) The parent edge, and the one that earns it: a per-site Z-drive meeting amplitude
        // damping, which is F113's own configuration reached through corollary 10's substitution.
        var driveA = new List<PauliTerm>();
        var dampB = new List<PauliTerm>();
        for (int l = 0; l < n; l++)
        {
            driveA.Add(PauliTerm.SingleSite(n, l, PauliLetter.Z, new Complex(0.25 + 0.5 * l, 0.0)));
            dampB.Add(PauliTerm.SingleSite(n, l, PauliLetter.Z, new Complex(0.125 * (l + 1), 0.0)));
        }
        cases.Add(new PhysicalGeneratorPolarityBreakCase(
            "F113 reduction (per-site Z-drive against damping)", n, driveA, dampB, PauliLetter.Z,
            "the configuration corollary 10 reduces to F113, and the case that earns the parent "
          + "edge; gate G7."));

        // (2) The positive control. Without one, a witness that reported zero everywhere would
        // look like a confirmation of balance rather than a broken instrument.
        var z0 = new[] { PauliTerm.SingleSite(n, 0, PauliLetter.Z, Complex.One) };
        cases.Add(new PhysicalGeneratorPolarityBreakCase(
            "positive control (A = B = Z₀)", n, z0, z0, PauliLetter.Z,
            $"the single shared string Z₀ is bit_b-odd with #Z = 1, so the value is exactly "
          + $"−4^(N+1) = {(-Math.Pow(4.0, n + 1)).ToString("G17", Inv)}; gate G8."));

        // (3) Corollary 1: a Hermitian H has B = 0 and is balanced whatever A is.
        cases.Add(new PhysicalGeneratorPolarityBreakCase(
            "Hermitian H (B = 0)", n, driveA, Array.Empty<PauliTerm>(), PauliLetter.Z,
            "corollary 1: with no anti-Hermitian part there is nothing for the odd content of A "
          + "to pair with, however odd it is."));

        // (4) Corollary 4: A and B both carry odd content, on strings that never meet.
        if (n >= 2)
        {
            var oddA = new[] { PauliTerm.SingleSite(n, 0, PauliLetter.Z, new Complex(1.5, 0.0)) };
            var oddB = new[] { PauliTerm.SingleSite(n, 1, PauliLetter.Z, new Complex(2.5, 0.0)) };
            cases.Add(new PhysicalGeneratorPolarityBreakCase(
                "disjoint odd supports", n, oddA, oddB, PauliLetter.Z,
                "corollary 4: both halves carry bit_b-odd content and the form still reads zero, "
              + "because the strings are different and the form is DIAGONAL. This is why "
              + "bit_b-odd content in A is necessary and NOT sufficient."));
        }

        // (5) The PT dimer, at its own exceptional point and then detuned. The single most
        // quotable output here: the quantity is blind to the PT-breaking threshold.
        if (n == 2)
        {
            var dimerA = new[]
            {
                new PauliTerm(new[] { PauliLetter.X, PauliLetter.X }, new Complex(0.5, 0.0)),
                new PauliTerm(new[] { PauliLetter.Y, PauliLetter.Y }, new Complex(0.5, 0.0)),
            };
            var gainLoss = new[]
            {
                PauliTerm.SingleSite(2, 0, PauliLetter.Z, new Complex(0.5, 0.0)),
                PauliTerm.SingleSite(2, 1, PauliLetter.Z, new Complex(-0.5, 0.0)),
            };
            cases.Add(new PhysicalGeneratorPolarityBreakCase(
                "PT dimer at its exceptional point (g = J = 1)", 2, dimerA, gainLoss, PauliLetter.Z,
                "corollary 7: exactly balanced at EVERY coupling including the EP, because the "
              + "bond letters XX and YY are bit_b-even and the gain-loss Z pair meets nothing. "
              + "The quantity is a detuning × gain-imbalance meter, not a PT meter; gate G11."));

            var detuned = dimerA.Concat(new[]
            {
                PauliTerm.SingleSite(2, 0, PauliLetter.Z, new Complex(0.4, 0.0)),
            }).ToArray();
            cases.Add(new PhysicalGeneratorPolarityBreakCase(
                "the same dimer with a 0.4·Z₀ detuning", 2, detuned, gainLoss, PauliLetter.Z,
                "add one Z to A and the same gain-loss profile now reads a nonzero value: the "
              + "detuning is what the meter sees, and the PT label decided nothing either way."));
        }

        // (6) The ONLY zero here that is a CANCELLATION rather than an empty sum, and therefore
        // the only one about which the two routes could disagree. Every other balanced case above
        // reaches zero because no string is shared, which the closed form settles without adding
        // anything; here two nonzero terms meet and cancel, so the matrix route has to reproduce a
        // subtraction rather than an absence. It is the registry entry's own worked example: a
        // uniform drive on a gain end and a loss end carries +ωγ against −ωγ.
        if (n >= 2)
        {
            var uniformDrive = new List<PauliTerm>();
            var gainLossProfile = new List<PauliTerm>();
            for (int l = 0; l < n; l++)
            {
                uniformDrive.Add(PauliTerm.SingleSite(n, l, PauliLetter.Z, new Complex(0.2, 0.0)));
                gainLossProfile.Add(PauliTerm.SingleSite(
                    n, l, PauliLetter.Z, new Complex(l % 2 == 0 ? 0.5 : -0.5, 0.0)));
            }
            if (n % 2 == 0)
            {
                cases.Add(new PhysicalGeneratorPolarityBreakCase(
                    "cancelling zero (uniform drive on a gain/loss profile)", n,
                    uniformDrive, gainLossProfile, PauliLetter.Z,
                    "the gain sites carry +ωγ and the loss sites −ωγ and the sum is exactly 0: a "
                  + "CANCELLATION, and the only zero here that the closed form reaches by "
                  + "subtraction rather than by an empty sum. Move the drive off uniformity and "
                  + "the same gain-loss profile reads nonzero, which is why no LABEL decides."));
            }
        }

        // (7) The third answer. Read under Π_X, the F113 configuration is not balanced, it is
        // SILENT: both ±i norms are individually and exactly zero.
        cases.Add(new PhysicalGeneratorPolarityBreakCase(
            "the F113 configuration read under Π_X", n, driveA, dampB, PauliLetter.X,
            "proof §(g), gate S8: every string here is bit_a-EVEN, so G_H sits in the +1 "
          + "eigenspace of Ad_{Π_X}² and there is no polarity content at all. Reading "
          + "\"balanced\" off this would confirm nothing."));

        return cases;
    }

    /// <summary>Run one configuration through both routes.</summary>
    public static PhysicalGeneratorPolarityBreakReading Read(PhysicalGeneratorPolarityBreakCase c)
    {
        if (c is null) throw new ArgumentNullException(nameof(c));

        var a = new PauliHamiltonian(c.N, c.ATerms).ToMatrix();
        var b = new PauliHamiltonian(c.N, c.BTerms).ToMatrix();
        var generator = NoJumpGenerator.Build(a, b);

        // σ = 0: G_H carries no dephasing dissipator, so the palindrome centre is not shifted.
        var coordinates = PolarityCoordinates.Decompose(generator, c.N, 0.0, c.DephaseLetter);

        double closedForm = PhysicalGeneratorPolarityBreak.PredictAsymmetry(
            c.ATerms, c.BTerms, c.N, c.DephaseLetter);

        double eps = Math.Pow(2.0, -52);
        return new PhysicalGeneratorPolarityBreakReading(
            Case: c,
            Measured: coordinates.Asymmetry,
            ClosedForm: closedForm,
            Residual: coordinates.Asymmetry - closedForm,
            NoiseFloor: eps * coordinates.MAntiNormSquared,
            MAntiNormSquared: coordinates.MAntiNormSquared,
            IsStructurallySilent: CarriesNoOddContent(c));
    }

    /// <summary>The EXACT letter-level silence test, decided before any float is looked at: if
    /// every string in A and in B is Π²-EVEN on this axis, then G_H lies in the +1 eigenspace of
    /// Ad_Π² and M_anti is the exact zero array, at every N and every coefficient. Proof §(g)
    /// states this in the bit_a case; it is the same argument on either axis.
    ///
    /// <para>A zero coefficient means the string is not in the operator at all, so it is skipped
    /// rather than read: a term that contributes nothing must not be able to flip an exact
    /// verdict, which is the trap <see cref="PolarityCoordinates.IsStructurallyDegenerate(IReadOnlyList{PauliTerm}, bool, PauliLetter)"/>
    /// records having walked into.</para></summary>
    public static bool CarriesNoOddContent(PhysicalGeneratorPolarityBreakCase c)
    {
        if (c is null) throw new ArgumentNullException(nameof(c));
        foreach (var term in c.ATerms.Concat(c.BTerms))
        {
            if (term.Coefficient == Complex.Zero) continue;
            if (PiOperator.SquaredEigenvalue(term.Letters, c.DephaseLetter) != +1) return false;
        }
        return true;
    }

    public string DisplayName => $"PhysicalGeneratorPolarityBreakWitness (F155 live, N={N})";

    public string Summary
    {
        get
        {
            var readings = Readings;
            var informative = readings.Where(r => r.IsInformative).ToArray();
            int agreeing = informative.Count(r => r.NoiseRatio <= AgreementSlack || r.IsBitExact);
            int exact = readings.Count(r => r.IsBitExact);
            double worst = informative
                .Where(r => !double.IsNaN(r.NoiseRatio))
                .Select(r => r.NoiseRatio)
                .DefaultIfEmpty(0.0)
                .Max();
            string verdicts = string.Join(", ", readings
                .GroupBy(r => r.Verdict)
                .OrderBy(g => g.Key, StringComparer.Ordinal)
                .Select(g => $"{g.Count()}×{g.Key}"));

            return $"N={N}: of {readings.Count} configurations, {informative.Length} are ones the "
                 + $"measured route could have refuted, and on those {agreeing} agree with the "
                 + $"closed form at the float noise ({exact} of all {readings.Count} bit-exactly), "
                 + $"worst |residual|/(eps·‖M_anti‖²) = {worst.ToString("G6", Inv)}; "
                 + $"verdicts {verdicts}.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                "the two routes",
                summary: "the closed form is exact arithmetic on Pauli coefficients and never "
                       + "builds a matrix; the measurement builds G_H = −i(A⊗I − I⊗Aᵀ) + "
                       + "(B⊗I + I⊗Bᵀ), takes the palindrome residual at σ = 0, and reads "
                       + "‖M₊‖² − ‖M₋‖² off it. Neither is a lookup of the other.",
                provenance: NodeProvenance.Live);

            foreach (var r in Readings)
            {
                yield return new InspectableNode(
                    r.Case.Name,
                    summary: $"{r.Verdict} on Π_{r.Case.DephaseLetter} at N={r.Case.N}: measured "
                           + $"{r.Measured.ToString("G17", Inv)} against closed form "
                           + $"{r.ClosedForm.ToString("G17", Inv)}, residual "
                           + $"{r.Residual.ToString("G6", Inv)}"
                           + (r.IsBitExact
                                ? " (bit-exact)"
                                : $" at |residual|/(eps·‖M_anti‖²) = {r.NoiseRatio.ToString("G6", Inv)}")
                           + (r.IsInformative
                                ? $", signal/noise {r.SignalToNoise.ToString("G6", Inv)}"
                                : ", and the measured route could NOT have refuted this reading "
                                + "(signal at or under its own float noise), so the agreement "
                                + "confirms nothing")
                           + $". {r.Case.What}",
                    provenance: NodeProvenance.Live);
            }
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
