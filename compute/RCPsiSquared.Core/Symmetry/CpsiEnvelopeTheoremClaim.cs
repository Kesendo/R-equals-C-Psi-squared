using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The CΨ envelope question for two qubits under local Z-dephasing (jump √γ Z_l on each
/// site), with CΨ = Tr(ρ²)·L₁(ρ)/3. Read literally, "for every two-qubit state and every H the successive
/// local maxima of CΨ are non-increasing" has counterexamples in every class, and its pointwise, absorbing
/// and local-control companions are false as well. What survives is a question about main peaks.
///
/// <list type="bullet">
///   <item><b>Local fields: the main peaks rise.</b> H = X⊗I + 0.37·I⊗Y from |01⟩. H is a sum of one-site
///         terms, so the state stays a product ρ₁⊗ρ₂, and a product state has P = P₁P₂ and
///         L₁ = (1 + ℓ₁)(1 + ℓ₂) − 1. At γ = 0 both qubits stay pure with ℓ₁ = |sin 2t| and
///         ℓ₂ = |sin 0.74t|, so CΨ(t) = [(1 + |sin 2t|)(1 + |sin 0.74t|) − 1]/3
///         (<see cref="LocalFieldCpsiAtGammaZero"/>). The two clocks are incommensurate, and the first two
///         maxima, smooth and strict, rise 0.734893 → 0.991268 (<see cref="LocalFieldFirstTwoMaxima"/>).
///         The propagator depends continuously on γ, so the rise survives weak damping; exact propagation
///         reads 0.731066 → 0.975722 at γ = 0.002 and 0.650531 → 0.688512 at γ = 0.05.</item>
///   <item><b>Number-conserving H: micro-maxima.</b> Damping moves the zeros of ρ₀₁,₁₀, which follows the
///         damped one-excitation clock, against the zeros of the (0,1)- and (1,2)-block coherences, which keep
///         the undamped clock (their dissipator is the uniform rate 2γ). Each zero is a V-corner of L₁; where
///         two such corners bracket a slope that changes sign, a micro-maximum can grow just above the trough,
///         below the next main peak. In a pure state the corners coincide at γ = 0 and split, which needs an
///         entry that passes through zero and is not generic. For ψ₀ = a|00⟩ + b|01⟩ + c|11⟩ under XXX the
///         trajectory is closed-form (<see cref="XxxPureCpsi"/>) with corners exactly at nπ/4 and nπ/ω,
///         ω = 2√(4 − γ²); for |0⟩ ⊗ R_y(π/4)|0⟩ at γ = 0.5 the micro-maximum 0.087079 lies 3.3·10⁻⁶ above both
///         corners and the next main peak is 0.098673 (<see cref="MicroMaximumExample"/>). In the mixed
///         (|i0⟩⟨i0| + |i+⟩⟨i+|)/2 under XXX at γ = 0.05 the corners sit 0.035 apart at γ = 0 and damping turns
///         the fixed one, a kink on a rising slope, into a trough: 0.234417 → 0.306486, while the main peaks
///         fall. These are examples of the mechanism, not a statement about all states.</item>
/// </list>
///
/// <para>What holds, each where it is proven: the named Bell+ channel formulas (F25's pure-Z closed
/// form is the direct typed parent), one-quarter as the algebraic bilinear maxval (the other parent),
/// and the conditional convergence implication: a continuous trajectory that converges to ρ* with
/// CΨ(ρ*) &lt; 1/4 and starts above 1/4 eventually stays below.</para>
///
/// <para>What stays open, and what this node's statement asks: under a number-conserving H, can a main peak,
/// the largest value of CΨ in a period T = 2π/Ω of the one-excitation block clock, ever rise above the one
/// before it? At γ = 0 U(T) is a phase on each excitation sector ([H, Z₁ + Z₂] = 0, Ω ≠ 0), so
/// CΨ(t + T) = CΨ(t) for every state and the main peaks are equal; for a pure state every maximum is equal, since CΨ = [(c + √p + √(s − p))² − 1]/3 is one-humped in
/// the one sinusoid p(t) = |ψ₀₁(t)|². A scan of 120 runs (pure and mixed) per period and 88 pure runs per
/// hump, at each of γ = 0.05, 0.2 and 0.5, finds no main peak that rises.</para>
///
/// <para>Tier <see cref="Tier.OpenQuestion"/>, because the statement is that question. The literal claim
/// the class is named for has exact counterexamples, so no Tier 1 statement remains to carry, and
/// <see cref="Tier.Retracted"/> would bury the question that is still asked. It has a tractable approach:
/// the γ = 0 flatness is exact (periodicity for every state, the lemma for pure ones), and the open step is
/// what damping does to a flat sequence of main peaks. The refutations are proven here as computing members
/// with their own tests, recomputed at inspect time; F25 and the quarter are Tier 1 parents.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> Part 5 +
/// <c>simulations/envelope_n2_rises.py</c> (exact propagation of every class, the main-peak scans) + the
/// finite N ≥ 3 rise atlas <c>compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs</c>.</para></summary>
public sealed class CpsiEnvelopeTheoremClaim : Claim
{
    /// <summary>The field on site 1 of the local-field counterexample, H = X⊗I + 0.37·I⊗Y.</summary>
    public const double LocalFieldX = 1.0;

    /// <summary>The field on site 2 of the local-field counterexample.</summary>
    public const double LocalFieldY = 0.37;

    /// <summary>The dephasing rate of the named micro-maximum example.</summary>
    public const double MicroMaximumGamma = 0.5;

    /// <summary>Direct parent: F25's named Bell+/local-Z closed form.</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    /// <summary>Direct parent: one-quarter as an algebraic bilinear maxval, not a dynamical absorber.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public CpsiEnvelopeTheoremClaim(
        F25CPsiBellPlusPi2Inheritance f25,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("Open question: for two qubits under a number-conserving H and local Z-dephasing, does the " +
               "largest value of CΨ = Tr(ρ²)·L₁/3 in each period of the one-excitation block clock, with the " +
               "windows anchored at a lowest point of the γ = 0 curve, never rise from one period to the next, " +
               "for every initial state? At γ = 0 CΨ repeats with that period for every state, and every " +
               "maximum of a pure state is equal. The literal statement that successive local maxima never " +
               "rise has counterexamples in every class: local fields raise the main peaks " +
               "(H = X⊗I + 0.37·I⊗Y from |01⟩, 0.734893 → 0.991268 at γ = 0), and under a number-conserving H " +
               "damping moves the zeros of ρ₀₁,₁₀ against those of the other coherences and can grow " +
               "micro-maxima just above a trough, below the next main peak, as a pure and a mixed example " +
               "show. Pointwise monotonicity, 1/4 as an absorbing boundary, and " +
               "protection below 1/4 against local unitaries are false. Exact: the named Bell+ channel " +
               "formulae, the algebraic quarter maxval, and the conditional convergence implication " +
               "(continuous ρ(t) → ρ*, CΨ(ρ*) < 1/4, start above ⇒ eventually stays below). Finite N/Q/K " +
               "atlas rows at N ≥ 3 are numerical readings only.",
               Tier.OpenQuestion,
               "docs/proofs/PROOF_MONOTONICITY_CPSI.md + " +
               "simulations/envelope_n2_rises.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/QuarterEnvelope.cs + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/EnvelopeBoundaryTests.cs + " +
               "experiments/ENVELOPE_RISE_BOUNDARY.md")
    {
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    /// <summary>Build with one shared algebraic-quarter parent across F25 and the direct edge.</summary>
    public static CpsiEnvelopeTheoremClaim Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        return new CpsiEnvelopeTheoremClaim(f25, quarter);
    }

    public static CpsiEnvelopeTheoremClaim Shared { get; } = Build();

    /// <summary>A local maximum of a CΨ curve: its time and its value.</summary>
    public readonly record struct LocalMaximum(double T, double Cpsi);

    /// <summary>The three rates at t = 0 of the pointwise example: purity, L₁ and CΨ.</summary>
    public readonly record struct CpsiRates(double Purity, double L1, double Cpsi);

    /// <summary>A pair of corners that damping has moved apart (their times and CΨ values), the
    /// micro-maximum between them, and the main peak that follows.</summary>
    public readonly record struct MicroMaximum(
        double FirstCornerT, double FirstCornerCpsi,
        double SecondCornerT, double SecondCornerCpsi,
        LocalMaximum Micro, LocalMaximum NextMainPeak);

    /// <summary>CΨ(t) at γ = 0 for |01⟩ under H = X⊗I + 0.37·I⊗Y: [(1 + |sin 2t|)(1 + |sin 0.74t|) − 1]/3.
    /// The state stays the pure product (cos t|0⟩ − i·sin t|1⟩) ⊗ (−sin 0.37t|0⟩ + cos 0.37t|1⟩); each
    /// factor's L₁ is |sin 2ht| for its field h, a product state has L₁ = (1 + ℓ₁)(1 + ℓ₂) − 1, and the
    /// purity is 1.</summary>
    public static double LocalFieldCpsiAtGammaZero(double t) =>
        ((1.0 + Math.Abs(Math.Sin(2.0 * LocalFieldX * t))) *
         (1.0 + Math.Abs(Math.Sin(2.0 * LocalFieldY * t))) - 1.0) / 3.0;

    /// <summary>The first two local maxima of <see cref="LocalFieldCpsiAtGammaZero"/>. At t = π/2 the factor
    /// |sin 2t| touches zero, a kink and a local minimum of the curve; the curve has one hump on [0, π/2]
    /// and one on [π/2, π], and golden-section search on each half locates it.</summary>
    public static (LocalMaximum First, LocalMaximum Second) LocalFieldFirstTwoMaxima() =>
        (GoldenSectionMaximum(LocalFieldCpsiAtGammaZero, 0.0, Math.PI / 2.0),
         GoldenSectionMaximum(LocalFieldCpsiAtGammaZero, Math.PI / 2.0, Math.PI));

    /// <summary>Exact CΨ(t) for ψ₀ = a|00⟩ + b|01⟩ + c|11⟩ (a, b, c ≥ 0, a² + b² + c² = 1) under
    /// H = XX + YY + ZZ with jumps √γ Z_l on both sites. |00⟩ and |11⟩ are eigenstates. The block
    /// (|01⟩, |10⟩) precesses at 4 about x with its coherence damped at 4γ, so with s = b² and
    /// ω = 2√(4 − γ²): r_z = s·e^{−2γt}[cos ωt + (2γ/ω) sin ωt], r_y = −(4s/ω)·e^{−2γt} sin ωt, r_x ≡ 0, and
    /// |ρ₀₁,₁₀| = |r_y|/2 vanishes at nπ/ω. The (0,1)- and (1,2)-block coherences are e^{−2γt} times their
    /// undamped values, of modulus ab|cos 2t|, ab|sin 2t|, cb|cos 2t|, cb|sin 2t|, vanishing at nπ/4; ρ₀₀,₁₁
    /// has modulus e^{−4γt}·ac. <c>simulations/envelope_n2_rises.py</c> checks the form against exact
    /// propagation (2.8·10⁻¹⁶).</summary>
    public static double XxxPureCpsi(double a, double b, double c, double gamma, double t)
    {
        double s = b * b;
        double omega = 2.0 * Math.Sqrt(4.0 - gamma * gamma);
        double e2 = Math.Exp(-2.0 * gamma * t);
        double rz = s * e2 * (Math.Cos(omega * t) + (2.0 * gamma / omega) * Math.Sin(omega * t));
        double ry = -(4.0 * s / omega) * e2 * Math.Sin(omega * t);
        double l1 = 2.0 * (e2 * b * (a + c) * (Math.Abs(Math.Cos(2.0 * t)) + Math.Abs(Math.Sin(2.0 * t)))
                           + e2 * e2 * a * c + Math.Abs(ry) / 2.0);
        double purity = Math.Pow(a, 4) + Math.Pow(c, 4) + (s * s + ry * ry + rz * rz) / 2.0
                        + 2.0 * e2 * e2 * s * (a * a + c * c) + 2.0 * Math.Pow(e2, 4) * a * a * c * c;
        return purity * l1 / 3.0;
    }

    /// <summary>CΨ(t) for the named micro-maximum example: |0⟩ ⊗ R_y(π/4)|0⟩ = cos(π/8)|00⟩ + sin(π/8)|01⟩
    /// under XXX at γ = <see cref="MicroMaximumGamma"/>.</summary>
    public static double MicroMaximumExampleCpsi(double t) =>
        XxxPureCpsi(Math.Cos(Math.PI / 8.0), Math.Sin(Math.PI / 8.0), 0.0, MicroMaximumGamma, t);

    /// <summary>The first corner pair of <see cref="MicroMaximumExampleCpsi"/>, split by damping: the (0,1)-coherence
    /// corner at π/4 and the block-coherence corner at π/ω, the micro-maximum between them (golden-section
    /// search on the open interval; it is one-humped there), and the main peak on (π/ω, π/2), whose next
    /// corner is at π/2.</summary>
    public static MicroMaximum MicroMaximumExample()
    {
        double omega = 2.0 * Math.Sqrt(4.0 - MicroMaximumGamma * MicroMaximumGamma);
        double first = Math.PI / 4.0, second = Math.PI / omega;
        return new MicroMaximum(
            first, MicroMaximumExampleCpsi(first),
            second, MicroMaximumExampleCpsi(second),
            GoldenSectionMaximum(MicroMaximumExampleCpsi, first, second),
            GoldenSectionMaximum(MicroMaximumExampleCpsi, second, Math.PI / 2.0));
    }

    /// <summary>Golden-section search for the maximum of a one-humped <paramref name="f"/> on
    /// [<paramref name="a"/>, <paramref name="b"/>]; 90 steps shrink the bracket by 0.618⁹⁰ ≈ 1.6·10⁻¹⁹.</summary>
    public static LocalMaximum GoldenSectionMaximum(Func<double, double> f, double a, double b)
    {
        double g = (Math.Sqrt(5.0) - 1.0) / 2.0;
        double c = b - g * (b - a), d = a + g * (b - a);
        double fc = f(c), fd = f(d);
        for (int i = 0; i < 90; i++)
        {
            if (fc > fd)
            {
                b = d; d = c; fd = fc;
                c = b - g * (b - a); fc = f(c);
            }
            else
            {
                a = c; c = d; fc = fd;
                d = a + g * (b - a); fd = f(d);
            }
        }
        double t = 0.5 * (a + b);
        return new LocalMaximum(t, f(t));
    }

    /// <summary>The pointwise example of <c>PROOF_MONOTONICITY_CPSI.md</c> Part 5, recomputed: the rates at
    /// t = 0 of purity, L₁ and CΨ for ρ₀ = ((I + X/2 + Z/2)/2) ⊗ |0⟩⟨0| under ρ' = −i[H, ρ] + (Z₁ρZ₁ − ρ)/4,
    /// with H = Y⊗I (<paramref name="withHamiltonian"/> true) or H = 0. The entries of ρ₀ and ρ' are small
    /// dyadic rationals, so P' and L₁' come out exactly: (−1/8, 3/4) with the local Hamiltonian and
    /// (−1/8, −1/4) under the dephasing alone. CΨ' = (P'·L₁ + P·L₁')/3 is one correctly rounded division by
    /// 3 of the exact dyadic 1/2 (or −1/4), so it equals the double 1.0/6.0 (or −1.0/12.0). The rise needs
    /// the Hamiltonian; local Z-dephasing by itself lowers purity and L₁ and so CΨ. A zero entry contributes
    /// its one-sided rate |ρ'_ij|, which is zero here: ρ₀ and ρ' share their support.</summary>
    public static CpsiRates PointwiseRiseRatesAtZero(bool withHamiltonian)
    {
        var i2 = PauliMatrix.Of(PauliLetter.I);
        var x = PauliMatrix.Of(PauliLetter.X);
        var y = PauliMatrix.Of(PauliLetter.Y);
        var z = PauliMatrix.Of(PauliLetter.Z);
        var ket0 = Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, 0 } });

        ComplexMatrix rho = ((i2 + x * 0.5 + z * 0.5) * 0.5).KroneckerProduct(ket0);
        ComplexMatrix h = withHamiltonian
            ? y.KroneckerProduct(i2)
            : Matrix<Complex>.Build.Dense(4, 4);
        ComplexMatrix z1 = z.KroneckerProduct(i2);
        ComplexMatrix rate = (h * rho - rho * h) * new Complex(0, -1) + (z1 * rho * z1 - rho) * 0.25;

        double purity = (rho * rho).Trace().Real;
        double purityRate = 2.0 * (rho * rate).Trace().Real;
        double l1 = 0.0, l1Rate = 0.0;
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
            {
                if (r == c) continue;
                double magnitude = rho[r, c].Magnitude;
                l1 += magnitude;
                l1Rate += magnitude > 0.0
                    ? (Complex.Conjugate(rho[r, c]) * rate[r, c]).Real / magnitude
                    : rate[r, c].Magnitude;
            }
        return new CpsiRates(purityRate, l1Rate, (purityRate * l1 + purity * l1Rate) / 3.0);
    }

    public override string DisplayName =>
        "CΨ envelope at N=2: can a main peak rise under a number-conserving H? (open; the literal statement fails)";

    public override string Summary
    {
        get
        {
            var (first, second) = LocalFieldFirstTwoMaxima();
            var micro = MicroMaximumExample();
            return "open: under a number-conserving H, can the largest value of CΨ in a period of the " +
                   "one-excitation block clock ever rise above the one before it? The literal successive-maxima " +
                   $"statement fails in every class: local fields raise the main peaks ({F6(first.Cpsi)} → " +
                   $"{F6(second.Cpsi)} at γ = 0, and the rise survives damping at γ = 0.002 and 0.05), and under a " +
                   "number-conserving H damping can grow micro-maxima just above a trough " +
                   $"({F6(micro.Micro.Cpsi)} → {F6(micro.NextMainPeak.Cpsi)} for |0⟩ ⊗ R_y(π/4)|0⟩ at γ = 0.5). " +
                   "Named Bell+ formulae, the algebraic quarter and the conditional convergence implication hold " +
                   $"({Tier.Label()})";
        }
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            var (first, second) = LocalFieldFirstTwoMaxima();
            var times = Enumerable.Range(0, 601).Select(k => k * 0.02).ToArray();
            yield return new InspectableNode("local fields on a product state: the main peaks rise",
                summary: "H = X⊗I + 0.37·I⊗Y from |01⟩ stays a product state; at γ = 0 " +
                         "CΨ(t) = [(1 + |sin 2t|)(1 + |sin 0.74t|) − 1]/3, first two maxima " +
                         $"{F6(first.Cpsi)} (t = {F4(first.T)}) → {F6(second.Cpsi)} (t = {F4(second.T)}), " +
                         "recomputed by golden-section search on the closed form. Continuity in γ carries " +
                         "the rise to small damping; exact propagation reads 0.731066 → 0.975722 at " +
                         "γ = 0.002 and 0.650531 → 0.688512 at γ = 0.05 (simulations/envelope_n2_rises.py).",
                payload: new InspectablePayload.Curve("CΨ(t) at γ = 0", times,
                    times.Select(LocalFieldCpsiAtGammaZero).ToArray(), "t", "CΨ"),
                provenance: NodeProvenance.Live);
            var micro = MicroMaximumExample();
            yield return new InspectableNode("number-conserving H: micro-maxima where damping moves one zero against another",
                summary: "cos(π/8)|00⟩ + sin(π/8)|01⟩ under XXX at γ = 0.5, closed form: the (0,1)-coherence " +
                         $"corner at π/4 (CΨ {F9(micro.FirstCornerCpsi)}) and the block-coherence corner at π/ω " +
                         $"= {F6(micro.SecondCornerT)} (CΨ {F9(micro.SecondCornerCpsi)}) enclose the micro-maximum " +
                         $"{F9(micro.Micro.Cpsi)} (t = {F6(micro.Micro.T)}); the next main peak is " +
                         $"{F7(micro.NextMainPeak.Cpsi)} (t = {F5(micro.NextMainPeak.T)}), a rising successive " +
                         "pair. At γ = 0 the two corners coincide and the pair is gone; a tuned pure state shows " +
                         "such a split at γ = 0.05 too, 7·10⁻⁹ high. In the mixed (|i0⟩⟨i0| + |i+⟩⟨i+|)/2 under XXX " +
                         "at γ = 0.05 the two corners, zeros of ρ₀₁,₁₀ and ρ₀₀,₁₀, sit 0.035 apart at γ = 0, and " +
                         "damping turns the fixed one, a kink on a rising slope, into a trough: 0.234417 → " +
                         "0.306486; all 13 of its rising pairs start at micro-maxima 3·10⁻⁶ to 5·10⁻⁵ above their " +
                         "troughs, while its main peaks fall. Examples of the mechanism, not a statement about " +
                         "all states.",
                provenance: NodeProvenance.Live);
            var withH = PointwiseRiseRatesAtZero(withHamiltonian: true);
            var dephasingOnly = PointwiseRiseRatesAtZero(withHamiltonian: false);
            yield return new InspectableNode("pointwise, absorbing and local-control statements fail",
                summary: "ρ₀ = ((I + X/2 + Z/2)/2) ⊗ |0⟩⟨0| (CΨ = 1/8) under H = Y⊗I with first-site Z " +
                         $"dephasing (Z₁ρZ₁ − ρ)/4: (P', L₁', CΨ')(0) = ({R(withH.Purity)}, {R(withH.L1)}, " +
                         $"{R(withH.Cpsi)}); with H = 0: ({R(dephasingOnly.Purity)}, {R(dephasingOnly.L1)}, " +
                         $"{R(dephasingOnly.Cpsi)}). A local Hadamard sends |00⟩ from CΨ = 0 to 1/3; when " +
                         "|0⟩ ⊗ (Bloch (sin 60°, 0, cos 60°)) has crossed 1/4 downward under Z-dephasing alone, " +
                         "the R_y that turns qubit 2's Bloch vector onto the equator (32.1°) lifts it back to " +
                         "0.2952; an active first-site Z pulse keeps CΨ = 1/8 yet turns CΨ'(0) from +1/6 to " +
                         "−1/3; a fixed local semigroup with H = 0 carries |00⟩ from 0 up through 1/4 toward " +
                         "its stationary |+0⟩ (CΨ = 1/3).",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the surviving conditional statement",
                summary: "continuity plus ρ(t) → ρ* with CΨ(ρ*) < 1/4 implies that any trajectory " +
                         "starting above 1/4 eventually stays below; locality or Markovianity alone does not.");
            yield return new InspectableNode("the open question: main peaks under a number-conserving H",
                summary: "the main peak of a period is the largest value of CΨ in one period T = 2π/Ω of the " +
                         "one-excitation block clock, the windows anchored at a lowest point of the γ = 0 curve. " +
                         "At γ = 0 ([H, Z₁ + Z₂] = 0, Ω ≠ 0) U(T) is a phase on each excitation sector, so " +
                         "CΨ(t + T) = CΨ(t) for every state and the main peaks are equal; for a pure state every " +
                         "maximum is equal, since CΨ = [(c + √p + √(s − p))² − 1]/3 is one-humped in the one " +
                         "sinusoid p(t) = |ψ₀₁(t)|². At N = 2 ZZ and a uniform field only rephase ρ_ij and are " +
                         "invisible to CΨ. Whether damping can ever make a main peak rise is not proven; 120 " +
                         "runs per period (pure and mixed) and 88 pure runs per hump, at each of γ = 0.05, 0.2 and 0.5, " +
                         "show no rise (simulations/envelope_n2_rises.py).");
            yield return F25;
            yield return Quarter;
        }
    }

    private static string F4(double v) => v.ToString("0.0000", CultureInfo.InvariantCulture);
    private static string F5(double v) => v.ToString("0.00000", CultureInfo.InvariantCulture);
    private static string F6(double v) => v.ToString("0.000000", CultureInfo.InvariantCulture);
    private static string F7(double v) => v.ToString("0.0000000", CultureInfo.InvariantCulture);
    private static string F9(double v) => v.ToString("0.000000000", CultureInfo.InvariantCulture);

    /// <summary>Print a rate of <see cref="PointwiseRiseRatesAtZero"/> as the small fraction it is.</summary>
    private static string R(double v)
    {
        foreach (int den in new[] { 1, 2, 3, 4, 6, 8, 12, 16, 24 })
        {
            double num = v * den;
            if (num == Math.Round(num))
                return den == 1
                    ? ((long)num).ToString(CultureInfo.InvariantCulture)
                    : $"{(num < 0 ? "−" : "")}{Math.Abs((long)num).ToString(CultureInfo.InvariantCulture)}/{den}";
        }
        return v.ToString("R", CultureInfo.InvariantCulture);
    }
}
