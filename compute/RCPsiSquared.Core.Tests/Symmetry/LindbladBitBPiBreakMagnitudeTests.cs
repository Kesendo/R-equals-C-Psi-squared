using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F113 LindbladBitBPiBreakMagnitude Claim metadata, static-helper, and
/// cross-validation tests. Covers:
///
/// <list type="bullet">
///   <item>Claim metadata (Z2Axis, BitATwin classification, Tier, Theorem text).</item>
///   <item><see cref="LindbladBitBPiBreakMagnitude.PredictAsymmetry"/> +
///         <see cref="LindbladBitBPiBreakMagnitude.PredictAsymmetryUniform"/> closed-form
///         arithmetic (uniform N=2/3/4 anchors, non-uniform sum, detailed-balance
///         cancellation, σ⁻ ↔ σ⁺ sign flip, empty / length-mismatch edge cases).</item>
///   <item>Cross-validation against the matrix-level F112 polarity decomposition
///         (the integration anchor): construct L via the same Lindblad pipeline that
///         <see cref="PolarityCoordinates"/> uses, compute the asymmetry numerically,
///         and verify it matches <see cref="LindbladBitBPiBreakMagnitude.PredictAsymmetry"/>
///         bit-exactly at N=2 (uniform Z-drive + σ⁻) and N=3 (non-uniform rates).</item>
/// </list>
///
/// <para>The cross-validation tests use σ⁻ only (γ_pump = 0); the σ⁺ branch is
/// covered by the helper-arithmetic sign-flip test. Uniform Z-dephasing is added to
/// the dissipator to exercise the "Z-dephasing contributes 0" scope claim from the
/// F113 statement.</para></summary>
public class LindbladBitBPiBreakMagnitudeTests
{
    private static LindbladBitBPiBreakMagnitude Make()
    {
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        return new LindbladBitBPiBreakMagnitude(
            new LindbladBitBPiBalance(
                new F108Part1Pi2EvenAlwaysPalindromic(part2),
                new LindbladBitAPiBalance(part2)));
    }

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_IsBitB() =>
        Assert.Equal(Z2Axis.BitB, Make().Z2Axis);

    [Fact]
    public void BitATwin_IsNull()
    {
        // F113 is intrinsically about Z-axis single-site drives crossed with σ⁻ / σ⁺
        // amplitude damping. The commutator algebra [Z, σ⁻] = −2·σ⁻ has no meaningful
        // bit_a-axis analog; no twin exists by design, not by gap. Matches F108 Part 3.
        Assert.Null(Make().BitATwin);
    }

    [Fact]
    public void BitATwinStatus_IsBitBSpecific() =>
        Assert.Equal(BitATwinClassification.BitBSpecific, Make().BitATwinStatus);

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // Bit-exact at N=2, 3, 4 via constructive parameter sweep + per-site
        // decomposition + sign / variant verification. Universal-N status remains
        // Tier1Candidate (open algebraic derivation of the (1/2)·4^N coefficient
        // from Π-eigenspace structure); not promoted to a separate Claim.
        Assert.Equal(Tier.Tier1Derived, Make().Tier);
    }

    [Fact]
    public void Theorem_Mentions_closed_form_and_Z_drive_and_amplitude_damping()
    {
        var theorem = Make().Theorem;
        Assert.Contains("(4^N / 2)", theorem);
        Assert.Contains("Z_l", theorem);
        Assert.Contains("σ⁻_l", theorem);
        Assert.Contains("σ⁺_l", theorem);
    }

    // ============================================================
    // PredictAsymmetryUniform: bit-exact anchors at N=2, 3, 4
    // ============================================================

    [Fact]
    public void PredictAsymmetryUniform_N2_GivesMinusSixteen()
    {
        // ω = γ_T1 = 1, γ_pump = 0: asymmetry = (2 / 2) · 4^2 · 1 · (0 − 1) = −16.
        // Standard physics convention: T1 cooling at positive ω gives negative
        // polarity asymmetry. Bit-exact match required (no floating noise from integer powers).
        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetryUniform(
            omega: 1.0, gammaT1: 1.0, gammaPump: 0.0, N: 2);
        Assert.Equal(-16.0, a);
    }

    [Fact]
    public void PredictAsymmetryUniform_N3_GivesMinusNinetySix()
    {
        // (3 / 2) · 4^3 · 1 · (0 − 1) = 1.5 · 64 · (−1) = −96.
        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetryUniform(
            omega: 1.0, gammaT1: 1.0, gammaPump: 0.0, N: 3);
        Assert.Equal(-96.0, a);
    }

    [Fact]
    public void PredictAsymmetryUniform_N4_GivesMinusFiveHundredTwelve()
    {
        // (4 / 2) · 4^4 · 1 · (0 − 1) = 2 · 256 · (−1) = −512.
        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetryUniform(
            omega: 1.0, gammaT1: 1.0, gammaPump: 0.0, N: 4);
        Assert.Equal(-512.0, a);
    }

    // ============================================================
    // PredictAsymmetry: per-site sum + variant tests
    // ============================================================

    [Fact]
    public void PredictAsymmetry_NonUniformRates_MatchesSumFormula()
    {
        // Reproduces f113_break_formula_derivation.py "Non-uniform rates" anchor
        // under standard physics convention (σ⁻ lowering): pure σ⁻ T1 cooling at
        // positive ω gives NEGATIVE asymmetry.
        // ω_l = (0.05, 0.1, 0.2), γ_T1,l = (0.001, 0.002, 0.003), γ_pump = 0, N=3.
        // (1/2)·4^3 · Σ ω_l · (0 − γ_T1,l) = 32 · (−0.00085) = −0.0272.
        var omegas = new[] { 0.05, 0.1, 0.2 };
        var t1 = new[] { 0.001, 0.002, 0.003 };
        var pump = new[] { 0.0, 0.0, 0.0 };

        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetry(omegas, t1, pump, N: 3);

        double expected = 0.5 * Math.Pow(4.0, 3) *
                          (0.05 * (0.0 - 0.001) + 0.1 * (0.0 - 0.002) + 0.2 * (0.0 - 0.003));
        Assert.Equal(expected, a, precision: 12);
        Assert.Equal(-0.0272, a, precision: 12);
    }

    [Fact]
    public void PredictAsymmetry_DetailedBalance_IsZero()
    {
        // γ_T1 = γ_pump per site ⇒ (γ_T1,l − γ_pump,l) = 0 per site ⇒ asymmetry = 0.
        // F113 scope statement: "detailed-balance σ⁻ + σ⁺ contributions cancel".
        var omegas = new[] { 0.13, 0.27, 0.5 };
        var t1 = new[] { 0.001, 0.004, 0.002 };
        var pump = new[] { 0.001, 0.004, 0.002 };

        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetry(omegas, t1, pump, N: 3);

        Assert.Equal(0.0, a);
    }

    [Fact]
    public void PredictAsymmetry_SigmaPlusFlipsSign()
    {
        // F113 sign convention: σ⁺ contributes opposite sign to σ⁻. Pure σ⁻ (γ_pump = 0)
        // vs pure σ⁺ (γ_T1 = 0) at same ω, same magnitude rate ⇒ asymmetry sign flips,
        // magnitude unchanged.
        var omegas = new[] { 0.13, 0.07 };
        var rate = new[] { 0.001, 0.002 };
        var zero = new[] { 0.0, 0.0 };

        double aSigmaMinus = LindbladBitBPiBreakMagnitude.PredictAsymmetry(
            omegas, gammaT1PerSite: rate, gammaPumpPerSite: zero, N: 2);
        double aSigmaPlus = LindbladBitBPiBreakMagnitude.PredictAsymmetry(
            omegas, gammaT1PerSite: zero, gammaPumpPerSite: rate, N: 2);

        Assert.Equal(-aSigmaMinus, aSigmaPlus, precision: 12);
        Assert.NotEqual(0.0, aSigmaMinus);
    }

    [Fact]
    public void PredictAsymmetry_EmptyArrays_IsZero()
    {
        // N=0 is a degenerate scope (no qubits, no sites). All three arrays must be
        // empty; the sum is 0 vacuously, the prefactor 4^0 / 2 = 0.5 multiplies 0 to 0.
        var empty = Array.Empty<double>();
        double a = LindbladBitBPiBreakMagnitude.PredictAsymmetry(empty, empty, empty, N: 0);
        Assert.Equal(0.0, a);
    }

    [Fact]
    public void PredictAsymmetry_LengthMismatch_Throws()
    {
        // omegas.Count != N ⇒ ArgumentException with message naming the bad parameter.
        var omegas3 = new[] { 0.1, 0.2, 0.3 };
        var rates2 = new[] { 0.001, 0.002 };

        Assert.Throws<ArgumentException>(() =>
            LindbladBitBPiBreakMagnitude.PredictAsymmetry(omegas3, rates2, rates2, N: 2));
        Assert.Throws<ArgumentException>(() =>
            LindbladBitBPiBreakMagnitude.PredictAsymmetry(rates2, omegas3, rates2, N: 2));
        Assert.Throws<ArgumentException>(() =>
            LindbladBitBPiBreakMagnitude.PredictAsymmetry(rates2, rates2, omegas3, N: 2));
    }

    // ============================================================
    // Cross-validation against PolarityCoordinates.Decompose
    // (the integration anchor: F113 prediction matches numerical
    // F112 asymmetry computed from the actual Liouvillian)
    // ============================================================

    [Fact]
    public void PredictAsymmetry_MatchesPolarityCoordinatesDecompose_AtN2_UniformZDrive()
    {
        // Welle 2 anchor: N=2, ω=0.13 per site, γ_T1=0.001 per site, γ_Z=0.005.
        // Standard physics convention (σ⁻ = lowering): T1 cooling at positive
        // ω gives NEGATIVE polarity asymmetry. F113-predicted = (1/2)·4^2 · 2 ·
        // 0.13 · (0 − 0.001) = 16 · 0.13 · (−0.001) = −2.08e-3.
        //
        // Cross-validation builds L via the same Lindblad pipeline that
        // PolarityCoordinates.Decompose uses (T1Dissipator + Z-dephasing,
        // standard physics σ⁻), computes the asymmetry numerically, and
        // verifies bit-exact sign + magnitude agreement with PredictAsymmetry.
        // This is the integration anchor pinning the C# closed form to the
        // C# Diagnostics-layer numerical pipeline.
        const int N = 2;
        const double omega = 0.13;
        const double gammaT1 = 0.001;
        const double gammaZ = 0.005;

        var zDriveTerms = new List<PauliTerm>
        {
            PauliTerm.SingleSite(N, 0, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
            PauliTerm.SingleSite(N, 1, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
        };

        double measured = ComputeAsymmetryNumerically(
            N, zDriveTerms,
            gammaZPerSite: new[] { gammaZ, gammaZ },
            gammaT1PerSite: new[] { gammaT1, gammaT1 });

        double predicted = LindbladBitBPiBreakMagnitude.PredictAsymmetry(
            omegasPerSite: new[] { omega, omega },
            gammaT1PerSite: new[] { gammaT1, gammaT1 },
            gammaPumpPerSite: new[] { 0.0, 0.0 },
            N: N);

        // Bit-exact closed-form anchor: 16 · 0.13 · (0 − 0.001) = −0.00208.
        Assert.Equal(-0.00208, predicted, precision: 12);

        // Sign + magnitude equality (cross-validation relative tolerance 1e-9;
        // the PolarityCoordinates pipeline is numerical via vec → Pauli transform).
        double relDiff = Math.Abs(measured - predicted) /
                         Math.Max(Math.Abs(predicted), 1e-15);
        Assert.True(relDiff < 1e-9,
            $"F113 predicted {predicted:E6} vs measured {measured:E6}, " +
            $"rel diff {relDiff:E3}");
    }

    [Fact]
    public void PredictAsymmetry_MatchesPolarityCoordinatesDecompose_AtN3_NonUniformRates()
    {
        // Same f113_break_formula_derivation.py "Non-uniform rates" anchor at N=3:
        // ω_l = (0.05, 0.1, 0.2), γ_T1,l = (0.001, 0.002, 0.003), γ_Z = 0 (pure
        // Z-drive + σ⁻ T1 isolation). Standard physics convention: σ⁻ T1 cooling
        // at positive ω gives NEGATIVE asymmetry = −0.0272.
        const int N = 3;
        var omegas = new[] { 0.05, 0.1, 0.2 };
        var t1 = new[] { 0.001, 0.002, 0.003 };
        var gammaZ = new[] { 0.0, 0.0, 0.0 };

        var zDriveTerms = new List<PauliTerm>();
        for (int l = 0; l < N; l++)
            zDriveTerms.Add(PauliTerm.SingleSite(
                N, l, PauliLetter.Z, coefficient: (Complex)(omegas[l] / 2.0)));

        double measured = ComputeAsymmetryNumerically(
            N, zDriveTerms, gammaZPerSite: gammaZ, gammaT1PerSite: t1);

        double predicted = LindbladBitBPiBreakMagnitude.PredictAsymmetry(
            omegasPerSite: omegas,
            gammaT1PerSite: t1,
            gammaPumpPerSite: new[] { 0.0, 0.0, 0.0 },
            N: N);

        Assert.Equal(-0.0272, predicted, precision: 12);

        double relDiff = Math.Abs(measured - predicted) /
                         Math.Max(Math.Abs(predicted), 1e-15);
        Assert.True(relDiff < 1e-9,
            $"F113 predicted {predicted:E6} vs measured {measured:E6}, " +
            $"rel diff {relDiff:E3}");
    }

    // ============================================================
    // Numerical asymmetry helper: builds L via the same Z-dephasing +
    // σ⁻ T1 Lindblad pipeline that PolarityCoordinates.Decompose uses
    // internally (T1Dissipator.Build + PalindromeResidual.Build), then
    // computes the ±i Π-eigenvalue projector norms on M_anti directly.
    // Avoids the Diagnostics-layer dependency (Core.Tests does not
    // reference Diagnostics); the σ-shift convention matches
    // PolarityCoordinates.Decompose: σ_shift = Σ_l γ_Z,l (γ_T1 is NOT
    // included; see PolarityCoordinates.Decompose's chain.SigmaGamma call
    // at compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs).
    // ============================================================

    private static double ComputeAsymmetryNumerically(
        int N,
        IReadOnlyList<PauliTerm> hTerms,
        IReadOnlyList<double> gammaZPerSite,
        IReadOnlyList<double> gammaT1PerSite)
    {
        // Build H from k-body Pauli terms.
        var H = new PauliHamiltonian(N, hTerms).ToMatrix();

        // Build L with Z-dephasing + σ⁻ amplitude damping (the T1Dissipator path).
        // T1Dissipator.Build expects both per-site lists at full length N.
        var L = T1Dissipator.Build(H, gammaZPerSite, gammaT1PerSite);

        // σ-shift convention matches PolarityCoordinates.Decompose (chain.SigmaGamma):
        // only the Z-dephasing rates enter the F1-palindrome σ-shift.
        double sigma = 0.0;
        for (int l = 0; l < N; l++) sigma += gammaZPerSite[l];

        // M = Π · L_pauli · Π⁻¹ + L_pauli + 2σ·I in the 4^N Pauli-string basis.
        // PalindromeResidual.Build performs the vec → Pauli transform internally,
        // matching PolarityCoordinates' construction bit-exactly.
        var M = PalindromeResidual.Build(L, N, sigma, PauliLetter.Z);

        var piOp = PiOperator.BuildFull(N, PauliLetter.Z);
        var piInv = piOp.ConjugateTranspose(); // Π is unitary signed permutation; Π⁻¹ = Π†.

        // Split M into M_sym (Π²-even) + M_anti (Π²-odd), then refine M_anti into
        // ±i Π-eigenvalue projectors (the F112 asymmetry = ‖M_+i‖² − ‖M_-i‖²).
        var piMpi = piOp * M * piInv;
        var mAnti = (M - piMpi) / 2.0;
        var piMAntiPiInv = piOp * mAnti * piInv;
        var mPlusHalf = (mAnti - Complex.ImaginaryOne * piMAntiPiInv) / 2.0;
        var mMinusHalf = (mAnti + Complex.ImaginaryOne * piMAntiPiInv) / 2.0;

        double normPlus = Math.Pow(mPlusHalf.FrobeniusNorm(), 2);
        double normMinus = Math.Pow(mMinusHalf.FrobeniusNorm(), 2);
        return normPlus - normMinus;
    }
}
