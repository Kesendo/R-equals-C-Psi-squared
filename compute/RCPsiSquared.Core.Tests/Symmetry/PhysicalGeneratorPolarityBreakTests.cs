using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F155 <see cref="PhysicalGeneratorPolarityBreak"/> Claim metadata and closed-form
/// arithmetic. The matrix-level cross-validation, where the closed form meets a directly built
/// generator, lives in <c>RCPsiSquared.Diagnostics.Tests</c> beside the witness, because it needs
/// the polarity decomposition.
///
/// <para>Every expectation here is written as a LITERAL derived by hand from the theorem, never
/// by calling the method under test in a different arrangement: a gate whose expectation comes
/// out of the thing it gates cannot fail.</para></summary>
public class PhysicalGeneratorPolarityBreakTests
{
    private static PhysicalGeneratorPolarityBreak Make()
    {
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        return new PhysicalGeneratorPolarityBreak(
            new LindbladBitBPiBreakMagnitude(
                new LindbladBitBPiBalance(
                    new F108Part1Pi2EvenAlwaysPalindromic(part2),
                    new LindbladBitAPiBalance(part2))));
    }

    private static PauliTerm Site(int n, int site, PauliLetter letter, double coefficient) =>
        PauliTerm.SingleSite(n, site, letter, new Complex(coefficient, 0.0));

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Tier_IsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Make().Tier);

    [Fact]
    public void Z2Axis_IsBitB() =>
        Assert.Equal(Z2Axis.BitB, Make().Z2Axis);

    [Fact]
    public void BitATwin_IsNull_BecauseNoBespokeTwinIsOwed() =>
        Assert.Null(Make().BitATwin);

    [Fact]
    public void BitATwinStatus_IsCoveredByHadamardDuality()
    {
        // The asymmetry is a Frobenius-norm statement on the 4^N Pauli space and Ad_{Q_zx} is
        // Frobenius-unitary with Q_zx·Π_Z·Q_zx⁻¹ = Π_X, which is case 1 of
        // docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md. Proof §(g) states the transport as the
        // identity asymmetry_X(A, B) = asymmetry_Z(hAh, hBh), gate S10 measures it.
        //
        // This does NOT contradict F113's BitBSpecific next door: that verdict is about a
        // physical CHANNEL (σ⁺/σ⁻ single out the Z energy ladder), and F155 makes no channel
        // assumption at all.
        Assert.Equal(BitATwinClassification.CoveredByHadamardDuality, Make().BitATwinStatus);
    }

    [Fact]
    public void Parent_IsTheF113BreakMagnitude() =>
        Assert.NotNull(Make().BreakMagnitude);

    [Fact]
    public void Anchor_NamesTheProofAndBothGateScripts()
    {
        string anchor = Make().Anchor;
        Assert.Contains("PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md", anchor);
        Assert.Contains("simulations/f155_polarity_break_bilinear.py", anchor);
        Assert.Contains("simulations/f155_dephase_siblings.py", anchor);
    }

    // ============================================================
    // The closed form
    // ============================================================

    [Fact]
    public void PredictAsymmetry_SingleZOnBothHalves_IsMinusFourToTheNPlusOne()
    {
        // A = B = Z₀ at N=3. The one shared string is Z⊗I⊗I: bit_b(Z) = 1 so it is in the
        // support, #Z = 1 so the weight is −1, and a·b = 1. Hence 4^4 · (−1) = −256 exactly.
        var a = new[] { Site(3, 0, PauliLetter.Z, 1.0) };
        var b = new[] { Site(3, 0, PauliLetter.Z, 1.0) };

        Assert.Equal(-256.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 3));
    }

    [Fact]
    public void PredictAsymmetry_TwoZLetters_WeightIsPlusOneAndSupportIsEven_SoItVanishes()
    {
        // Z⊗Z⊗I has bit_b = 0, which is OUTSIDE the support: an even number of bit_b letters
        // leaves the string Π²-even. The weight would have been (−1)² = +1, and it never
        // gets to apply. This separates the two halves of the rule, which a single-letter
        // fixture cannot.
        var zz = new[] { new PauliTerm(
            new[] { PauliLetter.Z, PauliLetter.Z, PauliLetter.I }, Complex.One) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(zz, zz, 3));
    }

    [Fact]
    public void PredictAsymmetry_YStringIsInTheSupportAndCarriesWeightPlusOne()
    {
        // Y₀ at N=2: bit_b(Y) = 1 so it is in the support, but #Z = 0, so the weight is +1 and
        // the sign is OPPOSITE to the Z case. Support and weight are two different letters'
        // business, and this is the fixture that shows it.
        var y = new[] { Site(2, 0, PauliLetter.Y, 1.0) };

        Assert.Equal(64.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(y, y, 2));
    }

    [Fact]
    public void PredictAsymmetry_IsDiagonal_DistinctStringsNeverPair()
    {
        // A carries Z₀ (in the support), B carries Z₁ (also in the support). They are DIFFERENT
        // strings, so the bilinear form pairs neither with the other and the sum is empty.
        var a = new[] { Site(2, 0, PauliLetter.Z, 1.0) };
        var b = new[] { Site(2, 1, PauliLetter.Z, 1.0) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
    }

    [Fact]
    public void PredictAsymmetry_BitBEvenContentIsSilentAgainstEveryB()
    {
        // X₀ has bit_b(X) = 0. Corollary 2: A entirely bit_b-even is balanced against every B,
        // however large the coefficients.
        var x = new[] { Site(2, 0, PauliLetter.X, 7.5) };
        var bBig = new[] { Site(2, 0, PauliLetter.X, 1000.0) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(x, bBig, 2));
    }

    [Fact]
    public void PredictAsymmetry_HermitianH_MeansBIsEmpty_AndTheResultIsExactlyZero()
    {
        // Corollary 1: B = 0 (a Hermitian H) is balanced whatever A is.
        var a = new[] { Site(2, 0, PauliLetter.Z, 3.0), Site(2, 1, PauliLetter.Y, -1.5) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(
            a, Array.Empty<PauliTerm>(), 2));
    }

    [Fact]
    public void PredictAsymmetry_BProportionalToIdentity_IsExactlyZero()
    {
        // Corollary 3: the identity string has bit_b = 0 and is outside the support.
        var a = new[] { Site(2, 0, PauliLetter.Z, 3.0) };
        var identity = new[] { new PauliTerm(
            new[] { PauliLetter.I, PauliLetter.I }, new Complex(2.5, 0.0)) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, identity, 2));
    }

    [Fact]
    public void PredictAsymmetry_RepeatedStringsAreSummed_NotOverwritten()
    {
        // A caller may pass an unreduced list. Two Z₀ terms at 1.0 and 2.0 are the single
        // string Z₀ at 3.0, so the value is 4^3 · (−1) · 3 · 1 = −192 exactly.
        var a = new[] { Site(2, 0, PauliLetter.Z, 1.0), Site(2, 0, PauliLetter.Z, 2.0) };
        var b = new[] { Site(2, 0, PauliLetter.Z, 1.0) };

        Assert.Equal(-192.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
    }

    [Fact]
    public void PredictAsymmetry_TheXConventionMovesTheSupportToBitA()
    {
        // Proof §(g): under Π_X the support is bit_a-odd and the weight counts X. Z₀ has
        // bit_a(Z) = 0, so the string that carries the whole Z-convention value is SILENT here,
        // both halves individually zero. This is the S8 shape at the level of the closed form.
        var z = new[] { Site(2, 0, PauliLetter.Z, 1.0) };
        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(
            z, z, 2, PauliLetter.X));

        // X₀ is bit_a-odd with #X = 1, so it carries weight −1 under Π_X, and it was the
        // silent one under Π_Z. The two conventions swap which string they can see.
        var x = new[] { Site(2, 0, PauliLetter.X, 1.0) };
        Assert.Equal(-64.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(
            x, x, 2, PauliLetter.X));
        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(x, x, 2));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void PredictAsymmetry_TheYConventionIsExactlyMinusTheZOne(int n)
    {
        // Proof §(g): Π_Y = Π_Z⁻¹ exactly, so inverting an order-4 map exchanges its ±i
        // eigenprojections and asymmetry_Y = −asymmetry_Z for EVERY superoperator, with no
        // hypothesis at all. Y and Z share the support (both graded by bit_b) and differ only in
        // which letter the weight counts, and on bit_b-odd strings (−1)^#Y = −(−1)^#Z. The
        // closed form must carry that identity, and the C# side had no gate on the Y branch at
        // all. NOT a novelty claim, and the distinction is one this repo has already paid for:
        // simulations/f112_klein_v4_cross_dephase_verify.py has been passing "Y" into the Python
        // polarity primitive for months (line 210), and asserting otherwise is the exact error
        // docs/CAUGHT_ERRORS.md records on 2026-08-19 about this same formula. What was missing
        // is a gate on THIS implementation, which is a smaller and checkable statement.
        var a = new List<PauliTerm>();
        var b = new List<PauliTerm>();
        for (int l = 0; l < n; l++)
        {
            a.Add(Site(n, l, PauliLetter.Z, 0.5 + l));
            a.Add(Site(n, l, PauliLetter.Y, 0.25 * (l + 1)));
            b.Add(Site(n, l, PauliLetter.Z, 0.125 * (l + 2)));
            b.Add(Site(n, l, PauliLetter.Y, -0.75 + l));
        }
        // A two-letter string, so the weight exponent is exercised past 1 in both conventions.
        var yz = new PauliLetter[n];
        for (int l = 0; l < n; l++) yz[l] = PauliLetter.I;
        yz[0] = PauliLetter.Y;
        yz[1] = PauliLetter.Z;
        a.Add(new PauliTerm(yz, new Complex(1.5, 0.0)));
        b.Add(new PauliTerm(yz, new Complex(-2.25, 0.0)));

        double z = PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, n, PauliLetter.Z);
        double y = PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, n, PauliLetter.Y);

        Assert.NotEqual(0.0, z);
        // Both routes are sums of the same dyadic products in the same order with flipped signs,
        // so this is an exact identity and there is nothing to tolerate.
        Assert.Equal(-z, y);
    }

    // ============================================================
    // The parent edge: F113 is the Z-diagonal special case
    // ============================================================

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void PredictAsymmetry_ReducesToF113_UnderTheSubstitutionCorollaryTenNames(int n)
    {
        // Corollary 10: a_{Z_l} = ω_l/2 and b_{Z_l} = (γ_T1,l − γ_pump,l)/4 turns F155's sum
        // into F113's (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l). The prefactors 4^(N+1) and 4^N/2
        // differ by exactly the factor 8 that substitution carries.
        var omegas = new double[n];
        var gammaT1 = new double[n];
        var gammaPump = new double[n];
        for (int l = 0; l < n; l++)
        {
            omegas[l] = 0.25 + 0.5 * l;
            gammaT1[l] = 0.125 * (l + 1);
            gammaPump[l] = 0.0625 * (l + 2);
        }

        var aTerms = new List<PauliTerm>();
        var bTerms = new List<PauliTerm>();
        for (int l = 0; l < n; l++)
        {
            aTerms.Add(Site(n, l, PauliLetter.Z, omegas[l] / 2.0));
            bTerms.Add(Site(n, l, PauliLetter.Z, (gammaT1[l] - gammaPump[l]) / 4.0));
        }

        double f113 = LindbladBitBPiBreakMagnitude.PredictAsymmetry(omegas, gammaT1, gammaPump, n);
        double f155 = PhysicalGeneratorPolarityBreak.PredictAsymmetry(aTerms, bTerms, n);

        // Both routes are exact sums of dyadic rationals here (the rates are chosen so), so the
        // two floats are equal bit for bit and there is nothing to tolerate.
        Assert.Equal(f113, f155);
    }

    [Fact]
    public void PredictAsymmetry_DetailedBalanceCancels_TheSameWayF113Does()
    {
        // γ_T1 = γ_pump site by site makes every b coefficient exactly zero.
        var a = new[] { Site(2, 0, PauliLetter.Z, 0.5), Site(2, 1, PauliLetter.Z, 0.25) };
        var b = new[] { Site(2, 0, PauliLetter.Z, 0.0), Site(2, 1, PauliLetter.Z, 0.0) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
    }

    [Fact]
    public void PredictAsymmetry_UniformDriveOnEveryGainLossSite_CancelsAcrossTheChain()
    {
        // The registry's own worked example: a gain end and a loss end with the SAME drive on
        // every site cancel, since one carries +ωγ and the other −ωγ. The label does not decide.
        var a = new[] { Site(2, 0, PauliLetter.Z, 0.2), Site(2, 1, PauliLetter.Z, 0.2) };
        var b = new[] { Site(2, 0, PauliLetter.Z, 0.5), Site(2, 1, PauliLetter.Z, -0.5) };

        Assert.Equal(0.0, PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
    }

    // ============================================================
    // Guards
    // ============================================================

    [Fact]
    public void PredictAsymmetry_RefusesAComplexCoefficient()
    {
        var a = new[] { PauliTerm.SingleSite(2, 0, PauliLetter.Z, new Complex(1.0, 0.5)) };
        var b = new[] { Site(2, 0, PauliLetter.Z, 1.0) };

        var ex = Assert.Throws<ArgumentException>(
            () => PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
        Assert.Contains("Hermitian", ex.Message);
    }

    [Fact]
    public void PredictAsymmetry_RefusesALetterCountMismatch()
    {
        var a = new[] { Site(3, 0, PauliLetter.Z, 1.0) };
        var b = new[] { Site(3, 0, PauliLetter.Z, 1.0) };

        Assert.Throws<ArgumentException>(
            () => PhysicalGeneratorPolarityBreak.PredictAsymmetry(a, b, 2));
    }

    [Fact]
    public void PredictAsymmetry_RefusesTheIdentityAsADephaseLetter() =>
        Assert.Throws<ArgumentException>(() => PhysicalGeneratorPolarityBreak.PredictAsymmetry(
            Array.Empty<PauliTerm>(), Array.Empty<PauliTerm>(), 2, PauliLetter.I));

    [Fact]
    public void PredictAsymmetry_RefusesNBelowOne() =>
        Assert.Throws<ArgumentOutOfRangeException>(
            () => PhysicalGeneratorPolarityBreak.PredictAsymmetry(
                Array.Empty<PauliTerm>(), Array.Empty<PauliTerm>(), 0));
}
