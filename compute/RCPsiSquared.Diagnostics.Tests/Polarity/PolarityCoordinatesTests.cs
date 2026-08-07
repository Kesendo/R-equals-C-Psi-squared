using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F81;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

public class PolarityCoordinatesTests
{
    private static ChainSystem MakeChain(int N = 3, double gammaZero = 0.05) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: gammaZero,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void Decompose_OrthogonalityInvariant_HoldsBitExact()
    {
        // ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖² bit-exactly.
        var chain = MakeChain();
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
        };
        var pol = PolarityCoordinates.Decompose(chain, terms);

        double sum = pol.MZeroNormSquared + pol.MPlusHalfNormSquared + pol.MMinusHalfNormSquared;
        Assert.Equal(pol.MNormSquared, sum, precision: 9);
        Assert.True(pol.OrthogonalityResidual < 1e-9,
            $"orthogonality residual {pol.OrthogonalityResidual:E3} exceeds 1e-9");
    }

    [Fact]
    public void MAntiNormSquared_EqualsDirectMAntiNorm_ToAFewUlp()
    {
        // Guards the identity RelativeAsymmetry's denominator rests on:
        // ‖M_plus_half‖² + ‖M_minus_half‖² = ‖M_anti‖². With A = M_anti, B = Π·A·Π⁻¹ the
        // halves are (A ∓ iB)/2, so the parallelogram law gives (‖A‖² + ‖B‖²)/2, and Π
        // unitary gives ‖B‖² = ‖A‖². NOT compared exactly, deliberately: the property sums two
        // Math.Pow(FrobeniusNorm(), 2) values (two roundings each) while the right-hand side
        // is a single entry-wise sum, so the two routes agree only to a few ulp except where
        // both land on 0.0 or a power of two. The identity is exact in exact arithmetic; the
        // float agreement is not, and an == here would be a claim the route cannot support.
        var chain = MakeChain();
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
        };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        var pi = PiDecomposition.Decompose(chain, terms);

        double directMAnti = pi.MAnti.Enumerate().Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary);
        Assert.Equal(directMAnti, pol.MAntiNormSquared, precision: 9);
        Assert.True(pol.MAntiNormSquared > 0.0,
            "XY+YZ has a Π²-odd part, so this fixture must NOT be degenerate; "
            + "a zero here would make the test vacuous");
    }

    [Theory]
    // Both rows are Π²-even, so M_anti = L_{H_odd} = 0 by F81 and there is no polarity
    // content to balance. They differ sharply in ‖M‖², which is the point of the pair.
    [InlineData("truly (Heisenberg): M itself vanishes by the F83 closed form", true)]
    [InlineData("non-truly Π²-even (YZ+ZY): ‖M‖² is large, M_anti still zero", false)]
    public void Decompose_Pi2EvenH_IsPolarityDegenerate_AndMNormIsNotTheDiscriminator(
        string because, bool trulyRow)
    {
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.5,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);
        var terms = trulyRow
            ? new[]
            {
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
                new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
            }
            : new[]
            {
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
                new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
            };
        var pol = PolarityCoordinates.Decompose(chain, terms);

        // Exact here for a different reason than the ulp note above: a vanishing M_anti makes
        // both half-norms exactly 0.0, and 0.0 survives every rounding route, so == is safe.
        Assert.True(pol.MAntiNormSquared == 0.0, $"{because}: got {pol.MAntiNormSquared:E3}");
        Assert.True(pol.IsPolarityDegenerate, because);
        // NaN, not 0.0, since 2026-08-07: with no polarity content the ratio is 0/0 and
        // there is no number to report. A 0.0 here read as "perfectly balanced", the one
        // reading this whole fixture pair exists to deny.
        Assert.True(double.IsNaN(pol.RelativeAsymmetry),
            $"{because}: degenerate ⇒ NO ratio, got {pol.RelativeAsymmetry:E3}");

        // The discriminating half. ‖M‖² is 0 on the truly row and 256 on the non-truly one,
        // yet BOTH are degenerate; so ‖M‖² cannot tell the two states of the balance apart
        // and is the wrong denominator. Under the retired max(‖M‖², 1e-15) the non-truly row
        // divided by 256 and looked like a genuine pass.
        if (trulyRow)
            Assert.True(pol.MNormSquared == 0.0, $"truly ⇒ M ≡ 0; got {pol.MNormSquared:E3}");
        else
            Assert.True(pol.MNormSquared > 1.0,
                $"non-truly Π²-even ⇒ ‖M‖² is large; got {pol.MNormSquared:E3}");
    }

    [Fact]
    public void Decompose_HermitianHHomogeneousC_AsymmetryIsZero()
    {
        // F112 typed scope: Hermitian H + bit_b-homogeneous c (Z-dephasing only) ⇒
        // ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly ⇒ Asymmetry = 0.
        var chain = MakeChain();
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
        };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.True(Math.Abs(pol.Asymmetry) < 1e-12,
            $"F112 in-scope asymmetry {pol.Asymmetry:E3} exceeds 1e-12");
    }

    [Fact]
    public void Decompose_MatchesF81_MsymEqualsMZero()
    {
        // Cross-check: M_zero == F81 M_sym bit-exactly, and M_plus_half + M_minus_half
        // == F81 M_anti bit-exactly.
        var chain = MakeChain();
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        var f81 = PiDecomposition.Decompose(chain, terms);
        var pol = PolarityCoordinates.Decompose(chain, terms);

        Assert.True((pol.MZero - f81.MSym).FrobeniusNorm() < 1e-12);
        var antiSum = pol.MPlusHalf + pol.MMinusHalf;
        Assert.True((antiSum - f81.MAnti).FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Decompose_HeisenbergZerosOut_M_IsZero()
    {
        // Truly H (XX+YY) + pure Z-dephasing: F1 holds exactly ⇒ M = 0 ⇒
        // all three polarity components are zero.
        var chain = MakeChain();
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.True(pol.MNormSquared < 1e-18,
            $"truly M norm² {pol.MNormSquared:E3} should vanish");
        Assert.True(pol.MZeroNormSquared < 1e-18);
        Assert.True(pol.MPlusHalfNormSquared < 1e-18);
        Assert.True(pol.MMinusHalfNormSquared < 1e-18);
        Assert.True(Math.Abs(pol.Asymmetry) < 1e-18);
    }

    /// <summary>The SUBSTANTIVE half of the bilinear sweep: Π²-ODD bonds, where F112's
    /// prediction has something to predict. Each of these carries a genuinely non-zero
    /// polarity content and the asymmetry vanishes on top of it, which is the confirmation.
    /// Split out from a nine-row theory on 2026-08-07: the other five rows are Π²-EVEN, so
    /// their M_anti is zero by construction and their "asymmetry ≈ 0" could not have failed.
    /// They are pinned below as degenerate instead, which is a different statement.</summary>
    [Theory]
    [InlineData("XY", PauliLetter.X, PauliLetter.Y)]
    [InlineData("YX", PauliLetter.Y, PauliLetter.X)]
    [InlineData("XZ", PauliLetter.X, PauliLetter.Z)]
    [InlineData("ZX", PauliLetter.Z, PauliLetter.X)]
    public void Decompose_F112Scope_AsymmetryIsZero_OnPi2OddBilinears(
        string label, PauliLetter a, PauliLetter b)
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(a, b) };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.False(IsSilent(chain, terms),
            $"{label} must carry polarity content for this assertion to mean anything");
        Assert.True(pol.MAntiNormSquared > 1.0,
            $"{label}: polarity content ‖M_anti‖² = {pol.MAntiNormSquared:E3} is not substantive");
        Assert.True(Math.Abs(pol.Asymmetry) < 1e-12,
            $"F112 in-scope asymmetry for {label} = {pol.Asymmetry:E3} exceeds 1e-12");
    }

    /// <summary>The Π²-EVEN bilinears. These are in F112's typed scope too, but they carry NO
    /// polarity content, so the right statement about them is that they are SILENT, not that
    /// they are balanced. Asserting "asymmetry ≈ 0" here is asserting 0 − 0 &lt; 1e-12, which no
    /// change to the physics could break. Note ZZ and YZ/ZY have a large ‖M‖²: this is not the
    /// M-vanishes case, it is the polarity-content-vanishes case, and only the second one is
    /// what the relative asymmetry divides by.</summary>
    [Theory]
    [InlineData("XX", PauliLetter.X, PauliLetter.X)]
    [InlineData("YY", PauliLetter.Y, PauliLetter.Y)]
    [InlineData("YZ", PauliLetter.Y, PauliLetter.Z)]
    [InlineData("ZY", PauliLetter.Z, PauliLetter.Y)]
    [InlineData("ZZ", PauliLetter.Z, PauliLetter.Z)]
    public void Decompose_Pi2EvenBilinears_AreSilent_NotBalanced(
        string label, PauliLetter a, PauliLetter b)
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(a, b) };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.True(IsSilent(chain, terms),
            $"{label} is Π²-even under Z-dephasing and must read as structurally silent");
        Assert.True(pol.MAntiNormSquared < 1e-18,
            $"{label}: ‖M_anti‖² = {pol.MAntiNormSquared:E3} should be structurally zero");
    }

    /// <summary>Structural silence via the production path, not a re-implementation: build the
    /// witness the same way StandardSet does and read its exact predicate.</summary>
    private static bool IsSilent(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> bondTerms) =>
        new LindbladBitBPiBalanceWitness(
            witnessName: "probe", chain: chain, bondTerms: bondTerms,
            gammaT1: null, expectedVerdict: "DEGENERATE").IsStructurallyDegenerate;

    [Fact]
    public void Decompose_F81ViolationMatches_PiDecomposition()
    {
        // Pass-through: PolarityCoordinates pipes through F81Violation from
        // PiDecomposition unchanged (bond-term overload).
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y) };
        var f81 = PiDecomposition.Decompose(chain, terms);
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.Equal(f81.F81Violation, pol.F81Violation, precision: 12);
    }

    [Fact]
    public void Decompose_KBodyOverload_MatchesBilinear_ForSameH()
    {
        // The k-body overload should reproduce the bilinear overload bit-exactly when
        // fed an H term list equivalent to expanding the bond terms across chain.Bonds.
        var chain = MakeChain();
        var bondTerms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
        };
        // Expand to k-body PauliTerm list using the chain's J coupling, exactly as
        // PauliHamiltonian.Bilinear does internally.
        var kBody = new List<PauliTerm>();
        foreach (var bond in chain.Bonds)
            foreach (var t in bondTerms)
                kBody.Add(PauliTerm.TwoSite(chain.N, bond.Site1, t.LetterA,
                                             bond.Site2, t.LetterB, (Complex)chain.J));

        var polBond = PolarityCoordinates.Decompose(chain, bondTerms);
        var polKBody = PolarityCoordinates.Decompose(chain, (IReadOnlyList<PauliTerm>)kBody);

        Assert.Equal(polBond.MNormSquared, polKBody.MNormSquared, precision: 9);
        Assert.Equal(polBond.MZeroNormSquared, polKBody.MZeroNormSquared, precision: 9);
        Assert.Equal(polBond.MPlusHalfNormSquared, polKBody.MPlusHalfNormSquared, precision: 9);
        Assert.Equal(polBond.MMinusHalfNormSquared, polKBody.MMinusHalfNormSquared, precision: 9);
        Assert.Equal(polBond.Asymmetry, polKBody.Asymmetry, precision: 9);
    }

    [Fact]
    public void Decompose_KBodyOverload_SingleSiteZDriveWithT1_BreaksBalance()
    {
        // Welle 2 structural counterexample (2026-05-26): H = ω·(Z₀+Z₁)/2 + σ⁻ T1 with
        // ω = 0.13, γ_T1 = 0.001 yields measurable F112 asymmetry. This is the witness
        // 5 of LindbladBitBPiBalanceWitness.StandardSet (synthetic isolation matching
        // the f95_angle_steering Tier-A dataset fits in
        // simulations/f112_hardware_lens_multi.py).
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.005,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);
        const double omega = 0.13;
        var terms = new List<PauliTerm>
        {
            PauliTerm.SingleSite(chain.N, 0, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
            PauliTerm.SingleSite(chain.N, 1, PauliLetter.Z, coefficient: (Complex)(omega / 2.0)),
        };
        var pol = PolarityCoordinates.Decompose(chain, terms, gammaT1: 0.001);

        // A single-site Z-drive is Π²-odd, so there IS polarity content here and the ratio
        // is a real contrast. Assert that first: without it a "break" could be read off a
        // vacuous 0/0, which is exactly the defect the ‖M‖² denominator used to hide.
        Assert.False(pol.IsPolarityDegenerate,
            "Z-drive is Π²-odd, so M_anti must be non-zero for the ratio to mean anything; "
            + $"got ‖M_anti‖² = {pol.MAntiNormSquared:E3}");

        // The quantity F113 predicts is the ABSOLUTE asymmetry, not the ratio:
        //   asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)
        // with the drive on both sites of the N=2 chain and no pump. Mind the factor 2:
        // F113's ω_l is the LARMOR rate and each term above carries ω/2 as its Pauli
        // coefficient, so ω_l = ω per site and Σ_l ω_l = 2ω, not ω.
        double predicted = (Math.Pow(4, chain.N) / 2.0) * (2.0 * omega) * (0.0 - 0.001);
        Assert.Equal(predicted, pol.Asymmetry, precision: 12);

        // The denominator, checked against an independent route rather than a literal:
        // ‖M_anti‖² recomputed entry-wise from PiDecomposition. This is the assertion that
        // discriminates the current definition from the retired max(‖M‖², 1e-15) one, which
        // on this fixture gives 3.845528e-3 against 7.692080e-3 (the ratio is gamma-independent,
        // so these coincide with the witness fixture's despite its different gamma_0). The two
        // degeneracy rows above cannot tell them apart (both denominators yield 0.0 there),
        // so without this line the guard would be decorative.
        // M_anti = M_plus_half + M_minus_half as MATRICES ((A−iB)/2 + (A+iB)/2 = A), normed
        // afterwards; independent of the sum-of-two-norms the property computes.
        double directMAnti = (pol.MPlusHalf + pol.MMinusHalf).Enumerate()
            .Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary);
        Assert.Equal(Math.Abs(pol.Asymmetry) / directMAnti, pol.RelativeAsymmetry, precision: 12);

        double rel = pol.RelativeAsymmetry;
        Assert.True(rel > 1e-6,
            $"Z-drive + T1 should break F112 balance; got rel asym = {rel:E3} (no break)");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Decompose_SilentAtNAboveTwo_IsDegenerateStructurally_NotByTheFloat(int n)
    {
        // Why this test exists: at N=2 a structurally silent M_anti lands on exact 0.0, so the
        // float test IsPolarityDegenerate was already right there, and until 2026-08-07 every
        // IsPolarityDegenerate assertion in this file sat at that one N, the one that cannot
        // break the claim. From N=3 up the same physics arrives as rounding (~1e-33), the
        // float test reads False, and before 2026-08-07 RelativeAsymmetry returned
        // |noise| / noise: a number with no meaning, free to exceed any threshold, on the
        // very hypothesis F112 proves balanced.
        //
        // N=2 is a row here rather than an assumption about a neighbouring file, and it is
        // held at the SAME γ as the others: the contrast has to be N and nothing else. The
        // old N=2 fixture in this file runs at γ = 0.5, a dyadic, which is an input that
        // cannot break an exactness claim (CLAUDE.md, "no rounding"); γ = 0.1 is not dyadic,
        // so the exact 0.0 asserted below at N=2 is a real property of the cancellation.
        var chain = MakeChain(N: n, gammaZero: 0.1);
        var heisenberg = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
        };
        var pol = PolarityCoordinates.Decompose(chain, heisenberg);

        // The premise, asserted rather than assumed, and it is the whole discriminator: the
        // float route sees it at N=2 and does NOT from N=3 up. If either half ever goes red
        // the test has stopped discriminating and must be re-sited, not relaxed.
        if (n == 2)
        {
            Assert.True(pol.MAntiNormSquared == 0.0,
                $"N=2: the cancellation is exact here, which is why the float test sufficed; "
                + $"got ‖M_anti‖² = {pol.MAntiNormSquared:E3}");
        }
        else
        {
            Assert.True(pol.MAntiNormSquared != 0.0 && pol.MAntiNormSquared < 1e-25,
                $"N={n}: the fixture must sit in the rounding regime the float test misses; "
                + $"got ‖M_anti‖² = {pol.MAntiNormSquared:E3}");
        }

        Assert.True(pol.IsStructurallySilent,
            $"N={n}: Heisenberg is Π²-even and Z-dephasing is homogeneous, so the exact "
            + "predicate must carry through from the term list");
        Assert.True(pol.IsPolarityDegenerate,
            $"N={n}: the structural verdict must govern, not the float ‖M_anti‖² == 0.0");
        Assert.True(double.IsNaN(pol.RelativeAsymmetry),
            $"N={n}: no polarity content ⇒ no ratio; got {pol.RelativeAsymmetry:E3}");
    }
}
