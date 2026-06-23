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

    [Theory]
    [InlineData("XX", PauliLetter.X, PauliLetter.X)]
    [InlineData("XY", PauliLetter.X, PauliLetter.Y)]
    [InlineData("YX", PauliLetter.Y, PauliLetter.X)]
    [InlineData("YY", PauliLetter.Y, PauliLetter.Y)]
    [InlineData("YZ", PauliLetter.Y, PauliLetter.Z)]
    [InlineData("ZY", PauliLetter.Z, PauliLetter.Y)]
    [InlineData("ZZ", PauliLetter.Z, PauliLetter.Z)]
    [InlineData("XZ", PauliLetter.X, PauliLetter.Z)]
    [InlineData("ZX", PauliLetter.Z, PauliLetter.X)]
    public void Decompose_F112Scope_AsymmetryIsZero_AcrossBilinearFamilies(
        string label, PauliLetter a, PauliLetter b)
    {
        // Every bilinear bond Hamiltonian with real coefficient is Hermitian; single-Pauli
        // Z-dephasing is bit_b-homogeneous (Z has bit_b=1, the only c_op). F112 predicts
        // Asymmetry = 0 bit-exact across all such H, irrespective of trichotomy class.
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(a, b) };
        var pol = PolarityCoordinates.Decompose(chain, terms);
        Assert.True(Math.Abs(pol.Asymmetry) < 1e-12,
            $"F112 in-scope asymmetry for {label} = {pol.Asymmetry:E3} exceeds 1e-12");
    }

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

        double rel = Math.Abs(pol.Asymmetry) / Math.Max(pol.MNormSquared, 1e-15);
        Assert.True(rel > 1e-6,
            $"Z-drive + T1 should break F112 balance; got rel asym = {rel:E3} (no break)");
    }
}
