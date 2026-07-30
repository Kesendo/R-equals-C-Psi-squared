using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F49;
using RCPsiSquared.Diagnostics.F81;

namespace RCPsiSquared.Diagnostics.Tests.F81;

public class PiDecompositionTests
{
    [Fact]
    public void TrulyHamiltonian_HasZero_M_And_F81_HoldsTrivially()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.MNormSquared < 1e-18);
        Assert.True(d.F81Violation < 1e-9);
    }

    [Fact]
    public void PurePi2OddHamiltonian_F81_HoldsBitExactly()
    {
        // F81 Tier 1: for pure Π²-odd H (no T1), M_anti = L_{H_odd} bit-exactly.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.F81Violation < 1e-9, $"F81 violation = {d.F81Violation:E3}");
    }

    [Fact]
    public void PurePi2EvenNonTruly_M_anti_IsZero_LHOdd_IsZero()
    {
        // YZ+ZY is Π²-even non-truly. H_odd = 0 (no Π²-odd terms).
        // F81 says M_anti = L_{H_odd} = 0 → M_anti = 0.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.LHOddNormSquared < 1e-18);
        Assert.True(d.MAntiNormSquared < 1e-9, $"M_anti norm² = {d.MAntiNormSquared:E3}");
        Assert.True(d.F81Violation < 1e-9);
    }

    [Fact]
    public void Decomposition_MSymPlusMAnti_EqualsM()
    {
        // M = M_sym + M_anti by construction.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        var sum = d.MSym + d.MAnti;
        var diff = sum - d.M;
        Assert.True(diff.FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Decomposition_NormSquared_Decomposes_Orthogonally()
    {
        // Frobenius orthogonality: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.Equal(d.MNormSquared, d.MSymNormSquared + d.MAntiNormSquared, 8);
    }

    [Fact]
    public void T1_BreaksF81_ViolationMatchesF82ClosedForm()
    {
        // Truly H + T1 → M_anti contains only the T1 contribution → F81 violation = ‖D_{T1, odd}‖_F.
        // F82 closed form: ‖D_{T1, odd}‖_F = √(Σγ_T1²) · 2^(N-1).
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        double[] gammaT1 = { 0.02, 0.02 };
        var d = PiDecomposition.Decompose(chain, terms, gammaT1);

        double expected = Math.Sqrt(gammaT1.Sum(g => g * g)) * Math.Pow(2, chain.N - 1);
        Assert.Equal(expected, d.F81Violation, precision: 9);
    }

    [Fact]
    public void PurePi2Odd_HasFiftyFiftySplit_PerProofStep8()
    {
        // F81 proof Step 8 / memory project_mirror_two_sides: pure Π²-odd H gives
        // ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 exactly. This locks the 50/50 prediction at
        // the Pythagoras level, beyond just orthogonality.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.Equal(d.MNormSquared / 2, d.MSymNormSquared, precision: 9);
        Assert.Equal(d.MNormSquared / 2, d.MAntiNormSquared, precision: 9);
    }

    [Fact]
    public void IdentityBearingOddTerm_EntersHOdd_SoTheIdentityHolds()
    {
        // (Z, I) is Π²-odd exactly like (X, Z): the class is a parity of the Y and Z counts,
        // which an identity leg does not change. Excluding such terms left H_odd short while
        // M_anti was not, so the F81 identity read as violated by ‖M_anti‖ = 16 at N = 3
        // against an empty operator. The Python twin carried the same defect.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.Z, PauliLetter.I) };
        var d = PiDecomposition.Decompose(chain, terms);

        // The discriminating pair: M_anti is built from M and Pi and reads 256.0 either way,
        // so it is the AGREEMENT with L_H_odd that moves (0 before, 256.0 now), and the
        // violation with it (16.0 before). The M_anti value is pinned only to keep the
        // agreement from being satisfied by both sides vanishing.
        Assert.True(d.F81Violation < 1e-9, $"F81 violated by {d.F81Violation}");
        Assert.Equal(256.0, d.MAntiNormSquared, precision: 9);
        Assert.Equal(d.MAntiNormSquared, d.LHOddNormSquared, precision: 9);
    }

    [Fact]
    public void IdentityBearingTrulyTerm_StaysOutOfHOdd()
    {
        // (X, I) has #Y = #Z = 0, so it is truly and contributes no antisymmetric part.
        // The counterpart of the test above: the widening must not swallow truly terms.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.I) };
        var d = PiDecomposition.Decompose(chain, terms);

        // Assert on L_H_odd, not on M_anti: M_anti is built from M and Pi before the odd-term
        // filter runs, so it is 0 here whatever the filter does and cannot witness a widening.
        // L_H_odd is the filter's own output, and it is what an over-widened filter would fill.
        Assert.True(d.LHOddNormSquared < 1e-15, $"L_H_odd = {d.LHOddNormSquared}");
        Assert.True(d.MAntiNormSquared < 1e-15, $"M_anti = {d.MAntiNormSquared}");
    }
}
