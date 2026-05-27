using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F81;
using RCPsiSquared.Diagnostics.Polarity;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.F81;

/// <summary>Welle 15 Task A: verifies the dephaseLetter parameter on
/// PiDecomposition.Decompose and the Klein-V₄ equivariance D · Π_Z · D = Π_Y
/// (Welle 12) that makes Y-dephase results bit-exact D-conjugates of Z-dephase
/// results. Also covers the PolarityCoordinates first overload threading.</summary>
public class PiDecompositionDephaseLetterTests
{
    [Fact]
    public void Default_DephaseLetter_Is_Z()
    {
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z) };

        var dDefault = PiDecomposition.Decompose(chain, terms);
        var dExplicitZ = PiDecomposition.Decompose(chain, terms, dephaseLetter: PauliLetter.Z);

        Assert.True((dDefault.M - dExplicitZ.M).FrobeniusNorm() < 1e-15);
        Assert.True((dDefault.MSym - dExplicitZ.MSym).FrobeniusNorm() < 1e-15);
        Assert.True((dDefault.MAnti - dExplicitZ.MAnti).FrobeniusNorm() < 1e-15);
    }

    [Fact]
    public void Y_Dephase_Result_Equals_D_Conjugate_Of_Z_Result_At_N_Equals_2()
    {
        // Welle 12 (PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N): D · Π_Z · D = Π_Y bit-exact.
        // Therefore M, M_sym, M_anti under Π_Y equal D · (Z-dephase counterpart) · D.
        // ZZ bond is bit_b-homogeneous (Z has bit_b = 1) and Π²-even (truly).
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z) };

        var dZ = PiDecomposition.Decompose(chain, terms, dephaseLetter: PauliLetter.Z);
        var dY = PiDecomposition.Decompose(chain, terms, dephaseLetter: PauliLetter.Y);
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(chain.N);

        var dConjMAnti = D * dZ.MAnti * D;
        var dConjMSym = D * dZ.MSym * D;
        Assert.True((dY.MAnti - dConjMAnti).FrobeniusNorm() < 1e-10,
            $"‖M_anti_Y − D·M_anti_Z·D‖_F = {(dY.MAnti - dConjMAnti).FrobeniusNorm():E3}");
        Assert.True((dY.MSym - dConjMSym).FrobeniusNorm() < 1e-10,
            $"‖M_sym_Y − D·M_sym_Z·D‖_F = {(dY.MSym - dConjMSym).FrobeniusNorm():E3}");
    }

    [Fact]
    public void Y_Dephase_Equivariance_At_N_Equals_3()
    {
        // Same identity at N=3 with two ZZ bonds. Slight tolerance relaxation for the
        // larger matrices (4^3 = 64).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z) };

        var dZ = PiDecomposition.Decompose(chain, terms, dephaseLetter: PauliLetter.Z);
        var dY = PiDecomposition.Decompose(chain, terms, dephaseLetter: PauliLetter.Y);
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(chain.N);

        Assert.True((dY.MAnti - D * dZ.MAnti * D).FrobeniusNorm() < 1e-9);
        Assert.True((dY.MSym - D * dZ.MSym * D).FrobeniusNorm() < 1e-9);
    }

    [Fact]
    public void PolarityCoordinates_First_Overload_Honors_DephaseLetter()
    {
        // ZZ bond, Hermitian H, bit_b-homogeneous c: F112 Z-dephase predicts Asymmetry = 0;
        // Klein-V₄ inherits Asymmetry = 0 on the Y-dephase via D-conjugation (Frobenius
        // norms are D-conjugation invariant).
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z) };

        var polZ = PolarityCoordinates.Decompose(chain, terms, dephaseLetter: PauliLetter.Z);
        var polY = PolarityCoordinates.Decompose(chain, terms, dephaseLetter: PauliLetter.Y);

        Assert.True(Math.Abs(polZ.Asymmetry) < 1e-12,
            $"Z-dephase asymmetry {polZ.Asymmetry:E3} should be zero (F112 Tier1Derived)");
        Assert.True(Math.Abs(polY.Asymmetry) < 1e-12,
            $"Y-dephase asymmetry {polY.Asymmetry:E3} should be zero (Klein-V₄ inheritance)");
    }

    [Fact]
    public void PolarityCoordinates_First_Overload_DephaseLetter_X_Asymmetry_Zero_For_Bit_A_Homogeneous_H_At_N_Equals_2()
    {
        // First C# test of the F112-X cross-dephase scope (Welle 13 Route 1).
        // XX bond is bit_a-homogeneous (X has bit_a = 1) and Hermitian; X-dephasing
        // gives bit_a-homogeneous c. Predicts Asymmetry = 0 bit-exact.
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.1);
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X) };

        var pol = PolarityCoordinates.Decompose(chain, terms, dephaseLetter: PauliLetter.X);
        Assert.True(pol.MNormSquared > 1e-6,
            "test design: M should be non-trivial so the asymmetry assertion is substantive");
        Assert.True(Math.Abs(pol.Asymmetry) < 1e-12,
            $"F112-X in-scope asymmetry {pol.Asymmetry:E3} exceeds 1e-12");
    }
}
