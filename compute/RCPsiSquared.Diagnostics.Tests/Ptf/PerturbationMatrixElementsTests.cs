using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.DZero;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

public class PerturbationMatrixElementsTests
{
    [Fact]
    public void MatrixElements_OnKernelBasis_HasExpectedShape()
    {
        // For the unperturbed Z-dephased XY chain L_A, the kernel has dim N+1 (sector projectors).
        // V_L = bond-perturbation. Matrix-element shape: (N+1) × (N+1).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        var vL = BondPerturbation.Build(chain.N, siteA: 0, siteB: 1, BondPerturbation.Kind.XY);

        var elems = PerturbationMatrixElements.Compute(sm.RightEigenvectors, sm.LeftCovectors, vL);
        Assert.Equal(sm.KernelDimension, elems.RowCount);
        Assert.Equal(sm.KernelDimension, elems.ColumnCount);
    }

    [Fact]
    public void MatrixElements_OnKernel_AreZero_ForBondPerturbationWithSameH()
    {
        // The kernel of L_A is sector-projector-spanned. The bond perturbation V_L = -i[H_bond, ·]
        // is a commutator; on the kernel modes (which commute with H = XY), V_L acts as zero.
        // So the kernel-block matrix elements ⟨W | V_L | M⟩ should be small.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        var vL = BondPerturbation.Build(chain.N, siteA: 0, siteB: 1, BondPerturbation.Kind.XY);

        var elems = PerturbationMatrixElements.Compute(sm.RightEigenvectors, sm.LeftCovectors, vL);
        // Sector projectors commute with the XY hopping (popcount-conserving) → matrix elements vanish.
        Assert.True(elems.FrobeniusNorm() < 1e-9,
            $"V_L acts as zero on kernel modes for XY-bond perturbation; got Frobenius {elems.FrobeniusNorm():E3}");
    }
}
