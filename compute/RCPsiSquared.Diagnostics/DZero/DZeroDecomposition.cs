using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.DZero;

/// <summary>Decompose ρ into d=0 (kernel of L) + d=2 (orthogonal) parts.
///
/// <para>ρ_d0 = lim_{t→∞} e^{Lt} ρ — the steady-state projection (what survives forever),
/// living in the kernel-of-L subspace built from <see cref="StationaryModes"/>.
/// ρ_d2 = ρ − ρ_d0 — the part that decoheres.</para>
///
/// <para>The d=2 framework reading: d²−2d=0 forces d=0 or d=2; d=2 is the qubit dimension,
/// d=0 is the substrate axis the qubit sits on. See <b>F4</b> in
/// <c>docs/ANALYTICAL_FORMULAS.md</c> (sector projectors as kernel modes) and
/// <c>docs/THE_POLARITY_LAYER.md</c> for the +0/0/−0 substrate reading.</para>
/// </summary>
public sealed record DZeroDecompositionResult(
    ComplexMatrix RhoD0,
    ComplexMatrix RhoD2,
    double D0Weight,
    double D2Norm,
    int KernelDimension);

public static class DZeroDecomposition
{
    /// <summary>Build the stationary modes fresh and decompose ρ.</summary>
    public static DZeroDecompositionResult Decompose(ComplexMatrix rho, ChainSystem chain, double tolerance = 1e-9) =>
        Decompose(rho, StationaryModes.Compute(chain, tolerance: tolerance));

    /// <summary>Decompose ρ using a precomputed <see cref="StationaryModesResult"/>. Avoids
    /// re-running the eigendecomposition of L when the caller already has one.</summary>
    public static DZeroDecompositionResult Decompose(ComplexMatrix rho, StationaryModesResult sm)
    {
        // Vec_F (column-major) flatten of rho
        int d = rho.RowCount;
        var rhoVec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                rhoVec[j * d + i] = rho[i, j];

        // Kernel projector P = R_kernel · W_kernel (biorthogonal)
        var pKernel = sm.RightEigenvectors * sm.LeftCovectors;
        var rhoD0Vec = pKernel * rhoVec;

        var rhoD0 = Matrix<Complex>.Build.Dense(d, d);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                rhoD0[i, j] = rhoD0Vec[j * d + i];
        rhoD0 = (rhoD0 + rhoD0.ConjugateTranspose()) / 2.0;

        var rhoD2 = rho - rhoD0;
        return new DZeroDecompositionResult(
            RhoD0: rhoD0,
            RhoD2: rhoD2,
            D0Weight: rhoD0.Trace().Real,
            D2Norm: rhoD2.FrobeniusNorm(),
            KernelDimension: sm.KernelDimension);
    }
}
