using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Shared helpers for F108 Part 1 and Part 2 closure tests. Both parts run
/// the same operator-level palindrome residual against the dephase-letter-matched
/// Π_5bilinear variant; this helper parameterizes that pipeline on the dephase
/// letter and the bilinear-set builder.</summary>
internal static class F108TestSupport
{
    public const double Gamma = 0.05;
    public const double ResidualTol = 1e-10;

    /// <summary>Build a chain of bilinears: place every <paramref name="bilinears"/>
    /// pair on every (b, b+1) bond.</summary>
    public static IReadOnlyList<PauliTerm> BuildChainFromBilinears(
        int N, IReadOnlyList<(PauliLetter, PauliLetter)> bilinears)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
            foreach (var (a, c) in bilinears)
                terms.Add(PauliTerm.TwoSite(N, b, a, b + 1, c, Complex.One));
        return terms;
    }

    /// <summary>Compute ‖Π_5bilinear · L · Π⁻¹ + L + 2σ·I‖_F in the Pauli-string basis,
    /// where L uses <paramref name="dephaseLetter"/> dephasing on every site and Π is
    /// the matching Π_5bilinear variant. Returns 0 (machine precision) when F108 Part
    /// 1 (Z-deph) or Part 2 (X-deph) holds.</summary>
    public static double ComputeOperatorResidual(
        int N, IReadOnlyList<PauliTerm> terms, PauliLetter dephaseLetter)
    {
        ComplexMatrix H;
        if (terms.Count == 0)
        {
            int d = 1 << N;
            H = Matrix<Complex>.Build.Dense(d, d);
        }
        else
        {
            H = new PauliHamiltonian(N, terms).ToMatrix();
        }

        var gammas = Enumerable.Repeat(Gamma, N).ToArray();
        var Lvec = PauliDephasingDissipator.Build(H, gammas, dephaseLetter);

        var transform = PauliBasis.VecToPauliBasisTransform(N);
        double invD = 1.0 / (1 << N);
        var Lpauli = (transform.ConjugateTranspose() * Lvec * transform) * invD;

        var pi = Pi5BilinearOperator.BuildFull(N, dephaseLetter);
        var piInv = pi.ConjugateTranspose();

        long d2 = 1L << (2 * N);
        double sigma = N * Gamma;
        var residual = pi * Lpauli * piInv + Lpauli +
            (Complex)(2.0 * sigma) * Matrix<Complex>.Build.DenseIdentity((int)d2);
        return residual.FrobeniusNorm();
    }
}
