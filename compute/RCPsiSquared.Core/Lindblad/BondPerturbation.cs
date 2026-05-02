using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>The "dynamics-of-dynamics" superoperator V_L^b = ∂L/∂J_b at a single bond.
///
/// Under L_total = L_A + δJ · V_L^b the Liouvillian L_A picks up a δJ-strength bond-(b, b+1)
/// Hamiltonian perturbation. V_L^b is the commutator superoperator with the bond's
/// unit-coupling Hamiltonian — perturbing L by adding a J-modulated bond is equivalent to
/// adding δJ·V_L^b.
///
/// Used by the PTF (Perspectival Time Field) workflow to compute first-order eigenvector
/// mixing of slow modes under a bond defect.
/// </summary>
public static class BondPerturbation
{
    public enum Kind
    {
        /// <summary>½(XX + YY) — XY hopping bilinear (popcount-conserving).</summary>
        XY,
        /// <summary>XX + YY + ZZ — Heisenberg.</summary>
        Heisenberg,
        /// <summary>XX only.</summary>
        XX,
        /// <summary>YY only.</summary>
        YY,
        /// <summary>ZZ only.</summary>
        ZZ,
    }

    public static ComplexMatrix Build(int N, int siteA, int siteB, Kind kind)
    {
        if (siteA == siteB || siteA < 0 || siteA >= N || siteB < 0 || siteB >= N)
            throw new ArgumentException($"bond ({siteA}, {siteB}) invalid for N={N}");

        var bonds = new[] { new Bond(siteA, siteB, 1.0) };
        var terms = TermsFor(kind);
        var Hpert = PauliHamiltonian.Bilinear(N, bonds, terms).ToMatrix();

        int d = 1 << N;
        var I = Matrix<Complex>.Build.DenseIdentity(d);
        return -Complex.ImaginaryOne * (Hpert.KroneckerProduct(I) - I.KroneckerProduct(Hpert.Transpose()));
    }

    private static IReadOnlyList<(PauliLetter, PauliLetter, Complex)> TermsFor(Kind kind) => kind switch
    {
        Kind.XY => new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, 0.5),
            (PauliLetter.Y, PauliLetter.Y, 0.5),
        },
        Kind.Heisenberg => new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, Complex.One),
            (PauliLetter.Y, PauliLetter.Y, Complex.One),
            (PauliLetter.Z, PauliLetter.Z, Complex.One),
        },
        Kind.XX => new (PauliLetter, PauliLetter, Complex)[] { (PauliLetter.X, PauliLetter.X, Complex.One) },
        Kind.YY => new (PauliLetter, PauliLetter, Complex)[] { (PauliLetter.Y, PauliLetter.Y, Complex.One) },
        Kind.ZZ => new (PauliLetter, PauliLetter, Complex)[] { (PauliLetter.Z, PauliLetter.Z, Complex.One) },
        _ => throw new ArgumentOutOfRangeException(nameof(kind)),
    };
}
