using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F81;

/// <summary>F81 Π-decomposition of M into symmetric/antisymmetric parts. Tier-1 result
/// (PROOF_F81_PI_CONJUGATION_OF_M.md):
///
/// <code>
///   Π · M · Π⁻¹ = M − 2 · L_{H_odd}      (pure Z-dephasing)
///
///   M_sym  = (M + Π·M·Π⁻¹) / 2
///   M_anti = (M − Π·M·Π⁻¹) / 2 = L_{H_odd}    (when F81 holds)
/// </code>
///
/// <para>For Π²-even-only H: L_{H_odd} = 0 and M is its own Π-conjugate. For Π²-odd H:
/// M_anti exactly equals the unitary commutator L_{H_odd} = −i[H_odd, ·]. M_sym and M_anti
/// are Frobenius-orthogonal: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².</para>
///
/// <para>For T1 amplitude damping, F81 is violated by the cooling/heating asymmetry:
/// f81_violation = √(Σ (γ_↓ − γ_↑)²) · 2^(N−1) per F84. With detailed balance
/// (γ_↓ = γ_↑) the violation vanishes even though both channels are active.</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F81 + F84 entries.</para>
/// </summary>
public sealed record PiDecompositionResult(
    ComplexMatrix M,
    ComplexMatrix MSym,
    ComplexMatrix MAnti,
    ComplexMatrix LHOdd,
    double F81Violation,
    double MNormSquared,
    double MSymNormSquared,
    double MAntiNormSquared,
    double LHOddNormSquared);

public static class PiDecomposition
{
    public static PiDecompositionResult Decompose(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms,
        IReadOnlyList<double>? gammaT1PerSite = null, double violationTolerance = 1e-7)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaZ = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();

        ComplexMatrix L;
        if (gammaT1PerSite is not null && gammaT1PerSite.Any(g => g != 0))
            L = T1Dissipator.Build(H, gammaZ, gammaT1PerSite);
        else
            L = PauliDephasingDissipator.BuildZ(H, gammaZ);

        var M = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, PauliLetter.Z);
        var pi = PiOperator.BuildFull(chain.N, PauliLetter.Z);
        var piInv = pi.ConjugateTranspose();
        var piMpi = pi * M * piInv;
        var mSym = (M + piMpi) / 2.0;
        var mAnti = (M - piMpi) / 2.0;

        var oddTerms = terms
            .Where(t => !t.IsTruly && t.LetterA != PauliLetter.I && t.LetterB != PauliLetter.I && t.Pi2Parity == 1)
            .ToList();

        ComplexMatrix lHOdd;
        if (oddTerms.Count > 0)
        {
            var hOdd = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, oddTerms.ToBilinearSpec(chain.J)).ToMatrix();
            lHOdd = CommutatorSuperoperatorInPauliBasis(hOdd, chain.N);
        }
        else
        {
            long d2 = 1L << (2 * chain.N);
            lHOdd = Matrix<Complex>.Build.Sparse((int)d2, (int)d2);
        }

        double f81Violation = (mAnti - lHOdd).FrobeniusNorm();
        return new PiDecompositionResult(
            M, mSym, mAnti, lHOdd,
            F81Violation: f81Violation,
            MNormSquared: Math.Pow(M.FrobeniusNorm(), 2),
            MSymNormSquared: Math.Pow(mSym.FrobeniusNorm(), 2),
            MAntiNormSquared: Math.Pow(mAnti.FrobeniusNorm(), 2),
            LHOddNormSquared: Math.Pow(lHOdd.FrobeniusNorm(), 2));
    }

    /// <summary>L_{H} = −i[H, ·] in the 4^N Pauli-string basis. Uses the cached vec_F transform
    /// from <see cref="PauliBasis.VecToPauliBasisTransform"/>.</summary>
    private static ComplexMatrix CommutatorSuperoperatorInPauliBasis(ComplexMatrix H, int N)
    {
        int d = 1 << N;
        var I = Matrix<Complex>.Build.DenseIdentity(d);
        var lVec = -Complex.ImaginaryOne * (H.KroneckerProduct(I) - I.KroneckerProduct(H.Transpose()));
        var transform = PauliBasis.VecToPauliBasisTransform(N);
        return (transform.ConjugateTranspose() * lVec * transform) / Math.Pow(2, N);
    }
}
