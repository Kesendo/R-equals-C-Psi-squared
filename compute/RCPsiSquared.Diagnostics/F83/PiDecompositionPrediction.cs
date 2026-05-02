using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F83;

/// <summary>F83 closed-form Π-decomposition of M: predict ‖M‖², ‖M_anti‖², ‖M_sym‖² and the
/// anti-fraction from the Hamiltonian's per-Π²-class Frobenius norms — without computing M.
///
/// <code>
///   ‖M‖²_F        = 4·‖H_odd‖²·2^N + 8·‖H_even_nontruly‖²·2^N
///   ‖M_anti‖²_F  = 2·‖H_odd‖²·2^N
///   ‖M_sym‖²_F   = 2·‖H_odd‖²·2^N + 8·‖H_even_nontruly‖²·2^N
///   anti-fraction = ‖M_anti‖² / ‖M‖² = 1 / (2 + 4·r)   with r = ‖H_even_nontruly‖²/‖H_odd‖²
/// </code>
///
/// <para>Closed-form anti-fractions: r=0 (pure Π²-odd) → 1/2 (F81 50/50); r=∞ (pure Π²-even
/// non-truly) → 0 (F81 100/0); r=1 → 1/6.</para>
///
/// <para>Hardware-confirmed at Marrakesh 2026-04-30 (Confirmations entry
/// "f83_pi2_class_signature_marrakesh"). See PROOF_F83_PI_DECOMPOSITION_RATIO and
/// PROOF_F85_KBODY_GENERALIZATION; docs/ANALYTICAL_FORMULAS.md F83 + F85 entries.</para>
/// </summary>
public sealed record PiDecompositionForecast(
    double MSquared,
    double MAntiSquared,
    double MSymSquared,
    double AntiFraction,
    double HOddSquared,
    double HEvenNonTrulySquared,
    double R);

public static class PiDecompositionPrediction
{
    public static PiDecompositionForecast Predict(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms)
    {
        var oddTerms = new List<PauliPairBondTerm>();
        var evenTerms = new List<PauliPairBondTerm>();
        foreach (var t in terms)
        {
            if (t.IsTruly) continue;
            if (t.LetterA == PauliLetter.I || t.LetterB == PauliLetter.I) continue; // F83 scope: 2-body, no identity
            (t.Pi2Parity == 1 ? oddTerms : evenTerms).Add(t);
        }

        double hOddSq = FrobeniusSquared(BuildOrZero(chain, oddTerms));
        double hEvenSq = FrobeniusSquared(BuildOrZero(chain, evenTerms));

        int dPow = 1 << chain.N;
        double mSq = 4 * hOddSq * dPow + 8 * hEvenSq * dPow;
        double mAntiSq = 2 * hOddSq * dPow;
        double mSymSq = mSq - mAntiSq;

        double r;
        double antiFraction;
        if (hOddSq < 1e-15)
        {
            r = double.PositiveInfinity;
            antiFraction = 0.0;
        }
        else
        {
            r = hEvenSq / hOddSq;
            antiFraction = mSq > 0 ? mAntiSq / mSq : 0.0;
        }

        return new PiDecompositionForecast(mSq, mAntiSq, mSymSq, antiFraction, hOddSq, hEvenSq, r);
    }

    /// <summary>Anti-fraction convenience wrapper: ‖M_anti‖² / ‖M‖² = 1 / (2 + 4·r).</summary>
    public static double AntiFraction(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms) =>
        Predict(chain, terms).AntiFraction;

    private static MathNet.Numerics.LinearAlgebra.Matrix<Complex>? BuildOrZero(ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> terms) =>
        terms.Count == 0 ? null : PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();

    private static double FrobeniusSquared(MathNet.Numerics.LinearAlgebra.Matrix<Complex>? m) =>
        m is null ? 0.0 : (m.ConjugateTranspose() * m).Trace().Real;
}
