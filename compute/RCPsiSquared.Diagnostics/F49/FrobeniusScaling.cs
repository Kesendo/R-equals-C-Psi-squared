using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F49;

/// <summary>F49 / F85: closed-form ‖M‖²_F = c_H · F(N, G) for the F1-palindrome residual,
/// generalized to k-body terms via Π²-class.
///
/// <para><b>F85 generalisation</b> (PROOF_F85_KBODY_GENERALIZATION.md):</para>
/// <code>
///   ‖M(L_Z)‖²    = Σ_k 4·c(k)·‖H_k‖²_F·2^N
///   ‖M(L_Z+T1)‖² = ‖M(L_Z)‖² + 4^(N−1) · [3·Σγ_T1² + 4·(Σγ_T1)²]
/// </code>
/// <para>where c(k) is the Π²-class factor:</para>
/// <list type="bullet">
///   <item>truly term (k-body criterion: #Y even AND #Z even) → c = 0 (M = 0 by Master Lemma)</item>
///   <item>Π²-odd non-truly → c = 1 (factor 4·2^N)</item>
///   <item>Π²-even non-truly → c = 2 (factor 8·2^N)</item>
/// </list>
/// <para>T1 contribution is Hamiltonian-independent, γ_Z-independent, M-orthogonal.
/// Verified at N=3..6 for arbitrary {γ_T1_l}.</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F49 + F85 entries.</para>
/// </summary>
public static class FrobeniusScaling
{
    /// <summary>‖M‖² = c_H · F(N, G) where F is the chain factor and c_H is the
    /// Hamiltonian-class constant supplied by the caller.</summary>
    public static double PredictNormSquared(ChainSystem chain, double cH, HamiltonianClass cls = HamiltonianClass.Main)
    {
        return cH * PalindromeResidualScaling.FactorChain(chain.N, cls);
    }

    /// <summary>‖M‖² from terms via the F85 Π²-class identity. Drops truly terms; sums
    /// per-Π²-class sub-Hamiltonian Frobenius norms with c-factor weights.</summary>
    public static double PredictNormSquaredFromTerms(ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> terms,
        IReadOnlyList<double>? gammaT1PerSite = null)
    {
        double zPart = 0;
        foreach (int targetParity in new[] { 1, 0 }) // 1 = Π²-odd (c=1), 0 = Π²-even non-truly (c=2)
        {
            var groupTerms = terms.Where(t => !t.IsTruly && t.Pi2Parity == targetParity).ToList();
            if (groupTerms.Count == 0) continue;
            var Hgroup = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, groupTerms.ToBilinearSpec(chain.J)).ToMatrix();
            double frobSq = (Hgroup.ConjugateTranspose() * Hgroup).Trace().Real;
            int cFactor = targetParity == 1 ? 1 : 2;
            zPart += 4 * (1 << chain.N) * cFactor * frobSq;
        }

        double t1Part = 0;
        if (gammaT1PerSite is not null)
        {
            if (gammaT1PerSite.Count != chain.N)
                throw new ArgumentException($"gamma_t1 list length {gammaT1PerSite.Count} != N {chain.N}");
            double sumG = gammaT1PerSite.Sum();
            double sumG2 = gammaT1PerSite.Sum(g => g * g);
            t1Part = Math.Pow(4, chain.N - 1) * (3 * sumG2 + 4 * sumG * sumG);
        }
        return zPart + t1Part;
    }
}
