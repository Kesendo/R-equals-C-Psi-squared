using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F87 Π²-class trichotomy: every Pauli-pair (or k-body) Hamiltonian on a chain
/// under single-letter dephasing is exactly one of {truly, soft, hard}.
///
/// <list type="bullet">
///   <item><b>truly</b>: ‖M‖ = 0 (palindrome residual vanishes; F1 holds bit-exactly)</item>
///   <item><b>soft</b>: ‖M‖ &gt; 0 BUT the L spectrum still pairs each λ with −λ−2σ within tolerance</item>
///   <item><b>hard</b>: spectrum pairing also fails, full Π-symmetry-broken</item>
/// </list>
///
/// Hardware-confirmed at Marrakesh 2026-04-26 (Δ(soft − truly) = −0.722 measured vs −0.723
/// Trotter-n3 prediction). See <see cref="Core.Confirmations.ConfirmationsRegistry"/> entry
/// "palindrome_trichotomy".
///
/// Dissipator-resonance law (verified 2026-05-01): F87-hardness lives in the Klein cell
/// matching the dephase letter's Klein index: Z = (0, 1), X = (1, 0), Y = (1, 1).
/// SU(2)-rotation-equivalent.
///
/// γ-universality of the verdict: the classification is evaluated at the chain's single
/// operating γ (<see cref="ChainSystem.GammaZero"/>, uniform across sites). Since 2026-06-10
/// a Hard verdict for a windowed diagonal-cell pair is γ-universal, hard at one γ is hard at
/// every γ &gt; 0 (the windowed all-γ theorem, <see cref="WindowedConverseAllGammaClaim"/>,
/// Pascal-Gram positivity F117, no residual), so the choice of γ does not affect the verdict
/// for that scope. Truly and Soft are γ-independent statements of the operator identity
/// (‖M‖ = 0) and the palindrome spectrum pairing by construction.
///
/// See docs/ANALYTICAL_FORMULAS.md F87 entry for the structural derivation.
/// </summary>
public static class PauliPairTrichotomy
{
    public static TrichotomyClass Classify(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms,
        double opTolerance = 1e-10, double spectrumTolerance = 1e-6,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (terms.Count == 0) return TrichotomyClass.Truly;

        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        var L = PauliDephasingDissipator.Build(H, gammaList, dephaseLetter);
        double sigma = chain.SigmaGamma;

        var M = PalindromeResidual.Build(L, chain.N, sigma, dephaseLetter);
        double opNorm = M.FrobeniusNorm();
        if (opNorm < opTolerance) return TrichotomyClass.Truly;

        var evd = L.Evd();
        var evals = evd.EigenValues.ToArray();
        if (SpectrumPairs(evals, sigma, spectrumTolerance))
            return TrichotomyClass.Soft;
        return TrichotomyClass.Hard;
    }

    /// <summary>k-body template overload (homogeneous k, no mixing): accepts a list
    /// of <see cref="PauliTerm"/> templates (each with Letters of length k ≤ N and a
    /// coefficient) and classifies the resulting chain Hamiltonian using the
    /// sliding-window k-body builder (<see cref="PauliKBodyChainExtensions.ChainKBody"/>).
    /// Parallels Python framework's <c>classify_pauli_pair</c> k-body dispatch path.
    ///
    /// <para>Mixed-body (k=2 + k≥3 in the same call) is NOT supported here; use
    /// only k=2 templates via this overload to mirror the existing
    /// <see cref="PauliPairBondTerm"/> overload, or only k≥3 templates for the
    /// F104 use case. Mixing requires two separate Classify calls and external
    /// H-sum; out of scope for F104.</para>
    ///
    /// <para>Pipeline identical to the k=2 overload: build H via ChainKBody,
    /// build L via dephasing dissipator, build M via palindrome residual, then
    /// truly-test (‖M‖_F &lt; opTolerance) followed by soft/hard test
    /// (greedy multiset eigenvalue pairing).</para>
    ///
    /// <para>Evaluated at the chain's single operating γ; for a windowed diagonal-cell
    /// pair a Hard verdict is γ-universal (hard at one γ is hard at every γ &gt; 0,
    /// <see cref="WindowedConverseAllGammaClaim"/>), so the choice of γ does not affect
    /// the verdict in that scope.</para></summary>
    public static TrichotomyClass Classify(ChainSystem chain, IReadOnlyList<PauliTerm> termTemplates,
        double opTolerance = 1e-10, double spectrumTolerance = 1e-6,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (termTemplates.Count == 0) return TrichotomyClass.Truly;

        var H = termTemplates.ChainKBody(chain.N);
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        var L = PauliDephasingDissipator.Build(H, gammaList, dephaseLetter);
        double sigma = chain.SigmaGamma;

        var M = PalindromeResidual.Build(L, chain.N, sigma, dephaseLetter);
        double opNorm = M.FrobeniusNorm();
        if (opNorm < opTolerance) return TrichotomyClass.Truly;

        var evd = L.Evd();
        var evals = evd.EigenValues.ToArray();
        if (SpectrumPairs(evals, sigma, spectrumTolerance))
            return TrichotomyClass.Soft;
        return TrichotomyClass.Hard;
    }

    /// <summary>Greedy pairing: for each unused λ, find the closest unused λ' to the target
    /// −λ−2σ. A self-pair (λ ≈ −σ, fixed point of the involution) consumes only λ.
    /// Returns true if every pairing distance is within <paramref name="tolerance"/>.</summary>
    private static bool SpectrumPairs(Complex[] evals, double sigma, double tolerance)
    {
        var used = new bool[evals.Length];
        double maxErr = 0;
        for (int i = 0; i < evals.Length; i++)
        {
            if (used[i]) continue;
            Complex target = -evals[i] - 2 * sigma;
            int bestJ = -1;
            double bestDist = double.PositiveInfinity;
            for (int j = 0; j < evals.Length; j++)
            {
                if (used[j]) continue;
                double dist = (evals[j] - target).Magnitude;
                if (dist < bestDist) { bestDist = dist; bestJ = j; }
            }
            if (bestJ < 0) return false;
            used[i] = true;
            if (bestJ != i) used[bestJ] = true;
            if (bestDist > maxErr) maxErr = bestDist;
        }
        return maxErr < tolerance;
    }
}
