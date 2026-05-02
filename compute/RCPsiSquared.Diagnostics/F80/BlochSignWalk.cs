using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F80;

/// <summary>F80 structural-identity prediction of M's spectrum for chain Π²-odd 2-body H.
///
/// Theorem F80 (verified bit-exact at N=3..7, see PROOF_F80_BLOCH_SIGNWALK.md):
/// <code>
///   Spec(M)_{nontrivial} = { ±2i · λ : λ ∈ Spec(H_non-truly) }
///   mult_M(±2i·λ) = mult_H_non-truly(λ) · 2^N
/// </code>
///
/// <para>The remaining bilinears must all be Π²-odd 2-body: P, Q ∈ {X, Y, Z} with
/// bit_b(P) + bit_b(Q) ≡ 1 mod 2 — i.e. one of (X,Y), (X,Z), (Y,X), (Z,X). Π²-even
/// non-truly bilinears (Y,Z and Z,Y) are richer and are NOT in F80's verified scope.</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F80 entry.</para>
/// </summary>
public static class BlochSignWalk
{
    /// <summary>F80 spectrum prediction. Returns a dictionary {imaginary eigenvalue → multiplicity}.
    /// If all input terms are truly, returns {0 → 4^N}. Bond coupling is taken from <c>chain.J</c>.</summary>
    public static IReadOnlyDictionary<double, int> PredictMSpectrumImaginaryParts(
        ChainSystem chain, IReadOnlyList<PauliPairBondTerm> terms)
    {
        var nonTruly = new List<PauliPairBondTerm>();
        foreach (var t in terms)
        {
            if (t.IsTruly) continue;
            if (t.LetterA == PauliLetter.I || t.LetterB == PauliLetter.I)
                throw new ArgumentException(
                    $"term ({t.LetterA},{t.LetterB}) contains identity (single-body); F80 covers chain Π²-odd 2-body only");
            if (t.Pi2Parity != 1)
                throw new ArgumentException(
                    $"term ({t.LetterA},{t.LetterB}) is Π²-even non-truly; F80 covers Π²-odd only");
            nonTruly.Add(t);
        }

        if (nonTruly.Count == 0)
            return new Dictionary<double, int> { { 0.0, 1 << (2 * chain.N) } };

        var Hnt = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, nonTruly.ToBilinearSpec(chain.J)).ToMatrix();
        var evd = Hnt.Evd();
        var evals = evd.EigenValues.Select(z => Math.Round(z.Real, 10)).ToArray();

        var counts = new Dictionary<double, int>();
        foreach (var e in evals)
            counts[e] = counts.TryGetValue(e, out var n) ? n + 1 : 1;

        // Translate to M-spectrum: imaginary part is 2·λ_H, multiplicity ×2^N.
        int dPow = 1 << chain.N;
        var result = new Dictionary<double, int>();
        foreach (var (ev, mult) in counts)
            result[2.0 * ev] = mult * dPow;
        return result;
    }
}
