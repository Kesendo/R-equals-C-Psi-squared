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
/// bit_b(P) + bit_b(Q) ≡ 1 mod 2, i.e. one of (X,Y), (X,Z), (Y,X), (Z,X). Π²-even
/// non-truly bilinears (Y,Z and Z,Y) are richer and are NOT in F80's verified scope.
/// The remaining bilinears must also share one Step-5 sign s: s = +1 for (X,Y) and (Y,X),
/// s = −1 for (X,Z) and (Z,X). A sum mixing the two signs has M(ρ) = −2i·H_L·ρ + 2i·ρ·H_R,
/// with H_L the s = +1 bonds and H_R the s = −1 bonds, whose spectrum is in general not
/// 2i·Spec(H); such a sum is outside the proven scope and is rejected, also where it happens to
/// match (two commuting terms at N=2). The predictor is scoped to the open chain, where F80's
/// sign walk is proven, so a non-chain <c>ChainSystem.Topology</c> is rejected.</para>
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
        if (chain.Topology != TopologyKind.Chain)
            throw new ArgumentException(
                $"chain topology is {chain.Topology}; this predictor is scoped to the open chain (TopologyKind.Chain), where F80's sign walk is proven");

        var nonTruly = new List<PauliPairBondTerm>();
        foreach (var t in terms)
        {
            if (t.IsTruly) continue;
            if (t.LetterA == PauliLetter.I || t.LetterB == PauliLetter.I)
                throw new ArgumentException(
                    $"term ({t.LetterA},{t.LetterB}) carries an identity leg; F80's sign-walk is stated for the identity-free Π²-odd pairs (XY, XZ, YX, ZX)");
            if (t.Pi2Parity != 1)
                throw new ArgumentException(
                    $"term ({t.LetterA},{t.LetterB}) is Π²-even non-truly; F80 covers Π²-odd only");
            nonTruly.Add(t);
        }

        // Step-5 sign of each Π²-odd pair: +1 when the pair carries Y, −1 when it carries Z.
        bool hasPlus = nonTruly.Any(t => t.LetterA == PauliLetter.Y || t.LetterB == PauliLetter.Y);
        bool hasMinus = nonTruly.Any(t => t.LetterA == PauliLetter.Z || t.LetterB == PauliLetter.Z);
        if (hasPlus && hasMinus)
            throw new ArgumentException(
                "the Π²-odd terms mix the Step-5 signs s = +1 ((X,Y), (Y,X)) and s = −1 ((X,Z), (Z,X)); " +
                "F80's identity Spec(M) = 2i·Spec(H) is proven for a same-sign sum, a mixed sum has " +
                "M(ρ) = −2i·H_L·ρ + 2i·ρ·H_R with H_L the s = +1 bonds and H_R the s = −1 bonds " +
                "(PROOF_F80_BLOCH_SIGNWALK), whose spectrum is in general not 2i·Spec(H), and is " +
                "outside the scope predicted here");

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
