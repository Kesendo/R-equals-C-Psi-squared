using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>How the F89 palindrome relates to the path-3 octic monodromy: the mirror SPLITS at the Galois
/// boundary. Surfaces <see cref="GaloisMonodromyWitness.MirrorMonodromy"/> and
/// <see cref="GaloisMonodromyWitness.PalindromeStrandPairing"/>.
///
/// <para>C-K (the base-space face passes): the q ↦ −q̄ reflection (from L(q)* = L(−q̄)) is an exact family
/// symmetry, so it automatically intertwines the monodromy; in the aligned labelling the induced strand
/// bijection comes out as the identity (σ_K = id), so each cluster EP carries the same braid as its
/// q ↦ −q̄ mirror.</para>
///
/// <para>C-T (the fibre face does not): the spectral fold λ ↦ −λ̄ − 8 (mirror about Re λ = −4, exact at the
/// real base q) induces a genuine strand involution σ_T (four fixed on the fold + two mirror-twin 2-cycles),
/// which is NON-CENTRAL: commuting with the full S_8 monodromy would force it central, but Z(S_8) = 1 and
/// σ_T ≠ id, so it does not commute with the braiding (not a loop-independent symmetry; conjugation by it is
/// still a nontrivial inner automorphism of S_8). Sibling of --root galoismonodromy and
/// --root branchpalindrome; reading reflections/ON_WHO_WATCHES_WHOM.md.</para></summary>
public sealed class MonodromyMirrorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public string DisplayName =>
        "F89 monodromy meets the mirror (live: q↦−q̄ intertwines the octic braiding; the Re=−4 fold is non-central, Z(S_8)=1)";

    public string Summary
    {
        get
        {
            var (specK, specT, reports, all) = GaloisMonodromyWitness.MirrorMonodromy();
            var (_, fixedPts, twoCyc, _, _) = GaloisMonodromyWitness.PalindromeStrandPairing();
            return $"the mirror splits at the Galois boundary. C-K: the q↦−q̄ reflection intertwines the " +
                   $"monodromy {(all ? "on each cluster EP (σ_K = id)" : "PARTIALLY (see children)")} " +
                   $"(spectral sanity to {specK.ToString("E1", Inv)}; forced by L(q)*=L(−q̄)). C-T: the Re=−4 fold " +
                   $"induces a non-central involution σ_T ({fixedPts} fixed on the fold, {twoCyc} mirror-twin " +
                   $"2-cycles; sanity to {specT.ToString("E1", Inv)}) that does NOT commute with the monodromy " +
                   $"(Z(S_8)=1), so it is not a loop-independent symmetry of the braiding (conjugation by it is " +
                   $"still an inner automorphism of S_8).";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (specK, specT, reports, all) = GaloisMonodromyWitness.MirrorMonodromy();
            var (sigmaT, fixedPts, twoCyc, braidInv, nBraids) = GaloisMonodromyWitness.PalindromeStrandPairing();

            yield return new InspectableNode("C-K: q↦−q̄ intertwines the monodromy (σ_K = identity, aligned labelling)",
                summary: $"spectral sanity conj(spec@+2)=spec@−2 to {specK.ToString("E2", Inv)}; the intertwining is " +
                         $"forced a priori by L(q)*=L(−q̄), and in the aligned labelling the induced strand bijection " +
                         $"is the identity, so each cluster EP carries the same braid as its q↦−q̄ mirror (τ(−q̄*) = " +
                         $"τ(q*)). Holds on {reports.Count(r => r.Intertwines)}/{reports.Count} near-real-axis cluster " +
                         $"EPs ({(all ? "all ✓" : "partial")}). The branch-locus palindrome lifted from the seams' " +
                         "positions to the braids they carry.");

            foreach (var r in reports)
                yield return new InspectableNode(
                    $"  EP q={r.Q.Real.ToString("0.000", Inv)}{r.Q.Imaginary.ToString("+0.000;-0.000", Inv)}i  (Re λ_EP={r.LambdaMidRe.ToString("0.00", Inv)})",
                    summary: $"τ(q*)=({string.Join(" ", r.TauPlus)})  τ(−q̄*)=({string.Join(" ", r.TauMinus)})  " +
                             $"intertwines={r.Intertwines}");

            yield return new InspectableNode("C-T: the Re=−4 fold is non-central, not a loop-independent braid symmetry (Z(S_8)=1)",
                summary: $"the real-q spectral fold λ↦−λ̄−8 (sanity to {specT.ToString("E2", Inv)}) induces the strand " +
                         $"involution σ_T (one-line) = [{string.Join(" ", sigmaT)}]: {fixedPts} fixed strands on the fold " +
                         $"(Re λ=−4), {twoCyc} mirror-twin 2-cycles. σ_T ≠ id, and a permutation commuting with the full " +
                         "S_8 monodromy would be central, but Z(S_8)=1, so σ_T does NOT commute with the braiding: the " +
                         "strand mirror-pairing is braided away around the seams. (Conjugation by σ_T is still an inner " +
                         "automorphism of S_8; the strictly stronger 'σ_T permutes the EP-transposition set' is " +
                         $"{(braidInv ? "invariant" : "numerically NOT invariant")} over {nBraids} cluster braids, but that " +
                         "is not what Z(S_8)=1 settles.) The from-below ground of 'who watches whom has no fixed answer'.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
