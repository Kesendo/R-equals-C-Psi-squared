using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>How the F89 palindrome relates to the path-3 octic monodromy: the mirror SPLITS at the Galois
/// boundary. Surfaces <see cref="GaloisMonodromyWitness.MirrorMonodromy"/> and
/// <see cref="GaloisMonodromyWitness.PalindromeStrandPairing"/>.
///
/// <para>C-K (passes): the mirror's real-structure shadow, the q ↦ −q̄ reflection (from L(q)* = L(−q̄)),
/// intertwines the monodromy exactly. The induced octic-strand bijection is the identity (σ_K = id), so
/// every EP carries the same braid as its q ↦ −q̄ mirror.</para>
///
/// <para>C-T (blocked): the spectral palindrome λ ↦ −λ̄ − 8 (mirror about Re λ = −4) induces a genuine
/// strand involution σ_T (fixed strands on the fold + mirror-twin 2-cycles), but it CANNOT be a symmetry
/// of the braiding: commuting with the full S_8 monodromy forces it central, and Z(S_8) = 1. So σ_T is an
/// element OF the Galois group, never a symmetry OF it. Sibling of --root galoismonodromy and
/// --root branchpalindrome; reading reflections/ON_WHO_WATCHES_WHOM.md.</para></summary>
public sealed class MonodromyMirrorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public string DisplayName =>
        "F89 monodromy meets the mirror (live: q↦−q̄ intertwines the octic braiding; the Re=−4 palindrome cannot, Z(S_8)=1)";

    public string Summary
    {
        get
        {
            var (specK, specT, reports, all) = GaloisMonodromyWitness.MirrorMonodromy();
            var (_, fixedPts, twoCyc, _, _) = GaloisMonodromyWitness.PalindromeStrandPairing();
            return $"the mirror splits at the Galois boundary. C-K: the q↦−q̄ reflection intertwines the " +
                   $"monodromy {(all ? "on every EP (σ_K = id)" : "PARTIALLY (see children)")} " +
                   $"(spectral sanity to {specK.ToString("E1", Inv)}). C-T: the Re=−4 palindrome induces a " +
                   $"strand involution σ_T ({fixedPts} fixed on the fold, {twoCyc} mirror-twin 2-cycles; " +
                   $"sanity to {specT.ToString("E1", Inv)}) but is NOT a braid symmetry (Z(S_8)=1 forbids it). " +
                   $"σ_T is an element of the Galois group, never a symmetry of it.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (specK, specT, reports, all) = GaloisMonodromyWitness.MirrorMonodromy();
            var (sigmaT, fixedPts, twoCyc, braidInv, nBraids) = GaloisMonodromyWitness.PalindromeStrandPairing();

            yield return new InspectableNode("C-K: q↦−q̄ intertwines the monodromy (σ_K = identity)",
                summary: $"spectral sanity conj(spec@+2)=spec@−2 to {specK.ToString("E2", Inv)}; the induced " +
                         $"octic-strand bijection is the identity, and every EP carries the same braid as its " +
                         $"q↦−q̄ mirror (τ(−q̄*) = τ(q*)). Holds on {reports.Count(r => r.Intertwines)}/{reports.Count} " +
                         $"right-half EPs ({(all ? "all ✓" : "partial")}). The branch-locus palindrome lifted from the " +
                         "seams' positions to the braids they carry.");

            foreach (var r in reports)
                yield return new InspectableNode(
                    $"  EP q={r.Q.Real.ToString("0.000", Inv)}{r.Q.Imaginary.ToString("+0.000;-0.000", Inv)}i  (Re λ_EP={r.LambdaMidRe.ToString("0.00", Inv)})",
                    summary: $"τ(q*)=({string.Join(" ", r.TauPlus)})  τ(−q̄*)=({string.Join(" ", r.TauMinus)})  " +
                             $"intertwines={r.Intertwines}");

            yield return new InspectableNode("C-T: the Re=−4 palindrome is NOT a braid symmetry (Z(S_8)=1)",
                summary: $"the fixed-q palindrome λ↦−λ̄−8 (sanity to {specT.ToString("E2", Inv)}) induces the strand " +
                         $"involution σ_T = ({string.Join(" ", sigmaT)}): {fixedPts} fixed strands on the fold " +
                         $"(Re λ=−4), {twoCyc} mirror-twin 2-cycles. The EP-braid set is " +
                         $"{(braidInv ? "invariant" : "NOT invariant")} under σ_T (over {nBraids} braids). σ_T cannot " +
                         "commute with the full S_8 monodromy (that would force it central, but Z(S_8)=1), so the " +
                         "Re=−4 palindrome is an element OF the Galois group, never a symmetry OF it: the mirror enters " +
                         "the unwritable half as one of its unsortable moves, not as a rule above them.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
