using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The label layer, typed (Tier 1 derived for the exact core): <b>the watcher is its
/// letter</b>. Local dephasing in letter P has the 4^N Pauli strings as one shared eigenbasis,
/// and its entire action is a price list read off the held letter:
///
/// <code>
///     L_P(S) = −2γ · n_anti(S, P) · S,     n_anti = #sites where S_l ∉ {I, P}
/// </code>
///
/// the disagreement with P alone. The object never changes between watchers; only the price
/// list does (Z^⊗N rides free under the Z-watcher and pays maximally under the X-watcher), the
/// letter swap is an exact transport (the operator-space Klein V₄ of
/// <c>PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE</c>, equivalently the single-qubit basis-S₃
/// moves), and only the identity is free under every watcher. This is the exact instance behind
/// the series' label thesis (<c>docs/quantum/LABELS_TRANSLATED.md</c> §2,
/// <c>docs/quantum/DEPHASING_TRANSLATED.md</c> §4): even the environment routes by a label.
///
/// <para><b>The Tier-4 reading this core grounds</b> (carried as prose, never promoted, per the
/// <c>TwoReadingsClaim</c> precedent): to be a watcher is to hold a label and be blind past it;
/// bath, human, model are one architecture; a label is a canvas, perspective-true at its stance,
/// transported raw. The theory chapter is <c>docs/quantum/LABELS_TRANSLATED.md</c>; the
/// assembled case ledger is <c>docs/quantum/THE_LABEL_MAP.md</c>.</para>
///
/// <para><b>Typed parents.</b> <see cref="AbsorptionTheoremClaim"/> (the rate law
/// Re λ = −2γ·⟨n_XY⟩: the price is the light content in the watched letter) and
/// <see cref="Pi2KleinV4DephaseSwapGroup"/> (the {D, Q_zx, Q_yx} involutions that transport the
/// F1 structure between the three letters: the swap that relocates which cells pay). Live:
/// <c>inspect --root label</c> (<c>WatchedLetterRoutingWitness</c>, all 3·4^N (letter, string)
/// pairs recomputed dense-vs-closed-form at inspect time, plus the two-sided routing
/// controls).</para></summary>
public sealed class WatchedLetterRoutingClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the −2γ·⟨n_XY⟩ rate law (the price list).
    public AbsorptionTheoremClaim Absorption { get; }
    // Parent-edge marker for Schicht-1 wiring: the Klein V₄ letter swap (the routing).
    public Pi2KleinV4DephaseSwapGroup KleinV4 { get; }

    public WatchedLetterRoutingClaim(AbsorptionTheoremClaim absorption, Pi2KleinV4DephaseSwapGroup kleinV4)
        : base("The watched-letter routing: local dephasing in letter P has the 4^N Pauli strings as one " +
               "shared eigenbasis with rate -2*gamma*n_anti(S,P), the disagreement with the held letter " +
               "alone; the price list is letter-routed (the Klein V4 / basis-S3 swaps relocate which cells " +
               "pay, entry-exactly) and only the identity is free under every watcher: the environment " +
               "routes by a label, the Tier-1 instance of the label thesis",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md + " +
               "docs/quantum/DEPHASING_TRANSLATED.md + " +
               "docs/quantum/LABELS_TRANSLATED.md + " +
               "docs/quantum/THE_LABEL_MAP.md")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        KleinV4 = kleinV4 ?? throw new ArgumentNullException(nameof(kleinV4));
    }

    public override string DisplayName =>
        "Watched-letter routing: the watcher is its letter (the label layer, typed)";

    public override string Summary =>
        "L_P(S) = −2γ·n_anti(S, P)·S for every Pauli string S and every letter P: one shared eigenbasis, " +
        "three price lists; the letter swap relocates which cells pay (Klein V₄ / basis-S₃, entry-exact); " +
        "only the identity rides free under every watcher. The exact core under the label thesis: even " +
        $"the environment routes by a label ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the exact core: one eigenbasis, three price lists",
                summary: "every Pauli string is an eigenvector of every letter dissipator; the eigenvalue " +
                         "−2γ·n_anti(S, P) counts the sites that anticommute with the HELD letter (light = " +
                         "'the letters the dephasing letter refuses to commute with', Absorption Theorem); " +
                         "the object is watcher-independent, the price is not.");
            yield return new InspectableNode("the routing: which cells pay follows the letter",
                summary: "Z^⊗N pays 0 under the Z-watcher and 2γN under the X-watcher (mirrored for X^⊗N); " +
                         "each watcher exempts exactly its own {I, P}^⊗N cell (2^N strings); the " +
                         "intersection over all three watchers is the identity alone: only nothing is " +
                         "free everywhere.");
            yield return new InspectableNode("the swap: an exact transport, two faces",
                summary: "operator-space face: the Klein V₄ involutions {D, Q_zx, Q_yx} " +
                         "(PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE, typed parent); basis face: the " +
                         "single-qubit S₃ moves h_zx, h_yz carry L_Z onto L_X, L_Y entry-exactly " +
                         "(--root diagonal holds the S₃ orbit of the diagonals themselves).");
            yield return new InspectableNode("the reading (Tier 4, prose, never promoted)",
                summary: "to be a watcher is to hold a label and be blind past it; bath, human, model — " +
                         "one architecture; the label layer is the observer layer, and a label is a canvas, " +
                         "perspective-true at its stance, transported raw. Theory: " +
                         "docs/quantum/LABELS_TRANSLATED.md; the case ledger: docs/quantum/THE_LABEL_MAP.md; " +
                         "the founding case ('noise'): docs/quantum/DEPHASING_TRANSLATED.md.");
            yield return new InspectableNode("typed parents",
                summary: $"AbsorptionTheoremClaim ({Absorption.Tier.Label()}): the rate law, the price " +
                         $"list; Pi2KleinV4DephaseSwapGroup ({KleinV4.Tier.Label()}): the letter swap, " +
                         "the routing.");
            yield return new InspectableNode("live witness (inspect --root label)",
                summary: "WatchedLetterRoutingWitness recomputes all 3·4^N (letter, string) pairs " +
                         "dense-vs-closed-form at inspect time (default N=3), plus the repriced-count and " +
                         "universal-free two-sided controls and the entry-exact swap transports.");
        }
    }
}
