using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The doormat for the Object Manager: the load-bearing house terms, two-three lines
/// each, written for a stranger with no physics PhD. One child node per term (DisplayName =
/// term, Summary = the plain-language definition). Compressed to CLI scale from
/// <c>docs/GLOSSARY.md</c> and the actual sources (the <see cref="Core.Knowledge.Tier"/> enum,
/// the <see cref="Core.Knowledge.Claim"/> base, the <see cref="Core.Confirmations.Confirmation"/>
/// record, the live witnesses, the open-arcs ledger, <c>docs/ANALYTICAL_FORMULAS.md</c>, and the
/// F87 trichotomy). Mounted as <c>--root glossary</c>; it needs no <c>--N</c>.</summary>
public static class GlossaryInspectableNode
{
    public static IInspectable Build()
    {
        var terms = new (string Term, string Definition)[]
        {
            ("claim",
                "A typed assertion the code can carry and check: a name, a tier label (how sure we are), " +
                "and an anchor (a pointer to the proof or formula entry it rests on). Inspect one with " +
                "--claim <ClassName>; the whole list lives under --root world."),

            ("tier",
                "The evidence scale on every claim. Tier 1 (derived) = proven analytically, bit-exact; " +
                "Tier 1 (candidate) = strong numerical witness, proof not yet closed; Tier 2 (empirical) = " +
                "a pattern seen across cases without proof; Tier 2 (hardware-verified) = predicted, then " +
                "measured on a real quantum chip; Open question = a known gap with a plan; Retracted = " +
                "made, then refuted by more data."),

            ("confirmation",
                "A prediction that a real quantum computer then confirmed. Each one records the machine " +
                "(e.g. ibm_marrakesh, ibm_kingston), the job id, the predicted value, and the measured " +
                "value, so you can check the match yourself. Browse them under --root world."),

            ("witness",
                "A claim's live lab. Instead of looking up a stored answer, a witness rebuilds the actual " +
                "matrix at inspect time and recomputes the evidence on the spot, then checks it against the " +
                "closed form. If the two meet, the claim stands in front of you, recomputed."),

            ("arc",
                "An entry in the open-arcs ledger: a line of work we started but did not finish. Each arc " +
                "records where it parked and the next concrete step; a leaf is either built (folded into a " +
                "claim) or retired (with the reason). It is the project's own to-do list. See --root arcs."),

            ("F-number",
                "An entry in docs/ANALYTICAL_FORMULAS.md, the formula registry running F1..F121. Each Fn is " +
                "one named result (a formula that replaces a matrix computation) with its own tier label. " +
                "F1 is the palindrome equation itself; later numbers refine and extend it."),

            ("palindrome / Π",
                "The central result. Under local Z-noise, a spin chain's list of decay rates reads the same " +
                "forwards and backwards: it is mirror-symmetric about the value −Σγ (minus the " +
                "total noise). Π (Pi) is the conjugation operator that proves it: it swaps every decay " +
                "rate with its mirror partner, so knowing one half gives the other exactly."),

            ("light / lens",
                "Two halves of a quantum state. The {X, Y} content oscillates and gets absorbed by the noise: " +
                "that part is the light. The {I, Z} content survives: that part is the lens. The decay rate " +
                "follows one law, Re(λ) = −2γ·⟨n_XY⟩: how fast a mode dies is " +
                "set by how much light it carries."),

            ("truly / soft / hard",
                "The F87 trichotomy: three ways a Hamiltonian's dephased spectrum can carry the mirror. " +
                "Truly = the term is inert, already palindromic for free. Soft = a second mirror (a chiral " +
                "sign flip) restores the palindrome. Hard = no such mirror exists and the palindrome is " +
                "genuinely broken."),

            ("Q and γ",
                "The two knobs. γ (gamma) is the dephasing rate, the clock's tick: how fast the " +
                "environment erases quantum coherence. Q = J/γ is the coupling J measured in clock " +
                "ticks, the only ratio an inside observer can see. K = γt is the accumulated dose " +
                "(elapsed time in ticks)."),

            ("--N and --n",
                "Two different sizes. --N is the number of qubits in the chain (the system size). --n is the " +
                "lower popcount of a coherence block: which (n, n+1) excitation-number block of that chain a " +
                "block-level root (fourmode, f86, c2*) works in. Roots that carry their own dimensions " +
                "(qudit, world, arcs, glossary) ignore --N entirely."),
        };

        var children = new IInspectable[terms.Length];
        for (int i = 0; i < terms.Length; i++)
            children[i] = new InspectableNode(displayName: terms[i].Term, summary: terms[i].Definition);

        return new InspectableNode(
            displayName: "glossary",
            summary: "the house language, two lines per term; start here, then --root world",
            children: children);
    }
}
