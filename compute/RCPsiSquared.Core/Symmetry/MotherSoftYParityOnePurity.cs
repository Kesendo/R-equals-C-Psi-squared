using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F109 (Tier1Derived): Under any single-letter dephase channel (Z, X, or Y),
/// every Pauli pair classified as soft and located in the Mother sector Klein (0, 0)
/// has shared y_par = 1. Sister to <see cref="TrulyYParityZeroPurity"/> (F107) on
/// the same y_par axis; together they pin the y_par signature of two of the four
/// soft + truly slots in the Mother sector across all dephase letters.
///
/// <para>Derivation chain (see PROOF_F109):</para>
/// <list type="number">
///   <item>Klein (0,0) constraint (#X+#Y even AND #Y+#Z even) forces all three (#X, #Y, #Z)
///         to share the same parity.</item>
///   <item>Per F107 per-dephase truly criteria + same-parity collapse: Klein (0,0) truly
///         under any dephase = all three even.</item>
///   <item>Klein (0,0) non-truly (Π²-even non-truly) = all three odd.</item>
///   <item>Klein (0,0) is Π²-EVEN under every dephase (bit_b=0 for Z/Y, bit_a=0 for X).</item>
///   <item><b>Empirically:</b> Π²-even non-truly pairs are SOFT (not hard) across
///         F103/F105/F106 (5346+ pairs verified, zero Π²-even hard observed).
///         Closed-form proof of this step is the open F108 Part 1 (block-restricted
///         palindrome lemma); F109 inherits this dependency.</item>
///   <item>Klein (0,0) soft term ⟹ #Y odd ⟹ y_par = 1; y_par-homogeneous pair: shared y_par = 1.</item>
/// </list>
///
/// <para>Empirical evidence: 1026 mother-soft classifications across F103 (3 × 21),
/// F105 (3 × 21), F106 (3 × 300), all y_par=1, zero y_par=0. F109 explains this
/// bit-exactly modulo the Step 5 dependency.</para>
///
/// <para>Cross-letter spot-check (Step 3 + algebra): Klein (0,0) non-truly k=3
/// terms have (#X, #Y, #Z) = (1, 1, 1) (only triple with all-odd and sum ≤ 3),
/// giving 3! = 6 XYZ-permutations. Unordered pairs with self: 6·7/2 = 21 (matches
/// F103/F105 (0, 21)). For k=4: 24 letter sequences (3 non-I + 1 I), pairs 24·25/2
/// = 300 (matches F106 (0, 300)).</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// sixth member of the YParity-axis Claim family. Proof:
/// <c>docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md</c>.</para></summary>
public sealed class MotherSoftYParityOnePurity : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    /// <summary>Per-dephase Klein (0,0) non-truly criterion (same for all three dephase
    /// letters after the Klein constraint + F107 truly-collapse). All three letter
    /// counts must be odd; only the (1, 1, 1) triple satisfies this at k ≤ 3.</summary>
    public string NonTrulyCriterion => "Klein (0,0) non-truly ⟺ #X, #Y, #Z all odd (forced by Klein same-parity + F107 truly = all even)";

    /// <summary>Step 5 dependency (open as F108 Part 1). Empirically verified but
    /// not yet closed-form derived.</summary>
    public string OpenDependency => "F108 Part 1: Π²-even non-truly ⟹ SOFT (not hard). Empirically verified across 5346+ pairs (F103/F105/F106), closed-form proof open (block-restricted palindrome lemma on Π²-eigenspace decomposition of L).";

    /// <summary>The theorem statement: mother (0,0) soft ⟹ shared y_par = 1.</summary>
    public string Theorem => "Klein (0,0) soft under any dephase D in {Z, X, Y} ⟹ pair y_par = 1; equivalently #Y(both terms) = 1 mod 2";

    /// <summary>Returns true iff the given Pauli term is a Klein (0, 0) non-truly
    /// term (= candidate to be classified soft under any dephase, per Step 3 of
    /// PROOF_F109). Per Step 6, if Step 5 holds for this term, it has y_par = 1.</summary>
    public static bool IsMotherNonTrulyCandidate(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        // Klein (0, 0): bit_a = 0, bit_b = 0.
        if ((term.TotalBitA & 1) != 0) return false;
        if (term.Pi2Parity != 0) return false;
        // Non-truly: at least one of #Y, #Z odd under Z-deph (any dephase by F107
        // collapse on Klein (0,0): equivalent to NOT all even). Per Step 3 of
        // PROOF_F109, Klein (0,0) non-truly ⟺ all three #X, #Y, #Z odd.
        return (term.Nx & 1) == 1 && (term.Ny & 1) == 1 && (term.Nz & 1) == 1;
    }

    /// <summary>Direct verification of F109 Step 6: if a Klein (0,0) term is a
    /// non-truly candidate (Step 3 hit), it must have y_par = 1. Returns true iff
    /// F109 holds for this input (true also when the term is not a non-truly
    /// candidate; F109 is conditional on Klein (0,0) non-truly classification).</summary>
    public static bool VerifyOnTerm(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (!IsMotherNonTrulyCandidate(term)) return true;
        return term.YParity == 1;
    }

    public MotherSoftYParityOnePurity()
        : base("F109 mother sector Klein (0,0) soft is y_par=1 pure (closed-form modulo F108 Part 1, all dephase letters)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F109 + " +
               "docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md")
    {
    }

    public override string DisplayName =>
        "F109 mother soft = y_par 1 pure (closed-form modulo F108 Part 1)";

    public override string Summary =>
        $"Theorem: {Theorem}. Steps 1-4 + 6 are closed-form via F107 + Klein same-parity collapse. " +
        $"Step 5 (Π²-even ⟹ soft) is empirically verified across 5346+ pairs but closed-form open. " +
        $"Empirical: 1026 mother-soft classifications (F103+F105+F106), zero y_par=0 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Non-truly criterion (Step 3)", summary: NonTrulyCriterion);
            yield return new InspectableNode("Open dependency (Step 5, F108 Part 1)", summary: OpenDependency);
            yield return new InspectableNode("Cross-letter spot-check",
                summary: "k=3: Klein (0,0) non-truly = 6 XYZ-perms ⟹ 21 unordered pairs (matches F103/F105 (0, 21) ×3). " +
                         "k=4: 24 sequences ⟹ 300 pairs (matches F106 (0, 300) ×3).");
            yield return new InspectableNode("Empirical evidence",
                summary: "F103 (N=4 k=3): mother soft (0, 21) ×3 dephase. F105 (N=5 k=3): same. " +
                         "F106 (N=4 k=4): (0, 300) ×3. Total: 1026 mother-soft, zero y_par=0.");
            yield return new InspectableNode("Sister claims on YParity axis",
                summary: "F107: truly ⟹ y_par=0 (closed-form). F109: mother soft ⟹ y_par=1 (closed-form mod F108 Part 1). " +
                         "Together pin two of the four trichotomy slots in Klein (0,0).");
            yield return new InspectableNode("Open siblings",
                summary: "F108 Part 1 (Π²-even-soft block-palindrome lemma) blocks F109's Tier1Derived upgrade to fully unconditional. " +
                         "F110: hard cells y_par-pure with Y-inversion remains the deeper open work.");
        }
    }
}
