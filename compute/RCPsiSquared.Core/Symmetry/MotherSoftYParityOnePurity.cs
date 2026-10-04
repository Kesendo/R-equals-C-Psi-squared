using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F109 (Tier1Derived): Under any single-letter dephase channel (Z, X, or Y),
/// every y_par-homogeneous Pauli pair classified as soft and located in the Mother sector
/// Klein (0, 0) has shared y_par = 1. Sister to <see cref="TrulyYParityZeroPurity"/> (F107) on
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
///   <item>Klein (0,0) non-truly pairs are SOFT (not hard) at every body count, coupling
///         and per-site rate: a Mother-sector string commutes with every letter string,
///         and the canonical mirror's flip-letter string (X^⊗N for Z- and Y-dephasing,
///         Z^⊗N for X-dephasing) anticommutes with every jump, so right multiplication by
///         it gives R·L·R⁻¹ = −L† − 2σ and the spectrum pairs about −σ (the colouring,
///         with F158's sufficiency step, PROOF_PALINDROME_TWO_END_COUNT §(e)). The F108
///         mirrors (<see cref="Pi5BilinearOperator"/>) reach only the Π²-D-even strings of
///         even weight, and an all-odd string has odd weight.</item>
///   <item>Klein (0,0) non-truly term ⟹ #Y odd ⟹ y_par = 1; a soft pair holds one, so a
///         y_par-homogeneous soft pair has shared y_par = 1.</item>
/// </list>
///
/// <para>Empirical evidence: 1026 mother-soft classifications across F103 (3 × 21),
/// F105 (3 × 21), F106 (3 × 300), all y_par=1, zero y_par=0. F109 explains this
/// bit-exactly.</para>
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

    /// <summary>Step 5 status: the Mother sector is coloured. Its strings commute with
    /// every letter string; the canonical mirror's lit flip-letter string reflects L,
    /// and F158's sufficiency step pairs the spectrum.</summary>
    public string Step5Status => "Step 5 by the colouring: a Mother-sector string commutes with X^⊗N, Y^⊗N and Z^⊗N; the canonical mirror's flip-letter string anticommutes with every jump, so R(ρ) = ρ·A^⊗N gives R·L·R⁻¹ = −L† − 2σ and the spectrum pairs about −σ (F158 §(e)), at every body count, coupling and per-site rate.";

    /// <summary>The theorem statement: a y_par-homogeneous mother (0,0) soft pair has shared y_par = 1.</summary>
    public string Theorem => "Klein (0,0) y_par-homogeneous pair soft under any dephase D in {Z, X, Y} ⟹ shared y_par = 1; equivalently #Y odd in both terms";

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

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F109
    /// pins the y_par=1 purity of mother Klein (0,0) soft cells across the Z₂³
    /// 8-cell decomposition; KleinEightCellClaim is the structural anchor.
    /// Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public MotherSoftYParityOnePurity(KleinEightCellClaim klein8)
        : base("F109 mother sector Klein (0,0) y_par-homogeneous soft pairs are y_par=1 pure across all three dephase letters (its non-truly pairs soft by the colouring, F158); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F109 + " +
               "docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md + " +
               "experiments/THE_PALINDROME_AS_A_COLOURING.md + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md + " +
               "simulations/f107_f110_route_gate.py (the exact gate of the route) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F109 y_par-homogeneous mother soft = y_par 1 pure (closed-form)";

    public override string Summary =>
        $"Theorem: {Theorem}. Steps 1-4 + 6 are closed-form via F107 + Klein same-parity collapse. " +
        $"Step 5 (Mother non-truly ⟹ soft) is the colouring of the Mother sector with F158, " +
        $"at every body count. " +
        $"Empirical: 1026 mother-soft classifications, zero y_par=0 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Non-truly criterion (Step 3)", summary: NonTrulyCriterion);
            yield return new InspectableNode("Step 5 status (the colouring, F158)", summary: Step5Status);
            yield return new InspectableNode("Cross-letter spot-check",
                summary: "k=3: Klein (0,0) non-truly = 6 XYZ-perms ⟹ 21 unordered pairs (matches F103/F105 (0, 21) ×3). " +
                         "k=4: 24 sequences ⟹ 300 pairs (matches F106 (0, 300) ×3).");
            yield return new InspectableNode("Empirical evidence",
                summary: "F103 (N=4 k=3): mother soft (0, 21) ×3 dephase. F105 (N=5 k=3): same. " +
                         "F106 (N=4 k=4): (0, 300) ×3. Total: 1026 mother-soft, zero y_par=0.");
            yield return new InspectableNode("Sister claims on YParity axis",
                summary: "F107: truly ⟹ y_par=0 (closed-form). F109: y_par-homogeneous mother soft ⟹ y_par=1 (closed-form). " +
                         "Together pin two of the four trichotomy slots in Klein (0,0).");
            yield return new InspectableNode("Open siblings",
                summary: "F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10, typed 2026-05-25): hard cells y_par-asymmetric " +
                         "with Y-inversion; with F107 (truly) and F109 (mother soft) it gives the y_par signature of the truly, " +
                         "Mother-soft and hard pairs. Closed-form 42:8/228:0 derivation Tier1Derived via F103 §6+§7 (2026-06-10).");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition anchor for the y_par axis F109 lives on.");
        }
    }
}
