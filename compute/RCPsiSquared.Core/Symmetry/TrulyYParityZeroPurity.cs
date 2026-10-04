using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F107 (Tier1Derived): F87 truly classification under any single-letter
/// dephase channel (Z, X, or Y) forces y_par = 0 on every term in a truly-classified
/// Pauli pair. The Z criterion follows from Π_Z = R·D at every body count
/// (PROOF_PI_FACTORS_AS_R_TIMES_D §4(d), (f)); the canonical Π_Y = Π_Z⁻¹ and Π_X, the Hadamard
/// conjugate of Π_Z (PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE), carry it to Y and X.
///
/// <para>Derivation chain (see PROOF_F107):</para>
/// <list type="number">
///   <item>Per dephase D, truly term ⟺ the two letters other than the canonical mirror's
///         flip letter (X for Z and Y, Z for X) both occur an even number of times.</item>
///   <item>Z-deph truly: #Y even AND #Z even. X-deph truly: #X even AND #Y even.
///         Y-deph truly: #Y even AND #Z even.</item>
///   <item>All three include #Y even. y_par = #Y mod 2. Hence truly ⟹ y_par = 0.</item>
///   <item>For y_par-homogeneous pairs: both terms have y_par = 0; pair y_par = 0.</item>
/// </list>
///
/// <para>Empirical evidence: 4524 truly classifications across F103 (N=4 k=3),
/// F105 (N=5 k=3), F106 (N=4 k=4); zero have y_par=1. F107 explains this bit-exactly.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.YParity"/>;
/// fifth member of the YParity-axis Claim family after F102 + F103/N4K3 + F105/N5K3 +
/// F106/N4K4. Proof: <c>docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md</c>.</para></summary>
public sealed class TrulyYParityZeroPurity : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.YParity;
    public Claim? BitATwin => null;

    /// <summary>Per-dephase truly criterion as discovered in PROOF_F107 Step 1.
    /// Always includes "#Y even" as a sub-condition.</summary>
    public string ZDephTrulyCriterion => "#Y even AND #Z even";
    public string XDephTrulyCriterion => "#X even AND #Y even";
    public string YDephTrulyCriterion => "#Y even AND #Z even (Π_Y = Π_Z⁻¹, so the same strings as under Z)";

    /// <summary>The theorem statement: truly classification under any dephase letter
    /// forces y_par = 0 on every Pauli term.</summary>
    public string Theorem => "For any dephase letter D in {Z, X, Y}: F87 truly criterion ⟹ #Y(term) even ⟹ y_par(term) = 0";

    /// <summary>Pair corollary: a y_par-homogeneous pair classified as truly has shared y_par = 0.</summary>
    public string PairCorollary => "Klein-homogeneous + y_par-homogeneous pair is truly ⟹ pair y_par = 0";

    /// <summary>Returns true iff the given Pauli term satisfies the truly criterion under
    /// the specified dephase letter. Direct implementation of PROOF_F107 Step 1 table.</summary>
    public static bool TrulyCriterionHolds(PauliTerm term, PauliLetter dephase)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        int nX = term.Nx;
        int nY = term.Ny;
        int nZ = term.Nz;
        return dephase switch
        {
            PauliLetter.Z => (nY & 1) == 0 && (nZ & 1) == 0,
            PauliLetter.X => (nX & 1) == 0 && (nY & 1) == 0,
            PauliLetter.Y => (nY & 1) == 0 && (nZ & 1) == 0,
            _ => throw new ArgumentException($"dephase must be X, Y, or Z; got {dephase}", nameof(dephase)),
        };
    }

    /// <summary>Direct verification of F107: if the truly criterion holds for the given
    /// (term, dephase) pair, then y_par(term) must be 0. Returns true iff F107 holds for
    /// this input (true also when criterion does not hold, as F107 is one-way truly⟹y_par=0).</summary>
    public static bool VerifyOnTerm(PauliTerm term, PauliLetter dephase)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (!TrulyCriterionHolds(term, dephase)) return true;
        return term.YParity == 0;
    }

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F107
    /// classifies pairs on the y_par axis of the Z₂³ 8-cell decomposition;
    /// KleinEightCellClaim is the structural anchor. Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public TrulyYParityZeroPurity(KleinEightCellClaim klein8)
        : base("F107 F87 truly classification forces y_par = 0 (closed-form from the canonical mirrors, all dephase letters); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F107 + " +
               "docs/proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md + " +
               "docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md + " +
               "docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md + " +
               "docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md + " +
               "simulations/f107_f110_route_gate.py (the exact gate of the route) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F107 truly forces y_par = 0 (closed-form, all dephase letters)";

    public override string Summary =>
        $"Theorem: {Theorem}. Per-dephase truly criteria all include '#Y even'. " +
        $"Empirical: 4524 truly classifications (F103+F105+F106), zero y_par=1 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Pair corollary", summary: PairCorollary);
            yield return new InspectableNode("Z-dephase truly criterion", summary: ZDephTrulyCriterion);
            yield return new InspectableNode("X-dephase truly criterion", summary: XDephTrulyCriterion);
            yield return new InspectableNode("Y-dephase truly criterion", summary: YDephTrulyCriterion);
            yield return new InspectableNode("Empirical evidence",
                summary: "F103 (N=4 k=3): 300 truly, 0 y_par=1. F105 (N=5 k=3): 300 truly, 0 y_par=1. " +
                         "F106 (N=4 k=4): 3924 truly, 0 y_par=1. Total: 4524 truly, 0 y_par=1.");
            yield return new InspectableNode("Sibling y_par-axis claims",
                summary: "F109 (MotherSoftYParityOnePurity, Tier1Derived): mother sector Klein (0,0) y_par-homogeneous soft pairs have y_par=1 " +
                         "across all dephase letters, its non-truly pairs soft by the colouring. " +
                         "F110 (HardCellYInversionPattern, Tier1Derived): hard Klein-homogeneous pairs only in the diagonal cell, " +
                         "y_par-asymmetric with Y-inversion. F111 (HardCellPureDTemplate, Tier1Derived): the pure-D template rule " +
                         "at k = N = 4. The soft pairs of the other cells are F103's (§3.3, §8).");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition (bit_a, bit_b, y_par) anchoring the y_par axis F107 lives on.");
        }
    }
}
