using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F61 closed form (Tier 1 proven, verified 69 configurations
/// N=2..7, 4 topologies, 4 γ profiles):
///
/// <code>
///   [L, Π²_X] = 0      exactly, for all N (isotropic Heisenberg/XY +
///                      Z-dephasing on any subset of sites, any graph,
///                      any site-dependent γ_k)
///
///   Π²_X = (−1)^n_XY   bit_a parity operator on Pauli strings
///
///   Consequences:
///     1. Every eigenmode has definite n_XY parity (even or odd)
///     2. Every SE density matrix has purely even n_XY content
///     3. No SE state can excite an odd-n_XY eigenmode (overlap = 0)
/// </code>
///
/// <para>This Π²_X operator is the global Z-string Z⊗N = ⊗_l Z_l: the (−1)^n_XY action
/// is conjugation by Z⊗N, typed as <c>ZGlobalMirror</c>. The operator identity
/// Z⊗N = Π²_X is the bit_a twin of F1² (X⊗N = Π²_Z), registered as a corollary of F61
/// in <c>docs/ANALYTICAL_FORMULAS.md</c>.</para>
///
/// <para>F61 is the bit_a Z₂ companion to F63's bit_b Z₂. The Pauli letter set
/// {I, X, Y, Z} factorises as <c>C₂ × C₂</c> indexed by <c>(bit_a, bit_b) =
/// (n_XY, w_YZ)</c>; the two Z₂ symmetries of L admitted by the d = 2 Pauli
/// algebra. Per F34/QUBIT_NECESSITY this is the maximal symmetry: no third
/// independent Z₂ classification exists.</para>
///
/// <list type="bullet">
///   <item><b>F61, bit_a parity (n_XY)</b>: L commutes with the n_XY parity
///         operator. Mechanism: Heisenberg bonds (XX, YY) flip exactly two
///         bit_a bits per term so n_XY parity is preserved; Z-dephasing
///         dissipator <c>Z ρ Z</c> has zero bit_a flips.</item>
///   <item><b>F63, bit_b parity (w_YZ)</b>: L commutes with the w_YZ parity
///         operator (<see cref="F63LCommutesPi2Pi2Inheritance"/>).</item>
///   <item><b>Combined</b>: 4-block decomposition <c>(bit_a, bit_b) ∈ {0,1}²</c>
///         with each block of dimension <c>4^(N − 1)</c>.</item>
/// </list>
///
/// <para>Pi2-Foundation anchors (parallel to F63):</para>
///
/// <list type="bullet">
///   <item><b>BlockCount = 4</b>: shared with F63; the C₂ × C₂ four-cell
///         decomposition. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1)
///         = <c>a_{−1} = 4</c> = d² for 1 qubit.</item>
///   <item><b>PerBlockDimension(N) = 4^(N − 1) = a_{3−2N}</b>: shared with
///         F63 / F38 / F39 / F1-T1.</item>
///   <item><b>Z₂SymmetryCount = 2</b> (joint with F63): the C₂ × C₂ is the
///         maximal Z₂ × Z₂ admitted by the d = 2 Pauli algebra (per F34).</item>
///   <item><b>SE-accessibility boundary</b>: single-excitation density
///         matrices have purely even n_XY → SE optimisers can only reach
///         even-n_XY modes. If a slower odd-n_XY eigenmode exists, its rate
///         is structurally beyond SE optimisation reach. This is F61's
///         engineering consequence; NOT directly Pi2-anchored, but inherits
///         from the bit_a Z₂ structure.</item>
/// </list>
///
/// <para>Scope: F61 holds for isotropic Heisenberg (XX+YY+ZZ) or XY (XX+YY),
/// any graph, any site-dependent γ_k. <b>Breaks for</b>: amplitude damping
/// (T1, note this is broader breakage than F63, which is preserved under
/// T1) and transverse fields (odd-n_XY terms in H). The asymmetry between
/// F61 and F63 break-conditions is the structural difference between bit_a
/// and bit_b parity.</para>
///
/// <para>Tier1Derived: F61 is Tier 1 proven (PROOF_PARITY_SELECTION_RULE).
/// Verified 69 configurations (chain N=2..7, ring N=2..6, star N=3..6,
/// complete N=3..4, 4 γ profiles + IBM-sacrifice at N=5 chain); the second
/// slow mode is SE-inaccessible in all 60 configs where it is distinct
/// (Frobenius ratio machine-zero). The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F61 +
/// <c>docs/proofs/PROOF_PARITY_SELECTION_RULE.md</c> +
/// <c>simulations/results/lens_survey/lens_survey_scaling.txt</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F63LCommutesPi2Pi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F61BitAParityPi2Inheritance : Claim, IZ2AxisClaim
{
    public F63LCommutesPi2Pi2Inheritance F63 { get; }
    public Pi2DyadicLadderClaim Ladder { get; }
    /// <summary>Joint with F63: the count of independent Z₂ symmetries L admits
    /// is <c>2</c>. F61 contributes bit_a (n_XY), F63 contributes bit_b (w_YZ).
    /// Per F34/QUBIT_NECESSITY the C₂ × C₂ is maximal.</summary>
    public int IndependentZ2SymmetryCount => F63.IndependentZ2SymmetryCount;

    /// <summary>The 4-block decomposition count: <c>4 = a_{−1}</c> on the dyadic
    /// ladder. Shared with F63; live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1).</summary>
    public double BlockCount => Ladder.Term(-1);

    /// <summary>Per-block dimension: <c>4^(N − 1) = a_{3−2N}</c> on the dyadic
    /// ladder. Shared with F63 / F38 / F39 / F1-T1; live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3 − 2N</c>).</summary>
    public double PerBlockDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F61 requires N ≥ 1.");
        return Ladder.Term(LadderIndexForPerBlock(N));
    }

    /// <summary>The dyadic-ladder index where the per-block dimension lands:
    /// <c>3 − 2N</c>. Same as F38, F39, F1-T1, F63.</summary>
    public int LadderIndexForPerBlock(int N) => 3 - 2 * N;

    /// <summary>True iff F61's per-block dimension agrees with F63's per-block
    /// dimension. Drift check on the bit_a Z₂ ↔ bit_b Z₂ joint structure.</summary>
    public bool PerBlockDimensionAgreesWithF63(int N) =>
        Math.Abs(PerBlockDimension(N) - F63.PerBlockDimension(N)) < 1e-12;

    /// <summary>True iff F61's BlockCount agrees with F63's BlockCount. Drift
    /// check on the joint C₂ × C₂ structure (both should anchor to a_{−1} =
    /// 4 on the ladder).</summary>
    public bool BlockCountAgreesWithF63() =>
        Math.Abs(BlockCount - F63.BlockCount) < 1e-12;

    /// <summary>The Z₂ axis F61 covers: bit_a (Π²_X = Z⊗N). F63 covers bit_b
    /// (Π²_Z = X⊗N). Together they generate the maximal Klein-Vierergruppe of
    /// independent Π² operators per F34/QUBIT_NECESSITY.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>F61 sits on the BitA axis; the BitATwin slot is always null
    /// (BitATwin is meaningful only for BitB Claims).</summary>
    public Claim? BitATwin => null;

    /// <summary>SE accessibility ceiling: every single-excitation density
    /// matrix has purely even n_XY content, so SE optimisers can ONLY reach
    /// even-n_XY eigenmodes. Returns <c>true</c> iff the supplied parity is
    /// even (= SE-accessible).</summary>
    public bool IsSeAccessible(int nXyParity) => nXyParity == 0;

    /// <summary>F61's break conditions, broader than F63's. Returns the set
    /// of operators/processes that break the bit_a Z₂ symmetry: amplitude
    /// damping (T1) and transverse fields (h_x, h_y). Note F63 is preserved
    /// under T1; the asymmetry IS the bit_a vs bit_b distinction.</summary>
    public IReadOnlyList<string> BreakConditions => new[]
    {
        "amplitude damping (T1; σ⁺/σ⁻ jump operators flip bit_a)",
        "transverse fields h_x · X_l (odd-n_XY single-site terms)",
        "transverse fields h_y · Y_l (odd-n_XY single-site terms)",
    };

    public F61BitAParityPi2Inheritance(
        F63LCommutesPi2Pi2Inheritance f63,
        Pi2DyadicLadderClaim ladder)
        : base("F61 [L, Π²_X] = 0 inherits from Pi2-Foundation: bit_a Z₂ companion to F63's bit_b Z₂",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F61 + " +
               "docs/proofs/PROOF_PARITY_SELECTION_RULE.md + " +
               "simulations/results/lens_survey/lens_survey_scaling.txt + " +
               "compute/RCPsiSquared.Core/Symmetry/F63LCommutesPi2Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        F63 = f63 ?? throw new ArgumentNullException(nameof(f63));
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F61 [L, Π²_X] = 0 (bit_a n_XY parity) as F63 sister-claim";

    public override string Summary =>
        $"[L, Π²_X] = 0 exactly all N (bit_a / n_XY parity); companion to F63 bit_b / w_YZ; together C₂ × C₂ maximal; " +
        $"SE accessibility ceiling: only even-n_XY modes reachable; breaks under T1 + transverse fields ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F61 closed form",
                summary: "[L, Π²_X] = 0 exactly all N; verified 69 configs (N=2..7, Chain/Star/Ring/Complete, 4 γ profiles); second slow mode SE-inaccessible in all 60 configs where it is distinct (Frobenius ratio machine-zero)");
            yield return new InspectableNode("F61 ↔ F63 sister reading",
                summary: "F61 = bit_a Z₂ (n_XY parity); F63 = bit_b Z₂ (w_YZ parity); together C₂ × C₂ maximal symmetry of d=2 Pauli algebra (per F34/QUBIT_NECESSITY)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "BlockCount = 4 = a_{-1} (shared with F63); PerBlockDim = 4^(N−1) = a_{3−2N} (shared); same shifts as F38, F39, F1-T1");
            yield return InspectableNode.RealScalar("BlockCount (= a_{-1} = 4)", BlockCount);
            yield return new InspectableNode("Z₂ axis",
                summary: $"F61 covers: bit_a (n_XY); F63 covers: bit_b (w_YZ); enum classification: {Z2Axis}");
            yield return new InspectableNode("SE accessibility corollary",
                summary: "every SE density matrix has purely even n_XY → SE optimisers can ONLY reach even-n_XY modes; if a slower odd-n_XY mode exists, its rate is structurally beyond SE optimisation reach");
            yield return new InspectableNode("F61 break conditions (broader than F63)",
                summary: "F61 breaks under T1 + transverse fields (h_x, h_y); F63 is preserved under T1, breaks only under transverse fields. Asymmetry IS the bit_a/bit_b distinction.");
            for (int N = 2; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"PerBlockDim = {PerBlockDimension(N):G6} (= a_{LadderIndexForPerBlock(N)}); " +
                             $"agrees with F63: {PerBlockDimensionAgreesWithF63(N)}");
            }
        }
    }
}
