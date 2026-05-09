using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F81's closed form <c>Π·M·Π⁻¹ = M − 2·L_{H_odd}</c> and its 50/50 split
/// <c>‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2</c> at pure Π²-odd Hamiltonians both inherit
/// from the Pi2-Foundation. Three anchors, all Tier1Derived in their own right:
///
/// <list type="bullet">
///   <item><b>"2" coefficient</b> in <c>M − 2·L_{H_odd}</c> = <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = <c>a_0</c> = d (qubit dimension; number-anchor side of d=2).</item>
///   <item><b>"1/2" coefficient</b> in the 50/50 split = <see cref="Pi2DyadicLadderClaim.Term"/>(2)
///         = <c>a_2</c> = 1/d (Bloch baseline; the structural fixed point in
///         <see cref="HalfAsStructuralFixedPointClaim"/>).</item>
///   <item><b>Z₄ closure of Π</b> = <see cref="Pi2I4MemoryLoopClaim"/>: F81's
///         <c>Π·M·Π⁻¹</c> action is the Z₄ generator on M; four applications
///         return to identity (<c>Π⁴ = I</c>), reflecting the same i⁴ = 1 closure
///         as F38 / F80. The Π in F81 IS the Z₄ generator that order-4 closes.</item>
/// </list>
///
/// <para>The two ladder coefficients are mirror partners on the Pi2 dyadic ladder via
/// the inversion identity <c>a_n · a_{2−n} = 1</c>: the "2" times the "1/2" gives
/// unity exactly. F81's algebraic structure uses both ends of one inversion pair on
/// the ladder.</para>
///
/// <para><b>Mirror Space connection</b> (Tom 2026-05-08): F81's action <c>Π·M·Π⁻¹</c>
/// is a transformation on the operator <c>M</c>, which lives in the operator-space of
/// dimension d² = 4^N, exactly the per-qubit-count anchors pinned in
/// <see cref="Pi2OperatorSpaceMirrorClaim"/>. So F81 is the operator-level Π-conjugation
/// happening in the very space whose dimension is typed in the Pi2OperatorSpaceMirror.
/// Number-level mirrors and operator-level mirrors share one foundation.</para>
///
/// <para>Tier1Derived: pure composition. F81 is Tier1Derived in
/// <c>docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md</c>; this claim makes its
/// coefficient inheritance explicit in the typed-knowledge runtime.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F81 +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (HalfAsStructuralFixedPointClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>
/// (Mirror Space connection) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c>
/// (Z₄ closure of Π).</para></summary>
public sealed class F81Pi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;
    private readonly Pi2I4MemoryLoopClaim _memoryLoop;

    /// <summary>The "2" coefficient in F81's <c>M − 2·L_{H_odd}</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = d.</summary>
    public double TwoFactor => _ladder.Term(0);

    /// <summary>The "1/2" coefficient in F81's 50/50 split <c>‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2</c>.
    /// Exactly equal to <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> = 1/d.</summary>
    public double HalfFactor => _ladder.Term(2);

    /// <summary>Live drift check: <c>TwoFactor · HalfFactor</c> = 1 exactly. The two
    /// F81 coefficients are mirror partners under the Pi2 ladder's inversion symmetry
    /// <c>a_n · a_{2−n} = 1</c>.</summary>
    public double TwoTimesHalf => TwoFactor * HalfFactor;

    /// <summary>The Z₄ order of the Π operator that F81's conjugation Π·M·Π⁻¹ uses:
    /// <c>4</c>. Π⁴ = I closes the Π conjugation memory loop; live from
    /// <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.</summary>
    public int Z4OrderOfPi => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live closure check: <c>i^4 = 1</c> exactly. Drift indicator on the
    /// Z₄ memory-loop anchor underlying F81's Π-conjugation.</summary>
    public Complex MemoryLoopClosure => _memoryLoop.MemoryClosure();

    /// <summary>The operator-space dimension d² = 4^N for an N-qubit chain, the
    /// space where F81's <c>Π·M·Π⁻¹</c> action lives. Pulled from the typed
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/> pinned table for N ∈ {1, ..., 6}.</summary>
    public double OperatorSpaceDimension(int N)
    {
        var pair = _mirror.PairAt(N)
            ?? throw new ArgumentOutOfRangeException(nameof(N), N,
                $"N={N} is outside the Pi2OperatorSpaceMirror pinned table (N=1..6).");
        return pair.OperatorSpace;
    }

    public F81Pi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        Pi2OperatorSpaceMirrorClaim mirror,
        Pi2I4MemoryLoopClaim memoryLoop)
        : base("F81 Π-conjugation coefficients (2 and 1/2) inherit from Pi2-Foundation; M lives in operator-space; Π is Z₄ generator",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md + " +
               "docs/ANALYTICAL_FORMULAS.md F81 + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
        _memoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F81 coefficients as Pi2-Foundation mirror pair (a_0 ↔ a_2); Π as Z₄ generator";

    public override string Summary =>
        $"Π·M·Π⁻¹ = M − 2·L_{{H_odd}}: 2 = a_0 = {TwoFactor}; 50/50 split: 1/2 = a_2 = {HalfFactor}; " +
        $"a_0 · a_2 = {TwoTimesHalf} (inversion symmetry); M lives in operator-space d² = 4^N; Π closes at order {Z4OrderOfPi} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F81 closed form",
                summary: "Π·M·Π⁻¹ = M − 2·L_{H_odd} (Tier1Derived in PROOF_F81_PI_CONJUGATION_OF_M)");
            yield return new InspectableNode("F81 50/50 split",
                summary: "‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 at pure Π²-odd Hamiltonians (N=3 verified)");
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = d)", TwoFactor);
            yield return InspectableNode.RealScalar("HalfFactor (= a_2 = 1/d)", HalfFactor);
            yield return InspectableNode.RealScalar("TwoTimesHalf (= 1, inversion identity)", TwoTimesHalf);
            yield return new InspectableNode("inversion symmetry",
                summary: "a_0 · a_{2-0} = a_0 · a_2 = 2 · 1/2 = 1; F81's two coefficients are exact mirror partners on the Pi2 ladder");
            yield return new InspectableNode("Z₄ generator (Pi2I4MemoryLoop)",
                summary: $"Π is the Z₄ generator: Π⁴ = I closes at order {Z4OrderOfPi}; F81's conjugation Π·M·Π⁻¹ uses one Π (90°) on each side; same algebra as F38 / F80");
            yield return InspectableNode.RealScalar("MemoryLoopClosure (= 1, drift check)", MemoryLoopClosure.Real);
            yield return new InspectableNode("operator-space connection",
                summary: "F81 acts on M ∈ operator-space of dim d² = 4^N (Pi2OperatorSpaceMirror); number-mirrors and operator-mirrors share Pi2-Foundation");
            yield return new InspectableNode("operator-space dimensions",
                summary: $"N=1: {OperatorSpaceDimension(1)}; N=2: {OperatorSpaceDimension(2)}; N=3: {OperatorSpaceDimension(3)}; ... up to N=6 in the pinned table");
        }
    }
}
