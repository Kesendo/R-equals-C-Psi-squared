using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F1's palindrome closed form <c>Π · L · Π⁻¹ = −L − 2Σγ · I</c> has its
/// "<c>2</c>" coefficient sitting on the Pi2-Foundation: it is exactly <c>a_0</c> on
/// the dyadic halving ladder, the qubit dimension d. The shift "2Σγ·I" is therefore
/// "<c>d · Σγ · I</c>", the Σγ accumulates the per-site dephasing rates, and the
/// "d" counts how the F1 derivation stacks two Π applications.
///
/// <para>The "+L" and "−L" form a Z₂ inversion in the Liouvillian operator algebra;
/// applying Π once gives <c>−L − 2σ·I</c>, applying it twice closes back to <c>+L</c>.
/// That two-step closure is the F1 palindrome reading: Π is involutive on the
/// algebra modulo the constant shift <c>2σ·I</c>, where the "2" is the qubit
/// dimension.</para>
///
/// <para><b>Z₄ memory-loop reading</b>: F1 is Layer 1 of the Pi2I4MemoryLoop ladder
/// (per <see cref="Pi2I4MemoryLoopClaim"/>'s own docstring). The "−1" in "−L" is
/// two 90° rotations summed: <c>i² = −1</c>, the Z₄ generator squared. So F1's
/// sign flip IS the Z₄ memory loop at half-period; F38's Π² = (−1)^w_YZ is the
/// same inversion at the Pauli-string level.</para>
///
/// <para>This is not a new identity, F1 is Tier1Derived in
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> and verified at N=8 (all 65,536 eigenvalues
/// paired, zero exceptions). What this claim makes typed is the <b>inheritance</b>: the
/// constant 2 in F1's closed form is not a free parameter; it is <c>a_0 = d</c> on
/// the Pi2 dyadic ladder, and the sign flip is i² on the Z₄ memory loop.</para>
///
/// <para>Tier1Derived: pure composition. Anchors:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F1 +
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c>.</para></summary>
public sealed class F1Pi2Inheritance : Claim, IF99AnchorBearing, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// typed via <see cref="BitATwinClaim"/> and exposed through the
    /// <see cref="IZ2AxisClaim.BitATwin"/> interface (wired 2026-05-25).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling: <see cref="F61BitAParityPi2Inheritance"/>
    /// (F1² on the Π²_X axis, i.e., Z⊗N as the bit_a parity operator) when
    /// constructed via the registry; <c>null</c> when constructed directly
    /// without the optional F61 ctor parameter (unit-test backward compat).
    /// Wired 2026-05-25 as a Schicht-1 closure of the previously unfilled
    /// twin slot. No cycle: F61 → F63 → F38 chain does not pass through F1.</summary>
    public Claim? BitATwin => BitATwinClaim;

    /// <summary>Filled when F1Pi2Inheritance was constructed with the F61 ctor
    /// parameter (the registry path); TrivialNotYetTyped when constructed
    /// without (the legacy unit-test path).</summary>
    public BitATwinClassification BitATwinStatus =>
        BitATwinClaim is not null
            ? BitATwinClassification.Filled
            : BitATwinClassification.TrivialNotYetTyped;
    /// <inheritdoc />
    /// <remarks>Parent role: feeds the F99 inheritance graph structurally
    /// but the claim's own value does not sit on the F86b α-axis.
    /// F1 palindrome shift is dimensional (the "2" in Π·L·Π⁻¹ = -L - 2Σγ·I
    /// is a_0 on the dyadic ladder = qubit dimension), not an F86b α value.</remarks>
    public F99AnchorRole F99Role => F99AnchorRole.Parent;

    /// <inheritdoc />
    public IReadOnlyList<double> F99AnchorValues { get; } = Array.Empty<double>();

    public Pi2DyadicLadderClaim Ladder { get; }
    public Pi2I4MemoryLoopClaim MemoryLoop { get; }

    /// <summary>F1 master palindrome identity, the typed parent. This claim's
    /// entire purpose is to lift F1's "2" coefficient and "−1" sign flip into the
    /// Pi2-Foundation lineage; before 2026-05-16 the F1PalindromeIdentity was
    /// referenced only by file-path string in the anchor field (a "ghost
    /// inheritance"). Added 2026-05-16 (Wave 6a) as a typed ctor parent: this
    /// is the single highest-leverage closure in the inheritance map sweep
    /// because F41, F43, F44, F68, F78, F79, F87 all take
    /// <see cref="F1Pi2Inheritance"/> as a typed parent; promoting F1 here
    /// brings them all into the F1 PalindromeIdentity lineage transitively
    /// without per-class ctor changes.
    ///
    /// <para><b>Bidirectional leverage</b>: the same hub property cuts both ways.
    /// Welle 7 (2026-05-25) wired <see cref="F61BitAParityPi2Inheritance"/> as F1's
    /// bit_a-twin ctor parent; that single edge propagates the entire BitA-axis
    /// chain (F61 → F63 → F38 → Pi2OperatorSpaceMirror → F88b → KleinFour →
    /// PolarityLayerOrigin → QubitDim → PolynomialFoundation) into all seven
    /// F-children listed above as transitive ancestors. The Welle 9 cleanup
    /// (2026-05-26 commit 79ec276) documented this by expanding eight test
    /// BuildBaseRegistry helpers that broke when their pre-Welle-7 minimal
    /// registries no longer satisfied F1's new dependency surface. F87 was the
    /// empirical addition to the prophesied list of children: F1's own
    /// docstring (written pre-F87) predicted F41/F43/F44/F68/F78/F79; the
    /// test-breakage swept up F87 as the seventh.</para></summary>
    public F1.F1PalindromeIdentity F1 { get; }

    /// <summary>F1's bit_a-twin: F61BitAParityPi2Inheritance (F1² on Π²_X axis,
    /// Z⊗N as bit_a parity operator). Wired 2026-05-25 to close the previously
    /// TrivialNotYetTyped BitA-twin slot. No cycle: F61's ctor chain (F61 → F63
    /// → F38) does not pass through F1. Nullable: legacy unit tests construct
    /// without F61 (BitATwinStatus stays TrivialNotYetTyped in that path); the
    /// registry-built F1Pi2Inheritance always has F61 wired (Filled). Renamed
    /// from F61BitATwin to BitATwinClaim in Welle 8 (2026-05-26) to align with
    /// the F38/F39/F63/X-Mirror BitB Claims' naming convention.</summary>
    public F61BitAParityPi2Inheritance? BitATwinClaim { get; }

    /// <summary>The "2" coefficient in F1's "−L − 2Σγ·I" closed form. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = d (qubit dimension).</summary>
    public double TwoFactor => Ladder.Term(0);

    /// <summary>The "−1" sign flip in F1's "−L" reading lives on the Z₄ memory
    /// loop: <c>i² = −1</c>, two 90° rotations summed (Pi2I4MemoryLoop Layer 1).
    /// Live computation through <see cref="Pi2I4MemoryLoopClaim.PowerOfI"/>(2).</summary>
    public Complex SignFlipFromZ4 => MemoryLoop.PowerOfI(2);

    public F1Pi2Inheritance(
        F1.F1PalindromeIdentity f1,
        Pi2DyadicLadderClaim ladder,
        Pi2I4MemoryLoopClaim memoryLoop,
        F61BitAParityPi2Inheritance? f61BitATwin = null)
        : base("F1 palindrome's 2 coefficient inherits from Pi2-Foundation (2 = a_0); sign flip = Z₄ i²; master identity = F1PalindromeIdentity; bit_a-twin = F61",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F1 + docs/proofs/MIRROR_SYMMETRY_PROOF.md + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (typed master parent) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs (typed bit_a-twin, registry-only)")
    {
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        MemoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
        BitATwinClaim = f61BitATwin;
    }

    public override string DisplayName =>
        "F1's 2 coefficient as Pi2-Foundation a_0; sign flip as Z₄ i²";

    public override string Summary =>
        $"Π·L·Π⁻¹ = −L − 2Σγ·I: 2 coefficient = a_0 = {TwoFactor}; −1 sign flip = i² on Z₄ memory loop ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F1 closed form",
                summary: "Π·L·Π⁻¹ = −L − 2Σγ·I (Tier1Derived in MIRROR_SYMMETRY_PROOF)");
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = d)", TwoFactor);
            yield return new InspectableNode("Z₄ Layer 1 reading",
                summary: "The −1 sign flip in '−L' = i² = two 90° rotations on Pi2I4MemoryLoop; F1 is named Layer 1 in Pi2I4MemoryLoop's own docstring");
            yield return InspectableNode.RealScalar("SignFlipFromZ4 (= i² = −1)", SignFlipFromZ4.Real);
            yield return new InspectableNode("Z₂ inversion",
                summary: "+L and −L − 2σ·I are Π-partners; two Π applications close back; the 2 is the qubit dimension d, not a free parameter");
        }
    }
}
