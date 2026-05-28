namespace RCPsiSquared.Core.Symmetry;

/// <summary>Structural classification of the BitA-twin slot for a Z2-axis Claim
/// in the cubic Z₂³ polarity architecture.
///
/// <para>For Claims on the BitB axis (Π²_Z = X⊗N), the slot can be in one of four
/// states: <see cref="Filled"/>, <see cref="TrivialNotYetTyped"/>,
/// <see cref="NeedsDerivation"/>, or <see cref="BitBSpecific"/>. For Claims on any
/// other axis (BitA, Klein2, YParity, Cubic3, NotApplicable), the slot is
/// <see cref="NotApplicableForThisAxis"/>.</para>
///
/// <para>This refines the bare "BitATwin is null" gap-detection into a structural
/// breakdown, turning Stage 2a (61 open BitA-twin slots) from a raw open-count
/// into a categorized inventory of "what kind of work each gap represents".
/// Per-Claim authors classify their own gap; <see cref="PolarityCubeMap"/> aggregates.</para></summary>
public enum BitATwinClassification
{
    /// <summary>The BitA twin Claim exists and is registered. Equivalent to the
    /// existing <c>BitATwin is not null</c> state.</summary>
    Filled = 0,

    /// <summary>The BitA twin is a mechanical letter-swap mirror of the BitB Claim
    /// (Z↔X, X↔Z, or similar SU(2) rotation), trivially derivable but not yet
    /// implemented as a separate typed Claim. Future low-cost work to fill.
    /// Example: F1Pi2Inheritance's twin is F61BitAParityPi2Inheritance — F61 IS
    /// typed but not wired as F1's BitATwin reference; the work is mechanical.</summary>
    TrivialNotYetTyped = 1,

    /// <summary>The BitA twin would be a substantive new structural finding
    /// requiring new derivation work. Not a mechanical mirror; the bit_a-axis
    /// statement involves different algebra than the bit_b-axis one. Future
    /// medium-to-high effort to fill, possibly with empirical anchor.</summary>
    NeedsDerivation = 2,

    /// <summary>The BitB Claim has no meaningful bit_a-axis twin. Either the
    /// algebraic content is intrinsically tied to Z-dephasing / bit_b structure
    /// (e.g., F1T1's amplitude-damping breaks bit_a per F61's break conditions),
    /// or the bit_a-axis analog would be a trivial rephrasing of an already-typed
    /// claim. No filling work required.</summary>
    BitBSpecific = 3,

    /// <summary>The Claim is NOT on the BitB axis (it's on BitA, Klein2, YParity,
    /// Cubic3, or NotApplicable). The BitA-twin slot semantics does not apply.</summary>
    NotApplicableForThisAxis = 4,

    /// <summary>The BitA twin exists as a corollary of the global Hadamard X↔Z duality
    /// (<c>docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md</c>): the claim reduces to Π/L
    /// spectrum, eigenspace, or operator-identity / Absorption-Theorem-popcount content,
    /// so its bit_a image holds by Q_zx / Hadamard conjugation. No bespoke typed twin
    /// Claim is owed (distinct from <see cref="Filled"/>, which has a typed twin, and
    /// from <see cref="NeedsDerivation"/>, which is genuinely open).</summary>
    CoveredByHadamardDuality = 5,
}
