using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Cubic-Z₂³-architecture inventory: aggregates the current set of
/// <see cref="IZ2AxisClaim"/> instances across the typed-knowledge layer and
/// exposes counts, ratios, and gap-lists. Inspectable via
/// <c>rcpsi inspect PolarityCubeMap</c>.
///
/// <para>The cubic structure has three independent Z₂ classifiers at the Pauli-term
/// level (bit_a, bit_b, Y-parity); only the first two surface as independent Π²
/// operators (per F34/QUBIT_NECESSITY). Y-parity becomes independent only at
/// k≥3-body terms; at k=2 the Z₂³ collapses to Klein-Vierergruppe Z₂².</para>
///
/// <para>The map's primary use is gap-detection: BitB Claims with null BitATwin
/// pointers indicate unfilled bit_a-twin slots (Stage 2a of the cubic-unpacking
/// arc).</para></summary>
public sealed class PolarityCubeMap : Claim
{
    public IReadOnlyList<IZ2AxisClaim> BitBClaims { get; }
    public IReadOnlyList<IZ2AxisClaim> BitAClaims { get; }
    public IReadOnlyList<IZ2AxisClaim> Klein2Claims { get; }
    public IReadOnlyList<IZ2AxisClaim> YParityClaims { get; }
    public IReadOnlyList<IZ2AxisClaim> Cubic3Claims { get; }
    public IReadOnlyList<IZ2AxisClaim> NotApplicableClaims { get; }

    public int TotalClaims => BitBClaims.Count + BitAClaims.Count + Klein2Claims.Count +
                              YParityClaims.Count + Cubic3Claims.Count + NotApplicableClaims.Count;

    public int OpenBitATwinSlots =>
        BitBClaims.Count(c => c.BitATwin is null);

    public int FilledBitATwinSlots =>
        BitBClaims.Count(c => c.BitATwin is not null);

    public double TwinCoverageRatio =>
        BitBClaims.Count == 0 ? 0.0 : (double)FilledBitATwinSlots / BitBClaims.Count;

    public IReadOnlyList<string> UnfilledTwinSlotNames =>
        BitBClaims.Where(c => c.BitATwin is null)
                  .Select(c => c.GetType().Name)
                  .OrderBy(n => n)
                  .ToList();

    /// <summary>BitB Claims whose BitA twin is mechanically derivable but not yet
    /// typed as a separate Claim. The 2026-05-25 F1↔F61 closure (followed by
    /// Welle 7's F38/F39/F63/X-Mirror twin wiring 2026-05-26) drained the
    /// historical examples; current count is typically 0 unless a unit-test
    /// constructs a BitB Claim directly (bypassing the registry path that wires
    /// the optional BitA-twin ctor parameter). Future low-cost work to fill.</summary>
    public int TrivialNotYetTypedTwinSlots =>
        BitBClaims.Count(c => c.BitATwinStatus == BitATwinClassification.TrivialNotYetTyped);

    /// <summary>BitB Claims whose BitA twin would require new substantive derivation
    /// (not a mechanical letter swap). Future medium-to-high effort to fill.</summary>
    public int NeedsDerivationTwinSlots =>
        BitBClaims.Count(c => c.BitATwinStatus == BitATwinClassification.NeedsDerivation);

    /// <summary>BitB Claims with no meaningful BitA twin (e.g., F1T1's amplitude
    /// damping intrinsically breaks bit_a symmetry per F61's BreakConditions, so
    /// no bit_a-axis analog exists). No filling work required.</summary>
    public int BitBSpecificTwinSlots =>
        BitBClaims.Count(c => c.BitATwinStatus == BitATwinClassification.BitBSpecific);

    /// <summary>BitB Claims whose BitA twin holds as a corollary of the global Hadamard
    /// X↔Z duality (PROOF_BIT_A_TWIN_VIA_HADAMARD.md) — Π/L spectrum, eigenspace, or
    /// operator-identity / Absorption-Theorem-popcount content. No bespoke typed twin owed.</summary>
    public int CoveredByHadamardDualityTwinSlots =>
        BitBClaims.Count(c => c.BitATwinStatus == BitATwinClassification.CoveredByHadamardDuality);

    /// <summary>Sub-list of unfilled BitA-twin slot names categorized as
    /// <see cref="BitATwinClassification.TrivialNotYetTyped"/>. These are the
    /// natural first batch for Stage 2a fill work.</summary>
    public IReadOnlyList<string> TrivialNotYetTypedTwinSlotNames =>
        BitBClaims.Where(c => c.BitATwinStatus == BitATwinClassification.TrivialNotYetTyped)
                  .Select(c => c.GetType().Name)
                  .OrderBy(n => n)
                  .ToList();

    public PolarityCubeMap(IReadOnlyList<IZ2AxisClaim> allClaims)
        : base("Polarity Cube Map: cubic Z₂³ architecture inventory across the Pi²-Inheritance claim set",
               // Tier2Empirical (not Tier1Derived) because the snapshot can include
               // Tier1Candidate parents (e.g., F1T1AmplitudeDampingPi2Inheritance) and
               // Tier2Empirical parents (e.g., Pi2KleinBilinearTable), and the builder
               // enforces parent.Tier >= child.Tier. The cubic Z₂³ architecture itself is
               // Tier1Derived (algebraic), but the inventory's strength is bounded by its
               // weakest contributor; demoting to Tier2Empirical makes the snapshot wirable
               // across the full Pi²-Inheritance claim set without artificially restricting
               // which Claims may participate.
               Tier.Tier2Empirical,
               "docs/PI2KB_INHERITANCE_MAP.md + " +
               "docs/SYMMETRY_FAMILY_INVENTORY.md")
    {
        if (allClaims is null) throw new ArgumentNullException(nameof(allClaims));

        BitBClaims = allClaims.Where(c => c.Z2Axis == Z2Axis.BitB).ToList();
        BitAClaims = allClaims.Where(c => c.Z2Axis == Z2Axis.BitA).ToList();
        Klein2Claims = allClaims.Where(c => c.Z2Axis == Z2Axis.Klein2).ToList();
        YParityClaims = allClaims.Where(c => c.Z2Axis == Z2Axis.YParity).ToList();
        Cubic3Claims = allClaims.Where(c => c.Z2Axis == Z2Axis.Cubic3).ToList();
        NotApplicableClaims = allClaims.Where(c => c.Z2Axis == Z2Axis.NotApplicable).ToList();
    }

    public override string DisplayName =>
        $"Polarity Cube Map ({TotalClaims} Claims; twin coverage {TwinCoverageRatio:P1})";

    public override string Summary =>
        $"BitB={BitBClaims.Count}, BitA={BitAClaims.Count}, Klein2={Klein2Claims.Count}, " +
        $"YParity={YParityClaims.Count}, Cubic3={Cubic3Claims.Count}; " +
        $"open BitA twin slots: {OpenBitATwinSlots} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Total Claims", TotalClaims);
            yield return InspectableNode.RealScalar("BitB count (F1² / Π²_Z axis)", BitBClaims.Count);
            yield return InspectableNode.RealScalar("BitA count (F61 / Π²_X axis)", BitAClaims.Count);
            yield return InspectableNode.RealScalar("Klein2 count (uses both Π² axes)", Klein2Claims.Count);
            yield return InspectableNode.RealScalar("YParity count (term-level k≥3)", YParityClaims.Count);
            yield return InspectableNode.RealScalar("Cubic3 count (full Z₂³)", Cubic3Claims.Count);
            yield return InspectableNode.RealScalar("Open BitA twin slots", OpenBitATwinSlots);
            yield return InspectableNode.RealScalar("Twin coverage ratio", TwinCoverageRatio);
            yield return new InspectableNode("BitA twin slot breakdown (Stage 2a)",
                summary: $"trivial-not-typed: {TrivialNotYetTypedTwinSlots}, " +
                         $"needs-derivation: {NeedsDerivationTwinSlots}, " +
                         $"bit_b-specific: {BitBSpecificTwinSlots}, " +
                         $"covered-by-Hadamard: {CoveredByHadamardDualityTwinSlots}");

            if (TrivialNotYetTypedTwinSlotNames.Count > 0)
            {
                yield return new InspectableNode("Trivial-not-yet-typed BitA twin slots (low-cost fill targets)",
                    summary: string.Join(", ", TrivialNotYetTypedTwinSlotNames.Take(10)) +
                             (TrivialNotYetTypedTwinSlotNames.Count > 10 ? $" ... (+{TrivialNotYetTypedTwinSlotNames.Count - 10} more)" : ""));
            }

            if (UnfilledTwinSlotNames.Count > 0)
            {
                yield return new InspectableNode("All unfilled BitA twin slots (top 10)",
                    summary: string.Join(", ", UnfilledTwinSlotNames.Take(10)) +
                             (UnfilledTwinSlotNames.Count > 10 ? $" ... (+{UnfilledTwinSlotNames.Count - 10} more)" : ""));
            }
        }
    }
}
