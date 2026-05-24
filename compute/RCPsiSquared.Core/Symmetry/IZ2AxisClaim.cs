using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Marker interface for Claims that carry a Z₂-axis classification within
/// the cubic Z₂³ polarity architecture. Implemented by every Pi²-Inheritance
/// Claim and by the Klein-cell decomposition Claims.
///
/// <para>The optional <see cref="BitATwin"/> pointer is meaningful only when
/// <see cref="Z2Axis"/> = <see cref="Symmetry.Z2Axis.BitB"/>: it points to the typed
/// bit_a sibling Claim (= the same theorem read on the Π²_X axis), or is
/// <c>null</c> when no such twin is currently typed (= unfilled twin slot,
/// a gap in the cubic-architecture coverage). Parent-edge inheritance (e.g.,
/// F62 listing F61 as a constructor parent) is NOT twinship; only a Claim
/// whose statement is the bit_a mirror of the BitB statement counts as a
/// twin.</para></summary>
public interface IZ2AxisClaim
{
    /// <summary>The Z₂-axis classification of this Claim within the cubic Z₂³
    /// polarity architecture. Required; no default.</summary>
    Z2Axis Z2Axis { get; }

    /// <summary>For <see cref="Symmetry.Z2Axis.BitB"/> Claims: the typed bit_a-twin
    /// sibling Claim if one exists, or <c>null</c> if the twin slot is unfilled.
    /// For all other axes: always <c>null</c>.</summary>
    Claim? BitATwin { get; }
}
