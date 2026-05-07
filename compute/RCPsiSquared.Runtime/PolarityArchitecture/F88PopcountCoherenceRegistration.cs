using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers <see cref="PopcountCoherenceClaim"/> with dual-parent inheritance:
/// at the operator level it descends from <see cref="KleinFourCellClaim"/> (F88's 4-cell
/// Π² decomposition is the operator-level root); at the state level it descends from
/// <see cref="PolarityLayerOriginClaim"/> (the +0/-0 polarity content of popcount-mirror
/// states is the state-level root). Both paths trace back to
/// <see cref="PolynomialFoundationClaim"/> via QubitDimensionalAnchorClaim.
///
/// <para>The Builder records both b.Get&lt;X&gt;() calls as separate Edges, and AncestorsOf
/// returns both paths in the transitive closure.</para></summary>
public static class F88PopcountCoherenceRegistration
{
    public static ClaimRegistryBuilder RegisterF88PopcountCoherence(this ClaimRegistryBuilder builder) =>
        builder.Register<PopcountCoherenceClaim>(b =>
        {
            _ = b.Get<KleinFourCellClaim>();        // operator-level parent
            _ = b.Get<PolarityLayerOriginClaim>();  // state-level parent
            return new PopcountCoherenceClaim();
        });
}
