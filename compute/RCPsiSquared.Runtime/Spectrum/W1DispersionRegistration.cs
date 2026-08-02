using RCPsiSquared.Core.Spectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Spectrum;

/// <summary>Registers <see cref="W1Dispersion"/>, the Tier1Derived closed-form dispersion of
/// the (0,1) coherence block (the span of the |0⟩⟨j| between the ferromagnet and the single
/// excitations) for a Heisenberg chain under uniform per-site Z-dephasing.
///
/// <para>Lives in its own <c>Runtime/Spectrum/</c> namespace rather than under
/// <c>F86Main/</c> because the block is not F86-specific: it stays exactly closed for any
/// XXZ anisotropy, and its uniform decay rate 2γ is Δ-blind, but the ω_k band itself needs
/// |Δ| = 1 (Δ = 1 supplies the Laplacian's degree diagonal; Δ = −1 gives the SIGNLESS
/// Laplacian, cospectral with it on a bipartite graph, so the open chain returns the same
/// band; Δ = 0, 0.5, 2 all fail; see <see cref="W1Dispersion"/>). F86's JW track is one consumer (the
/// block's decay rate 2γ anchors <c>JwBondQPeakUnified</c>'s NEW-NEW Lorentzian width); other
/// F-theorems with chain-Hamiltonian content can register edges to <see cref="W1Dispersion"/>
/// without inheriting an F86 dependency.</para>
///
/// <para>Parameterised by (N, J, γ): the registry holds one concrete dispersion. Multiple
/// parameter points require multiple registries (mirrors the <see cref="F86Main"/> pattern
/// for parameterised laws like <c>TPeakLaw(γ)</c> and <c>QEpLaw(g_eff)</c>).</para></summary>
public static class W1DispersionRegistration
{
    public static ClaimRegistryBuilder RegisterW1Dispersion(
        this ClaimRegistryBuilder builder,
        int N,
        double J,
        double gammaZero) =>
        builder.Register<W1Dispersion>(_ => new W1Dispersion(N, J, gammaZero));
}
