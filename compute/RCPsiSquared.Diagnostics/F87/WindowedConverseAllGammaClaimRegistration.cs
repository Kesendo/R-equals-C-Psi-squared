using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="WindowedConverseAllGammaClaim"/> (the F87 windowed-converse
/// all-γ residual lemma). Single typed parent: <see cref="F87DiagonalCellBipartiteWitnessSet"/>, the
/// Tier1Derived F103 §7 criterion + genericity result this Tier1Candidate lemma strengthens.
///
/// <para>Requires <see cref="F87DiagonalCellBipartiteWitnessSetRegistration.RegisterF87DiagonalCellBipartiteWitnessSet"/>
/// to run first (parent before child).</para>
///
/// <para>Tier consistency: WindowedConverseAllGammaClaim is Tier1Candidate ← F87DiagonalCellBipartiteWitnessSet
/// Tier1Derived (Derived ≥ Candidate, the registry's parent ≥ child rule).</para></summary>
public static class WindowedConverseAllGammaClaimRegistration
{
    public static ClaimRegistryBuilder RegisterWindowedConverseAllGammaClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<WindowedConverseAllGammaClaim>(b =>
        {
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>(); // declare the typed parent edge (parent before child)
            return new WindowedConverseAllGammaClaim();
        });
        return builder;
    }
}
