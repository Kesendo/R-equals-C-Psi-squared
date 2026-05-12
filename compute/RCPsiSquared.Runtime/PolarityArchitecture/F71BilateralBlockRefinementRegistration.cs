using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F71BilateralBlockRefinement"/>: Z₂ × Z₂ refinement
/// of the joint-popcount sectors via the two independent F71 mirrors (F71_col on the bra side,
/// F71_row on the ket side). Has TWO typed parents — <see cref="JointPopcountSectors"/>
/// (the γ-blind block-diagonal foundation that the bilateral split refines further) and
/// <see cref="F71MirrorBlockRefinement"/> (the diagonal Z₂ within this Z₂ × Z₂; F71-even =
/// (+,+) ⊕ (−,−), F71-odd = (+,−) ⊕ (−,+); the bilateral refinement is strictly finer).</summary>
public static class F71BilateralBlockRefinementRegistration
{
    public static ClaimRegistryBuilder RegisterF71BilateralBlockRefinement(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F71BilateralBlockRefinement>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            var f71Diagonal = b.Get<F71MirrorBlockRefinement>();
            return new F71BilateralBlockRefinement(sectors, f71Diagonal);
        });
}
