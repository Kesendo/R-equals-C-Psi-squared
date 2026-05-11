using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="InhomogeneousGammaF71BreakingWitness"/>: the
/// γ-asymmetric F71-breaking witness. Has TWO typed parents — <see cref="JointPopcountSectors"/>
/// (γ-blind: joint-popcount block-diagonality independent of γ_l) and
/// <see cref="F71MirrorBlockRefinement"/> (γ-symmetric: exact only when γ_l = γ_{N-1-l}).
/// Together they predict that F71-refinement off-block Frobenius scales with the F71
/// asymmetry norm of the γ-distribution.</summary>
public static class InhomogeneousGammaF71BreakingWitnessRegistration
{
    public static ClaimRegistryBuilder RegisterInhomogeneousGammaF71BreakingWitness(
        this ClaimRegistryBuilder builder) =>
        builder.Register<InhomogeneousGammaF71BreakingWitness>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            var f71 = b.Get<F71MirrorBlockRefinement>();
            return new InhomogeneousGammaF71BreakingWitness(sectors, f71);
        });
}
