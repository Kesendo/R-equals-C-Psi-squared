using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F71MirrorBlockRefinement"/>: Z₂ refinement of the
/// joint-popcount sectors via the chain spatial-mirror operator P_F71. Has exactly one typed
/// parent — <see cref="JointPopcountSectors"/> — whose joint-popcount block-diagonal
/// structure F71 refines further.</summary>
public static class F71MirrorBlockRefinementRegistration
{
    public static ClaimRegistryBuilder RegisterF71MirrorBlockRefinement(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F71MirrorBlockRefinement>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            return new F71MirrorBlockRefinement(sectors);
        });
}
