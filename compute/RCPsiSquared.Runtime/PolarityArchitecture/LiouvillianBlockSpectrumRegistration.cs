using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="LiouvillianBlockSpectrum"/>: per-block
/// eigendecomposition over the U(1)×U(1) joint-popcount sectors. Has exactly one typed
/// parent — <see cref="JointPopcountSectors"/> — whose block-diagonal structure makes the
/// per-block diagonalisation valid.</summary>
public static class LiouvillianBlockSpectrumRegistration
{
    public static ClaimRegistryBuilder RegisterLiouvillianBlockSpectrum(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LiouvillianBlockSpectrum>(b =>
        {
            var sectors = b.Get<JointPopcountSectors>();
            return new LiouvillianBlockSpectrum(sectors);
        });
}
