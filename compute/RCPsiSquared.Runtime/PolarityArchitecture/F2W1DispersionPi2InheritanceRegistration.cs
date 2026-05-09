using RCPsiSquared.Core.Spectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F2W1DispersionPi2Inheritance"/>:
/// F2 w=1 Liouvillian dispersion ω_k = 4J·(1−cos(πk/N)). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (BandwidthPrefactor) and <c>a_0 = 2</c> (HoppingFactor).</item>
///   <item><see cref="W1Dispersion"/>: typed Spectrum primitive (Tier1Derived
///         in D10) with the same closed form. F2 inherits from W1Dispersion as
///         the F-formula registry-formal entry; W1Dispersion is the typed
///         primitive shared by the F86 JW track and other consumers. F2 → W1Dispersion
///         closes a long-standing gap (W1Dispersion had 0 F-formula descendants).</item>
/// </list>
///
/// <para>Tier consistency: F2 is Tier 1 proven D10; both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="W1DispersionRegistration.RegisterW1Dispersion"/>(N, J, γ).</para></summary>
public static class F2W1DispersionPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF2W1DispersionPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F2W1DispersionPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var w1 = b.Get<W1Dispersion>();
            return new F2W1DispersionPi2Inheritance(ladder, w1);
        });
}
