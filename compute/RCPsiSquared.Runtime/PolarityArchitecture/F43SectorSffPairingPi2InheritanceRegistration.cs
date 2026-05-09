using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F43SectorSffPairingPi2Inheritance"/>:
/// F43 sector-level palindromic pairing K_freq(w) = K_freq(N−w). Two typed
/// parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (XorRateCoefficient in 2·N·γ XOR-sector decay rate; same anchor
///         as F1 TwoFactor, F50 DecayRateFactor, F66 UpperPoleCoefficient).</item>
///   <item><see cref="F1Pi2Inheritance"/>: F43 is the sector-level reading
///         of F1's palindromic identity. F1's Π conjugation acts on Pauli
///         strings; F43 reads its action at the sector-weight partition
///         (w → N−w identity at the spectral-statistics level).</item>
/// </list>
///
/// <para>F43 completes a four-claim F1-family of palindromic readings:</para>
/// <list type="bullet">
///   <item>F44: pair RATIO (eigenvalue level, artanh closed form)</item>
///   <item>F68: pair SUM (eigenvalue level, α_b + α_p = 2γ₀)</item>
///   <item>F41: pair DIFFERENCE (time domain, period t_Pi)</item>
///   <item>F43: SECTOR-WEIGHT pairing (operator-space partition, K_freq(w) = K_freq(N−w))</item>
/// </list>
///
/// <para>Tier consistency: F43 is Tier 1 proven D09 (one-line consequence of
/// Π's commutation with the sector projector at the spectral level); valid
/// for Heisenberg chain, Z-dephasing, all N. Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F43SectorSffPairingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF43SectorSffPairingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F43SectorSffPairingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F43SectorSffPairingPi2Inheritance(ladder, f1);
        });
}
