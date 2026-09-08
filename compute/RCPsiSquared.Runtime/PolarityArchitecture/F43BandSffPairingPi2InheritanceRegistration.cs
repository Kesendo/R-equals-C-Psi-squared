using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F43BandSffPairingPi2Inheritance"/>.
///
/// <para>The Pi2 ladder supplies the coefficient 2 in the rate reflection
/// <c>d ↦ 2Nγ−d</c>. F1 supplies the multiplicity-preserving spectral
/// palindrome. F43 reads that exact identity on reflected decay-rate bands;
/// it does not introduce invariant fixed-XY-weight eigenvalue sectors. The
/// N+1 endpoint multiplicity is scoped to strictly positive uniform dephasing;
/// at γ=0 the endpoints collapse into the larger commutator kernel.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F43BandSffPairingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF43BandSffPairingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F43BandSffPairingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F43BandSffPairingPi2Inheritance(ladder, f1);
        });
}
