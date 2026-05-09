using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F4StationaryModeCountPi2Inheritance"/>:
/// F4's Clebsch-Gordan stationary-mode count Stat(N) = Σ_J m(J)·(2J+1)².
/// One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (IrrepDimensionCoefficient in 2J+1; same anchor as F1 TwoFactor,
///         F50 DegeneracyFactor, F66 UpperPoleCoefficient, F43
///         XorRateCoefficient).</item>
/// </list>
///
/// <para>Tier consistency: F4 is Tier 1 via Schur-Weyl SU(2)-irrep
/// decomposition; valid for Heisenberg Hamiltonian, Σγ = 0, all N.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F4StationaryModeCountPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF4StationaryModeCountPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F4StationaryModeCountPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F4StationaryModeCountPi2Inheritance(ladder);
        });
}
