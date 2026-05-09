using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F82T1AmplitudeDampingPi2Inheritance"/>:
/// F82's <c>Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{T1, odd}</c> as F81-corollary
/// with T1 amplitude damping correction. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (coefficient) and <c>a_{2−N} = 2^(N−1)</c> (scaling factor).</item>
///   <item><see cref="F81Pi2Inheritance"/>: mother claim. F82 reduces to
///         F81 exactly at γ_T1 = 0; the typed mother-corollary edge.
///         Pattern parallel to F77 → F75 mother claim wired today.</item>
/// </list>
///
/// <para>Tier consistency: F82 is Tier 1 proven; verified bit-exact at
/// N=2..5 (5e-16 residual). All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// the F81 chain (Pi2OperatorSpaceMirror + Pi2I4MemoryLoop +
/// F88PopcountCoherence + F88StaticDyadicAnchor + F81Pi2Inheritance).</para></summary>
public static class F82T1AmplitudeDampingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF82T1AmplitudeDampingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F82T1AmplitudeDampingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f81 = b.Get<F81Pi2Inheritance>();
            return new F82T1AmplitudeDampingPi2Inheritance(ladder, f81);
        });
}
