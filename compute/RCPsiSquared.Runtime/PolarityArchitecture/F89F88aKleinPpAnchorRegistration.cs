using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89F88aKleinPpAnchor"/>: the F89-to-F88a/F88b
/// bridge claim that F89's bond Hamiltonian (XX, YY) sits in F88a's Klein Pp cell and F89's
/// ρ_cc initial state is the F88b popcount-(1, 2) coherence pair. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="KleinFourCellClaim"/>: the F88a two-axis Π² decomposition; F89's bond
///         bilinears XX and YY are Pp-resident (Π²_Z = Π²_X = +1).</item>
///   <item><see cref="F89TopologyOrbitClosure"/>: F89's bond-graph orbit closure.</item>
/// </list>
///
/// <para>Tier consistency: all three Tier1Derived.</para>
///
/// <para>Requires upstream registrations: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers KleinFourCellClaim) + <c>RegisterF89TopologyOrbitClosure</c>.</para></summary>
public static class F89F88aKleinPpAnchorRegistration
{
    public static ClaimRegistryBuilder RegisterF89F88aKleinPpAnchor(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89F88aKleinPpAnchor>(b =>
        {
            var klein = b.Get<KleinFourCellClaim>();
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89F88aKleinPpAnchor(klein, f89);
        });
}
