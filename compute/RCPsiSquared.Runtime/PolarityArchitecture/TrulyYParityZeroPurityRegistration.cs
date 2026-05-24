using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="TrulyYParityZeroPurity"/> (F107).
///
/// <para>Standalone Claim: no constructor parents (fifth YParity-axis member after
/// F102 + F103 N4K3 + F105 N5K3 + F106 N4K4). The first DERIVED-not-EMPIRICAL Claim
/// in the cubic Z₂³ architecture: F103/F105/F106 carry frozen empirical counts; F107
/// carries a closed-form theorem statement with per-dephase truly criteria and a
/// verification helper. Registered into the typed-knowledge graph so it is visible to
/// <see cref="PolarityCubeMap"/> aggregation and the inspector.</para></summary>
public static class TrulyYParityZeroPurityRegistration
{
    public static ClaimRegistryBuilder RegisterTrulyYParityZeroPurity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TrulyYParityZeroPurity>(_ => new TrulyYParityZeroPurity());
}
