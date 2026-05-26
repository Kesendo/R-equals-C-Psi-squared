using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="TrulyYParityZeroPurity"/> (F107).
///
/// <para>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F107 is the
/// first DERIVED-not-EMPIRICAL Claim in the cubic Z₂³ architecture; the parent
/// edge makes its anchor on the Z₂³ 8-cell decomposition explicit. Wired
/// 2026-05-26.</para></summary>
public static class TrulyYParityZeroPurityRegistration
{
    public static ClaimRegistryBuilder RegisterTrulyYParityZeroPurity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TrulyYParityZeroPurity>(b =>
            new TrulyYParityZeroPurity(b.Get<KleinEightCellClaim>()));
}
