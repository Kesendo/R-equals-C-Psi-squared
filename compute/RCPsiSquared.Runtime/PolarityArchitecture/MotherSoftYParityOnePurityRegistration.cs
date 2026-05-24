using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="MotherSoftYParityOnePurity"/> (F109).
///
/// <para>Standalone Claim: no constructor parents. Sister to
/// <see cref="TrulyYParityZeroPurity"/> (F107) on the same y_par axis; together
/// they pin the y_par signature of truly + mother-soft cells across all dephase
/// letters. Sixth YParity-axis member.</para></summary>
public static class MotherSoftYParityOnePurityRegistration
{
    public static ClaimRegistryBuilder RegisterMotherSoftYParityOnePurity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MotherSoftYParityOnePurity>(_ => new MotherSoftYParityOnePurity());
}
