using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F39DetPiBitAInheritance"/>:
/// the BitA twin of <see cref="F39DetPiPi2Inheritance"/> on the Π_X axis.
///
/// <para>Standalone Claim (no ctor parents on the BitA side; the typed twin edge
/// from F39 BitB is wired via F39's optional ctor parameter). Must be registered
/// BEFORE <see cref="F39DetPiPi2Inheritance"/> so the BitB Claim's
/// <c>b.Get&lt;F39DetPiBitAInheritance&gt;()</c> resolves.</para></summary>
public static class F39DetPiBitAInheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF39DetPiBitAInheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F39DetPiBitAInheritance>(_ => new F39DetPiBitAInheritance());
}
