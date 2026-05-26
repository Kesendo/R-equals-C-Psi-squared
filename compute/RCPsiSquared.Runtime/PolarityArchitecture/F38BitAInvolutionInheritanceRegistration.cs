using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F38BitAInvolutionInheritance"/>:
/// the BitA twin of <see cref="F38Pi2InvolutionPi2Inheritance"/> on the Π²_X axis.
///
/// <para>Standalone Claim (no ctor parents on the BitA side; the typed twin edge
/// from F38 BitB is wired via F38's optional ctor parameter). Must be registered
/// BEFORE <see cref="F38Pi2InvolutionPi2Inheritance"/> so the BitB Claim's
/// <c>b.Get&lt;F38BitAInvolutionInheritance&gt;()</c> resolves. Matches the F108
/// Part 1 ↔ Part 2 pattern.</para></summary>
public static class F38BitAInvolutionInheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF38BitAInvolutionInheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F38BitAInvolutionInheritance>(_ => new F38BitAInvolutionInheritance());
}
