using RCPsiSquared.Core.F1;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1PalindromeIdentity"/>: the Liouvillian
/// palindrome identity Π·L·Π⁻¹ = −L − 2Σγ·I is the root F-chain claim, parameter-free
/// and Tier1Derived. Registered as a stand-alone edge because
/// <see cref="Pi2FamilyRegistration"/> covers only the Pi2-Foundation trunk
/// (PolynomialFoundation → QubitDimensionalAnchor → polarity / klein / half / quarter /
/// argmax-maxval primitives); F1 is the F-chain entry point and is consumed as a typed
/// parent by F80 and F81 (and transitively by their downstream F82 / F84 / ... chains).
///
/// <para>Required by: <see cref="F80FactorPi2InheritanceRegistration"/>,
/// <see cref="F81Pi2InheritanceRegistration"/>, and transitively by their downstream
/// Registrations (F82, F84, ...).</para></summary>
public static class F1PalindromeIdentityRegistration
{
    public static ClaimRegistryBuilder RegisterF1PalindromeIdentity(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1PalindromeIdentity>(_ => new F1PalindromeIdentity());
}
