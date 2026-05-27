using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2KleinV4DephaseSwapGroup"/> — the Klein-V₄
/// subgroup of unitary involutions on the 4^N Pauli basis realizing the dephase-letter
/// Klein V₄ {I, Z↔Y, Z↔X, Y↔X} on the F1 palindrome family {Π_Z, Π_X, Π_Y}.
///
/// <para>Standalone Tier1Derived primitive: no ctor parents. The Claim provides the
/// operator-space lift of the Klein V₄ on dephase letters and is consumed downstream
/// by F1-family transfer arguments (any F1 identity proven under one dephase letter
/// transfers to the other two via unitary conjugation by the appropriate element of
/// {I, D, Q_zx, H}). Pattern follows <see cref="Pi2DyadicLadderRegistration"/>
/// (lambda-built, no <c>b.Get&lt;X&gt;()</c> edges).</para>
///
/// <para>Welle 12 Task 3 (2026-05-27). Anchored on
/// <c>docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md</c> (Welle 12 Task 1) and
/// <c>docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c> (Welle 12 Task 2).
/// </para></summary>
public static class Pi2KleinV4DephaseSwapGroupRegistration
{
    public static ClaimRegistryBuilder RegisterPi2KleinV4DephaseSwapGroup(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2KleinV4DephaseSwapGroup>(_ => new Pi2KleinV4DephaseSwapGroup());
}
