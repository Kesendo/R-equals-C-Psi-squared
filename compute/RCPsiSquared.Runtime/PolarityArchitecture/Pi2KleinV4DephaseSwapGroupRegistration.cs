using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2KleinV4DephaseSwapGroup"/> — the Klein-V₄
/// subgroup of unitary involutions on the 4^N Pauli basis realizing the dephase-letter
/// Klein V₄ {I, Z↔Y, Z↔X, Y↔X} on the F1 palindrome family {Π_Z, Π_X, Π_Y}.
///
/// <para>Standalone Tier1Derived primitive: no ctor parents. The Claim provides the
/// operator-space lift of the Klein V₄ on dephase letters and is consumed downstream
/// by F1-family transfer arguments. Two transfer routes exist: Route 1 (per-axis
/// structural re-run, e.g. F112 Welle 13) transfers any identity depending only on
/// F38 Π_d² eigenvalue + Pauli-support disjointness between all three dephase letters.
/// Route 2 (Hadamard transport, only the {I, Q_zx} subgroup) lifts to a Hilbert-space
/// unitary U_H^⊗N and transports the Lindblad-form L itself; D and Q_yx (= H) are
/// operator-space-only and do not lift, so they only transport Π_d and Π_d²-graded
/// quantities, not L. Pattern follows <see cref="Pi2DyadicLadderRegistration"/>
/// (lambda-built, no <c>b.Get&lt;X&gt;()</c> edges).</para>
///
/// <para>Welle 12 Task 3 (2026-05-27). Anchored on
/// <c>docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md</c> (Welle 12 Task 1) and
/// <c>docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c> (Welle 12 Task 2).
/// The Welle 13 follow-up <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>
/// makes the per-axis transfer (Route 1) explicit for F112 and clarifies the
/// Hilbert-lift caveat on D and H (Route 2 partial).</para></summary>
public static class Pi2KleinV4DephaseSwapGroupRegistration
{
    public static ClaimRegistryBuilder RegisterPi2KleinV4DephaseSwapGroup(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2KleinV4DephaseSwapGroup>(_ => new Pi2KleinV4DephaseSwapGroup());
}
