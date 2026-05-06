namespace RCPsiSquared.Core.Knowledge;

/// <summary>Optional interface for Claims with pinned witness tables. The
/// <see cref="DriftReport"/> returned by <see cref="Verify"/> says whether the pinned
/// baseline matches the live recomputation. The Runtime's <c>DriftCheckSession</c> iterates
/// every Claim implementing this interface and collects the reports.
///
/// <para><b>Implementer's contract: Verify must compare against an INDEPENDENT source of
/// truth.</b> A common failure mode is for <c>Verify</c> to recompute the same closed-form
/// expression that the Claim's own property already uses; this is structurally tautological
/// and cannot detect drift. A correct <c>Verify</c> either: (a) builds the underlying object
/// (e.g. the Lindbladian) from primitives and computes the witness from scratch, or (b)
/// reads a pinned table populated from a different code path (e.g. hardware confirmation,
/// historical log) and compares against the live property. The first F73 implementation in
/// <c>PalindromeResidualScalingClaim</c> is reflexive-only and is documented as a known
/// limitation; replacing it with a Lindbladian-based check is a follow-up.</para></summary>
public interface IDriftCheckable
{
    DriftReport Verify();
}
