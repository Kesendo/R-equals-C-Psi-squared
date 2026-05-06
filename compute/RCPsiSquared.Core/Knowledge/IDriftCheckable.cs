namespace RCPsiSquared.Core.Knowledge;

/// <summary>Optional interface for Claims with pinned witness tables. The
/// <see cref="DriftReport"/> returned by <see cref="Verify"/> says whether the pinned
/// baseline matches the live recomputation. The Runtime's <c>DriftCheckSession</c> iterates
/// every Claim implementing this interface and collects the reports.</summary>
public interface IDriftCheckable
{
    DriftReport Verify();
}
