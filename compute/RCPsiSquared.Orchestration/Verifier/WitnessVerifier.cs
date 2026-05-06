using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Verifier;

/// <summary>Layer 3 consumer 2 (verifier). Wraps <see cref="DriftCheckSession"/> and
/// exposes a structured <see cref="VerifyResult"/> with totals, drifted count, and the
/// individual reports. Used as an xunit assertion in CI: drift count = 0 means the pinned
/// witness tables are still in sync with their live derivations.</summary>
public sealed class WitnessVerifier
{
    private readonly DriftCheckSession _session;

    public WitnessVerifier(ClaimRegistry registry)
    {
        _session = new DriftCheckSession(registry);
    }

    public VerifyResult VerifyAll()
    {
        var reports = _session.VerifyAll();
        var drifted = reports.Count(r => r.IsDrift);
        return new VerifyResult(reports.Count, drifted, reports);
    }

    public DriftReport Verify<T>() where T : Claim => _session.Verify<T>();
}
