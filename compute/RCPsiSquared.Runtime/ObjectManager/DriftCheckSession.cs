using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>On-demand drift validation. Iterates every Claim that implements
/// <see cref="IDriftCheckable"/> and collects its <see cref="DriftReport"/>. Does not
/// throw; the caller decides whether drift is fatal.</summary>
public sealed class DriftCheckSession
{
    private readonly ClaimRegistry _registry;

    public DriftCheckSession(ClaimRegistry registry)
    {
        _registry = registry;
    }

    public IReadOnlyList<DriftReport> VerifyAll() =>
        _registry.All()
            .OfType<IDriftCheckable>()
            .Select(c => c.Verify())
            .ToList();

    public DriftReport Verify<T>() where T : Claim
    {
        var claim = _registry.Get<T>();
        if (claim is not IDriftCheckable checkable)
            throw new InvalidOperationException(
                $"Claim {typeof(T).Name} does not implement IDriftCheckable.");
        return checkable.Verify();
    }
}
