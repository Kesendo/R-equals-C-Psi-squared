namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Thrown by <see cref="ClaimRegistryBuilder.Build"/> when registered claims violate
/// a structural rule (Tier inheritance, missing parent, cycle, duplicate, anchor file missing).
/// Carries a <see cref="Rule"/> string so xunit asserts can match deterministically.</summary>
public sealed class InvariantViolationException : Exception
{
    public string Rule { get; }
    public IReadOnlyList<Type> Path { get; }
    public Type? OffendingClaim { get; }
    public string Hint { get; }

    public InvariantViolationException(string rule, string message, string hint,
        Type? offendingClaim = null, IReadOnlyList<Type>? path = null)
        : base(message)
    {
        Rule = rule;
        Path = path ?? Array.Empty<Type>();
        OffendingClaim = offendingClaim;
        Hint = hint;
    }
}
