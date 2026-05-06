namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Internal sentinel: thrown by the recording <see cref="IBuilderContext"/> when
/// a factory calls <c>b.Get&lt;X&gt;()</c> and X is registered but not yet resolved. The
/// builder catches this, defers the current factory, and continues the resolution loop.</summary>
internal sealed class DependencyNotYetAvailableException : Exception
{
    public Type RequestedType { get; }
    public DependencyNotYetAvailableException(Type requestedType)
        : base($"Dependency {requestedType.Name} is registered but not yet resolved.")
    {
        RequestedType = requestedType;
    }
}
