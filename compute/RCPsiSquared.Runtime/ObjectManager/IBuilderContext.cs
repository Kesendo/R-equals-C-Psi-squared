using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>The argument passed to factory lambdas inside <see cref="ClaimRegistryBuilder.Register{T}(Func{IBuilderContext, T})"/>.
/// Each call to <see cref="Get{T}"/> records an inheritance edge from T to the current
/// resolving Claim. Used at <see cref="ClaimRegistryBuilder.Build"/> time and made unavailable
/// after the registry is built.</summary>
public interface IBuilderContext
{
    T Get<T>() where T : Claim;
}
