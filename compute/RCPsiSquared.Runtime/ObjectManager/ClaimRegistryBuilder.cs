using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

public sealed class ClaimRegistryBuilder
{
    private readonly Dictionary<Type, Func<IBuilderContext, Claim>> _factories = new();

    public ClaimRegistryBuilder Register<T>(Func<IBuilderContext, T> factory) where T : Claim
    {
        _factories[typeof(T)] = ctx => factory(ctx);
        return this;
    }

    public ClaimRegistry Build()
    {
        // Stub: empty registry only. Real resolution lands in Task 5.
        if (_factories.Count > 0)
            throw new NotImplementedException("Resolution algorithm lands in Task 5.");
        return new ClaimRegistry(
            new Dictionary<Type, Claim>(),
            Array.Empty<Edge>(),
            Array.Empty<Type>());
    }
}
