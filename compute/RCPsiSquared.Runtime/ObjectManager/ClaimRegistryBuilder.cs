using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>γ-style builder. Each <see cref="Register{T}"/> call records a factory lambda;
/// the lambda's <c>b.Get&lt;X&gt;()</c> calls double as edge declarations. <see cref="Build"/>
/// performs topological resolution by deferred construction: factories that touch unresolved
/// dependencies throw <see cref="DependencyNotYetAvailableException"/>; the builder defers
/// them and retries until every Claim is resolved or no progress is made (cycle).</summary>
public sealed class ClaimRegistryBuilder
{
    private readonly Dictionary<Type, Func<IBuilderContext, Claim>> _factories = new();
    private bool _built;

    /// <summary>Records a factory for <typeparamref name="T"/>. The factory receives an
    /// <see cref="IBuilderContext"/> whose <c>Get&lt;X&gt;()</c> calls are the sole mechanism
    /// for declaring that T depends on X.</summary>
    public ClaimRegistryBuilder Register<T>(Func<IBuilderContext, T> factory) where T : Claim
    {
        if (_built)
            throw new InvalidOperationException("ClaimRegistry already built; the builder is single-shot.");
        _factories[typeof(T)] = ctx => factory(ctx);
        return this;
    }

    /// <summary>Resolves all registered factories in topological order by deferred construction.
    /// Throws <see cref="InvariantViolationException"/> (Rule = "Cycle") when a circular
    /// dependency prevents forward progress.</summary>
    public ClaimRegistry Build()
    {
        if (_built)
            throw new InvalidOperationException("ClaimRegistry already built; use the returned registry.");
        _built = true;

        var resolved = new Dictionary<Type, Claim>();
        var edges = new List<Edge>();
        var order = new List<Type>();
        var pending = new HashSet<Type>(_factories.Keys);

        while (pending.Count > 0)
        {
            int progress = 0;
            // Snapshot to allow modification of pending while iterating.
            foreach (var type in pending.ToList())
            {
                Type currentlyResolving = type;
                var ctx = new RecordingBuilderContext(resolved, currentlyResolving, edges);
                try
                {
                    var claim = _factories[type](ctx);
                    resolved[type] = claim;
                    order.Add(type);
                    pending.Remove(type);
                    progress++;
                }
                catch (DependencyNotYetAvailableException)
                {
                    // Defer this factory; the missing dependency may resolve in a later pass.
                }
            }

            if (progress == 0)
                throw new InvariantViolationException(
                    rule: "Cycle",
                    message: "Dependency cycle detected among: " + string.Join(", ", pending.Select(t => t.Name)),
                    hint: "Restructure: split a Claim, or invert one edge.",
                    path: pending.ToArray());
        }

        return new ClaimRegistry(resolved, edges, order);
    }

    /// <summary>The recording context: returns resolved instances and records an edge,
    /// or throws <see cref="DependencyNotYetAvailableException"/> if the dependency is
    /// registered but not yet resolved.</summary>
    private sealed class RecordingBuilderContext : IBuilderContext
    {
        private readonly Dictionary<Type, Claim> _resolved;
        private readonly Type _resolvingType;
        private readonly List<Edge> _edges;

        public RecordingBuilderContext(
            Dictionary<Type, Claim> resolved, Type resolvingType, List<Edge> edges)
        {
            _resolved = resolved;
            _resolvingType = resolvingType;
            _edges = edges;
        }

        public T Get<T>() where T : Claim
        {
            if (_resolved.TryGetValue(typeof(T), out var claim))
            {
                _edges.Add(new Edge(
                    Parent: typeof(T),
                    Child: _resolvingType,
                    Reason: $"{_resolvingType.Name} requested {typeof(T).Name} via b.Get<{typeof(T).Name}>()"));
                return (T)claim;
            }
            throw new DependencyNotYetAvailableException(typeof(T));
        }
    }
}
