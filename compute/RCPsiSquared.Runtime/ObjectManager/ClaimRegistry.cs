using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Immutable result of a successful <see cref="ClaimRegistryBuilder.Build"/>. Holds
/// every registered Claim by exact runtime type, the typed Edges between them, and the
/// topological order in which they were resolved. Cross-query helpers (<see cref="AllOfTier"/>,
/// <see cref="AncestorsOf{T}"/>, etc.) are implemented in Task 12.</summary>
public sealed class ClaimRegistry
{
    private readonly Dictionary<Type, Claim> _claims;
    private readonly IReadOnlyList<Edge> _edges;
    private readonly IReadOnlyList<Type> _topologicalOrder;

    internal ClaimRegistry(
        Dictionary<Type, Claim> claims,
        IReadOnlyList<Edge> edges,
        IReadOnlyList<Type> topologicalOrder)
    {
        _claims = claims;
        _edges = edges;
        _topologicalOrder = topologicalOrder;
    }

    public T Get<T>() where T : Claim
    {
        if (!_claims.TryGetValue(typeof(T), out var claim))
            throw new KeyNotFoundException($"Claim {typeof(T).Name} is not registered.");
        return (T)claim;
    }

    public Claim Get(Type claimType)
    {
        if (!_claims.TryGetValue(claimType, out var claim))
            throw new KeyNotFoundException($"Claim {claimType.Name} is not registered.");
        return claim;
    }

    public bool TryGet<T>(out T claim) where T : Claim
    {
        if (_claims.TryGetValue(typeof(T), out var c))
        {
            claim = (T)c;
            return true;
        }
        claim = default!;
        return false;
    }

    public bool Contains<T>() where T : Claim => _claims.ContainsKey(typeof(T));

    public IEnumerable<Claim> All() => _claims.Values;

    public IReadOnlyList<Type> TopologicalOrder => _topologicalOrder;

    public IEnumerable<Edge> AllEdges() => _edges;

    public IReadOnlyList<Claim> AllOfTier(Tier t) =>
        _claims.Values.Where(c => c.Tier == t).ToList();

    public IReadOnlyList<Claim> AncestorsOf<T>() where T : Claim =>
        AncestorsOf(typeof(T));

    public IReadOnlyList<Claim> AncestorsOf(Type childType)
    {
        var seen = new HashSet<Type>();
        var queue = new Queue<Type>();
        foreach (var e in _edges.Where(e => e.Child == childType))
            queue.Enqueue(e.Parent);
        while (queue.Count > 0)
        {
            var t = queue.Dequeue();
            if (!seen.Add(t)) continue;
            foreach (var e in _edges.Where(e => e.Child == t))
                queue.Enqueue(e.Parent);
        }
        return seen.Select(t => _claims[t]).ToList();
    }

    public IReadOnlyList<Claim> DescendantsOf<T>() where T : Claim =>
        DescendantsOf(typeof(T));

    public IReadOnlyList<Claim> DescendantsOf(Type parentType)
    {
        var seen = new HashSet<Type>();
        var queue = new Queue<Type>();
        foreach (var e in _edges.Where(e => e.Parent == parentType))
            queue.Enqueue(e.Child);
        while (queue.Count > 0)
        {
            var t = queue.Dequeue();
            if (!seen.Add(t)) continue;
            foreach (var e in _edges.Where(e => e.Parent == t))
                queue.Enqueue(e.Child);
        }
        return seen.Select(t => _claims[t]).ToList();
    }

    public IReadOnlyList<Edge> EdgesFrom<T>() where T : Claim =>
        _edges.Where(e => e.Parent == typeof(T)).ToList();

    public IReadOnlyList<Edge> EdgesInto<T>() where T : Claim =>
        _edges.Where(e => e.Child == typeof(T)).ToList();

    public IReadOnlyList<OpenQuestion> OpenQuestionsBlockedBy<T>() where T : Claim =>
        DescendantsOf<T>().OfType<OpenQuestion>().ToList();
}
