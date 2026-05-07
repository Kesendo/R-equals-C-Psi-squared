using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Immutable result of a successful <see cref="ClaimRegistryBuilder.Build"/>. Holds
/// every registered Claim by exact runtime type, the typed Edges between them, and the
/// topological order in which they were resolved. Exposes lookup, iteration, and cross-query
/// helpers (<see cref="AllOfTier"/>, <see cref="AncestorsOf{T}"/>, <see cref="DescendantsOf{T}"/>,
/// <see cref="EdgesFrom{T}"/>, <see cref="EdgesInto{T}"/>, <see cref="OpenQuestionsBlockedBy{T}"/>).</summary>
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

    public int Count => _claims.Count;

    public IReadOnlyList<Type> TopologicalOrder => _topologicalOrder;

    public IEnumerable<Edge> AllEdges() => _edges;

    public IReadOnlyList<Claim> AllOfTier(Tier t) =>
        _claims.Values.Where(c => c.Tier == t).ToList();

    public IReadOnlyList<Claim> AncestorsOf<T>() where T : Claim => AncestorsOf(typeof(T));

    public IReadOnlyList<Claim> AncestorsOf(Type childType) =>
        Walk(seed: childType, edgeMatches: e => e.Child, edgeOther: e => e.Parent);

    public IReadOnlyList<Claim> DescendantsOf<T>() where T : Claim => DescendantsOf(typeof(T));

    public IReadOnlyList<Claim> DescendantsOf(Type parentType) =>
        Walk(seed: parentType, edgeMatches: e => e.Parent, edgeOther: e => e.Child);

    public IReadOnlyList<Edge> EdgesFrom<T>() where T : Claim =>
        _edges.Where(e => e.Parent == typeof(T)).ToList();

    public IReadOnlyList<Edge> EdgesInto<T>() where T : Claim =>
        _edges.Where(e => e.Child == typeof(T)).ToList();

    public IReadOnlyList<OpenQuestion> OpenQuestionsBlockedBy<T>() where T : Claim =>
        DescendantsOf<T>().OfType<OpenQuestion>().ToList();

    /// <summary>Generic BFS over recorded edges. <paramref name="edgeMatches"/> picks the
    /// endpoint that must equal the current frontier node; <paramref name="edgeOther"/>
    /// picks the endpoint to enqueue. Used both upward (Ancestors) and downward (Descendants)
    /// by swapping the two selectors.</summary>
    private IReadOnlyList<Claim> Walk(
        Type seed,
        Func<Edge, Type> edgeMatches,
        Func<Edge, Type> edgeOther)
    {
        var seen = new HashSet<Type>();
        var queue = new Queue<Type>();
        foreach (var e in _edges)
            if (edgeMatches(e) == seed) queue.Enqueue(edgeOther(e));
        while (queue.Count > 0)
        {
            var t = queue.Dequeue();
            if (!seen.Add(t)) continue;
            foreach (var e in _edges)
                if (edgeMatches(e) == t) queue.Enqueue(edgeOther(e));
        }
        return seen.Select(t => _claims[t]).ToList();
    }
}
