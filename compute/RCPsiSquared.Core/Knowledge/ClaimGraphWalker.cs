using System.Collections.Concurrent;
using System.Reflection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Reflection-based traversal of the typed-Claim inheritance graph.
/// Walks any <see cref="Claim"/>'s public properties whose declared type is
/// <see cref="Claim"/> or a subclass, recursively, collecting every reachable
/// claim into a set.
///
/// <para>Every public property whose declared type derives from <see cref="Claim"/>
/// is an executable parent edge. The walker follows exactly those edges; a
/// same number, similar shape, filename, or prose cross-reference creates no
/// ancestry. Parentless claims therefore return only themselves. Focused
/// regressions cover a genuine application-to-formula edge and a separate
/// parentless formula.</para>
///
/// <para>Implementation: BFS over typed Claim properties. Each property
/// that returns a Claim subclass is followed; null returns are skipped.
/// Cycles are prevented by tracking visited references.</para>
/// </summary>
public static class ClaimGraphWalker
{
    /// <summary>Walk the parent graph starting from <paramref name="root"/>,
    /// returning every reachable <see cref="Claim"/> (including <paramref name="root"/>
    /// itself) in BFS order.</summary>
    public static IReadOnlyList<Claim> WalkReachable(Claim root)
    {
        if (root is null) throw new ArgumentNullException(nameof(root));
        var visited = new HashSet<Claim>(ReferenceEqualityComparer.Instance);
        var ordered = new List<Claim>();
        var queue = new Queue<Claim>();
        queue.Enqueue(root);
        visited.Add(root);

        while (queue.Count > 0)
        {
            var c = queue.Dequeue();
            ordered.Add(c);
            foreach (var parent in TypedClaimProperties(c))
            {
                if (visited.Add(parent)) queue.Enqueue(parent);
            }
        }
        return ordered;
    }

    /// <summary>Filter <see cref="WalkReachable"/> down to claims that
    /// implement a given marker interface <typeparamref name="TInterface"/>
    /// (e.g. IF99AnchorBearing). The combination "walk from root + filter
    /// by interface" recovers an anchor map for any anchor scheme from any
    /// graph entry point.</summary>
    public static IReadOnlyList<TInterface> ReachableImplementing<TInterface>(Claim root)
        where TInterface : class
        => WalkReachable(root).OfType<TInterface>().ToList();

    // Property list per Type never changes; cache to avoid re-scanning on
    // every BFS visit. Test runs that walk the graph from multiple roots
    // visit the same claim types repeatedly.
    private static readonly ConcurrentDictionary<Type, PropertyInfo[]> _propCache = new();

    private static IEnumerable<Claim> TypedClaimProperties(Claim c)
    {
        PropertyInfo[] props = _propCache.GetOrAdd(c.GetType(), t =>
            t.GetProperties(BindingFlags.Instance | BindingFlags.Public)
             .Where(p => typeof(Claim).IsAssignableFrom(p.PropertyType)
                         && p.GetIndexParameters().Length == 0)
             .ToArray());
        foreach (var prop in props)
        {
            Claim? value;
            // Bare catch is intentional: some Claim subclass property getters
            // throw on partially-constructed claims (e.g. when a typed parent
            // is being walked during inheritance graph construction). The
            // walker treats inaccessible parents as not-yet-reachable rather
            // than failing the whole traversal.
            try { value = prop.GetValue(c) as Claim; }
            catch { continue; }
            if (value is not null) yield return value;
        }
    }
}
