using System.IO;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Locate the repository root by walking up from the current process directory until
/// a directory containing <c>CLAUDE.md</c> is found. Returns <c>null</c> if no such ancestor
/// exists (e.g. the test process is running in an isolated temp directory).</summary>
file static class RepoRootLocator
{
    // Lazy<T> guarantees thread-safe single-execution and observably-published result.
    // Replaces a hand-rolled (_cached, _searched) pair that had a race window where
    // _searched was set true before _cached was assigned, causing parallel xunit runs
    // to observe null and skip repo-root resolution.
    private static readonly Lazy<string?> _root = new(SearchUpForRepoRoot, isThreadSafe: true);

    public static string? Find() => _root.Value;

    private static string? SearchUpForRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (File.Exists(Path.Combine(dir.FullName, "CLAUDE.md")))
                return dir.FullName;
            dir = dir.Parent;
        }
        return null;
    }
}

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
        if (_factories.ContainsKey(typeof(T)))
            throw new InvariantViolationException(
                rule: "DuplicateRegistration",
                message: $"Type {typeof(T).Name} registered twice. The builder is single-pass; remove one of the Register calls.",
                hint: $"Search the registration code for '.Register<{typeof(T).Name}>'.",
                offendingClaim: typeof(T));
        _factories[typeof(T)] = ctx => factory(ctx);
        return this;
    }

    /// <summary>Resolves all registered factories in topological order by deferred construction.
    /// Throws <see cref="InvariantViolationException"/> (Rule = "MissingParent") when a factory
    /// calls <c>b.Get&lt;T&gt;()</c> for an unregistered type, or (Rule = "Cycle") when a
    /// circular dependency prevents forward progress.</summary>
    public ClaimRegistry Build()
    {
        if (_built)
            throw new InvalidOperationException("ClaimRegistry already built; use the returned registry.");
        _built = true;

        var resolved = new Dictionary<Type, Claim>();
        var edges = new List<Edge>();
        var order = new List<Type>();
        var pending = new HashSet<Type>(_factories.Keys);
        var registered = new HashSet<Type>(_factories.Keys);

        while (pending.Count > 0)
        {
            int progress = 0;
            foreach (var type in pending.ToList())
            {
                Type currentlyResolving = type;
                // Tentative edge buffer: edges accumulate per attempt and commit only on
                // success. Without this, a deferred factory that resolved some parents
                // before throwing on a later one would re-add those edges every retry pass.
                var attemptEdges = new List<Edge>();
                var ctx = new RecordingBuilderContext(resolved, registered, currentlyResolving, attemptEdges);
                try
                {
                    var claim = _factories[type](ctx);
                    resolved[type] = claim;
                    order.Add(type);
                    pending.Remove(type);
                    edges.AddRange(attemptEdges);
                    progress++;
                }
                catch (DependencyNotYetAvailableException)
                {
                    // defer; tentative edges discarded
                }
            }

            if (progress == 0)
                throw new InvariantViolationException(
                    rule: "Cycle",
                    message: "Dependency cycle detected among: " + string.Join(", ", pending.Select(t => t.Name)),
                    hint: "Restructure: split a Claim, or invert one edge.",
                    path: pending.ToArray());
        }

        // Tier inheritance: parent must be at least as strong as child.
        foreach (var edge in edges)
        {
            var parent = resolved[edge.Parent];
            var child = resolved[edge.Child];
            if (!TierStrength.IsAtLeastAsStrong(parent.Tier, child.Tier))
                throw new InvariantViolationException(
                    rule: "TierInheritance",
                    message: $"Tier inheritance violation: {edge.Child.Name} ({child.Tier}) ← {edge.Parent.Name} ({parent.Tier}). " +
                             $"parent.Tier must be at least as strong as child.Tier.",
                    hint: $"Either downgrade {edge.Child.Name} or strengthen the foundation under {edge.Parent.Name}.",
                    offendingClaim: edge.Child,
                    path: new[] { edge.Parent, edge.Child });
        }

        // Anchor-file existence (best-effort heuristic: any anchor token containing ".md").
        // Anchor strings often pack multiple references separated by " + " or " / "; we split on
        // those, take tokens with .md, and check each individual file path.
        foreach (var (type, claim) in resolved)
        {
            if (string.IsNullOrWhiteSpace(claim.Anchor)) continue;

            var tokens = claim.Anchor.Split(new[] { " + ", " / ", ", " }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var raw in tokens)
            {
                var token = raw.Trim();
                if (!token.Contains(".md", StringComparison.Ordinal)) continue;

                // Strip trailing markers before checking: ' ' / '#' / '(' for section
                // labels, ':' for line-number suffixes like 'docs/X.md:251' used in some
                // Pi2 claim anchors.
                var cleanedTokenEnd = token.IndexOfAny(new[] { ' ', '#', '(', ':' });
                var path = cleanedTokenEnd > 0 ? token[..cleanedTokenEnd] : token;
                if (string.IsNullOrEmpty(path)) continue;

                // Resolve from repo root when available; fall back to process CWD.
                var repoRoot = RepoRootLocator.Find();
                var resolved2 = repoRoot != null
                    ? Path.Combine(repoRoot, path)
                    : path;
                if (!File.Exists(resolved2))
                    throw new InvariantViolationException(
                        rule: "AnchorFileMissing",
                        message: $"Claim {type.Name} anchors at '{path}', file not found relative to working directory.",
                        hint: $"Either write the anchor file at '{path}' or correct the Anchor on {type.Name}.",
                        offendingClaim: type);
            }
        }

        return new ClaimRegistry(resolved, edges, order);
    }

    /// <summary>The recording context: returns resolved instances and records an edge,
    /// throws <see cref="InvariantViolationException"/> (Rule = "MissingParent") if the
    /// requested type was never registered, or throws <see cref="DependencyNotYetAvailableException"/>
    /// if the type is registered but not yet resolved (signals the builder to defer).</summary>
    private sealed class RecordingBuilderContext : IBuilderContext
    {
        private readonly Dictionary<Type, Claim> _resolved;
        private readonly HashSet<Type> _registered;
        private readonly Type _resolvingType;
        private readonly List<Edge> _edges;

        public RecordingBuilderContext(
            Dictionary<Type, Claim> resolved,
            HashSet<Type> registered,
            Type resolvingType,
            List<Edge> edges)
        {
            _resolved = resolved;
            _registered = registered;
            _resolvingType = resolvingType;
            _edges = edges;
        }

        public T Get<T>() where T : Claim
        {
            if (!_registered.Contains(typeof(T)))
                throw new InvariantViolationException(
                    rule: "MissingParent",
                    message: $"Missing dependency: {_resolvingType.Name} requested b.Get<{typeof(T).Name}>(), but {typeof(T).Name} is not registered.",
                    hint: $"Register {typeof(T).Name} via builder.Register<{typeof(T).Name}>(...).",
                    offendingClaim: _resolvingType);

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
