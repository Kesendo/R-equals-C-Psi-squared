using System.IO;
using System.Reflection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.Tests.AnchorAudit;

/// <summary>Reflection-based anchor audit: walks every <see cref="Claim"/> subclass in
/// <c>RCPsiSquared.Core</c> and <c>RCPsiSquared.Diagnostics</c>, reads its
/// <see cref="Claim.Anchor"/>, applies the same tokenisation heuristic the Runtime's
/// <c>ClaimRegistryBuilder</c> uses (split on " + ", " / ", ", "; take tokens with .md;
/// strip trailing markers via IndexOfAny(' ', '#', '(', ':')), and asserts each .md path
/// exists relative to the repo root. Equivalent to running the AnchorFileMissing invariant
/// against every typed Claim in both production assemblies at once, without registering
/// them individually.</summary>
public class CoreAnchorAuditTests
{
    [Fact]
    public void AllCoreClaims_HaveResolvableAnchorFiles() =>
        AuditAssembly(typeof(Claim).Assembly, "Core", minProbed: 10);

    [Fact]
    public void AllDiagnosticsClaims_HaveResolvableAnchorFiles()
    {
        // Load Diagnostics by referencing a known type from that assembly. We use
        // reflection to find one without taking a hard compile-time dependency on a
        // specific type, in case Diagnostics' surface evolves.
        var diagnosticsAssembly = AppDomain.CurrentDomain.GetAssemblies()
            .FirstOrDefault(a => a.GetName().Name == "RCPsiSquared.Diagnostics")
            ?? Assembly.Load("RCPsiSquared.Diagnostics");
        AuditAssembly(diagnosticsAssembly, "Diagnostics", minProbed: 1);
    }

    private static void AuditAssembly(Assembly assembly, string label, int minProbed)
    {
        var claimTypes = assembly.GetTypes()
            .Where(t => !t.IsAbstract && typeof(Claim).IsAssignableFrom(t))
            .ToList();

        var repoRoot = FindRepoRoot();
        Assert.NotNull(repoRoot);

        var failures = new List<string>();
        int probed = 0;

        foreach (var type in claimTypes)
        {
            // Skip claims that cannot be instantiated parameterlessly; their anchor is
            // determined by the constructor and may not be probeable without parameters.
            // The audit catches the bulk; parameter-dependent claims are exercised by
            // the per-family registration tests.
            var ctor = type.GetConstructor(Type.EmptyTypes);
            if (ctor is null) continue;

            Claim claim;
            try { claim = (Claim)ctor.Invoke(null); }
            catch { continue; }
            probed++;

            var anchor = claim.Anchor;
            if (string.IsNullOrWhiteSpace(anchor)) continue;

            var tokens = anchor.Split(new[] { " + ", " / ", ", " }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var raw in tokens)
            {
                var token = raw.Trim();
                if (!token.Contains(".md", StringComparison.Ordinal)) continue;

                var stripIdx = token.IndexOfAny(new[] { ' ', '#', '(', ':' });
                var path = stripIdx > 0 ? token[..stripIdx] : token;
                if (string.IsNullOrEmpty(path)) continue;

                var fullPath = Path.Combine(repoRoot, path);
                if (!File.Exists(fullPath))
                    failures.Add($"{type.Name} anchor token '{path}' not found at '{fullPath}'");
            }
        }

        // Sanity: ensure we actually exercised a non-trivial number of claims; otherwise
        // a vacuous pass would slip through if reflection broke.
        Assert.True(probed >= minProbed,
            $"{label} audit probed only {probed} parameterless Claims; expected >= {minProbed}");

        Assert.True(failures.Count == 0,
            $"Stale anchor files in {label} ({failures.Count}):\n  - " +
            string.Join("\n  - ", failures));
    }

    private static string? FindRepoRoot()
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
