using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.AnchorAudit;

/// <summary>Shared reflection over the typed-knowledge <see cref="Claim"/> population. Both
/// audit tests in this folder, <see cref="CoreAnchorAuditTests"/> and
/// <see cref="RegistryWiringAuditTests"/>, walk the same two production assemblies with the
/// same "concrete Claim subclass" predicate; this is the single source of both.</summary>
internal static class ClaimReflection
{
    /// <summary>The two production assemblies that declare <see cref="Claim"/> subclasses:
    /// RCPsiSquared.Core (located via <see cref="Claim"/>) and RCPsiSquared.Diagnostics
    /// (located via <see cref="KnowledgeRegistryFactory"/>).</summary>
    public static IReadOnlyList<Assembly> ClaimAssemblies { get; } = new[]
    {
        typeof(Claim).Assembly,
        typeof(KnowledgeRegistryFactory).Assembly,
    };

    /// <summary>Every concrete (non-abstract) <see cref="Claim"/> subclass declared in
    /// <paramref name="assembly"/>.</summary>
    public static IEnumerable<Type> ConcreteClaimTypesIn(Assembly assembly) =>
        assembly.GetTypes().Where(t => !t.IsAbstract && typeof(Claim).IsAssignableFrom(t));

    /// <summary>Every concrete <see cref="Claim"/> subclass across <see cref="ClaimAssemblies"/>.</summary>
    public static List<Type> AllConcreteClaimTypes() =>
        ClaimAssemblies.SelectMany(ConcreteClaimTypesIn).ToList();
}
