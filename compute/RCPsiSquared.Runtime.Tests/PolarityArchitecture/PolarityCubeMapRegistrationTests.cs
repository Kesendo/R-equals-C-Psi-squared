using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>End-to-end smoke + drift-guard for <see cref="PolarityCubeMap"/>'s registration.
///
/// <para>The 60+ b.Get&lt;T&gt;() dependencies in
/// <c>PolarityCubeMapRegistration.RegisterPolarityCubeMap</c> are otherwise unverified at
/// the registry level; the drift-guard catches the failure mode "new BitB Pi²-Inheritance
/// Claim added in Core but forgotten in the registration list".</para>
///
/// <para><b>Diagnostics-assembly exclusion</b>: the reflection scan filters to the Core
/// assembly only. F87Pi2Inheritance lives in RCPsiSquared.Diagnostics; Runtime does not
/// reference Diagnostics, so the Runtime-side PolarityCubeMap snapshot does not include
/// it, and the drift-guard count must agree with that scope. The parallel
/// architecture-enforcement test for Diagnostics IZ2AxisClaim coverage lives in
/// RCPsiSquared.Diagnostics.Tests/F87/F87Pi2InheritanceZ2AxisTests.cs.</para></summary>
public class PolarityCubeMapRegistrationTests
{
    private static readonly Assembly CoreAssembly = typeof(Claim).Assembly;

    /// <summary>Concrete Core-assembly IZ2AxisClaim types intentionally NOT yet wired into
    /// <see cref="KnowledgeRegistryFactory.BuildDefault"/>. Mirrors the deferred-claim
    /// subset on <c>RegistryWiringAuditTests.DeferredClaims</c> that also implements
    /// IZ2AxisClaim; the drift-guard test below allows these by name. When a deferred
    /// IZ2AxisClaim is wired into the factory, remove it from both lists.</summary>
    private static readonly IReadOnlySet<string> DeferredIZ2AxisClaims =
        new HashSet<string>(StringComparer.Ordinal)
        {
            // Pi2-Inheritance claim with no Register* call yet; tracked on
            // RegistryWiringAuditTests.DeferredClaims (ungrouped section).
            "XGlobalEigenstateMirrorPi2Inheritance",
        };

    [Fact]
    public void PolarityCubeMap_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        Assert.True(registry.Contains<PolarityCubeMap>());
        Assert.NotNull(registry.Get<PolarityCubeMap>());
    }

    [Fact]
    public void PolarityCubeMap_TotalClaims_MatchesIZ2AxisClaimInventory()
    {
        // Drift-guard: every concrete Core-assembly IZ2AxisClaim must appear in the
        // PolarityCubeMap registration unless explicitly deferred via
        // DeferredIZ2AxisClaims. If a future BitB Claim is added but the Register-list
        // in PolarityCubeMapRegistration.cs is not updated, the counts diverge and this
        // test fails with the missing type names.
        //
        // Scope = Core only: F87Pi2Inheritance lives in Diagnostics; Runtime does not
        // reference Diagnostics so the Runtime-side snapshot cannot see it. See the
        // class-level docstring for the parallel Diagnostics-side enforcement test.
        var coreIZ2AxisClaimTypes = CoreAssembly.GetTypes()
            .Where(t => t.IsClass && !t.IsAbstract)
            .Where(t => typeof(IZ2AxisClaim).IsAssignableFrom(t))
            .ToHashSet();

        Assert.NotEmpty(coreIZ2AxisClaimTypes);

        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        var registeredIZ2AxisTypes = registry.All()
            .Where(c => c is IZ2AxisClaim)
            .Select(c => c.GetType())
            .Where(t => t.Assembly == CoreAssembly)
            .ToHashSet();

        var missingFromRegistry = coreIZ2AxisClaimTypes
            .Except(registeredIZ2AxisTypes)
            .Select(t => t.Name)
            .Where(n => !DeferredIZ2AxisClaims.Contains(n))
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();

        Assert.True(missingFromRegistry.Count == 0,
            "Core-side IZ2AxisClaim types not registered in BuildDefault (drift-guard for "
            + "PolarityCubeMapRegistration.cs); wire them in or add to DeferredIZ2AxisClaims: "
            + string.Join(", ", missingFromRegistry));

        // Stale-allowlist check: every DeferredIZ2AxisClaims entry must match a real,
        // unregistered Core IZ2AxisClaim. Mirrors RegistryWiringAuditTests' rot-resistant
        // allowlist pattern.
        var registeredNames = registeredIZ2AxisTypes.Select(t => t.Name).ToHashSet(StringComparer.Ordinal);
        var allCoreNames = coreIZ2AxisClaimTypes.Select(t => t.Name).ToHashSet(StringComparer.Ordinal);

        var staleDeferred = DeferredIZ2AxisClaims
            .Where(registeredNames.Contains)
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();
        Assert.True(staleDeferred.Count == 0,
            "DeferredIZ2AxisClaims entries now registered; remove from allowlist: "
            + string.Join(", ", staleDeferred));

        var phantomDeferred = DeferredIZ2AxisClaims
            .Where(n => !allCoreNames.Contains(n))
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();
        Assert.True(phantomDeferred.Count == 0,
            "DeferredIZ2AxisClaims entries match no Core IZ2AxisClaim type (typo/rename): "
            + string.Join(", ", phantomDeferred));

        // The cube map's TotalClaims should at least cover every wired Core IZ2AxisClaim;
        // it may be larger because Diagnostics-side claims could also be registered.
        var expectedRegistered = coreIZ2AxisClaimTypes.Count - DeferredIZ2AxisClaims.Count;
        Assert.True(cubeMap.TotalClaims >= expectedRegistered,
            $"PolarityCubeMap.TotalClaims = {cubeMap.TotalClaims} is less than the "
            + $"wired Core IZ2AxisClaim count = {expectedRegistered}; the registration "
            + "list in PolarityCubeMapRegistration.cs is missing entries.");
    }

    [Fact]
    public void PolarityCubeMap_BitAClaims_ContainsExactlyF61()
    {
        // F61BitAParityPi2Inheritance is currently the only typed BitA Claim
        // (Π²_X = Z⊗N axis). If a second BitA Claim is added in the future, this
        // assertion intentionally fails so the new entry is reviewed for whether it
        // is a genuine bit_a sibling.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        Assert.Single(cubeMap.BitAClaims);
        Assert.Equal(nameof(F61BitAParityPi2Inheritance),
            cubeMap.BitAClaims[0].GetType().Name);
    }
}
