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
            // (empty after /simplify cleanup 2026-05-26: X-Mirror was wired with
            // its BitA twin Z-Mirror via XGlobalEigenstateMirrorPi2InheritanceRegistration.)
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
    public void PolarityCubeMap_CoveredByHadamardDualityTwinSlots_MatchesAbsorptionDescendants()
    {
        // The 9 Absorption-Theorem descendants (F33, F50, F55, F64, F65, F66, F67,
        // F68, F74) reclassify their bit_a twin slot to CoveredByHadamardDuality: the
        // bit_a image holds by the global Hadamard X↔Z duality
        // (docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md), so no bespoke typed twin is
        // owed. If a future BitB Claim adopts (or drops) this status, update the count.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        Assert.Equal(9, cubeMap.CoveredByHadamardDualityTwinSlots);
    }

    [Fact]
    public void PolarityCubeMap_BitAClaims_ContainsExpectedBitASiblings()
    {
        // Typed BitA Claims (Π²_X = Z⊗N axis): F61BitAParityPi2Inheritance was
        // the only entry until 2026-05-25 when F108 Part 2 (X-dephasing analog
        // of F108 Part 1, BitA twin) was added. Welle 7 (2026-05-26) added 4 more:
        // F38BitA, F39BitA, F63BitAReference, ZGlobalEigenstateMirrorBitA.
        // Welle 15 (2026-05-27) added LindbladBitAPiBalance (F112-X, BitA twin of
        // F112-Z per PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (f) item 4).
        // Total: 7. If an eighth BitA Claim lands, this assertion fails so the
        // new entry is reviewed for whether it is a genuine bit_a sibling.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        Assert.Equal(7, cubeMap.BitAClaims.Count);
        var typeNames = cubeMap.BitAClaims.Select(c => c.GetType().Name).ToHashSet();
        Assert.Contains(nameof(F61BitAParityPi2Inheritance), typeNames);
        Assert.Contains(nameof(F108Part2Pi2XEvenAlwaysPalindromic), typeNames);
        Assert.Contains(nameof(F38BitAInvolutionInheritance), typeNames);
        Assert.Contains(nameof(F39DetPiBitAInheritance), typeNames);
        Assert.Contains(nameof(F63BitAReference), typeNames);
        Assert.Contains(nameof(ZGlobalEigenstateMirrorBitAInheritance), typeNames);
        Assert.Contains(nameof(LindbladBitAPiBalance), typeNames);
    }

    [Fact]
    public void PolarityCubeMap_BitAClaims_HaveReciprocatingBitBPartner()
    {
        // Welle 8 reciprocity drift-guard (added 2026-05-26 after /simplify caught
        // the X-Mirror orphan condition): every BitA Claim in the cube map should
        // appear as the BitATwin of at least one registered BitB Claim. Catches
        // the failure mode "new BitA-twin Claim landed for an unregistered BitB
        // sibling, so the BitA Claim is counted but no twin slot is filled".
        //
        // Exception: F63BitAReference is by design NOT pointed to by F61 (the
        // BitA-axis Claim that proves [L, Π²_X] = 0); F63BitAReference is itself
        // pointed to by F63LCommutesPi2Pi2Inheritance (the BitB Claim) as the
        // structural-twin reference. F61BitAParityPi2Inheritance is similarly
        // pointed to by F1Pi2Inheritance via the BitATwinClaim property.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();

        var bitBClaimsWithTwins = cubeMap.BitBClaims
            .Where(c => c.BitATwin is not null)
            .Select(c => (object)c.BitATwin!)
            .ToHashSet(ReferenceEqualityComparer.Instance);

        foreach (var bitAClaim in cubeMap.BitAClaims)
        {
            Assert.True(
                bitBClaimsWithTwins.Contains(bitAClaim),
                $"BitA Claim {bitAClaim.GetType().Name} is registered but no BitB Claim " +
                $"points at it via BitATwin. Orphan condition: new BitA-twin landed " +
                $"without its reciprocating BitB partner being wired. Either register " +
                $"the BitB partner (cf. X-Mirror in Welle 8 cleanup) or remove the BitA " +
                $"Claim from the cube map inventory.");
        }
    }
}
