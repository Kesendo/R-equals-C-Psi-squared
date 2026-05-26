using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.AnchorAudit;

/// <summary>Reflection-based wiring audit: walks every concrete <see cref="Claim"/> subclass
/// in RCPsiSquared.Core and RCPsiSquared.Diagnostics and asserts each is either registered in
/// <see cref="KnowledgeRegistryFactory.BuildDefault"/> or listed on the explicit
/// <see cref="DeferredClaims"/> allowlist.
///
/// <para>Catches the "new Claim subclass added, forgot to wire it into the registry" failure
/// mode, which the per-anchor <c>CoreAnchorAuditTests</c> cannot see: a Claim with a valid
/// anchor file still counts as unwired if no <c>Register*</c> call ever references it.</para>
///
/// <para>An entry on the allowlist names the task tracking its wiring. As those tasks land,
/// the claim moves off the allowlist; <see cref="DeferredClaims_AreUnregisteredAndReal"/>
/// fails if an allowlist entry is actually registered (stale) or matches no Claim type (typo
/// or renamed class), so the allowlist cannot rot.</para></summary>
public class RegistryWiringAuditTests
{
    /// <summary>Concrete Claim subclasses intentionally NOT yet wired into
    /// <see cref="KnowledgeRegistryFactory.BuildDefault"/>. Each entry is the Claim's type
    /// name; the grouping comments name the task that will wire it.</summary>
    private static readonly IReadOnlySet<string> DeferredClaims =
        new HashSet<string>(StringComparer.Ordinal)
        {
            // Concrete Claim subclasses not yet wired into
            // KnowledgeRegistryFactory.BuildDefault(), grouped by the task that wires each
            // cluster. As a task lands, delete its lines here; DeferredClaims_AreUnregistered-
            // AndReal fails on any entry that becomes registered, so this list cannot rot.

            // Group D structural (task #5 triage): heavy parameterised block-spectrum /
            // Jordan-Wigner / Prosen compute primitives. Each has a private ctor + a public
            // Build(N, sector, gamma, ...) factory taking data, not Claim parents; the two
            // KleinFourGroupSelfPaired* even require even N (incompatible with the registry's
            // N=5 default). Invoked per-analysis, not registry singletons. The one genuine
            // bridge claim in this cluster, F89F88aKleinPpAnchor, was wired into the factory.
            "JwBondClusterPairAffinity",
            "JwBondQPeakPrediction",
            "JwBondQPeakUnified",
            "JwDispersionDProjection",
            "JwSlaterPairArnoldiEig",
            "JwSlaterPairBasis",
            "JwSlaterPairF1PalindromeProbe",
            "JwSlaterPairLProjection",
            "JwSlaterPairShiftInvertArnoldi",
            "JwSlaterPairSparseLBuilder",
            "KleinFourGroupSelfPairedRefinement",
            "KleinFourGroupSelfPairedSparseLBuilder",
            "LiouvillianSectorSweep",
            "OneSidedSectorClosedForm",

            // Group B structural (task #6 triage): F86 Sammelbecken. Almost all are
            // per-CoherenceBlock parameterised diagnostics (private ctor + Build(block)
            // factory) or types instantiated several times as children of F86KnowledgeBase,
            // not registry singletons. C2DirectionAFalsificationProbe and C2DirectionC-
            // FalsificationProbe additionally run a ~5 s N=5..8 empirical scan in their
            // Build(), too heavy for the BuildDefault() hot path. HardwareConfirmationClaim
            // is a generic per-Confirmation adapter whose anchor is a runtime value, not a
            // singleton. The 4 genuine standalone singleton claims in this cluster
            // (C2BareDoubledPtfClosedForm, F86HwhmClosedFormClaim, IbmBlockCpsiHardwareTable,
            // PolarityPairQPeakDecompositionClaim) were wired into the factory.
            "BlockCpsiTrajectory",
            "C2BlockCpsiQScan",
            "C2BlockJwDecomposition",
            "C2BondCoupling",
            "C2BondKModeProfile",
            "C2BondLQPeakScan",
            "C2BondLQProjection",
            "C2BondLensComparison",
            "C2DirectionAFalsificationProbe",
            "C2DirectionCFalsificationProbe",
            "C2EffectiveSpectrum",
            "C2FullBlockEigenAnatomy",
            "C2FullBlockPairAnatomy",
            "C2FullBlockSigmaAnatomy",
            "C2HwhmRatio",
            "C2InterChannelAnalytical",
            "C2InterChannelProjector",
            "C2KShape",
            "C2SvdBlockProjectedMagnitude",
            "C2UniversalShapeDerivation",
            "HardwareConfirmationClaim",
            "OrbitKTrend",
            "PerBondQPeakWitnessTable",
            "PerF71OrbitKModeTable",
            "PerF71OrbitKTable",
            "PerF71OrbitLQPeakTable",
            "ShapeFunctionWitnesses",
            "SigmaZeroCommutatorNormClaim",
            "TwoLevelEpModel",

            // Ungrouped (triaged by the owning task when reached): ChiralKClaim is a chiral
            // K-symmetry foundation; XGlobalEigenstateMirrorPi2Inheritance is a Pi2-inheritance
            // claim with no Register* call yet; OpenQuestion is a generic open-question holder
            // (confirm whether it is meant to be a registered singleton at all).
            "ChiralKClaim",
            "OpenQuestion",
            "XGlobalEigenstateMirrorPi2Inheritance",

            // Structural (stay deferred by design, not a wiring gap): F87CanonicalWitness is a
            // parameterised witness type instantiated as the children of the registered
            // F87StandardWitnessSet, never a standalone registry singleton.
            "F87CanonicalWitness",

            // Structural (stay deferred by design, not a wiring gap):
            // LindbladBitBPiBalanceWitness is a parameterised witness type (chain N, H terms,
            // optional γ_T1, expected verdict). The StandardSet factory requires N=2 (witness
            // 5 is a hard-coded 2-site Z-drive); the registry's default N=5 Heisenberg chain
            // is incompatible. Discoverable via the static StandardSet(chain) factory and the
            // 10 LindbladBitBPiBalanceWitnessTests in RCPsiSquared.Diagnostics.Tests/Polarity/,
            // never a standalone registry singleton.
            "LindbladBitBPiBalanceWitness",
        };

    [Fact]
    public void ClaimTypeNames_AreUniqueAcrossCoreAndDiagnostics()
    {
        var duplicates = ClaimReflection.AllConcreteClaimTypes()
            .GroupBy(t => t.Name)
            .Where(g => g.Count() > 1)
            .Select(g => g.Key)
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();

        Assert.True(duplicates.Count == 0,
            "Claim type names must be unique across Core + Diagnostics for the wiring audit "
            + $"to key on the simple name; found {duplicates.Count} collision(s):\n  - "
            + string.Join("\n  - ", duplicates));
    }

    [Fact]
    public void EveryClaim_IsRegistered_OrExplicitlyDeferred()
    {
        var registered = RegisteredClaimNames();

        var unwired = ClaimReflection.AllConcreteClaimTypes()
            .Select(t => t.Name)
            .Where(n => !registered.Contains(n) && !DeferredClaims.Contains(n))
            .Distinct()
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();

        Assert.True(unwired.Count == 0,
            $"{unwired.Count} Claim(s) are neither registered in "
            + "KnowledgeRegistryFactory.BuildDefault() nor on the DeferredClaims allowlist. "
            + "Wire each into the factory, or add it to DeferredClaims with the task that "
            + "tracks it:\n  - " + string.Join("\n  - ", unwired));
    }

    [Fact]
    public void DeferredClaims_AreUnregisteredAndReal()
    {
        var registered = RegisteredClaimNames();
        var allClaimNames = ClaimReflection.AllConcreteClaimTypes()
            .Select(t => t.Name)
            .ToHashSet(StringComparer.Ordinal);

        var stale = DeferredClaims
            .Where(registered.Contains)
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();
        Assert.True(stale.Count == 0,
            $"{stale.Count} DeferredClaims entry(ies) are now registered; remove them from "
            + "the allowlist:\n  - " + string.Join("\n  - ", stale));

        var phantom = DeferredClaims
            .Where(n => !allClaimNames.Contains(n))
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();
        Assert.True(phantom.Count == 0,
            $"{phantom.Count} DeferredClaims entry(ies) match no Claim type (typo or renamed "
            + "class):\n  - " + string.Join("\n  - ", phantom));
    }

    private static HashSet<string> RegisteredClaimNames() =>
        KnowledgeRegistryFactory.BuildDefault()
            .All()
            .Select(c => c.GetType().Name)
            .ToHashSet(StringComparer.Ordinal);
}
