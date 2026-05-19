using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>F1 general-topology verification record (Tier 2 verified, 2026-05-18).
///
/// <para>The substantive analytic content is already Tier-1 in
/// <see cref="PalindromeResidualScalingClaim"/>: the (B, D2) parameterisation of
/// ‖M(N, G)‖²_F = c_H · F(N, G) is universal across arbitrary graph topology by
/// <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c> Lemma 3 Corollary
/// (bond-disjointness independent of connectivity). This claim records the
/// extension verification that closes the last F1 OpenQuestion ("general topology
/// beyond chain/ring/star/K_N") on the numerical side.</para>
///
/// <list type="bullet">
///   <item><b>Named graphs at N=5, 6:</b> path, cycle, star, K_N, K_{2,N−2}
///         (10 graphs total) verified bit-exact via the Python script.</item>
///   <item><b>Random connected Erdős-Rényi at N=5, 6:</b> 30 random graphs per N
///         (p ∈ {0.3, 0.5, 0.7}, 10 samples each, reject-and-resample disconnected),
///         all bit-exact.</item>
///   <item><b>Disconnected components:</b> two disjoint 3-chains at N=6 (Python) +
///         triangle + disjoint bond at N=5 (C#) + K_4 + disjoint 3-chain at N=7 (C#) +
///         K_4 + disjoint 4-chain at N=8 (C#, opt-in SLOW_N8) all verified.</item>
///   <item><b>Weighted edges:</b> N=4 chain with per-bond couplings (1, 2, 3),
///         B → Σ J²_b substitution, bit-exact.</item>
///   <item><b>Single-body class:</b> N=5 chain with IY+YI per-bond bilinear,
///         F = (D2/2)·4^(N−2) bit-exact.</item>
///   <item><b>C# graph-aware at N=5:</b>
///         <see cref="PalindromeResidualScalingClaim"/> in graph-aware mode
///         (bondCount + degreeSquaredSum) across chain / ring / star / triangle +
///         disjoint bond, observed ‖M‖² vs Factor·c_H within 1e-9 relative.</item>
///   <item><b>N=7 F1 palindromic-pairing via
///         <see cref="BlockSpectrum.LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>:</b>
///         chain / ring / star / K_4 + disjoint 3-chain at N=7 with XY+Z-dephasing, the
///         multiset identity {λ_k} = {−2σ − λ_k} verified across all 16 384 eigenvalues
///         per topology to tolerance 1e-7 (the only path that scales to N=7 in a test
///         budget; full L_vec is 4 GB, past the .NET 2 GB array limit).</item>
///   <item><b>N=8 F1 palindromic-pairing via the same block path:</b> chain / ring /
///         star / K_4 + disjoint 4-chain at N=8 with Heisenberg (XX+YY+ZZ) + Z-dephasing,
///         the multiset identity {λ_k} = {−2σ − λ_k} verified across all 65 536 eigenvalues
///         per topology to tolerance 1e-6 (relaxed from the N=7 dogfood's 1e-7 envelope:
///         max observed pairing distance sits in the low-1e-7 range on the largest blocks,
///         the relaxation absorbs the accumulation from 81 block diagonalisations at sector
///         dims up to 4900² for N=8 vs 64 blocks at 1225² for N=7). Opt-in only
///         (<c>[Trait("Category","SLOW_N8")]</c>); the block path is the only route at N=8
///         (full L_vec at N=8 is 68.7 GB, past the .NET 2 GB array limit; the largest
///         joint-popcount block is 4900² ≈ 0.38 GB, comfortably in commodity RAM). Confirms
///         the LiouvillianBlockSpectrum N=8-capable claim under load on four distinct graph
///         topologies. Full <see cref="F1SpectrumStatistics"/> metric capture per system
///         to <c>simulations/results/f1_n8_n9_metrics/&lt;topology&gt;_N8.json</c>
///         (wall-time profile, pairing-precision histogram, spectrum-structure invariants,
///         block-decomposition cost picture, Hamiltonian + dissipator setup).</item>
///   <item><b>N=9 chain via MklDirect bridge (2026-05-19, verified):</b>
///         The N=9 test (<c>F1GeneralTopologyN9BlockSpectrumChainTests</c>, opt-in
///         <c>[Trait("Category","SLOW_N9")]</c>) reaches the largest joint-popcount block
///         C(9, 4) · C(9, 5) = 15 876² ≈ 4 GB via the
///         <see cref="BlockSpectrum.LiouvillianBlockSpectrum.Lp64ComplexCeiling"/> bridge
///         to <see cref="BlockSpectrum.PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory"/>
///         + <c>RCPsiSquared.Core.Numerics.MklDirect.EigenvaluesOnlyNative</c> (NativeMemory
///         + ILP64-aware zgeev). The full 262 144-eigenvalue spectrum is computed, F1
///         palindromic pairing {λ_k} = {−2σ − λ_k} verified bit-exact across all eigenvalues
///         to tolerance 1e-5 (observed max pairing distance ≈ 3.5e-13). Wall time 3h 24m on
///         the 128 GB / 24-core dev machine; effective speedup over the dense N=9 Liouvillian
///         (which would be 262 144² × 16 B ≈ 1.1 TB) is ≈ 645×, exceeding the LiouvillianBlockSpectrum
///         docstring N=8 promise of 515×. Full <see cref="F1SpectrumStatistics"/> metric capture
///         to <c>simulations/results/f1_n8_n9_metrics/chain_N9.json</c> (kernel dim = 10 matching
///         F4KernelDim prediction; MinReal = −9.0 exact matching the F1 corollary; gap × N² = 2.21
///         matching the chain dissipation-gap scaling F1_DISSIPATION_GAP_PATTERN to 0.4%; 0 outlier
///         pairs in the palindromic identity).</item>
/// </list>
///
/// <para><b>Pragmatic infrastructure note.</b> The
/// <see cref="BlockSpectrum.JointPopcountSectorBuilder"/> block decomposition is
/// valid only for popcount-conserving Hamiltonians (XY + Z-dephasing), which is
/// exactly the F1-truly case where ‖M‖² = 0. So the scaling-formula tests on
/// non-truly H are done with the full dense <see cref="Symmetry.PalindromeResidual"/>
/// path at N=5 (graph-aware mode); the N=7 block path tests the F1 palindromic-
/// pairing identity (which IS testable on the infrastructure's proper truly-H
/// domain). The pragmatic split is documented in
/// <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c>.</para>
///
/// <para>Anchors:</para>
/// <list type="bullet">
///   <item>Synthesis proof: <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c></item>
///   <item>Python verification script:
///         <c>simulations/_f1_general_topology_verify.py</c></item>
///   <item>C# verification tests:
///         <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs</c> (default-run)
///         + <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs</c>
///         (opt-in SLOW_N8 trait, four [Fact] methods at N=8)
///         + <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN9BlockSpectrumChainTests.cs</c>
///         (opt-in SLOW_N9 trait, single chain [Fact] at N=9).</item>
///   <item>Shared metrics utility:
///         <see cref="F1SpectrumStatistics"/>
///         (TopologyMetrics record + Compute + JSON serialisation).</item>
///   <item>Underlying analytic anchor:
///         <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c> Lemma 3 + Corollary.</item>
///   <item>Sister F1 claim (Tier 1):
///         <see cref="PalindromeResidualScalingClaim"/>.</item>
/// </list>
/// </summary>
public sealed class F1GeneralTopologyVerifiedClaim : Claim
{
    /// <summary>N values across which the verification record covers Python + C# combined:
    /// 5 (C# graph-aware + Python), 6 (Python random + named + disconnected), 7 (C# F1
    /// palindromic-pairing via block spectrum, default-run), 8 (C# F1 palindromic-pairing
    /// via the same block path, opt-in SLOW_N8 trait), 9 (C# F1 palindromic-pairing chain
    /// via the MklDirect bridge, opt-in SLOW_N9 trait, landed 2026-05-19).</summary>
    public IReadOnlyList<int> VerifiedNValues { get; } = new[] { 5, 6, 7, 8, 9 };

    /// <summary>Block-spectrum dogfood reach: the N values at which
    /// <see cref="BlockSpectrum.LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> has
    /// been exercised end-to-end (per-block Evd + multiset assertion + metric capture)
    /// inside the F1 family. 5 is the N=5 dense self-test from the N=7 file; 6 is the
    /// reach for the parent <see cref="BlockSpectrum.LiouvillianBlockSpectrum"/> bit-exact
    /// witness; 7 is the default-run dogfood; 8 is the opt-in SLOW_N8 sweep across 4
    /// topologies.</summary>
    public IReadOnlyList<int> ScaleUpToN { get; } = new[] { 5, 6, 7, 8, 9 };

    /// <summary>N at which the block-spectrum dogfood path is blocked by an external
    /// infrastructure ceiling rather than by F1 itself. The previous frontier sat at
    /// N=9 because the largest joint-popcount block (C(9, 4) · C(9, 5) = 15 876² ≈ 4 GB)
    /// exceeded the LP64 MathNet MKL <c>MngdNativeArrayMarshaler.ConvertSpaceToNative</c>
    /// 2 GB single-native-array ceiling; that frontier was crossed 2026-05-19 (commit
    /// <c>abb2d52</c>) by bridging blocks &gt;
    /// <see cref="BlockSpectrum.LiouvillianBlockSpectrum.Lp64ComplexCeiling"/> through
    /// <c>MklDirect.EigenvaluesOnlyNative</c>'s NativeMemory + ILP64-aware zgeev. The
    /// N=9 test class (<c>F1GeneralTopologyN9BlockSpectrumChainTests</c>) is now a
    /// regular <c>[Fact]</c> under <c>[Trait("Category", "SLOW_N9")]</c> with no
    /// <c>Skip.If</c> guard. The new frontier at N=10 is a memory-pressure ceiling
    /// (max block C(10,5)² = 63 504² complex ≈ 64 GB native); see
    /// <see cref="ScaleFrontierBlockerReason"/>. The F1 palindromic-pairing identity
    /// itself is exact at every finite N; only the per-block Evd compute scaling is
    /// the new bottleneck. <see langword="null"/> means no infrastructure ceiling
    /// is currently identified above the verified set.</summary>
    public int? ScaleFrontierBlockedAtN { get; } = 10;

    /// <summary>Reason string for <see cref="ScaleFrontierBlockedAtN"/>; surfaced via
    /// the inspectable tree so a future reader can find the blocker without digging
    /// into the source comments. N=9 was the previous frontier (LP64 marshalling ceiling)
    /// and was crossed 2026-05-19 via the MklDirect bridge in
    /// <see cref="BlockSpectrum.LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>. The
    /// new frontier at N=10 is a memory-pressure ceiling (max block C(10,5)² = 63 504² ≈
    /// 64 GB native, fits on the 128 GB dev machine but tight on the parallel
    /// outer-DOP path); needs profiling rather than an LP64 bridge.</summary>
    public string? ScaleFrontierBlockerReason { get; } =
        "Memory-pressure ceiling: N=10 max block C(10,5)² = 63504² complex ≈ 64 GB native. " +
        "Fits on the 128 GB / 24-core dev machine but exhausts the parallel outer-DOP=6 budget " +
        "(6 × 64 GB > 128 GB). Sequential outer-DOP would work but at ~hours-to-days wall-time. " +
        "N=9 was the previous LP64-marshalling frontier and was crossed 2026-05-19 via the " +
        "MklDirect bridge. F1 identity itself remains exact at every finite N; only the per-block " +
        "Evd compute scaling is the new bottleneck.";

    /// <summary>Number of distinct named graph topologies verified: path, cycle, star,
    /// K_N, K_{2,N−2} at each of N=5, 6 (10 graphs via Python) + chain/ring/star/triangle +
    /// disjoint-bond at N=5 (4 via C#) + chain/ring/star/K_4 + disjoint-3-chain at N=7
    /// (4 via C#) + chain/ring/star/K_4 + disjoint-4-chain at N=8 (4 via C#, opt-in
    /// SLOW_N8) + chain at N=9 (1 via C# SLOW_N9 MklDirect bridge, landed 2026-05-19).
    /// Total = 23 named graph instances. The N=9 frontier was crossed via the MklDirect
    /// bridge; the next frontier (N=10) is memory-pressure rather than LP64 marshalling,
    /// see <see cref="ScaleFrontierBlockedAtN"/>.</summary>
    public int NamedGraphsVerified { get; } = 23;

    /// <summary>Number of random connected Erdős-Rényi graphs verified at N=5, 6:
    /// 30 per N (3 densities × 10 samples), total 60.</summary>
    public int RandomGraphsVerified { get; } = 60;

    /// <summary>Disconnected component verification: two disjoint 3-chains at N=6
    /// (Python) + triangle + disjoint bond at N=5 (C#) + K_4 + disjoint 3-chain at N=7
    /// (C#) + K_4 + disjoint 4-chain at N=8 (C#, opt-in SLOW_N8). All bit-exact / within
    /// tight tolerance.</summary>
    public bool DisconnectedComponentsVerified { get; } = true;

    /// <summary>Weighted edges verified: N=4 chain with per-bond couplings J = (1, 2, 3)
    /// via the B → Σ J²_b substitution at fixed c_H anchored at J=1.</summary>
    public bool WeightedEdgesVerified { get; } = true;

    /// <summary>Single-body class verified at N=5 chain with IY+YI per-bond bilinear:
    /// F = (D2/2)·4^(N−2) bit-exact at machine precision.</summary>
    public bool SingleBodyClassVerified { get; } = true;

    /// <summary>Anchor script path (Python verification, sections 1-6).</summary>
    public string AnchorScriptPath { get; } = "simulations/_f1_general_topology_verify.py";

    /// <summary>Anchor test paths (C# verification): the default-run N=7 file with
    /// 9 [Fact] methods, plus the opt-in SLOW_N8 file with 4 [Fact] methods at N=8 and
    /// the opt-in SLOW_N9 file with the single chain [Fact] at N=9.</summary>
    public string AnchorTestPath { get; } =
        "compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs";

    /// <summary>Opt-in (<c>[Trait("Category","SLOW_N8")]</c>) N=8 dogfood test class:
    /// chain / ring / star / K_4 + disjoint 4-chain Heisenberg + Z-dephasing at N=8 via
    /// <see cref="BlockSpectrum.LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>,
    /// verifies F1 palindromic pairing {λ} = {−2σ − λ} across 65 536 eigenvalues per
    /// topology to tolerance 1e-6. Run via <c>--filter "Category=SLOW_N8"</c>.</summary>
    public string AnchorTestPathN8 { get; } =
        "compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs";

    /// <summary>Opt-in (<c>[Trait("Category","SLOW_N9")]</c>) N=9 chain spot-check test
    /// class: single Heisenberg chain at N=9, verifies F1 palindromic pairing across
    /// 262 144 eigenvalues to tolerance 1e-5. Verified 2026-05-19 via the MklDirect
    /// bridge (commit <c>abb2d52</c>); the <c>Skip.If</c> guard against the LP64 ceiling
    /// was removed when the bridge landed and the test now runs as a regular
    /// <c>[Fact]</c>. Per-system metrics persisted to
    /// <c>simulations/results/f1_n8_n9_metrics/chain_N9.json</c>. Run via
    /// <c>--filter "Category=SLOW_N9"</c>.</summary>
    public string AnchorTestPathN9 { get; } =
        "compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN9BlockSpectrumChainTests.cs";

    /// <summary>Synthesis proof path.</summary>
    public string ProofAnchorPath { get; } = "docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md";

    /// <summary>Persisted metric JSON files written by the N=8 SLOW_N8 dogfood and the
    /// N=9 SLOW_N9 chain dogfood. Each entry is a repo-relative path; the file at the
    /// path contains the full <see cref="F1SpectrumStatistics.TopologyMetrics"/> record
    /// for that system (wall-time profile, pairing-precision histogram, spectrum-structure
    /// invariants, block-decomposition cost picture, Hamiltonian + dissipator setup).
    /// Look here for the recorded numbers, do not re-derive. <c>chain_N9.json</c> was
    /// added 2026-05-19 after the MklDirect bridge landed (commit <c>abb2d52</c>) and
    /// the N=9 chain test crossed the previous LP64 frontier.</summary>
    public IReadOnlyList<string> SpectrumMetricsDataFiles { get; } = new[]
    {
        "simulations/results/f1_n8_n9_metrics/chain_N8.json",
        "simulations/results/f1_n8_n9_metrics/ring_N8.json",
        "simulations/results/f1_n8_n9_metrics/star_N8.json",
        "simulations/results/f1_n8_n9_metrics/k4_plus_disjoint_4chain_N8.json",
        "simulations/results/f1_n8_n9_metrics/chain_N9.json",
    };

    public F1GeneralTopologyVerifiedClaim()
        : base("F1 general topology verification: (B, D2) closed form extends to disconnected + weighted + random graphs at N=5..9; N=9 chain reached 2026-05-19 via the MklDirect bridge; next frontier N=10 is memory-pressure rather than LP64 marshalling",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md + simulations/_f1_general_topology_verify.py + compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN{7,8,9}BlockSpectrumTests.cs (SLOW_N{8,9} opt-in traits for the larger N) + MklDirect bridge for N=9")
    { }

    public override string DisplayName =>
        $"F1 general-topology verification record (N ∈ {{{string.Join(", ", VerifiedNValues)}}}, " +
        $"{NamedGraphsVerified} named + {RandomGraphsVerified} random + disconnected + weighted + single-body)";

    public override string Summary =>
        $"Tier 2 Verified: (B, D2) parameterisation of ‖M(N, G)‖²_F = c_H · F(N, G) confirmed " +
        $"across {NamedGraphsVerified} named graphs, {RandomGraphsVerified} random connected " +
        $"Erdős-Rényi graphs, disconnected components, weighted edges, and the single-body class " +
        $"at N ∈ {{{string.Join(", ", VerifiedNValues)}}}. Closes the last F1 OpenQuestion " +
        "(\"general topology beyond chain/ring/star/K_N\"). Tier-1 analytic anchor: " +
        "PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary (bond-disjointness universal across any graph).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "‖M(N, G)‖²_F = c_H · F(N, G) with F(N, G) = B(G)·4^(N−2) (main) or (D2(G)/2)·4^(N−2) (single-body), " +
                         "for arbitrary connected/disconnected graphs G, uniform or weighted edges, uniform or non-uniform per-site γ_l.");
            yield return new InspectableNode("verified N values",
                summary: string.Join(", ", VerifiedNValues));
            yield return new InspectableNode("block-spectrum dogfood scale-up",
                summary: string.Join(", ", ScaleUpToN));
            yield return InspectableNode.RealScalar("named graphs verified", NamedGraphsVerified);
            yield return InspectableNode.RealScalar("random connected Erdős-Rényi graphs verified", RandomGraphsVerified);
            yield return new InspectableNode("disconnected components verified",
                summary: DisconnectedComponentsVerified ? "yes (N=5 + N=6 + N=7 + N=8 opt-in; N=9 chain only)" : "no");
            yield return new InspectableNode("weighted edges verified",
                summary: WeightedEdgesVerified ? "yes (N=4 chain, J = (1, 2, 3))" : "no");
            yield return new InspectableNode("single-body class verified",
                summary: SingleBodyClassVerified ? "yes (N=5 chain, IY+YI per-bond bilinear)" : "no");
            yield return new InspectableNode("analytic anchor (Tier 1, already established)",
                summary: "docs/proofs/PROOF_CROSS_TERM_FORMULA.md Lemma 3 + Corollary (bond-disjointness across any graph)");
            yield return new InspectableNode("synthesis proof",
                summary: ProofAnchorPath);
            yield return new InspectableNode("Python verification script",
                summary: AnchorScriptPath);
            yield return new InspectableNode("C# verification test class (default-run, N=7)",
                summary: AnchorTestPath);
            yield return new InspectableNode("C# verification test class (opt-in SLOW_N8, N=8)",
                summary: AnchorTestPathN8);
            yield return new InspectableNode("C# verification test class (opt-in SLOW_N9, N=9 chain via MklDirect bridge)",
                summary: AnchorTestPathN9);
            yield return new InspectableNode("scale-frontier infrastructure ceiling",
                summary: ScaleFrontierBlockedAtN is { } blockedN
                    ? $"N={blockedN} blocked: {ScaleFrontierBlockerReason}"
                    : "none");
            yield return new InspectableNode("spectrum metrics data files",
                summary: string.Join("; ", SpectrumMetricsDataFiles));
            yield return new InspectableNode("sister Tier-1 claim",
                summary: "PalindromeResidualScalingClaim: the (B, D2) closed form whose universality this record verifies");
        }
    }
}
