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
///         triangle + disjoint bond at N=5 (C#) + K_4 + disjoint 3-chain at N=7 (C#)
///         all verified.</item>
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
///         <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs</c></item>
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
    /// palindromic-pairing via block spectrum).</summary>
    public IReadOnlyList<int> VerifiedNValues { get; } = new[] { 5, 6, 7 };

    /// <summary>Number of distinct named graph topologies verified: path, cycle, star,
    /// K_N, K_{2,N−2} at each of N=5, 6 (10 graphs via Python) + chain/ring/star/triangle +
    /// disjoint-bond at N=5 (4 via C#) + chain/ring/star/K_4 + disjoint-3-chain at N=7
    /// (4 via C#). Total = 18 named graph instances.</summary>
    public int NamedGraphsVerified { get; } = 18;

    /// <summary>Number of random connected Erdős-Rényi graphs verified at N=5, 6:
    /// 30 per N (3 densities × 10 samples), total 60.</summary>
    public int RandomGraphsVerified { get; } = 60;

    /// <summary>Disconnected component verification: two disjoint 3-chains at N=6
    /// (Python) + triangle + disjoint bond at N=5 (C#) + K_4 + disjoint 3-chain at N=7
    /// (C#). All bit-exact / within tight tolerance.</summary>
    public bool DisconnectedComponentsVerified { get; } = true;

    /// <summary>Weighted edges verified: N=4 chain with per-bond couplings J = (1, 2, 3)
    /// via the B → Σ J²_b substitution at fixed c_H anchored at J=1.</summary>
    public bool WeightedEdgesVerified { get; } = true;

    /// <summary>Single-body class verified at N=5 chain with IY+YI per-bond bilinear:
    /// F = (D2/2)·4^(N−2) bit-exact at machine precision.</summary>
    public bool SingleBodyClassVerified { get; } = true;

    /// <summary>Anchor script path (Python verification, sections 1-6).</summary>
    public string AnchorScriptPath { get; } = "simulations/_f1_general_topology_verify.py";

    /// <summary>Anchor test path (C# verification, 9 [Fact] methods).</summary>
    public string AnchorTestPath { get; } =
        "compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs";

    /// <summary>Synthesis proof path.</summary>
    public string ProofAnchorPath { get; } = "docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md";

    public F1GeneralTopologyVerifiedClaim()
        : base("F1 general topology verification: (B, D2) closed form extends to disconnected + weighted + random graphs at N=5..7",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md + simulations/_f1_general_topology_verify.py + compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs")
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
            yield return InspectableNode.RealScalar("named graphs verified", NamedGraphsVerified);
            yield return InspectableNode.RealScalar("random connected Erdős-Rényi graphs verified", RandomGraphsVerified);
            yield return new InspectableNode("disconnected components verified",
                summary: DisconnectedComponentsVerified ? "yes (N=5 + N=6 + N=7)" : "no");
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
            yield return new InspectableNode("C# verification test class",
                summary: AnchorTestPath);
            yield return new InspectableNode("sister Tier-1 claim",
                summary: "PalindromeResidualScalingClaim — the (B, D2) closed form whose universality this record verifies");
        }
    }
}
