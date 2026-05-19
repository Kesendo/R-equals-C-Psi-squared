using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F4 disconnected-graph extension: kernel dimension of the Heisenberg
/// Liouvillian L_H factorises across connected components of the underlying graph.
///
/// <code>
///     dim ker L_H(G)  =  Π_c (|c| + 1)
/// </code>
///
/// where the product runs over connected components c of the graph G and |c| is
/// the number of qubits in component c. This is the bond-disjoint generalisation
/// of the connected-graph behaviour: for any single connected component of size n,
/// popcount conservation under XX+YY+ZZ gives n+1 distinct popcount labels
/// {0, 1, 2, ..., n}, each labelling exactly one kernel mode of L_H restricted to
/// the component. For a disconnected graph, the Liouvillian splits as
/// L_H(G) = L_H(G_1) ⊕ L_H(G_2) ⊕ ... and the joint kernel is the tensor product
/// of per-component kernels, giving the multiplicative formula above.
///
/// <para><b>Empirical anchors (bit-exact):</b> 4 anchors from the N=8 F1 SLOW_N8 sweep 2026-05-18 plus 1 anchor from the N=9 chain MklDirect-bridge run 2026-05-19.</para>
///
/// <list type="table">
///   <item>chain N=8, 1 component, sizes {8}: predicted 9, observed 9.</item>
///   <item>ring N=8, 1 component, sizes {8}: predicted 9, observed 9.</item>
///   <item>star N=8, 1 component, sizes {8}: predicted 9, observed 9.</item>
///   <item>K_4 + disjoint 4-chain N=8, 2 components, sizes {4, 4}: predicted 5·5 = 25, observed 25.</item>
///   <item>chain N=9, 1 component, sizes {9}: predicted 10, observed 10 (via the MklDirect bridge 2026-05-19, the first cross-N anchor verifying the connected-case formula at the new scale frontier).</item>
/// </list>
///
/// <para>The connected-graph case (kernel dim = N+1 for one component of size N)
/// is the Heisenberg specialisation of <see cref="F4StationaryModeCountPi2Inheritance"/>
/// restricted to the popcount-conserving (joint-J = 0 invariant) sector: each total-spin
/// J multiplet contributes one popcount-eigenvalue, the J=N/2 fully-symmetric multiplet
/// covers popcount {0, 1, ..., N}. The new content here is the multiplicative
/// factorisation across components, which the four N=8 data anchors verify and which
/// follows analytically from per-component closure (DEGENERACY_PALINDROME Result 2)
/// plus the standard tensor-sum kernel factorisation (see Tier 1 derived paragraph
/// below).</para>
///
/// <para><b>Tier 1 derived (2026-05-19).</b> The connected-case upper bound
/// <c>dim ker L_H(G_c) ≤ |c|+1</c> is closed by
/// <c>experiments/DEGENERACY_PALINDROME.md</c> Result 2 (magnetization conservation:
/// the only operators satisfying both <c>[H, Q] = 0</c> and <c>D(Q) = 0</c> are the
/// |c|+1 diagonal functions of total S_z, i.e. the identity plus the |c| popcount
/// projectors P_n; nothing else can be added to exhaust the kernel). Combined with
/// the popcount-projector lower-bound construction (Section 1 of
/// <c>PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md</c>) this gives the per-component
/// equality <c>dim ker L_H(G_c) = |c|+1</c>. The multi-component product follows
/// from per-component closure plus standard tensor-sum factorisation of L_H across
/// disjoint components (Section 2 of the proof file): for
/// <c>L_H(G) = L_H(G_1) ⊗ I + I ⊗ L_H(G_2)</c> the kernel factorises as
/// <c>ker L_H(G) = ker L_H(G_1) ⊗ ker L_H(G_2)</c> (standard Lyapunov/Sylvester
/// kernel result). Together: <c>dim ker L_H(G) = Π_c (|c|+1)</c>. The 4 N=8 anchors
/// remain bit-exact corroboration of the now-closed derivation (kernel dim is an
/// integer, no rounding tolerance). The per-weight ker breakdown in
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> Appendix (2026-05-17) corroborates
/// the boundary upper-bound across chain / ring / star / K_n at N=3..5.</para>
///
/// <para><b>What this claim is NOT.</b> It is not F4 itself (the Clebsch-Gordan
/// closed form for Stat(N) = Σ_J m(J)·(2J+1)² which counts ALL stationary modes,
/// not just kernel-popcount-labelled ones, on a connected Heisenberg chain).
/// F4 sees the full (2J+1)² multiplet weight per total spin; this claim sees
/// only the popcount label and how it factorises across components. The two
/// agree at the connected case for the trivial popcount-only sector (J=N/2:
/// m = 1, contributes 1 kernel mode per popcount value → N+1 total).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md</c> (primary) +
/// <c>experiments/DEGENERACY_PALINDROME.md</c> (Result 2: <c>d_real(0) = N+1</c>, proven via magnetization conservation; this is the connected-case upper-bound closure) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F4StationaryModeCountPi2Inheritance.cs</c> (parent F4 Clebsch-Gordan) +
/// <c>compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs</c> (sister Tier 2 verification record / data record) +
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> (Appendix 2026-05-17: per-weight ker breakdown corroborates the boundary upper-bound across chain/ring/star/K_n at N=3..5) +
/// the four <c>simulations/results/f1_n8_n9_metrics/&lt;topology&gt;_N8.json</c> data files.</para>
///
/// <para><b>Sister sharpenings from the same 2026-05-18 / 2026-05-19 N=8 + Python anchor sweep.</b> Three companion
/// documents in <c>experiments/</c> and <c>hypotheses/</c> capture the other structural findings the same data
/// surfaced, each extending an existing April 2026 framework rather than introducing new physics:
/// <list type="bullet">
///   <item><c>experiments/STAR_CONFOCAL_LIMIT.md</c> (Tier 1 candidate, 2026-05-19): star topology saturates
///         <c>Im(λ) = σ = N·γ</c> exactly when J = 2γ, via SU(2) / Schur-Weyl hub-spoke decomposition.
///         The point-focus limit of the OPTICAL_CAVITY_ANALYSIS framework. 5 anchors at N=3..6 + N=8.</item>
///   <item><c>hypotheses/F1_DISSIPATION_GAP_PATTERN.md</c> (Tier 3, extended 2026-05-19): cross-topology
///         cross-N gap data. Chain shows clean <c>gap × N² ≈ 2.20</c> for N ≥ 4 (5 anchors); ring and star
///         follow different scaling laws. Bond count alone does not predict the gap.</item>
///   <item><c>hypotheses/STAR_SPECTRUM_COMPACTNESS.md</c> (Tier 3, partially resolved 2026-05-19): Reading 1
///         (MaxImag = σ cap) closed by STAR_CONFOCAL_LIMIT; Reading 2 (S_(N−1) irrep multiplicity accounting
///         for the 30× distinct-binned ratio) still open.</item>
/// </list>
/// These four documents together (this claim + the three sister readings) form the May 2026 sharpening of the
/// April 2026 spectrum-structure framework (DEGENERACY_PALINDROME + WEIGHT2_KERNEL + BURES_DEGENERACY +
/// OPTICAL_CAVITY_ANALYSIS).</para></summary>
public sealed class F4KernelDimensionByComponentsClaim : Claim
{
    /// <summary>Repo-relative paths of the four N=8 SLOW_N8 sweep JSON files whose
    /// <c>KernelDimension</c> fields anchor this claim. Look here for the raw numbers;
    /// do not re-derive.</summary>
    public IReadOnlyList<string> AnchorDataFiles { get; } = new[]
    {
        "simulations/results/f1_n8_n9_metrics/chain_N8.json",
        "simulations/results/f1_n8_n9_metrics/ring_N8.json",
        "simulations/results/f1_n8_n9_metrics/star_N8.json",
        "simulations/results/f1_n8_n9_metrics/k4_plus_disjoint_4chain_N8.json",
        "simulations/results/f1_n8_n9_metrics/chain_N9.json",
    };

    /// <summary>Predict the kernel dimension of the Heisenberg Liouvillian for a graph
    /// with the given connected-component sizes. Returns <c>Π_c (|c| + 1)</c> over the
    /// component sizes. For a single connected graph of size N, pass <c>[N]</c> and the
    /// result is N+1 (the F4 popcount-sector count).
    ///
    /// <para>Throws <see cref="ArgumentNullException"/> if <paramref name="componentSizes"/>
    /// is null; <see cref="ArgumentException"/> if the list is empty or any size is non-positive
    /// (every component must have ≥ 1 qubit).</para></summary>
    public int Predict(IReadOnlyList<int> componentSizes)
    {
        if (componentSizes is null) throw new ArgumentNullException(nameof(componentSizes));
        if (componentSizes.Count == 0)
            throw new ArgumentException("componentSizes must contain at least one component.", nameof(componentSizes));

        long product = 1;
        for (int i = 0; i < componentSizes.Count; i++)
        {
            int size = componentSizes[i];
            if (size <= 0)
                throw new ArgumentException(
                    $"Each component size must be ≥ 1; got componentSizes[{i}] = {size}.",
                    nameof(componentSizes));
            product *= (long)(size + 1);
            if (product > int.MaxValue)
                throw new OverflowException(
                    $"Predicted kernel dimension overflows int32 for componentSizes = [{string.Join(", ", componentSizes)}].");
        }
        return (int)product;
    }

    /// <summary>The four N=8 SLOW_N8 anchor rows surfaced as a typed table for the
    /// inspectable tree. Each tuple is (topology label, componentSizes, predicted, observed).
    /// Observed values are taken from the JSON <c>KernelDimension</c> field of the four
    /// data files listed in <see cref="AnchorDataFiles"/>.</summary>
    public IReadOnlyList<(string Topology, IReadOnlyList<int> ComponentSizes, int Predicted, int Observed)> EmpiricalAnchorsN8 { get; } =
        new (string, IReadOnlyList<int>, int, int)[]
        {
            ("chain N=8",                       new[] { 8 },    9,  9),
            ("ring N=8",                        new[] { 8 },    9,  9),
            ("star N=8",                        new[] { 8 },    9,  9),
            ("K_4 + disjoint 4-chain N=8",      new[] { 4, 4 }, 25, 25),
            ("chain N=9 (MklDirect bridge)",    new[] { 9 },    10, 10),
        };

    public F4KernelDimensionByComponentsClaim()
        : base("F4 disconnected-graph extension: dim ker L_H = Π_c (|c|+1) over connected components c; bit-exact at N=8 across 4 topologies (chain/ring/star/K_4+disjoint-4-chain) + N=9 chain via the MklDirect bridge 2026-05-19",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md (primary) + " +
               "experiments/DEGENERACY_PALINDROME.md (Result 2: d_real(0) = N+1, proven via magnetization conservation) + " +
               "compute/RCPsiSquared.Core/Symmetry/F4StationaryModeCountPi2Inheritance.cs (parent F4 Clebsch-Gordan) + " +
               "compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs (sister Tier 2 verification record / data record) + " +
               "docs/proofs/PROOF_WEIGHT1_DEGENERACY.md (Appendix 2026-05-17: per-weight ker breakdown corroborates the boundary upper-bound) + " +
               "simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json")
    { }

    public override string DisplayName =>
        "F4 disconnected-graph extension: dim ker L_H = Π_c (|c|+1)";

    public override string Summary =>
        $"dim ker L_H = Π_c (|c|+1) over connected components c; verified bit-exact at N=8 on chain/ring/star (1 component, dim=9) and K_4 + disjoint 4-chain (2 components, dim=5·5=25), plus N=9 chain (1 component, dim=10) via the MklDirect bridge 2026-05-19; {Tier.Label()}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "For graph G with connected components c, dim ker L_H(G) = Π_c (|c|+1). Specialises to N+1 for any single connected component of size N (matches F4 popcount-sector count).");

            yield return new InspectableNode("Tier 1 derived",
                summary: "Connected-case upper bound dim ker ≤ |c|+1 closed by DEGENERACY_PALINDROME Result 2 (magnetization conservation: identity + |c| popcount projectors exhaust the kernel); multi-component product follows from tensor-sum kernel factorisation L_H(G) = L_H(G_1) ⊗ I + I ⊗ L_H(G_2). 5 bit-exact anchors corroborate (4 at N=8 + 1 at N=9 chain). See PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md § \"Upper-bound closure (resolved 2026-05-18)\".");

            yield return new InspectableNode("relationship to F4 (parent)",
                summary: "F4 itself gives Stat(N) = Σ_J m(J)·(2J+1)² (full multiplet weight on a connected Heisenberg chain); this claim restricts to the popcount label and extends multiplicatively across disconnected components. F4 + this claim cover both axes (multiplet structure within a component, factorisation across components).");

            foreach (var (topology, sizes, predicted, observed) in EmpiricalAnchorsN8)
                yield return new InspectableNode($"N=8 anchor: {topology}",
                    summary: $"components = [{string.Join(", ", sizes)}], Π(|c|+1) = {predicted}, observed kernel dim = {observed} (bit-exact)");

            yield return new InspectableNode("anchor data files",
                summary: string.Join("; ", AnchorDataFiles));
        }
    }
}
