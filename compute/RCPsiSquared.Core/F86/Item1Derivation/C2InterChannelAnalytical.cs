using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): SVD-top inter-channel vectors |u_0⟩ ∈ HD=1 subspace
/// and |v_0⟩ ∈ HD=3 subspace of <c>V_inter = P_HD1† · M_H_total · P_HD3</c>. These are
/// the EP-partner modes of the 4-mode picture (PROOF_F86_QPEAK Item 1).
///
/// <para><b>Tier outcome: Tier2Verified.</b> The constructor probes whether a closed-form
/// ansatz lands at machine-precision overlap (≥ 1−1e-10) against the numerical
/// <see cref="InterChannelSvd"/>. If it does, <see cref="Tier"/> would be <c>Tier1Derived</c>;
/// otherwise the U0/V0 are the numerical SVD vectors (Tier2Verified) and
/// <see cref="PendingDerivationNote"/> records the ansätze tried so a future session
/// can pick up. Both paths satisfy the test contract: the analytical U0/V0 overlap the
/// numerical reference at machine precision (trivially so for Tier2Verified).</para>
///
/// <para><b>Algebraic obstruction observed at this session.</b> The top singular value
/// σ_0 of V_inter is <i>exactly degenerate</i> at even N (degeneracy = 2 at N=6 and N=8,
/// degeneracy = 1 at N=5 and N=7). The chain-mirror R splits the 2D top eigenspace into
/// one R-even and one R-odd direction with the SAME σ_0. Inside that 2D subspace the SVD
/// "top vector" is library-dependent (MathNet's SVD picks one canonical direction; numpy
/// picks another); a closed-form ansatz can only match the numerical reference at
/// machine precision if it accidentally coincides with MathNet's tiebreaker. This means
/// any honest Tier1Derived must EITHER (a) lift the test contract from "single vector
/// overlap" to "projector onto 2D top eigenspace overlap", OR (b) reverse-engineer the
/// MathNet SVD tiebreaker. Both are real work, neither was achieved in the A3 time-box.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2), Item 5
/// (σ_0 → 2√2 conjecture).</para>
/// </summary>
public sealed class C2InterChannelAnalytical : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>Top singular vector |u_0⟩ in the HD=1 subspace, lifted to the full block basis
    /// (length = Mtot = |P_n| · |P_{n+1}|, unit-normalised).</summary>
    public ComplexVector U0 { get; }

    /// <summary>Top singular vector |v_0⟩ in the HD=3 subspace, lifted to the full block basis
    /// (unit-normalised).</summary>
    public ComplexVector V0 { get; }

    /// <summary>Top singular value σ_0 of V_inter. Empirically approaches 2√2 as N grows
    /// (Item 5 conjecture); bit-exact 2√2 at N=7 to ~10⁻¹⁵.</summary>
    public double Sigma0 { get; }

    /// <summary>True iff a closed-form ansatz reached overlap ≥ 1−1e-10 against the numerical
    /// <see cref="InterChannelSvd"/>; false if U0/V0 are the numerical fallback.</summary>
    public bool IsAnalyticallyDerived { get; }

    /// <summary>Non-null iff <see cref="Tier"/> is <c>Tier2Verified</c>: a concrete summary of
    /// the ansätze that were tried, why they failed, and the most promising next directions.
    /// Visible in the inspection tree so a future session can pick up where this one left off.</summary>
    public string? PendingDerivationNote { get; }

    /// <summary>Public factory: validates c=2, runs <see cref="Resolve"/>, then constructs
    /// the instance with the Tier and vectors that <see cref="Resolve"/> returned. Use this
    /// instead of a public constructor so that <see cref="Resolution.Tier"/> is the single
    /// source of truth for the Claim's Tier — there is no second hard-coded value to keep
    /// in sync. A future session promoting to Tier1Derived only needs to return a different
    /// Resolution from <see cref="Resolve"/>; the base() call below picks it up automatically.</summary>
    public static C2InterChannelAnalytical Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2InterChannelAnalytical applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var resolved = Resolve(block);
        return new C2InterChannelAnalytical(block, resolved);
    }

    /// <summary>Private constructor: <see cref="Resolution.Tier"/> is the single source of
    /// truth for the Claim's Tier. All Tier/vector data flows from one Resolution instance,
    /// so Tier and IsAnalyticallyDerived cannot drift internally.</summary>
    private C2InterChannelAnalytical(CoherenceBlock block, Resolution resolved)
        : base("c=2 inter-channel SVD-top analytical",
               resolved.Tier,
               Item1Anchors.Root)
    {
        Block = block;
        U0 = resolved.U0;
        V0 = resolved.V0;
        Sigma0 = resolved.Sigma0;
        IsAnalyticallyDerived = resolved.IsAnalyticallyDerived;
        PendingDerivationNote = resolved.PendingDerivationNote;
    }

    /// <summary>Resolution record carrying both the Tier verdict and the resolved vectors.
    /// Returned by <see cref="Resolve"/> so the constructor can populate properties from a
    /// single SVD call; if a future session promotes to Tier1Derived, the new ansatz code
    /// returns a Resolution with <c>Tier.Tier1Derived</c> and analytical U0/V0.</summary>
    private readonly record struct Resolution(
        Tier Tier,
        ComplexVector U0,
        ComplexVector V0,
        double Sigma0,
        bool IsAnalyticallyDerived,
        string? PendingDerivationNote);

    private static Resolution Resolve(CoherenceBlock block)
    {
        // Numerical reference (always required, both as fallback and as ansatz oracle).
        var numerical = InterChannelSvd.Build(block, hd1: 1, hd2: 3);

        // === Tier1Derived ansatz attempts go here ===
        // None of the simple closed-form ansätze tried in the A3 session reached the
        // 1e-10 overlap bar. See PendingDerivationNote below for the directions tried.
        // Future sessions can attempt new ansätze and, on success, return
        //   new Resolution(Tier.Tier1Derived, u_ans, v_ans, sigma_ans, true, null);

        // For now: Tier2Verified — adopt the numerical SVD vectors as the canonical U0/V0.
        return new Resolution(
            Tier: Tier.Tier2Verified,
            U0: numerical.U0InFullBlock,
            V0: numerical.V0InFullBlock,
            Sigma0: numerical.Sigma0,
            IsAnalyticallyDerived: false,
            PendingDerivationNote: BuildPendingDerivationNote(block, numerical));
    }

    private static string BuildPendingDerivationNote(CoherenceBlock block, InterChannelSvd numerical)
    {
        int N = block.N;
        // Inspect top-eigenvalue degeneracy of V V† (governs whether |u_0⟩ direction is unique).
        var pHd1 = HdSubspaceProjector.Build(block, hammingDistance: 1);
        var pHd3 = HdSubspaceProjector.Build(block, hammingDistance: 3);
        var vInter = pHd1.ConjugateTranspose() * block.Decomposition.MhTotal * pHd3;
        var vvh = vInter * vInter.ConjugateTranspose();
        var evd = vvh.Evd();
        // Eigenvalues are real and non-negative; sort descending.
        var sorted = evd.EigenValues.Select(z => z.Real).OrderByDescending(x => x).ToArray();
        double top = sorted[0];
        int degen = sorted.Count(x => Math.Abs(x - top) / Math.Max(top, 1e-15) < 1e-8);

        return $"A3 time-box hit Tier2Verified at N={N}; σ_0={numerical.Sigma0:G6}, " +
               $"top-eigenvalue degeneracy of V V† = {degen}.\n" +
               "Algebraic obstruction: σ_0 is EXACTLY degenerate at even N (deg=2 at N=6, N=8). The\n" +
               "chain-mirror R splits the 2D top eigenspace into one R-even and one R-odd direction\n" +
               "with the same σ_0; inside that subspace, the 'top SVD vector' direction is\n" +
               "library-dependent (MathNet vs numpy disagree). Any closed-form ansatz can only land\n" +
               "machine-precision overlap (1e-10) against `numerical.U0InFullBlock` if it coincides\n" +
               "with MathNet's tiebreaker — which is implementation detail, not physics.\n" +
               "\n" +
               "Ansätze tried (overlaps at N=5, 6, 7, 8 in this session):\n" +
               "  - ψ_1(s)·ψ_1(e) symmetric: |overlap| ~ 0.016, 0.0007, 0.007, 0.001\n" +
               "  - ψ_k(s)·ψ_k(e) for k=1..5: best k=4 at N=7 → 0.246; not converging.\n" +
               "  - ψ_a(s)·ψ_b(e) cross-modes for (a,b) in 1..3 squared: peak 0.21 at (2,2) N=5; no scaling.\n" +
               "  - separable rank-1 of u_matrix u[s,e]: u_matrix is rank >= 2 across N (rank 2 at\n" +
               "    N=5, 7; rank 3 at N=6, 8). Single rank-1 ansatz cannot match.\n" +
               "\n" +
               "Promising next directions for a future session:\n" +
               "  (a) Lift the test contract: compare PROJECTOR onto top-2D eigenspace\n" +
               "      (U_top · U_top†) instead of single direction. The R-even and R-odd basis\n" +
               "      vectors of that subspace ARE library-independent.\n" +
               "  (b) Diagonalise V V† directly in the (s, e) bond-pair basis: matrix entries are\n" +
               "      sums of bond-flip overlaps, ψ_k-mode-weighted. Numerical evidence: at N=7\n" +
               "      (where deg=1) the R-even u_0 has support only on (s, e) with s+e odd or\n" +
               "      e in middle sites — the structure suggests a ψ_2-like mode coupling.\n" +
               "  (c) Bogoliubov / Jordan-Wigner approach: V_inter is the reduction of M_H_total\n" +
               "      to 1↔2 excitations; in the JW fermion picture this should diagonalise via\n" +
               "      free-fermion modes. Item 5 in PROOF_F86_QPEAK.md notes this approach for\n" +
               "      σ_0 → 2√(2(c−1)) but it has not been carried out for the singular vectors.\n" +
               "  (d) For c=2 specifically: V_inter has structure tied to bond-decomposition of\n" +
               "      M_H_per_bond[b], so |u_0⟩ may decompose as a sum over bonds with sine-mode\n" +
               "      weights. The B1/B2/B3 cross-coupling tasks may produce the right ansatz as\n" +
               "      a side-effect.";
    }

    public override string DisplayName =>
        $"c=2 SVD-top inter-channel analytical (N={Block.N})";

    public override string Summary =>
        $"σ_0={Sigma0:G6}, IsAnalyticallyDerived={IsAnalyticallyDerived} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("σ_0", Sigma0);
            yield return InspectableNode.RealScalar("σ_0 / (2√2)", Sigma0 / (2.0 * Math.Sqrt(2.0)));
            yield return new InspectableNode(
                "IsAnalyticallyDerived",
                summary: IsAnalyticallyDerived ? "true (Tier1Derived)" : "false (Tier2Verified, numerical fallback)");
            yield return InspectableNode.RealScalar("|U0| L2 norm", U0.L2Norm());
            yield return InspectableNode.RealScalar("|V0| L2 norm", V0.L2Norm());
            if (PendingDerivationNote is not null)
                yield return new InspectableNode("PendingDerivationNote", summary: PendingDerivationNote);
        }
    }
}
