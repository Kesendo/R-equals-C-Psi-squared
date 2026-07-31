using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The isotropic Heisenberg star under UNIFORM Z-dephasing has the
/// imaginary spectral spread
///
/// <code>
///     Im_max(star, N, J)  =  J · N / 2          for all N ≥ 3, J ≥ 0,
/// </code>
///
/// independent of γ. Equivalently <c>Im/σ = Q/2</c> with <c>σ = N·γ</c>.
///
/// <para><b>What is star-specific and what is not.</b> Three statements have to be
/// kept apart. (1) <c>ΔE_max(H_star) = J·N/2</c> is star-specific and exact.
/// (2) The SATURATION <c>max|Im λ_L| = ΔE_max(H)</c> is NOT star-specific: it
/// holds for every isotropic Heisenberg graph (chain, ring, complete, an
/// asymmetric |Aut|=1 graph, disconnected), so saturating carries no star
/// information; Section 6 of the proof derives it in general.
/// (3) J·N/2 is the MINIMUM of ΔE_max over connected graphs, with the star the
/// unique minimiser, verified exhaustively at N = 4, 5, 6.</para>
///
/// <para><b>Why this works.</b> The star Hamiltonian factors through the total
/// leaf spin <c>S_L = Σ_{k=1}^{N-1} S_k</c>: <c>H_star = J · S_0 · S_L</c>. Its
/// Casimir spectrum within each fixed-S_L ≥ 1/2 sector is two-fold
/// (S_tot = S_L ± 1/2), with energy gap <c>ΔE(S_L) = J·(S_L + 1/2)</c>. Maximum
/// at S_L = (N-1)/2: <c>ΔE_max = J·N/2</c>.
///
/// The upper bound is a numerical-range argument: D is self-adjoint in the
/// Hilbert-Schmidt product because the jump operators Z_l are HERMITIAN, so
/// <c>Im λ = ⟨v,Kv⟩/i</c> lies in the convex hull of the spectrum of
/// <c>K = -i[H,·]</c>. Hermiticity is the load-bearing hypothesis, not the
/// informal "adds only real decay": for H = 0 and the non-Hermitian jump
/// c = I + iY the spectrum is {0, 0, -2±2i}, breaking the bound outright.
///
/// The realisers are <c>|β_k⟩⟨ferro|</c>, where |ferro⟩ is a fully polarised
/// computational state and |β_k⟩ the E_min state at Hamming rung k = 1..N-1;
/// uniform dephasing acts on these as the SCALAR -2γk, giving
/// <c>λ = -2γk + i·J·N/2</c> (the conjugate mode |ferro⟩⟨β_k| carries -i), giving
/// 2(N-1) distinct eigenvalues of multiplicity 2, i.e. 4(N-1) with multiplicity. A generic
/// <c>|Ψ_+⟩⟨Ψ_−|</c> between the extremal multiplets is NOT an eigenoperator
/// (worst residual 0.94/1.15/1.20 at N=3/4/5); the extremal block is not
/// D-invariant. See <c>docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md</c>.</para>
///
/// <para><b>Scope.</b> Uniform γ is required: with site-dependent γ_l the
/// saturation fails (N=4, γ=(2.0,0.5,0.5,0.5) gives 1.5927 against 2.0), though
/// the upper bound survives. The fully polarised extreme is required: XY and
/// transverse-field Ising chains do not saturate under the same dephasing.</para>
///
/// <para><b>Empirical anchors: 25 distinct (N, Q) combinations over 29 runs, agreeing
/// to a worst relative deviation of 1.98e-14. These are machine-precision agreements,
/// not bit-exact: no stored anchor equals J·N/2 exactly in floating point.</b></para>
///
/// <list type="bullet">
///   <item>Q-sweep at γ₀ = 0.05, 24 anchors (Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5}
///         × N ∈ {3, 4, 5, 6}): all Im/σ = Q/2 to machine precision (worst rel.
///         err 1.98e-14, at N=6). Five of the six Q values are canonical anchors
///         from <c>docs/Q_REGIME_ANCHORS.md</c>; Q = 0.5 is not in that table.
///         Files: <c>simulations/results/q_sweep_anchor/star_N{3..6}_Q*.json</c>.</item>
///   <item>N=8 anchor at the Marrakesh convention (γ=0.5, J=1, Q=2): Im_max =
///         4.000000000000002, σ = 4. From the SLOW_N8 sweep (commit 89f725e),
///         file <c>star_N8.json</c>.</item>
///   <item>Python re-runs at γ=0.5/J=1 for N=3..6: <c>star_N{3..6}_python.json</c>
///         (the same Q=2 column as the sweep, and the same code path, so these
///         confirm reproducibility rather than cross-implementation agreement).</item>
/// </list>
///
/// <para><b>Sister bound (same Casimir technique, N=4-specific).</b> The
/// <see cref="RingN4DihedralLockClaim"/> covers the analogous Im_max = (3/4)·J·N
/// = 3·J saturation on the 4-cycle = K_{2,2}, where the sublattice split is
/// |A| = |B| = 2 instead of star's |A| = 1, |B| = N-1. The two claims share
/// the same proof skeleton: bipartite split → all-pairs bonding →
/// H = J·S_A·S_B Casimir form → maximum-S_tot ferromagnet eigenmode
/// realises the bound. See PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md
/// "Which topologies admit the elementary Casimir derivation" for the joint
/// characterisation. Note that skeleton yields the CLOSED FORM; the saturation
/// itself is common to every Heisenberg topology.</para>
///
/// <para><b>Cavity-picture reading.</b> This claim is the typed-statement form of
/// <c>experiments/STAR_CONFOCAL_LIMIT.md</c> (the point-focus limit of the
/// optical-cavity framework). The Marrakesh-convention reading "Im/σ = 1 when
/// J = 2γ" is the Q = 2 row of the Im/σ = Q/2 lock; the γ-independent statement
/// makes clear the closed form is a property of star geometry, not of any
/// particular (J, γ) point. The SATURATION, by contrast, is common to every
/// Heisenberg topology and is not the star's (see the node below).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md</c> (primary derivation) +
/// <c>experiments/STAR_CONFOCAL_LIMIT.md</c> (the cavity-picture sibling reading) +
/// <c>compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs</c> (sister bound via same Casimir technique) +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> (parent F1) +
/// 24 <c>simulations/results/q_sweep_anchor/star_N{3..6}_Q*.json</c> + <c>simulations/results/f1_n8_n9_metrics/star_N8.json</c>.</para>
///
/// <para><b>Companion sister claims from the same 2026-05-19 May-sprint
/// promotion:</b> <see cref="RingN4DihedralLockClaim"/> (N=4-specific Im-max bound)
/// and <see cref="F4KernelDimensionByComponentsClaim"/> (kernel-dim factorisation
/// across components). All three landed Tier 1 derived 2026-05-19 from the same
/// SLOW_N8 + N=9-bridge + Q-sweep sprint.</para></summary>
public sealed class StarImMaxBoundClaim : Claim
{
    /// <summary>The dimensionless saturation prefactor: Im_max = Coefficient · J · N
    /// = (1/2) · J · N. Equivalently Im/σ = Coefficient · Q = (1/2) · Q.</summary>
    public double Coefficient { get; } = 1.0 / 2.0;

    /// <summary>The minimum N at which the star topology is well-defined (a hub
    /// plus at least two leaves). For N=2 the "star" is just a single bond
    /// (chain N=2), and for N=3 the star = chain N=3 by graph isomorphism;
    /// the saturation result holds for any N ≥ 3, including the N=3 isomorphic
    /// case where chain and ring also saturate (triangle = K_3 symmetry).</summary>
    public int MinN { get; } = 3;

    /// <summary>Repo-relative paths of the Q-sweep + Marrakesh-convention JSON
    /// files whose <c>MaxImag</c> fields anchor this claim. Look here for the raw
    /// numbers; do not re-derive.</summary>
    public IReadOnlyList<string> AnchorDataFiles { get; } = new[]
    {
        // 24 Q-sweep anchors at γ₀=0.05 (Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} × N ∈ {3,4,5,6}).
        "simulations/results/q_sweep_anchor/star_N3_Q0.5000.json",
        "simulations/results/q_sweep_anchor/star_N3_Q1.0000.json",
        "simulations/results/q_sweep_anchor/star_N3_Q1.5000.json",
        "simulations/results/q_sweep_anchor/star_N3_Q1.7321.json",
        "simulations/results/q_sweep_anchor/star_N3_Q2.0000.json",
        "simulations/results/q_sweep_anchor/star_N3_Q2.5000.json",
        "simulations/results/q_sweep_anchor/star_N4_Q0.5000.json",
        "simulations/results/q_sweep_anchor/star_N4_Q1.0000.json",
        "simulations/results/q_sweep_anchor/star_N4_Q1.5000.json",
        "simulations/results/q_sweep_anchor/star_N4_Q1.7321.json",
        "simulations/results/q_sweep_anchor/star_N4_Q2.0000.json",
        "simulations/results/q_sweep_anchor/star_N4_Q2.5000.json",
        "simulations/results/q_sweep_anchor/star_N5_Q0.5000.json",
        "simulations/results/q_sweep_anchor/star_N5_Q1.0000.json",
        "simulations/results/q_sweep_anchor/star_N5_Q1.5000.json",
        "simulations/results/q_sweep_anchor/star_N5_Q1.7321.json",
        "simulations/results/q_sweep_anchor/star_N5_Q2.0000.json",
        "simulations/results/q_sweep_anchor/star_N5_Q2.5000.json",
        "simulations/results/q_sweep_anchor/star_N6_Q0.5000.json",
        "simulations/results/q_sweep_anchor/star_N6_Q1.0000.json",
        "simulations/results/q_sweep_anchor/star_N6_Q1.5000.json",
        "simulations/results/q_sweep_anchor/star_N6_Q1.7321.json",
        "simulations/results/q_sweep_anchor/star_N6_Q2.0000.json",
        "simulations/results/q_sweep_anchor/star_N6_Q2.5000.json",
        // N=8 Marrakesh-convention anchor (Q=2 at γ=0.5, J=1) from SLOW_N8 sweep.
        "simulations/results/f1_n8_n9_metrics/star_N8.json",
        // Python re-runs at γ=0.5, J=1 (Q=2) for N=3..6. Same Q-column as the
        // q_sweep_anchor Q=2 row above AND the same code path: f1_q_sweep_anchor.py
        // imports its `run` from f1_topology_heisenberg_small_n_anchor.py, which
        // builds L with framework.lindblad.lindbladian_z_dephasing. So these four
        // confirm reproducibility, not cross-implementation agreement.
        "simulations/results/f1_n8_n9_metrics/star_N3_python.json",
        "simulations/results/f1_n8_n9_metrics/star_N4_python.json",
        "simulations/results/f1_n8_n9_metrics/star_N5_python.json",
        "simulations/results/f1_n8_n9_metrics/star_N6_python.json",
    };

    /// <summary>Predict the maximum imaginary spectral spread of the
    /// Heisenberg + Z-dephasing Liouvillian on the N-site star for a given
    /// Hamiltonian coupling J. Returns <c>(1/2) · J · N</c>.
    ///
    /// <para>Throws <see cref="ArgumentOutOfRangeException"/> if <paramref name="N"/>
    /// is below <see cref="MinN"/>; <see cref="ArgumentException"/> if
    /// <paramref name="J"/> is not finite or is negative.</para></summary>
    public double Predict(int N, double J)
    {
        if (N < MinN)
            throw new ArgumentOutOfRangeException(nameof(N),
                $"Star saturation derivation requires N ≥ {MinN}; got N = {N}.");
        CasimirBoundClaimHelpers.RequireFiniteNonNegative(J, nameof(J));
        return Coefficient * J * N;
    }

    /// <summary>Predict the dimensionless ratio Im_max / σ where σ = N·γ for given
    /// Q = J/γ. Returns <c>(1/2) · Q</c> (independent of N, by Q-universality).</summary>
    public double PredictImOverSigma(double Q)
    {
        CasimirBoundClaimHelpers.RequireFiniteNonNegative(Q, nameof(Q));
        return Coefficient * Q;
    }

    /// <summary>The 24 Q-sweep anchor rows surfaced as a typed table for the
    /// inspectable tree (4 N-values × 6 Q-values). Each tuple is (Q label, N, Q value,
    /// predicted Im/σ). These are the PREDICTIONS of <see cref="PredictImOverSigma"/>,
    /// not measurements: this class does not read the JSON files, and an earlier
    /// version that reported an "observed" column was populating it from the
    /// prediction, so it could never disagree. The measured agreement is recorded
    /// once, as <see cref="WorstAnchorRelativeDeviation"/>; the raw numbers live in
    /// the files listed in <see cref="AnchorDataFiles"/> and are re-checked by
    /// <c>simulations/star_saturation_gate.py</c>. Q labels and Q values come from
    /// <see cref="CasimirBoundClaimHelpers.QSweepAnchorLabels"/> (shared with
    /// <see cref="RingN4DihedralLockClaim.EmpiricalAnchors"/>).</summary>
    public IReadOnlyList<(string QLabel, int N, double Q, double PredictedImOverSigma)> EmpiricalAnchors { get; } =
        BuildEmpiricalAnchors();

    /// <summary>The worst relative deviation between a stored anchor's <c>MaxImag</c>
    /// and the predicted J·N/2, over all 29 runs (24 Q-sweep + N=8 + 4 Python
    /// re-runs). Measured 1.98e-14, largest at N=6. Non-zero at every anchor, which
    /// is why this claim says "machine precision" and not "bit-exact".</summary>
    public double WorstAnchorRelativeDeviation { get; } = 1.98e-14;

    private static (string, int, double, double)[] BuildEmpiricalAnchors()
    {
        var labels = CasimirBoundClaimHelpers.QSweepAnchorLabels;
        var rows = new (string, int, double, double)[labels.Count * 4];
        int i = 0;
        for (int N = 3; N <= 6; N++)
            foreach (var (label, q) in labels)
                rows[i++] = (label, N, q, 0.5 * q);
        return rows;
    }

    public StarImMaxBoundClaim()
        : base("Im_max(star, N, J) = J·N/2 at uniform γ, independent of γ; closed-form via the H_star = J·S_0·S_L hub-leaf Casimir, realised by the |β_k⟩⟨ferro| coherences; J·N/2 is the MINIMUM of ΔE_max over connected graphs (star the unique minimiser, exhaustive at N ≤ 6), while the saturation itself is common to every Heisenberg topology; 25 distinct (N, Q) anchors over 29 runs from the 2026-05-19 Q-sweep and the SLOW_N8 sweep, worst relative deviation 1.98e-14",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md (primary derivation) + " +
               "experiments/STAR_CONFOCAL_LIMIT.md (the cavity-picture sibling reading) + " +
               "compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs (sister bound via same Casimir technique) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (parent F1) + " +
               "simulations/results/q_sweep_anchor/star_N{3..6}_Q*.json (24 anchors) + " +
               "simulations/results/f1_n8_n9_metrics/star_N8.json (1 anchor at N=8, Q=2)")
    { }

    public override string DisplayName =>
        "Star Im-max: Im_max(star, N, J) = J·N/2, the minimal Heisenberg spread";

    public override string Summary =>
        $"Im_max(star, N, J) = (1/2)·J·N at uniform γ, independent of γ; closed-form via SU(2)/Schur-Weyl hub-leaf Casimir on H_star = J·S_0·S_L; J·N/2 is the minimum of ΔE_max over connected graphs (star the unique minimiser, exhaustive N ≤ 6), the saturation itself being common to every Heisenberg topology; verified at 24 Q-sweep anchors (Q ∈ {{0.5, 1.0, 1.5, √3, 2.0, 2.5}} × N ∈ {{3, 4, 5, 6}}) plus N=8 at Q=2, worst relative deviation 1.98e-14; {Tier.Label()}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "Im_max(star, N, J) = (1/2)·J·N for all N ≥ 3, J ≥ 0, under UNIFORM Z-dephasing; equivalently Im/σ = (1/2)·Q. Independent of γ.");

            yield return new InspectableNode("Tier 1 derived",
                summary: "Closed form via star Casimir factorisation: H_star = J·S_0·S_L = (J/2)·(S²_tot − 3/4 − S²_L); within a fixed S_L ≥ 1/2 sector the hub couples to S_L ± 1/2 with gap J·(S_L + 1/2); maximum at S_L = (N-1)/2 gives ΔE_max = J·N/2. Upper bound by numerical range: D is Hilbert-Schmidt self-adjoint because the jump operators Z_l are Hermitian, so Im λ lies in the convex hull of the spectrum of -i[H,·]. Realised by |β_k⟩⟨ferro| (fully polarised extreme against the E_min state at Hamming rung k), on which uniform dephasing acts as the scalar -2γk, giving λ = -2γk + i·J·N/2 (conjugate mode -i), i.e. 4(N-1) eigenvalues with multiplicity. See PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md.");

            yield return new InspectableNode("what is NOT star-specific",
                summary: "The saturation max|Im λ_L| = ΔE_max(H) holds for every isotropic Heisenberg topology tested (chain, ring, complete, T-shape, disconnected), so it carries no star information. Star-specific are the closed form ΔE_max(H_star) = J·N/2 and the fact that this is the MINIMUM over connected graphs: exhaustive search at N = 4, 5, 6 (38, 728, 26704 connected graphs) gives minimum exactly J·N/2, attained by the N stars and nothing else.");

            yield return new InspectableNode("scope fences",
                summary: "Uniform γ required: site-dependent γ_l breaks the saturation (N=4, γ=(2.0,0.5,0.5,0.5) gives 1.5927 vs 2.0), though the upper bound survives. Ferromagnetic extreme required: XY and transverse-field Ising chains do not saturate. Hermitian jumps required, or the BOUND fails: H = 0 with c = I + iY has spectrum {0, 0, -2±2i} against a bound of 0. Gate: simulations/star_saturation_gate.py, 129 checks.");

            yield return new InspectableNode("relationship to ring N=4",
                summary: "Same skeleton for the CLOSED FORM (bipartite split → all-pairs bonding → H = J·S_A·S_B Casimir). Star uses sublattice sizes |A|=1, |B|=N-1; Ring N=4 uses |A|=|B|=2. Star gives (1/2)·J·N for any N; Ring N=4 gives (3/4)·J·N N=4-specifically (the bipartite-complete coincidence C_N = K_{N/2,N/2} fails for N > 4). What the skeleton yields is the gap value, not the saturation, which needs no bipartite structure.");

            yield return new InspectableNode("Q-sweep anchor count",
                summary: "24 anchors at γ₀=0.05 from the 2026-05-19 Q-sweep (Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5} × N ∈ {3, 4, 5, 6}). All Im/σ = Q/2 to machine precision, worst relative deviation 1.98e-14.");

            yield return new InspectableNode("N=8 Marrakesh anchor",
                summary: "Im_max = 4.000000000000002, σ = N·γ = 4, so Im/σ = 1.0 to machine precision (relative deviation 4.4e-16). From the SLOW_N8 sweep 2026-05-18, star_N8.json. Equivalent statement: Im_max = J·N/2 = 4 at J=1, N=8.");

            yield return new InspectableNode("anchor data files",
                summary: $"{AnchorDataFiles.Count} JSON files anchor this claim; see AnchorDataFiles property for full list.");
        }
    }
}
