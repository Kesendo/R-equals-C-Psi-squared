using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 (popcount-n, popcount-(n+1)) block CΨ trajectory — empirical test of
/// Question B from <c>project_rcpsi_to_f86_open_questions</c>: does the c-stratum
/// coherence block of an N-qubit chain inherit the 1/4 Mandelbrot-cardioid boundary
/// of R=CΨ²? Generalised to any chromaticity c ≥ 2 (lifted from c=2-only on 2026-05-07).
///
/// <para><b>Initial state (algebraic anchor):</b> the maximally-coherent pure
/// superposition |ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2, with |D_k⟩ the popcount-k Dicke state.
/// Its (popcount-n, popcount-(n+1)) coherence block has uniform amplitudes
/// 1/(2·√(C(N,n)·C(N,n+1))) at every entry (Cauchy-Schwarz saturating). By direct
/// counting: C_block(0) = Σ |ρ_ab|² = M_block/(4·M_block) = 1/4 EXACTLY for any N
/// AND any chromaticity c. **Chromaticity-universal.** With the max-coherence-
/// saturation normalization Ψ_block = ℓ₁ / ℓ₁_max, Ψ_block(0) = 1. Hence CΨ_block(0)
/// = 1/4 by algebraic identity, sitting exactly on the R=CΨ² Mandelbrot boundary at
/// any c.</para>
///
/// <para><b>Closed-form trajectory under pure dephasing:</b>
/// <code>
/// C_block(t) = (1/4) · Σ_{k=0..c-1} (M_{HD=2k+1} / M_block) · exp(−4γ·(2k+1)·t)
/// </code>
/// where M_{HD=2k+1} = C(N, n−k)·C(N−n+k, k)·C(N−n, k+1) is the count of basis pairs
/// at Hamming distance 2k+1, and M_block = C(N,n)·C(N,n+1). For c=2, n=1: yields
/// (1/(2N))·exp(−4γt) + ((N−2)/(4N))·exp(−12γt). For c=3, n=2: three exponentials.
/// Trajectory is dephasing-dominated under uniform J (channel-uniform initial sits in
/// the H-kernel via F73 sum-rule).</para>
///
/// <para><b>Computational form.</b> The block-Liouvillian L_block ∈ ℂ^(M×M) with
/// M = M_block is built once via <c>BlockLDecomposition.AssembleUniform</c>. The
/// initial Liouville-space vector ρ_0 ∈ ℂ^M has every entry equal to 1/(2·√M).
/// Time evolution uses eigendecomposition of L_block (computed once); each time step
/// is a matrix-vector product.</para>
///
/// <para><b>Tier:</b> Tier1Derived for the algebraic identity CΨ_block(0) = 1/4 at any
/// (N, c, n). Tier2Verified for the closed-form trajectory (verified numerically;
/// proof of channel-uniform = H-kernel structure is the F73 sum-rule, established).
/// Anchor: <c>docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md</c> Layer 2 (subsystem
/// crossing) generalised to coherence-blocks at arbitrary chromaticity.</para>
/// </summary>
public sealed class BlockCpsiTrajectory : Claim
{
    public CoherenceBlock Block { get; }
    public double Q { get; }
    public IReadOnlyList<double> TimeGrid { get; }
    public IReadOnlyList<double> CBlockTrajectory { get; }
    public IReadOnlyList<double> PsiBlockTrajectory { get; }
    public IReadOnlyList<double> CPsiBlockTrajectory { get; }

    public double CPsiInitial => CPsiBlockTrajectory[0];
    public double CPsiFinal => CPsiBlockTrajectory[^1];

    public BlockCpsiTrajectory(
        CoherenceBlock block, double q,
        IReadOnlyList<double> timeGrid,
        IReadOnlyList<double> cTrajectory,
        IReadOnlyList<double> psiTrajectory,
        IReadOnlyList<double> cPsiTrajectory)
        : base("c=2 block CΨ trajectory (Layer 2 inheritance test)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md Layer 2")
    {
        Block = block;
        Q = q;
        TimeGrid = timeGrid;
        CBlockTrajectory = cTrajectory;
        PsiBlockTrajectory = psiTrajectory;
        CPsiBlockTrajectory = cPsiTrajectory;
    }

    /// <summary>Builds the trajectory at fixed uniform Q across the given time grid.
    /// Equivalent to <see cref="BuildPerBond"/> with <c>perBondCouplings = [Q·γ₀, Q·γ₀, ...]</c>.
    /// The block must be c=2.</summary>
    public static BlockCpsiTrajectory Build(
        CoherenceBlock block, double q, IReadOnlyList<double> timeGrid)
    {
        var bondCouplings = Enumerable.Repeat(q * block.GammaZero, block.NumBonds).ToArray();
        return BuildPerBond(block, q, bondCouplings, timeGrid);
    }

    /// <summary>Builds the trajectory at non-uniform per-bond couplings J_b. The
    /// reference Q is recorded in <see cref="Q"/> (the uniform-J value the per-bond
    /// configuration is closest to). Engages Q-dependence: at uniform J the
    /// channel-uniform Dicke probe lives entirely in the H-kernel (per F73 sum-rule),
    /// so the trajectory is Q-independent. Non-uniform J breaks the F73 sum-rule —
    /// per-bond <c>V_b[α, j]</c> cross-block entries survive — and the trajectory
    /// becomes Q-dependent.</summary>
    public static BlockCpsiTrajectory BuildPerBond(
        CoherenceBlock block, double q,
        IReadOnlyList<double> bondCouplings,
        IReadOnlyList<double> timeGrid)
    {
        if (block.C < 2)
            throw new ArgumentException(
                $"BlockCpsiTrajectory requires chromaticity c ≥ 2; got c={block.C} (N={block.N}, n={block.LowerPopcount}). " +
                "c=1 stratum has only HD=1 channel (no inter-channel coupling) — trivial single-exponential decay.",
                nameof(block));
        if (timeGrid is null) throw new ArgumentNullException(nameof(timeGrid));
        if (timeGrid.Count == 0)
            throw new ArgumentException("timeGrid must be non-empty", nameof(timeGrid));

        int M = block.Basis.MTotal;
        double sqrtM = Math.Sqrt(M);

        // L_block from per-bond couplings (handles uniform and non-uniform alike).
        var L = block.Decomposition.AssembleAt(bondCouplings);

        // Initial Liouville-space vector ρ_0: uniform value 1/(2·√M) at every block entry.
        // This represents the (popcount-1, popcount-2) coherence content of the pure-state
        // |ψ⟩ = (|D_1⟩ + |D_2⟩)/√2; algebraically C_block(0) = M · |ρ_0|² = M / (4M) = 1/4.
        double initVal = 1.0 / (2.0 * sqrtM);
        var rho0 = ComplexVector.Build.Dense(M, _ => new Complex(initVal, 0.0));

        // Eigendecomposition for repeated time-stepping: exp(L·t)·v = R·diag(exp(λ_i·t))·R⁻¹·v.
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var lambda = evd.EigenValues;
        // Pre-compute coefficients c_i = (R⁻¹·ρ_0)_i.
        var c0 = Rinv * rho0;

        // Max ℓ₁ (Cauchy-Schwarz bound under our diagonal-occupation): for the initial
        // diagonal occupation pattern (popcount-1 sector at 1/(2N), popcount-2 sector at
        // 1/(2·C(N,2)) per state), |ρ_ab|_max = 1/(2·√(N·C(N,2))) = 1/(2·√M). The
        // saturating ℓ₁ across all M entries is M · 1/(2·√M) = √M/2. We use this as
        // the Ψ normalisation: Ψ_block = ℓ₁ / (√M/2). Initial Ψ_block(0) = 1 by saturation.
        double ell1Max = sqrtM / 2.0;

        var cTraj = new List<double>(timeGrid.Count);
        var psiTraj = new List<double>(timeGrid.Count);
        var cPsiTraj = new List<double>(timeGrid.Count);

        var expDiag = ComplexVector.Build.Dense(M);
        foreach (double t in timeGrid)
        {
            for (int i = 0; i < M; i++) expDiag[i] = Complex.Exp(lambda[i] * t) * c0[i];
            var rhoT = R * expDiag;

            // C_block = Σ |ρ_ab|² over block (block-purity content).
            // Ψ_block = ℓ₁(ρ_block) / ℓ₁_max where ℓ₁_max = √M/2 (Cauchy-Schwarz
            // saturating bound under initial diagonal occupation).
            double C = 0.0;
            double ell1 = 0.0;
            for (int i = 0; i < M; i++)
            {
                double mag = rhoT[i].Magnitude;
                C += mag * mag;
                ell1 += mag;
            }
            double psi = ell1 / ell1Max;
            cTraj.Add(C);
            psiTraj.Add(psi);
            cPsiTraj.Add(C * psi);
        }

        return new BlockCpsiTrajectory(block, q, timeGrid.ToArray(),
            cTraj, psiTraj, cPsiTraj);
    }

    /// <summary>First time at which CΨ_block falls strictly below the given threshold.
    /// Returns <c>null</c> if it never does within the time grid.</summary>
    public double? FirstCrossingTime(double threshold = 0.25)
    {
        for (int i = 0; i < CPsiBlockTrajectory.Count; i++)
            if (CPsiBlockTrajectory[i] < threshold)
                return TimeGrid[i];
        return null;
    }

    public override string DisplayName =>
        $"BlockCpsiTrajectory (N={Block.N}, Q={Q:F3}, M={Block.Basis.MTotal})";

    public override string Summary =>
        $"CΨ_block: {CPsiInitial:F4} → {CPsiFinal:F4} at Q={Q:F3}, " +
        $"first cross 1/4 at t={FirstCrossingTime():F3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", Block.N);
            yield return InspectableNode.RealScalar("M (block dim)", Block.Basis.MTotal);
            yield return InspectableNode.RealScalar("Q", Q, "F4");
            yield return InspectableNode.RealScalar("γ₀", Block.GammaZero, "F4");
            yield return InspectableNode.RealScalar("CΨ_block(t=0)", CPsiInitial, "F6");
            yield return InspectableNode.RealScalar("CΨ_block(t=last)", CPsiFinal, "F6");
            yield return InspectableNode.RealScalar("t_peak (= 1/(4γ₀))", 1.0 / (4.0 * Block.GammaZero), "F4");
            double? tCross = FirstCrossingTime();
            yield return new InspectableNode("first crossing 1/4",
                summary: tCross.HasValue ? $"t={tCross.Value:F4}" : "never (within grid)");

            // Trajectory points as flat children so the inspect tree shows the data.
            var rows = new List<IInspectable>();
            for (int i = 0; i < TimeGrid.Count; i++)
            {
                rows.Add(new InspectableNode(
                    $"t={TimeGrid[i]:F3}",
                    summary: $"C={CBlockTrajectory[i]:F4}, Ψ={PsiBlockTrajectory[i]:F4}, CΨ={CPsiBlockTrajectory[i]:F4}"));
            }
            yield return InspectableNode.Group("trajectory (t, C, Ψ, CΨ)", rows.ToArray());
        }
    }
}
