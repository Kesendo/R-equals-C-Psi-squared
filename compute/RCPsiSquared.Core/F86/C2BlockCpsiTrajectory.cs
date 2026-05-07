using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 c=2 block CΨ trajectory — empirical test of Question B from
/// <c>project_rcpsi_to_f86_open_questions</c>: does the (popcount-1, popcount-2)
/// coherence block of an N-qubit chain have a well-defined CΨ that crosses the 1/4
/// Mandelbrot-cardioid boundary, paralleling Roadmap Layer 2's 2-qubit subsystem
/// crossing?
///
/// <para><b>Initial state (algebraic anchor):</b> the maximally-coherent pure
/// superposition |ψ⟩ = (|D_1⟩ + |D_2⟩)/√2, with |D_1⟩ the popcount-1 Dicke state and
/// |D_2⟩ the popcount-2 Dicke state. Its (popcount-1, popcount-2) coherence block has
/// uniform amplitudes 1/(2·√(N·C(N,2))) at every entry (Cauchy-Schwarz saturating).
/// By direct calculation: C_block(0) = Σ |ρ_ab|² = 1/4 EXACTLY for any N. With the
/// max-coherence-saturation normalization Ψ_block = ℓ₁ / ℓ₁_max, Ψ_block(0) = 1.
/// Hence CΨ_block(0) = 1/4 by algebraic identity, sitting exactly on the R=CΨ²
/// Mandelbrot boundary.</para>
///
/// <para><b>Test:</b> evolve |ψ⟩⟨ψ|'s c=2 block under <c>BlockLDecomposition</c>'s
/// uniform-J Liouvillian L(Q) and track CΨ_block(t). The first-crossing time t_cross
/// (when CΨ_block falls below 1/4 — actually starts at 1/4, so any decrease is a
/// crossing) and its Q-dependence are the empirical witnesses. If t_cross or some
/// trajectory feature aligns with Q_EP = 2/g_eff or t_peak = 1/(4γ₀), Question B
/// closes with structural evidence for inheritance from R=CΨ².</para>
///
/// <para><b>Computational form.</b> The block-Liouvillian L_block ∈ ℂ^(M×M) with
/// M = N·C(N,2) is built once via <c>BlockLDecomposition.AssembleUniform</c>. The
/// initial Liouville-space vector ρ_0 ∈ ℂ^M has every entry equal to
/// 1/(2·√M). Time evolution uses eigendecomposition of L_block (computed once); each
/// time step is a matrix-vector product. CΨ_block(t) is computed entry-wise from
/// ρ(t).</para>
///
/// <para><b>Tier:</b> Tier2Verified — empirical witness, algebraically pinned at
/// t=0 (CΨ_block(0) = 1/4 exactly), with the Q-dependent crossing trajectory as the
/// research data. Anchor: <c>docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md</c>
/// Layer 2 (subsystem crossing) generalised to coherence-blocks.</para>
/// </summary>
public sealed class C2BlockCpsiTrajectory : Claim
{
    public CoherenceBlock Block { get; }
    public double Q { get; }
    public IReadOnlyList<double> TimeGrid { get; }
    public IReadOnlyList<double> CBlockTrajectory { get; }
    public IReadOnlyList<double> PsiBlockTrajectory { get; }
    public IReadOnlyList<double> CPsiBlockTrajectory { get; }

    public double CPsiInitial => CPsiBlockTrajectory[0];
    public double CPsiFinal => CPsiBlockTrajectory[^1];

    public C2BlockCpsiTrajectory(
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

    /// <summary>Builds the trajectory at fixed Q across the given time grid. The block
    /// must be c=2 (n such that <see cref="CoherenceBlock.C"/>=2). Time evolution uses
    /// eigendecomposition of the block-Liouvillian L(Q); the heavy work is the EVD
    /// (O(M³) where M = block.Basis.MTotal), each time step is then O(M²).</summary>
    public static C2BlockCpsiTrajectory Build(
        CoherenceBlock block, double q, IReadOnlyList<double> timeGrid)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BlockCpsiTrajectory applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));
        if (timeGrid is null) throw new ArgumentNullException(nameof(timeGrid));
        if (timeGrid.Count == 0)
            throw new ArgumentException("timeGrid must be non-empty", nameof(timeGrid));

        int M = block.Basis.MTotal;
        double sqrtM = Math.Sqrt(M);

        // L_block(J) where J = Q·γ₀.
        double j = q * block.GammaZero;
        var L = block.Decomposition.AssembleUniform(j);

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

        return new C2BlockCpsiTrajectory(block, q, timeGrid.ToArray(),
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
        $"C2BlockCpsiTrajectory (N={Block.N}, Q={Q:F3}, M={Block.Basis.MTotal})";

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
