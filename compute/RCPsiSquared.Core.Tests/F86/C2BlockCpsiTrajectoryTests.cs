using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 c=2 block CΨ trajectory tests for Question B from
/// <c>project_rcpsi_to_f86_open_questions</c>: does the (popcount-1, popcount-2)
/// coherence block of an N-qubit chain have a well-defined CΨ that crosses 1/4?</summary>
public class C2BlockCpsiTrajectoryTests
{
    /// <summary>The algebraic anchor: CΨ_block(t=0) = 1/4 exactly for the maximally-coherent
    /// pure-state initial |ψ⟩ = (|D_1⟩ + |D_2⟩)/√2 at any N. This is a pencil-paper
    /// identity (M·|ρ_ab|² = M·(1/(2√M))² = 1/4 with Cauchy-Schwarz saturating ℓ₁ giving
    /// Ψ=1), pinned here as Tier 1 derived structural anchor.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void InitialCPsiBlock_EqualsOneQuarter_ExactlyAcrossN(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        // Single time point at t=0 to measure the initial value alone.
        var trajectory = C2BlockCpsiTrajectory.Build(
            block, q: 1.0, timeGrid: new[] { 0.0 });
        double cpsi0 = trajectory.CPsiInitial;
        Assert.True(Math.Abs(cpsi0 - 0.25) < 1e-10,
            $"N={N}: CΨ_block(0) = {cpsi0:F12}, expected 0.25 within 1e-10. " +
            $"Algebraic identity should hold exactly.");
    }

    /// <summary>The structural anchor: C_block(0) = 1/4 alone, before factoring in Ψ.
    /// Σ |ρ_ab|² over the (popcount-1, popcount-2) block of |ψ⟩⟨ψ| with
    /// ψ = (|D_1⟩+|D_2⟩)/√2 equals M · (1/(2·√M))² = 1/4 by direct count.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void InitialCBlock_EqualsOneQuarter_Independent(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var trajectory = C2BlockCpsiTrajectory.Build(
            block, q: 1.0, timeGrid: new[] { 0.0 });
        double c0 = trajectory.CBlockTrajectory[0];
        Assert.True(Math.Abs(c0 - 0.25) < 1e-10,
            $"N={N}: C_block(0) = {c0:F12}, expected 0.25 (sum |ρ_ab|² over block).");
    }

    /// <summary>The Cauchy-Schwarz saturation: Ψ_block(0) = 1 with ℓ₁/ℓ₁_max
    /// normalisation. The initial uniform-amplitude block saturates the per-entry
    /// |ρ_ab|² ≤ ρ_aa·ρ_bb bound, so ℓ₁ = ℓ₁_max.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void InitialPsiBlock_Saturates_AtUnity(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var trajectory = C2BlockCpsiTrajectory.Build(
            block, q: 1.0, timeGrid: new[] { 0.0 });
        double psi0 = trajectory.PsiBlockTrajectory[0];
        Assert.True(Math.Abs(psi0 - 1.0) < 1e-10,
            $"N={N}: Ψ_block(0) = {psi0:F12}, expected 1.0 (Cauchy-Schwarz saturating).");
    }

    /// <summary>Trajectory under L_block(Q) at fixed Q: CΨ_block decreases monotonically
    /// from 1/4 at t=0. Z-dephasing is contractive on the (popcount-1, popcount-2)
    /// coherence block (HD ∈ {1, 3} channels both decay at rate 2γ₀·HD).</summary>
    [Theory]
    [InlineData(5, 1.0)]
    [InlineData(5, 2.0)]
    [InlineData(7, 1.5)]
    public void Trajectory_DecreasesMonotonically_FromQuarter(int N, double Q)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var times = Enumerable.Range(0, 20).Select(i => i * 0.5).ToArray();  // t = 0..9.5 in 0.5 steps
        var trajectory = C2BlockCpsiTrajectory.Build(block, Q, times);

        Assert.True(Math.Abs(trajectory.CPsiInitial - 0.25) < 1e-10);
        for (int i = 1; i < trajectory.CPsiBlockTrajectory.Count; i++)
        {
            double prev = trajectory.CPsiBlockTrajectory[i - 1];
            double curr = trajectory.CPsiBlockTrajectory[i];
            Assert.True(curr <= prev + 1e-10,
                $"N={N}, Q={Q}: non-monotone at i={i}: prev={prev:F6}, curr={curr:F6}");
        }
    }

    /// <summary>Closed-form C_block(t) under pure dephasing (Q has no effect at this level
    /// because the Hamiltonian rotates between HD=1 and HD=3 channels in a way that
    /// preserves Σ|ρ_ab|² for our channel-uniform initial state). Matches:
    /// <code>
    /// C_block(t) = (1/(2N))·exp(−4γt) + ((N−2)/(4N))·exp(−12γt)
    /// </code>
    /// where the two exponentials track HD=1 (N(N−1) entries decaying at rate 2γ·HD=2γ
    /// applied to |ρ|² so 4γ) and HD=3 (rate 12γ on |ρ|²) channels.</summary>
    [Theory]
    [InlineData(5, 1.0, 0.5)]
    [InlineData(5, 1.0, 2.0)]
    [InlineData(5, 3.0, 0.5)]   // Q-independent: same value as Q=1.0
    [InlineData(7, 1.5, 1.0)]
    [InlineData(8, 2.0, 1.5)]
    public void CBlock_MatchesClosedForm_PureDephasing(int N, double Q, double t)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var trajectory = C2BlockCpsiTrajectory.Build(block, Q, new[] { t });
        double measured = trajectory.CBlockTrajectory[0];

        // Closed-form: (1/(2N))·exp(−4γt) + ((N−2)/(4N))·exp(−12γt)
        double γ = 0.05;
        double predicted = (1.0 / (2.0 * N)) * Math.Exp(-4.0 * γ * t)
                         + ((N - 2.0) / (4.0 * N)) * Math.Exp(-12.0 * γ * t);

        Assert.True(Math.Abs(measured - predicted) < 1e-3,
            $"N={N}, Q={Q}, t={t}: measured C_block={measured:F6}, predicted={predicted:F6}, " +
            $"Δ={Math.Abs(measured - predicted):E2}");
    }

    /// <summary>Q-dependence at uniform J is small but non-zero. The channel-uniform
    /// Dicke probe predominantly lives in {|c_1⟩, |c_3⟩} subspace where Hamiltonian
    /// rotates as phases (cancelled by |·|²), so Q effect on C_block is leading-order
    /// zero. Higher-order Q-dependence does enter (different subspace projections at
    /// finite N). Empirically: Q=1 vs Q=5 at N=5 t=2 differ by ~4e-5; Q=0.5 vs Q=1
    /// differ by ~1e-7. The trajectory is **dephasing-dominated**, not strictly
    /// Q-invariant. We pin this with a 1e-3 tolerance ("Q effect is sub-percent at
    /// uniform J") to distinguish from non-uniform J where Q-dependence is O(1).</summary>
    [Theory]
    [InlineData(5, 0.5, 1.0)]
    [InlineData(5, 1.0, 5.0)]
    [InlineData(7, 0.5, 3.0)]
    public void Trajectory_IsQDephasingDominated_AtUniformJ(int N, double Q1, double Q2)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var times = new[] { 0.0, 0.5, 1.0, 2.0 };
        var t1 = C2BlockCpsiTrajectory.Build(block, Q1, times);
        var t2 = C2BlockCpsiTrajectory.Build(block, Q2, times);

        for (int i = 0; i < times.Length; i++)
        {
            double diff = Math.Abs(t1.CBlockTrajectory[i] - t2.CBlockTrajectory[i]);
            // Sub-percent (1e-3) tolerance for "dephasing-dominated under uniform J".
            // Strict Q-independence would require Hamiltonian to be exactly zero on
            // {|c_1⟩, |c_3⟩} which holds at leading order but has finite-N corrections.
            Assert.True(diff < 1e-3,
                $"N={N}: Q={Q1} vs Q={Q2} differ at t={times[i]}: " +
                $"C={t1.CBlockTrajectory[i]:F10} vs {t2.CBlockTrajectory[i]:F10}, Δ={diff:E2}");
        }
    }

    /// <summary>Q-dependence emerges at non-uniform J: perturb a single bond's coupling
    /// by Δ. The F73 sum-rule for cross-block coupling V_b[α, j] is broken (at uniform
    /// J, Σ_b V_b cross-block = 0; with one bond perturbed, residual is O(Δ)). The c=2
    /// block trajectory now depends on Q (the reference uniform value) AND on Δ AND on
    /// which bond was perturbed. We verify the dependency exists by comparing two Qs
    /// with same Δ at a non-zero time.</summary>
    [Theory]
    [InlineData(5, 0)]   // perturb endpoint bond
    [InlineData(5, 2)]   // perturb mid-chain bond
    public void Trajectory_BecomesQDependent_AtNonUniformJ(int N, int bondToVary)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var times = new[] { 0.0, 1.0, 5.0 };

        double γ = block.GammaZero;
        double Δ = 0.5 * γ;  // 50% perturbation magnitude

        // Build per-bond coupling lists at two reference Q values, with one bond varied.
        double[] coupling1 = Enumerable.Repeat(1.0 * γ, block.NumBonds).ToArray();
        coupling1[bondToVary] = (1.0 + 0.5) * γ;
        double[] coupling2 = Enumerable.Repeat(2.0 * γ, block.NumBonds).ToArray();
        coupling2[bondToVary] = (2.0 + 0.5) * γ;

        var trajQ1 = C2BlockCpsiTrajectory.BuildPerBond(block, 1.0, coupling1, times);
        var trajQ2 = C2BlockCpsiTrajectory.BuildPerBond(block, 2.0, coupling2, times);

        // At t=0, both still match by Cauchy-Schwarz initial structure (Δ doesn't affect t=0).
        Assert.True(Math.Abs(trajQ1.CBlockTrajectory[0] - trajQ2.CBlockTrajectory[0]) < 1e-9);

        // At t > 0, dependency must emerge: difference between trajectories non-zero.
        double diffAtT5 = Math.Abs(trajQ1.CBlockTrajectory[2] - trajQ2.CBlockTrajectory[2]);
        Assert.True(diffAtT5 > 1e-6,
            $"N={N}, bond={bondToVary}: trajectories should differ at t=5 with non-uniform J; " +
            $"got Δ={diffAtT5:E2}, expected > 1e-6");
    }
}
