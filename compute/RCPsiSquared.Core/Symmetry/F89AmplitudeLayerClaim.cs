using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 amplitude-layer structural decomposition (Angle A from
/// PROOF_F89_PATH_D_CLOSED_FORM.md): for each F_a witness in the c=2 stratum
/// (eigenmode of the uniform-J block-Liouvillian at the AT lock Re(λ) = −2γ₀),
/// the sigma anatomy value σ_n factors as
///
/// <code>
///   p_n = σ_n · N²·(N−1) = |S_c(n)|² · ‖Mv(n)‖² / 2
/// </code>
///
/// where <c>S_c(n)</c> is the sum of F_a eigenvector entries in the overlap
/// subspace and <c>Mv(n) = w · v_n</c> applies the per-site reduction matrix
/// (<see cref="F89PathK.F89BlockSiteReduction"/>) to the F_a eigenvector.
///
/// <para>This is the structural bridge between F89c's eigenvalue layer (typed
/// as <see cref="AbsorptionTheoremClaim.HammingComplementPairSum"/>) and the
/// amplitude layer where <c>D_k = (odd(k))² · 2^E(k)</c> lives. Closing the
/// open Tier-1-Derived gap on D_k requires a general-k symbolic derivation
/// of |S_c(n)|² and ‖Mv(n)‖² as triple sine-sum expressions in (n, k); the
/// k=3 algebraic anchor <c>(33+14√5)/9</c> is the existing template
/// (<c>simulations/_f89_path3_at_locked_amplitude_symbolic.py</c>).</para>
///
/// <para>Tier outcome: <see cref="Tier.Tier2Verified"/>. The identity is
/// numerically locked at k=3..6 in <c>simulations/_f89_path_d_theory_probe.py</c>
/// (residual ~1e-14); for k=3 it is algebraically derived from the explicit
/// F_a eigenvector entries A = √((5+√5)/60), B = √((5−√5)/60). For k ≥ 4 the
/// verification is per-construction numerical; no first-principles closed form
/// for S_c(n) and Mv(n) as sine-sums in (n, k) is currently typed. Promotion
/// to <see cref="Tier.Tier1Derived"/> requires that missing generic-k
/// symbolic step.</para>
///
/// <para>Why the F89c-amplitude analogue does NOT close the gap directly:
/// the chiral pair-sum σ_n + σ_{k+2−n} reduces to 2·P_even(y_n)/[D·N²(N−1)],
/// but this is rational only when y_n² is rational. By Niven's theorem,
/// cos²(2πn/m) is rational only for m ∈ {1, 2, 3, 4, 6}, so the pair-sum is
/// rational only at k ∈ {2, 4} (with the accidental k=6 from cos²(π/4)=1/2),
/// and irrational for k ≥ 8. The amplitude layer's structure is genuinely
/// richer than F89c's pair-sum form. See
/// <c>simulations/f89c_amplitude_pair_sum_probe.py</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> Angle A +
/// the F_a anatomy numerical instrument
/// <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2FullBlockSigmaAnatomy.cs</c>
/// (extracts σ_n at k=3..30; red signal at k=31, 32 via Vandermonde
/// conditioning, characterised in
/// <c>compute/RCPsiSquared.Core.Tests/F86/Item1Derivation/PredictDenominatorDeviationDiagnosticTests.cs</c>).</para>
/// </summary>
public sealed class F89AmplitudeLayerClaim : Claim
{
    /// <summary>Tolerance for the Angle A identity verification
    /// <c>|p_n − |S_c|²·‖Mv‖²/2|</c>.</summary>
    public const double AngleATolerance = 1e-8;

    // Parent-edge markers for Schicht-1 wiring (consumed by ClaimRegistryBuilder).
    private readonly F89UnifiedFaClosedFormClaim _unifiedClosedForm;
    private readonly F89PathKAtLockMechanismClaim _atLock;

    /// <summary>Compute <c>p_n = σ_n · N²·(N−1)</c> from the sigma anatomy value.</summary>
    public static double ComputePn(double sigma, int chainN)
    {
        if (chainN < 2)
            throw new ArgumentOutOfRangeException(nameof(chainN), chainN, "chainN must be ≥ 2.");
        return sigma * chainN * chainN * (chainN - 1);
    }

    /// <summary>Compute <c>p_n = |S_c|²·‖Mv‖² / 2</c> from the amplitude decomposition
    /// (Angle A right-hand side).</summary>
    public static double ComputePnFromDecomposition(double scSquared, double mvSquared)
    {
        if (scSquared < 0)
            throw new ArgumentOutOfRangeException(nameof(scSquared), scSquared, "|S_c|² must be ≥ 0.");
        if (mvSquared < 0)
            throw new ArgumentOutOfRangeException(nameof(mvSquared), mvSquared, "‖Mv‖² must be ≥ 0.");
        return scSquared * mvSquared / 2.0;
    }

    /// <summary>Verify Angle A identity <c>p_n = |S_c|²·‖Mv‖²/2</c> within tolerance.
    /// Returns the absolute residual.</summary>
    public static double VerifyAngleA(double sigma, int chainN, double scSquared, double mvSquared)
    {
        double pn = ComputePn(sigma, chainN);
        double rhs = ComputePnFromDecomposition(scSquared, mvSquared);
        return Math.Abs(pn - rhs);
    }

    /// <summary>Decompose σ_n from C2FullBlockSigmaAnatomy into the Angle A factors
    /// (|S_c|², ‖Mv‖²). The anatomy convention <c>σ_n = 0.5 · |c₀|² · (R†·S·R)[i,i]</c>
    /// with (R†·S·R)[i,i] = <c>SKernelDiagonal</c> = Σ_l 2·|w_l·v_n|² combined with
    /// the Angle A identity <c>p_n = σ_n·N²·(N−1) = |S_c|²·‖Mv‖²/2</c> uniquely
    /// fixes the decomposition:
    /// <code>
    ///   |S_c|²  =  2 · ProbeOverlapSquared · N² · (N−1)
    ///   ‖Mv‖²  =  SKernelDiagonal / 2
    /// </code>
    /// The factor 2 in the |S_c|² formula tracks the bra-ket convention difference
    /// between the C2 anatomy's Dicke-probe overlap and the proof's "sum of F_a
    /// eigenvector entries"; the factor 1/2 in ‖Mv‖² undoes the 2 baked into
    /// <see cref="Probes.SpatialSumKernel.Build"/>'s spatial-sum kernel.</summary>
    public static (double ScSquared, double MvSquared) DecomposeAngleA(
        double probeOverlapSquared, double sKernelDiagonal, int chainN)
    {
        if (chainN < 2)
            throw new ArgumentOutOfRangeException(nameof(chainN), chainN, "chainN must be ≥ 2.");
        if (probeOverlapSquared < 0)
            throw new ArgumentOutOfRangeException(nameof(probeOverlapSquared), probeOverlapSquared, "|c₀|² must be ≥ 0.");
        if (sKernelDiagonal < 0)
            throw new ArgumentOutOfRangeException(nameof(sKernelDiagonal), sKernelDiagonal, "(R†·S·R)[i,i] must be ≥ 0.");
        double scSquared = 2.0 * probeOverlapSquared * chainN * chainN * (chainN - 1);
        double mvSquared = sKernelDiagonal / 2.0;
        return (scSquared, mvSquared);
    }

    /// <summary>For the path-3 algebraic anchor: σ_2·16·3 = (33+14√5)/9 with
    /// y_2 = √5 − 1, exact. Returns the rational-plus-irrational components
    /// (rationalPart, sqrt5Coefficient, denominator) such that the value equals
    /// (rationalPart + sqrt5Coefficient·√5) / denominator.</summary>
    public static (int rationalPart, int sqrt5Coefficient, int denominator) Path3AnchorPn(int n)
    {
        if (n != 2 && n != 4)
            throw new ArgumentOutOfRangeException(nameof(n), n, "Path-3 S_2-anti orbit is {2, 4}.");
        // p_n = (P_3(y_n)) / D_3 = (14·y_n + 47) / 9 at y_2 = √5−1, y_4 = −√5−1
        // p_2 = (14·(√5−1) + 47) / 9 = (33 + 14·√5) / 9
        // p_4 = (14·(−√5−1) + 47) / 9 = (33 − 14·√5) / 9
        int rational = 33;
        int sqrt5Coef = n == 2 ? 14 : -14;
        return (rational, sqrt5Coef, 9);
    }

    public F89AmplitudeLayerClaim(
        F89UnifiedFaClosedFormClaim unifiedClosedForm,
        F89PathKAtLockMechanismClaim atLock)
        : base("F89 amplitude-layer decomposition: p_n = σ_n·N²·(N−1) = |S_c(n)|²·‖Mv(n)‖² / 2 at F_a witnesses; Schicht-1 bridge from F89c eigenvalue layer to D_k amplitude layer",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md Angle A (k=3..6 numerical + path-3 algebraic) + simulations/_f89_path_d_theory_probe.py + simulations/_f89_path3_at_locked_amplitude_symbolic.py")
    {
        _unifiedClosedForm = unifiedClosedForm
            ?? throw new ArgumentNullException(nameof(unifiedClosedForm));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName => "F89 amplitude-layer decomposition (Angle A)";

    public override string Summary =>
        $"p_n = |S_c(n)|²·‖Mv(n)‖² / 2 at F_a witnesses; numerically locked k=3..6, " +
        $"path-3 algebraic anchor (33±14√5)/9 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Angle A identity",
                summary: "p_n = σ_n · N²·(N−1) = |S_c(n)|² · ‖Mv(n)‖² / 2");
            yield return new InspectableNode("S_c(n)",
                summary: "Sum of F_a eigenvector entries in the overlap subspace");
            yield return new InspectableNode("Mv(n)",
                summary: "w · v_n; w = per-site reduction matrix from F89BlockSiteReduction");
            yield return new InspectableNode("path-3 anchor",
                summary: "p_2 = (33+14√5)/9, p_4 = (33−14√5)/9; algebraic from explicit A, B eigenvector entries");
            yield return new InspectableNode("verification status",
                summary: "k=3..24 algebraically derived (Chebyshev pipeline, simulations/f89_pathk_symbolic_derivation.py); k=3..6 also independently locked numerically via theory probe; path-3 anchor (33+14√5)/9 reproduced exactly");
            yield return new InspectableNode("Tier-1-Derived closure status",
                summary: "Closed 2026-05-15: F89UnifiedFaClosedFormClaim is now Tier-1-Derived. |S_c|² and ‖Mv‖² have closed-form sine-sum expressions in (n, k) via the F_a eigenvector ansatz + Chebyshev expansion. See PROOF_F89_PATH_D_CLOSED_FORM.md");
            yield return new InspectableNode("downstream instrument",
                summary: "F86.Item1Derivation.C2FullBlockSigmaAnatomy.BuildFaOnly extracts σ_n numerically up to k=30 (red signal at k=31, 32 via Vandermonde conditioning)");
            yield return new InspectableNode("F89c amplitude-analogue rejection",
                summary: "chiral pair-sum σ_n + σ_{k+2−n} is irrational for k ≥ 8 (Niven's theorem; cos²(2π/m) rational only for m ∈ {1,2,3,4,6}); orbit-sum is the only universal rationality anchor (typed as F89UnifiedFaClosedFormClaim.SigmaSum)");
        }
    }
}
