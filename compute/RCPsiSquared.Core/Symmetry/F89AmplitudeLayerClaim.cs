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
/// to the F_a eigenvector.
///
/// <para>This is the structural bridge between F89c's eigenvalue layer (typed
/// as <see cref="AbsorptionTheoremClaim.HammingComplementPairSum"/>) and the
/// amplitude layer where <c>D_k = (odd(k))² · 2^E(k)</c> lives. From the F_a
/// eigenvector ansatz <c>v_n[(i, (j, l))] = sign(i − other) · ψ_n(other) / √k</c>
/// with <c>ψ_n(j) = √(2/(k+2)) · sin(πn(j+1)/(k+2))</c>, both factors reduce
/// to explicit sine-sum closed forms in (n, k):</para>
///
/// <code>
///   |S_c(n)|²  =  (2 / (k·m)) · (Σ_{o=0..k} sin((o+1)·πn/m) · (k − 2o))²
///   ‖Mv(n)‖²  =  (2 / (k·m)) · Σ_{l=0..k} sin²((l+1)·πn/m) · (k − 2l)²
///   m = k + 2
/// </code>
///
/// <para>Tier outcome: <see cref="Tier.Tier1Derived"/> (promoted 2026-05-16).
/// The general-k closed forms above are exposed as
/// <see cref="ComputeScSquaredClosedForm"/> and <see cref="ComputeMvSquaredClosedForm"/>,
/// each independently derived from the F_a eigenvector ansatz (no Vandermonde
/// fit, no per-k tabulation). Angle A then matches the proof identity
/// p_n = P_k(y_n) / D_k bit-exactly at every (k, n) sampled in
/// <c>F89AmplitudeLayerClaimTests</c>. The k=3 algebraic anchor
/// <c>(33±14√5)/9</c> is the path-3 instance of the same formula
/// (<c>simulations/_f89_path3_at_locked_amplitude_symbolic.py</c>); promotion
/// rests on the same Chebyshev-expansion + orbit-polynomial-reduction pipeline
/// that closed <see cref="F89UnifiedFaClosedFormClaim"/> on 2026-05-15
/// (<see cref="F89PathPolynomialPipeline"/>).</para>
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

    /// <summary>Closed-form |S_c(n)|² for path-k at Bloch index n, evaluated as
    /// (2/(k·m)) · (Σ_{o=0..k} sin((o+1)·πn/m)·(k−2o))² with m = k+2. Direct
    /// reduction from the F_a eigenvector ansatz; no Vandermonde fit, no
    /// per-k tabulation. Valid for any k ≥ 3 and any n in the S_2-anti
    /// Bloch orbit {2, 4, ..., 2·⌊(k+1)/2⌋}.</summary>
    public static double ComputeScSquaredClosedForm(int k, int n)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 3.");
        int m = k + 2;
        double thetaN = Math.PI * n / m;
        double sum = 0;
        for (int o = 0; o <= k; o++)
            sum += Math.Sin((o + 1) * thetaN) * (k - 2 * o);
        return 2.0 * sum * sum / (k * m);
    }

    /// <summary>Closed-form ‖Mv(n)‖² for path-k at Bloch index n, evaluated as
    /// (2/(k·m)) · Σ_{l=0..k} sin²((l+1)·πn/m)·(k−2l)² with m = k+2. Direct
    /// reduction from the F_a eigenvector ansatz via per-site reduction matrix
    /// w; no Vandermonde fit, no per-k tabulation. Valid for any k ≥ 3 and any
    /// n in the S_2-anti Bloch orbit.</summary>
    public static double ComputeMvSquaredClosedForm(int k, int n)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 3.");
        int m = k + 2;
        double thetaN = Math.PI * n / m;
        double sum = 0;
        for (int l = 0; l <= k; l++)
        {
            double s = Math.Sin((l + 1) * thetaN);
            int weight = k - 2 * l;
            sum += s * s * weight * weight;
        }
        return 2.0 * sum / (k * m);
    }

    public F89AmplitudeLayerClaim(
        F89UnifiedFaClosedFormClaim unifiedClosedForm,
        F89PathKAtLockMechanismClaim atLock)
        : base("F89 amplitude-layer decomposition: p_n = σ_n·N²·(N−1) = |S_c(n)|²·‖Mv(n)‖² / 2 at F_a witnesses; closed-form |S_c|² and ‖Mv‖² as sine sums in (n, k) via the F_a eigenvector ansatz; Schicht-1 bridge from F89c eigenvalue layer to D_k amplitude layer",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md Angle A; |S_c|² and ‖Mv‖² as explicit sine sums (ComputeScSquaredClosedForm, ComputeMvSquaredClosedForm); bit-exact verified against F89UnifiedFaClosedFormClaim.Sigma in F89AmplitudeLayerClaimTests; path-3 algebraic anchor (33+14√5)/9; same Chebyshev-expansion + orbit-polynomial-reduction pipeline that closed F89UnifiedFaClosedFormClaim 2026-05-15 (F89PathPolynomialPipeline)")
    {
        _unifiedClosedForm = unifiedClosedForm
            ?? throw new ArgumentNullException(nameof(unifiedClosedForm));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName => "F89 amplitude-layer decomposition (Angle A)";

    public override string Summary =>
        $"p_n = |S_c(n)|²·‖Mv(n)‖² / 2 at F_a witnesses; both factors closed-form sine sums in (n, k) " +
        $"via the F_a eigenvector ansatz; path-3 algebraic anchor (33±14√5)/9 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Angle A identity",
                summary: "p_n = σ_n · N²·(N−1) = |S_c(n)|² · ‖Mv(n)‖² / 2");
            yield return new InspectableNode("|S_c(n)|² closed form",
                summary: "(2/(k·m)) · (Σ_{o=0..k} sin((o+1)·πn/m)·(k−2o))², m = k+2; ComputeScSquaredClosedForm");
            yield return new InspectableNode("‖Mv(n)‖² closed form",
                summary: "(2/(k·m)) · Σ_{l=0..k} sin²((l+1)·πn/m)·(k−2l)², m = k+2; ComputeMvSquaredClosedForm");
            yield return new InspectableNode("path-3 anchor",
                summary: "p_2 = (33+14√5)/9, p_4 = (33−14√5)/9; the closed forms above reproduce these algebraically");
            yield return new InspectableNode("verification status",
                summary: "Bit-exact match against F89UnifiedFaClosedFormClaim.Sigma at every (k, n) sampled in F89AmplitudeLayerClaimTests; closed-form factors independently evaluable at arbitrary k ≥ 3");
            yield return new InspectableNode("Tier-1-Derived closure",
                summary: "Promoted 2026-05-16. Both |S_c|² and ‖Mv|² individually closed-form sine sums in (n, k); same Chebyshev pipeline that closed F89UnifiedFaClosedFormClaim 2026-05-15. See PROOF_F89_PATH_D_CLOSED_FORM.md Angle A.");
            yield return new InspectableNode("downstream instrument (legacy numerical)",
                summary: "F86.Item1Derivation.C2FullBlockSigmaAnatomy.BuildFaOnly extracts σ_n numerically up to k=30; superseded as authoritative source by the analytical closed forms above and by F89PathPolynomialPipeline for D_k itself");
            yield return new InspectableNode("F89c amplitude-analogue rejection",
                summary: "chiral pair-sum σ_n + σ_{k+2−n} is irrational for k ≥ 8 (Niven's theorem; cos²(2π/m) rational only for m ∈ {1,2,3,4,6}); orbit-sum is the only universal rationality anchor (typed as F89UnifiedFaClosedFormClaim.SigmaSum)");
        }
    }
}
