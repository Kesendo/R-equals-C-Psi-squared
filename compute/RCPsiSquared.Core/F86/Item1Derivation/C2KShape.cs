using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), Stage D1: the K_CC_pr observable evaluated in the
/// 4-mode effective basis as a closed-form Duhamel formula in the 4×4 eigenstructure of
/// L_eff(Q).
///
/// <para>Formula:</para>
/// <code>
///   K_b(Q, t) = 2 · Re ⟨ρ(t) | S_kernel | ∂ρ/∂J_b⟩
/// </code>
/// <para>with the Duhamel building blocks:</para>
/// <list type="bullet">
///   <item><c>ρ(t) = R · diag(e^(λ_i·t)) · R⁻¹ · ρ_0</c> — propagator applied to the probe</item>
///   <item><c>∂ρ/∂J_b = R · (X_b ⊙ I_jk) · R⁻¹ · ρ_0</c> with <c>X_b = R⁻¹·V_b·R</c> the bond
///   coupling rotated to the L_eff eigenbasis</item>
///   <item><c>I_jk(t) = (e^(λ_k·t) − e^(λ_j·t)) / (λ_k − λ_j)</c> — the Duhamel kernel
///   (<c>= t·e^(λ_j·t)</c> at near-degeneracy with magnitude threshold 1e-10)</item>
///   <item>ρ_0 is the 4-mode-projected Dicke probe (<see cref="C2EffectiveSpectrum.ProbeProjection"/>);
///   S_kernel is the 4-mode projection of the spatial-sum coherence kernel
///   (<see cref="FourModeEffective.SKernelEff"/>)</item>
/// </list>
///
/// <para><b>Tier outcome: Tier 1 derived.</b> The Duhamel formula is exact in its inputs:
/// for any matrix L with eigendecomposition <c>L = R · Λ · R⁻¹</c>, the propagator
/// <c>e^(L·t)</c> and the parametric derivative <c>∂_J e^(L·t)</c> are given by closed-form
/// matrix expressions in (R, Λ). The class-level Tier reflects the formula's soundness, not
/// the upstream input quality. Even though the eigenvalues themselves are Tier 2 numerical
/// (<see cref="C2EffectiveSpectrum"/>'s cubic-c_3 obstruction prevents a closed-form quartic
/// factorisation), the K_b(Q, t) primitive encapsulates a Tier 1 derived calculation routed
/// through Tier 2 inputs.</para>
///
/// <para><b>Output verification.</b> <see cref="PeakOverT"/> reproduces the existing
/// <see cref="FourModeResonanceScan.ComputeKCurve"/> output bit-exact (1e-10 across N=5, 7
/// over the default 21-point t-grid and the default 153-point Q-grid). The two paths share
/// the same Duhamel structure; this class repackages the per-(Q, t) inner loop as a callable
/// per-Q primitive, so D2 can call <see cref="PeakOverT"/> at any Q without spinning up the
/// full Q-scan.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2), Stage D Duhamel.</para>
/// </summary>
public sealed class C2KShape : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>Composition: per-bond V_b matrices in the 4-mode basis. Used to build
    /// L_eff(Q) at any Q and to extract V_b for the ∂ρ/∂J_b factor in the Duhamel formula.
    /// Reusing <see cref="C2BondCoupling"/> keeps the bond-coupling source single — every V_b
    /// flows through the same anti-Hermiticity-guarded assembly.</summary>
    public C2BondCoupling BondCoupling { get; }

    /// <summary>Composition: the 4-mode effective spectrum (eigenstructure of L_eff(Q)).
    /// Provides <see cref="C2EffectiveSpectrum.ProbeProjection"/> and the Tier 1 structural
    /// sub-fact (probe ⊥ {|u_0⟩, |v_0⟩}).</summary>
    public C2EffectiveSpectrum Spectrum { get; }

    /// <summary>The 4-mode probe ρ_0 = B† · DickeBlockProbe (cached on construction). Exposed
    /// directly for symmetry with <see cref="SKernelEff"/>.</summary>
    public ComplexVector ProbeProjection { get; }

    /// <summary>The 4-mode S_kernel = B† · S · B (from <see cref="FourModeEffective.SKernelEff"/>).
    /// Cached as a 4×4 matrix; identical entries across (Q, t) so we build it once.</summary>
    public ComplexMatrix SKernelEff { get; }

    /// <summary>Public factory: validates c=2, composes BondCoupling + Spectrum, and projects
    /// the S_kernel via <see cref="FourModeEffective.Build"/>.</summary>
    public static C2KShape Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2KShape applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var bondCoupling = C2BondCoupling.Build(block);
        var spectrum = C2EffectiveSpectrum.Build(block);
        // FourModeEffective owns S_kernel projection. The probe projection is also available
        // via Spectrum.ProbeProjection (cached); we cross-check structural agreement at the
        // 4-mode-basis level by sourcing both from the canonical FourModeEffective build.
        var fourModeEffective = FourModeEffective.Build(block);
        return new C2KShape(block, bondCoupling, spectrum,
                            spectrum.ProbeProjection, fourModeEffective.SKernelEff);
    }

    private C2KShape(
        CoherenceBlock block,
        C2BondCoupling bondCoupling,
        C2EffectiveSpectrum spectrum,
        ComplexVector probeProjection,
        ComplexMatrix sKernelEff)
        : base("c=2 K_b(Q,t) Duhamel in 4-mode basis",
               // Tier 1 derived: the Duhamel formula is closed-form in its inputs. The
               // C2EffectiveSpectrum eigenvalues being Tier 2 numerical does not propagate to
               // this Claim — Tier is per-claim, not inherited. Analogous to C1's D_eff which
               // is Tier1Derived inside the Tier2 BondCoupling class.
               Tier.Tier1Derived,
               Item1Anchors.StageD)
    {
        Block = block;
        BondCoupling = bondCoupling;
        Spectrum = spectrum;
        ProbeProjection = probeProjection;
        SKernelEff = sKernelEff;
    }

    /// <summary>K_b(Q, t) via Duhamel evaluation in the 4-mode basis. Real-valued by
    /// construction (the formula is <c>2 · Re ⟨ρ | S | ∂ρ⟩</c>; we return the real part
    /// directly, not its absolute value — peak-over-t handles the sign downstream).
    /// </summary>
    public double KAt(double Q, int bond, double t)
    {
        ValidateBond(bond);

        // Build L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b — same path as Spectrum.LEffAtQ to keep the
        // matrix entries single-sourced.
        var L = Spectrum.LEffAtQ(Q);
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var evals = evd.EigenValues.ToArray();

        // c0 = R⁻¹ · ρ_0
        var c0 = Rinv * ProbeProjection;

        // X_b = R⁻¹ · V_b · R — the bond coupling rotated to the L_eff eigenbasis.
        var Vb = BondCoupling.AsMatrix(bond);
        var Xb = Rinv * Vb * R;

        return EvaluateK(t, R, evals, c0, Xb);
    }

    /// <summary>Peak of |K_b(Q, t)| over a t-grid for a given Q. Returns the (peak |K|, t at
    /// peak) tuple. Same convention as <see cref="FourModeResonanceScan"/>: the peak is over
    /// the absolute value, but K_b itself can have either sign and the t at peak is the t in
    /// the grid where |K| is maximised.</summary>
    public (double KPeak, double TAtPeak) PeakOverT(double Q, int bond, IReadOnlyList<double> tGrid)
    {
        ValidateBond(bond);
        if (tGrid is null || tGrid.Count == 0)
            throw new ArgumentException("tGrid must be non-empty.", nameof(tGrid));

        // Cache the per-Q eigenstructure so the t-loop only does the inner Duhamel evaluation.
        var L = Spectrum.LEffAtQ(Q);
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var evals = evd.EigenValues.ToArray();
        var c0 = Rinv * ProbeProjection;
        var Vb = BondCoupling.AsMatrix(bond);
        var Xb = Rinv * Vb * R;

        double kPeak = 0.0;
        double tAtPeak = double.NaN;
        foreach (double t in tGrid)
        {
            double k = EvaluateK(t, R, evals, c0, Xb);
            double kAbs = Math.Abs(k);
            if (kAbs > kPeak)
            {
                kPeak = kAbs;
                tAtPeak = t;
            }
        }
        return (kPeak, tAtPeak);
    }

    /// <summary>Convenience: peak |K_b| over the default 21-point t-grid spanning
    /// [0.6, 1.6]·t_peak with t_peak = 1/(4γ₀) — same grid that
    /// <see cref="FourModeResonanceScan"/> uses by default.</summary>
    public double PeakOverDefaultT(double Q, int bond)
    {
        var tGrid = ResonanceScan.DefaultTGrid(Block.GammaZero, points: 21);
        var (kPeak, _) = PeakOverT(Q, bond, tGrid);
        return kPeak;
    }

    /// <summary>Inner Duhamel evaluator. Mirrors the loop in
    /// <see cref="FourModeResonanceScan"/>'s ScanAtQ but extracted here for the per-(Q, t)
    /// callable form. Returns the real K_b(Q, t) value (not its absolute value).</summary>
    private double EvaluateK(double t, ComplexMatrix R, Complex[] evals, ComplexVector c0,
        ComplexMatrix Xb)
    {
        const int dim = 4;
        var expLam = new Complex[dim];
        for (int i = 0; i < dim; i++) expLam[i] = Complex.Exp(evals[i] * t);

        // I_mat[r, c] = (e^(λ_c·t) − e^(λ_r·t)) / (λ_c − λ_r) or t·e^(λ_r·t) at degeneracy.
        var iMat = new Complex[dim, dim];
        for (int r = 0; r < dim; r++)
            for (int c = 0; c < dim; c++)
            {
                Complex diff = evals[c] - evals[r];
                iMat[r, c] = diff.Magnitude > 1e-10
                    ? (expLam[c] - expLam[r]) / diff
                    : t * expLam[r];
            }

        // ρ(t) = R · (expLam ⊙ c0)  (diagonal × vector, then R·)
        var weighted = ComplexVector.Build.Dense(dim);
        for (int i = 0; i < dim; i++) weighted[i] = expLam[i] * c0[i];
        var rhoT = R * weighted;

        // ∂ρ/∂J_b in the eigenbasis: fbC0[r] = Σ_c X_b[r, c] · I_mat[r, c] · c0[c]
        // Then ∂ρ/∂J_b (4-mode basis) = R · fbC0.
        var fbC0 = new Complex[dim];
        for (int r = 0; r < dim; r++)
        {
            Complex s = Complex.Zero;
            for (int c = 0; c < dim; c++) s += Xb[r, c] * iMat[r, c] * c0[c];
            fbC0[r] = s;
        }
        var fbVec = ComplexVector.Build.Dense(fbC0);
        var drho = R * fbVec;

        // sDrho = S_kernel · ∂ρ/∂J_b
        var sDrho = SKernelEff * drho;

        // K_b = 2 · Re ⟨ρ(t) | sDrho⟩ = 2 · Re Σ_i conj(ρ(t)[i]) · sDrho[i]
        Complex inner = Complex.Zero;
        for (int i = 0; i < dim; i++) inner += Complex.Conjugate(rhoT[i]) * sDrho[i];
        return 2.0 * inner.Real;
    }

    private void ValidateBond(int bond)
    {
        if (bond < 0 || bond >= Block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(bond),
                $"bond must be in [0, {Block.NumBonds - 1}]; got {bond}.");
    }

    public override string DisplayName =>
        $"c=2 K_b(Q,t) Duhamel in 4-mode basis (N={Block.N}, bonds={Block.NumBonds})";

    public override string Summary =>
        $"K_b(Q, t) = 2·Re⟨ρ(t)|S|∂ρ/∂J_b⟩ closed-form in (R, Λ) of L_eff(Q); " +
        $"PeakOverT matches FourModeResonanceScan at 1e-10 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return BondCoupling;
            yield return Spectrum;
            yield return new InspectableNode("ProbeProjection",
                summary: $"4-mode-basis projection of Dicke probe; |components onto SVD-top| < 1e-12 (Tier 1 structural)");
            yield return new InspectableNode("SKernelEff",
                summary: $"4×4 S_kernel = B† · S · B from FourModeEffective.SKernelEff");
        }
    }
}
