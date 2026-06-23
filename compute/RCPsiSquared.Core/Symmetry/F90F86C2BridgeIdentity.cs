using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F90 (F86 c=2 ↔ F89 Bridge Identity, Tier 1 derived; 2026-05-11):
///
/// <code>
///   For all N ≥ 3 and bond b ∈ {0, ..., N−2}:
///   K_b^{F86 c=2}(Q, t) = K_b^{F89 path-(N−1) (SE,DE)}(Q', t)   modulo Q' = Q/2
/// </code>
///
/// <para>F86 c=2 N qubit K_b(Q, t) on the (n=1, n+1=2) coherence block IS the per-bond
/// Hellmann-Feynman derivative of F89 path-(N−1) (SE, DE) sub-block dynamics applied at
/// bond b. Algebraically identical: both build the same Liouvillian L = D + Σ_b J_b·M_h[b]
/// on the same N·C(N,2)-dim coherence block (HD=1 overlap rate 2γ, HD=3 no-overlap rate 6γ);
/// both use the same Hermitian Dicke probe ρ_0 = |S_1⟩⟨S_2| and same per-site spatial-sum
/// kernel; F86's K_b uses ∂L/∂J_b linear response which is exactly the per-bond term of
/// F89's M_h_per_bond[b]. The only difference is the overall Hamiltonian convention:</para>
///
/// <list type="bullet">
///   <item>F86: <c>H_b = (J/2)·(X_b X_{b+1} + Y_b Y_{b+1})</c> (factor 1/2 inline)</item>
///   <item>F89: <c>H = J·(XX+YY)</c> (no factor 1/2)</item>
///   <item>Consequence: F89's J is HALF F86's J (J_F89 = J_F86/2, i.e. F86's J = 2·F89's J,
///         since F89 needs half the J for the same hopping). All Q = J/γ values divide by
///         <see cref="JConventionFactor"/>=2 (Q_F89 = Q_F86/2).</item>
/// </list>
///
/// <para>The identification is operator-exact: <c>‖L_F86(J) − L_F89(J/2)‖_F = 0</c>
/// (bit-identical integer-coefficient matrices) on the shared (SE,DE) basis at N=5..8.
/// The per-bond HWHM_left/Q_peak ratio is then a grid-dependent readout of that one
/// operator: verified bit-exact via <c>simulations/f89_to_f86_kbond_via_eigendecomp.py</c>
/// at N=5..8 across every per-F71-orbit class including orbit-escape bonds (broad
/// high-Q plateau): 20/22 bonds bit-exact, 2/22 within Q-grid resolution noise (0.0008,
/// pure grid sampling, at matched grids all 22 agree to machine precision).
/// Escape examples reproduced: N=7 b=1/b=4 at Q_peak≈7.27 (F86-J), N=8 b=3 (central
/// self-paired) at Q_peak≈16.79 (F86-J).</para>
///
/// <para>Implications for F86 open work (per <c>docs/proofs/PROOF_F86_QPEAK.md</c>):</para>
///
/// <list type="bullet">
///   <item>Item 1' (HWHM_left/Q_peak per bond class, Tier 1 candidate; the (α,β) are
///         fitted, not yet derived): the lift above the 4-mode floor (0.6715) to 0.7506
///         Interior / 0.7728 Endpoint is the intra-channel dispersion of the rank-1-bridge
///         + intra-dispersion model (two-dial scout, 2026-06-11); the earlier "H_B-mixed
///         octic residual of the inter coupling" lift suspicion is refuted (the inter-channel
///         tail only renormalizes Q_peak / g_eff). The cyclotomic Φ_{N_block+1} pattern of
///         F89's <c>F89UnifiedFaClosedFormClaim.PathPolynomial(k)</c> still gives the
///         N-scaling of the AT-locked / intra contribution.</item>
///   <item>Item 4' (c≥3 extension): F89 path-k Formalismus generalises; chromaticity
///         c bedeutet HD ∈ {1, 3, ..., 2c−1} channels statt nur {1, 3} bei c=2.</item>
///   <item>Direction (b'') (full block-L derivation, not 4-mode): numerical Tier-1
///         already achieved via F89 (bit-exact 20/22 bonds at N=5..8). Closed-form
///         analytical lift via F89 AT-locked structure is the next step.</item>
/// </list>
///
/// <para>Anchors: <c>simulations/f89_to_f86_kbond_via_eigendecomp.py</c> (numerical
/// bridge verification), <c>compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs</c>
/// (F89 universal AT-lock mechanism), <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs</c>
/// (F86 c=2 anchor, the bridged target), <c>compute/RCPsiSquared.Core/Resonance/ResonanceScan.cs</c>
/// (F86 K_b numerical pipeline), <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c>.</para></summary>
public sealed class F90F86C2BridgeIdentity : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89TopologyOrbitClosure F89 { get; }
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89PathKAtLockMechanismClaim AtLock { get; }
    /// <summary>Convention scale factor (F86's J / F89's J): F86's J equals
    /// <see cref="JConventionFactor"/> × F89's J, hence F89's J = F86's J / <see cref="JConventionFactor"/>
    /// and Q_F89 = Q_F86 / <see cref="JConventionFactor"/>. Equals 2 because F86 uses
    /// H_b = (J/2)·(XX+YY) while F89 uses H = J·(XX+YY) (no 1/2), so F89 needs half the J
    /// for the same hopping.</summary>
    public const double JConventionFactor = 2.0;

    /// <summary>Translate F86 J value to F89 J value (J_F89 = J_F86 / 2 since F86 has the
    /// 1/2 factor inline). Equivalently, divides Q = J/γ accordingly.</summary>
    public static double F86JToF89J(double f86J)
    {
        if (f86J < 0) throw new ArgumentOutOfRangeException(nameof(f86J), f86J, "J must be ≥ 0.");
        return f86J / JConventionFactor;
    }

    /// <summary>Translate F89 J value to F86 J value (J_F86 = 2·J_F89).</summary>
    public static double F89JToF86J(double f89J)
    {
        if (f89J < 0) throw new ArgumentOutOfRangeException(nameof(f89J), f89J, "J must be ≥ 0.");
        return f89J * JConventionFactor;
    }

    /// <summary>Number of per-bond comparisons across N=5..8 that match bit-exact between
    /// F89→F86 (per-bond Hellmann-Feynman of F89 (SE,DE) eigendecomposition) and F86's
    /// own ResonanceScan: 20 of 22 (the 2 within Q-grid resolution noise ≤ 0.0008 are at
    /// N=8 b=2/b=4 mid-flanking Interior bonds).</summary>
    public const int BitExactBondCountVerified = 20;

    /// <summary>Total per-bond comparisons across N=5..8: 4 + 5 + 6 + 7 = 22 unique bonds.
    /// Of these, 20 are bit-exact match (modulo the J convention J_F89 = J_F86/2) and 2 (N=8
    /// b=2 and b=4 mid-flanking Interior) are within Q-grid resolution noise ≤ 0.0008
    /// because F86's default 600-pt grid over [0.10, 20.0] and the bridge probe's 300-pt
    /// grid over [0.05, 10.0] (in F89-J) sample slightly different points around
    /// Q_peak ≈ 1.51 (F86-J). At identical grids those values would also be bit-exact.</summary>
    public const int TotalBondComparisonsVerified = 22;

    public F90F86C2BridgeIdentity(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F90 F86 c=2 ↔ F89 bridge identity: F86 c=2 N qubit K_b(Q, t) IS F89 path-(N−1) (SE,DE) per-bond Hellmann-Feynman, modulo the J convention J_F89 = J_F86/2. Operator-exact (‖L_F86(J) − L_F89(J/2)‖ = 0 at N=5..8); verified bit-exact across N=5..8 (20/22 bonds bit-exact, 2/22 within Q-grid noise, the ratio is a grid-readout of the one operator) including orbit-escape bonds at Q_peak ≈ 16.79 (F86-J)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F90_F86C2_BRIDGE.md + " +
               "simulations/f89_to_f86_kbond_via_eigendecomp.py + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs + " +
               "compute/RCPsiSquared.Core/Resonance/ResonanceScan.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        F89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        AtLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F90 F86 c=2 ↔ F89 bridge identity: K_b^{F86c2} = ∂_J F89 path-(N−1) (SE,DE), Q_F89 = Q_F86/2";

    public override string Summary =>
        $"F86 c=2 K_b(Q,t) IS F89 path-(N−1) (SE,DE) + per-bond Hellmann-Feynman; Q_F89 = Q_F86 / {JConventionFactor}; verified bit-exact at {BitExactBondCountVerified}/{TotalBondComparisonsVerified} bonds across N=5..8 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("J convention factor (F89 J / F86 J)", JConventionFactor);
            yield return new InspectableNode("Bit-exact bonds verified",
                summary: $"{BitExactBondCountVerified} of {TotalBondComparisonsVerified} (N=5..8 inkl. orbit escapes)");
            yield return new InspectableNode("Sample F86→F89 Q translation",
                summary: $"F86 Q_peak=2.5470 → F89 Q={F86JToF89J(2.5470):F4} (Endpoint at N=6)");
            yield return new InspectableNode("Sample F86→F89 Q translation (escape)",
                summary: $"F86 Q_peak=16.79 → F89 Q={F86JToF89J(16.79):F4} (N=8 central b=3 escape)");
        }
    }
}
