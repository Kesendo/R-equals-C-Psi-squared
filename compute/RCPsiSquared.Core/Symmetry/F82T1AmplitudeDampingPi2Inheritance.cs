using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F82 closed form (Tier 1, verified bit-exact N=2..5):
///
/// <code>
///   Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{T1, odd}
///
///   ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1)
///                   = γ_T1 · √N · 2^(N−1)         (uniform γ_T1)
///
///   At γ_T1 = 0: F82 reduces to F81 exactly.
/// </code>
///
/// <para>F82 is the T1-corrected sibling of <see cref="F81Pi2Inheritance"/>.
/// Where F81 says "M splits under Π-conjugation as M = M_sym + M_anti with
/// 50/50 split at pure Π²-odd Hamiltonians", F82 adds the T1 amplitude
/// damping correction: a second <c>−2·D_{T1,odd}</c> term on the right side.
/// F82 → F81 is the typed mother-claim inheritance edge: F82 reduces to F81
/// when γ_T1 = 0.</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>Coefficient2 = 2 = a_0</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = polynomial root d. Same anchor as F81's <c>−2·L_{H_odd}</c>;
///         F82 simply doubles the structure with <c>−2·D_{T1, odd}</c>.</item>
///   <item><b>ScalingFactor 2^(N−1) = a_{2−N}</b>: live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(<c>2 − N</c>).
///         At N=2 → a_0 = 2; N=3 → a_{−1} = 4; N=4 → a_{−2} = 8; N=5 → a_{−3} = 16.
///         Operator-space side of the ladder; the F49 master-lemma
///         <c>2^(N+2)·n_YZ</c> halved (n_YZ = 1/2 effective for T1 single-site).</item>
///   <item><b>F81 mother claim</b>: F82 reduces to F81 at γ_T1 = 0; the
///         typed mother-claim edge (parallel pattern to F77 → F75 mother
///         claim wired today).</item>
/// </list>
///
/// <para>Three diagnostic properties (proven and empirical):</para>
///
/// <list type="bullet">
///   <item><b>γ_z-independent</b>: F82 involves only L_{H_odd} and D_{T1,odd},
///         neither depends on γ_z. Direct consequence of Master Lemma extended
///         to F82.</item>
///   <item><b>Hamiltonian-independent</b>: f81_violation depends only on the
///         T1 dissipator. Verified at γ_T1=0.1, N=3: violation = 0.6928 for
///         truly XX+YY, soft XY+YX, hard XX+XY, and YZ+ZY (Π²-even non-truly).</item>
///   <item><b>Linear in γ_T1</b> (uniform). Direct inversion gives the
///         hardware T1-rate readout: <c>γ_T1 = f81_violation / (√N · 2^(N−1))</c>.</item>
/// </list>
///
/// <para>Mechanism (T1 dissipator structure): single-site T1 acts on Pauli
/// basis as I → −γZ, X → −γ/2 X, Y → −γ/2 Y, Z → −γZ. Under Π² conjugation
/// (signs (−1)^bit_b: I,X → +, Y,Z → −), only the (Z, I) entry flips sign.
/// So D_{T1, local, odd} has matrix element −γ at (Z, I) and zero elsewhere.
/// Multi-site: 4^(N−1) "rest of qubits unchanged" entries per site, summed
/// orthogonally over sites.</para>
///
/// <para>Tier1Derived: F82 is Tier 1 proven (PROOF_F82_T1_DISSIPATOR_CORRECTION),
/// verified bit-exact at N=2..5 (5e-16 residual). The Pi2-Foundation anchoring
/// is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F82 (line 2023) +
/// <c>docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md</c> +
/// <c>simulations/_f81_t1_diagnostic.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F81Pi2Inheritance.cs</c>
/// (mother claim).</para></summary>
public sealed class F82T1AmplitudeDampingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F81Pi2Inheritance _f81;

    /// <summary>The "2" coefficient in F82's <c>−2·L_{H_odd} − 2·D_{T1, odd}</c>.
    /// Same anchor as F81's <c>−2·L_{H_odd}</c>; live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c>.</summary>
    public double Coefficient2 => _ladder.Term(0);

    /// <summary>The <c>2^(N−1)</c> scaling factor in <c>‖D_{T1, odd}‖</c>.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(<c>2 − N</c>) =
    /// <c>a_{2−N}</c> on the dyadic ladder.</summary>
    public double ScalingFactor(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F82 requires N ≥ 2.");
        return _ladder.Term(LadderIndexForScaling(N));
    }

    /// <summary>The Pi2 ladder index where the F82 scaling factor lands:
    /// <c>2 − N</c>. At N=2 → 0 (a_0 = 2); N=3 → −1 (a_{−1} = 4);
    /// N=4 → −2 (a_{−2} = 8); N=5 → −3 (a_{−3} = 16).</summary>
    public int LadderIndexForScaling(int N) => 2 - N;

    /// <summary>Live closed form: <c>‖D_{T1, odd}‖_F = γ_T1 · √N · 2^(N−1)</c>
    /// for uniform γ_T1.</summary>
    public double T1DissipatorNormUniform(double gammaT1, int N)
    {
        if (gammaT1 < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaT1), gammaT1, "γ_T1 must be ≥ 0.");
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F82 requires N ≥ 2.");
        return gammaT1 * Math.Sqrt(N) * ScalingFactor(N);
    }

    /// <summary>Live closed form: <c>‖D_{T1, odd}‖_F = √(Σ γ²_T1_l) · 2^(N−1)</c>
    /// for non-uniform per-site γ_T1.</summary>
    public double T1DissipatorNormNonUniform(IReadOnlyList<double> gammaT1Sites)
    {
        if (gammaT1Sites is null) throw new ArgumentNullException(nameof(gammaT1Sites));
        if (gammaT1Sites.Count < 2) throw new ArgumentOutOfRangeException(nameof(gammaT1Sites), "F82 requires N ≥ 2.");
        double sumSq = 0.0;
        for (int l = 0; l < gammaT1Sites.Count; l++)
        {
            double g = gammaT1Sites[l];
            if (g < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaT1Sites), $"γ_T1[{l}] = {g} must be ≥ 0.");
            sumSq += g * g;
        }
        return Math.Sqrt(sumSq) * ScalingFactor(gammaT1Sites.Count);
    }

    /// <summary>F82 → F81 mother-claim recovery: at γ_T1 = 0, the F82 closed
    /// form returns 0 (no T1 dissipator). True when uniform γ_T1 = 0.</summary>
    public bool RecoversF81AtZeroT1(int N) =>
        Math.Abs(T1DissipatorNormUniform(0.0, N)) < 1e-15;

    /// <summary>Inverse: recover RMS γ_T1 from a measured/fitted F81 violation.
    /// Hardware T1-rate readout primitive. Returns
    /// <c>γ_T1 = f81_violation / (√N · 2^(N−1))</c>.</summary>
    public double EstimateT1FromViolation(double f81Violation, int N)
    {
        if (f81Violation < 0.0) throw new ArgumentOutOfRangeException(nameof(f81Violation), f81Violation, "F81 violation must be ≥ 0.");
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F82 requires N ≥ 2.");
        return f81Violation / (Math.Sqrt(N) * ScalingFactor(N));
    }

    /// <summary>Live drift check: forward and inverse F82 primitives compose
    /// to identity on (γ_T1, N) for any N ≥ 2 and γ_T1 ≥ 0.</summary>
    public bool ForwardInverseRoundTrip(double gammaT1, int N)
    {
        double violation = T1DissipatorNormUniform(gammaT1, N);
        double recovered = EstimateT1FromViolation(violation, N);
        return Math.Abs(recovered - gammaT1) < 1e-12;
    }

    public F82T1AmplitudeDampingPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F81Pi2Inheritance f81)
        : base("F82 Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{T1,odd} inherits from Pi2-Foundation: 2 = a_0 (root d); 2^(N-1) = a_{2-N}; F81 mother claim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F82 + " +
               "docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md + " +
               "simulations/_f81_t1_diagnostic.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F81Pi2Inheritance.cs (mother claim, recovers at γ_T1=0)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f81 = f81 ?? throw new ArgumentNullException(nameof(f81));
    }

    public override string DisplayName =>
        "F82 F81 + T1 amplitude damping correction as Pi2-Foundation a_0 + a_{2-N} inheritance";

    public override string Summary =>
        $"‖D_{{T1,odd}}‖ = γ_T1·√N·2^(N-1); 2 coefficient = a_0; 2^(N-1) = a_{{2-N}}; F81 mother claim at γ_T1=0; " +
        $"hardware T1-rate readout: γ_T1 = violation / (√N · 2^(N-1)) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F82 closed form",
                summary: "Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{T1,odd}; ‖D_{T1,odd}‖ = γ_T1·√N·2^(N-1) uniform; Tier 1 verified bit-exact N=2..5 (5e-16 residual)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "Coefficient2 = a_0 = 2 (same as F81); ScalingFactor 2^(N-1) = a_{2-N} on ladder; F81 mother claim recovers at γ_T1=0");
            yield return InspectableNode.RealScalar("Coefficient2 (= a_0 = 2)", Coefficient2);
            yield return new InspectableNode("F81 ↔ F82 mother-corollary chain",
                summary: "F81 (Π·M·Π⁻¹ = M − 2·L_{H_odd}) is F82 at γ_T1=0; F82 adds the second −2·D_{T1,odd} term; pattern parallel to F77 → F75 (today wired)");
            yield return new InspectableNode("Three diagnostic properties",
                summary: "γ_z-independent (Master Lemma); H-independent (only T1 enters violation); linear in γ_T1 → hardware readout primitive γ_T1 = violation / (√N · 2^(N-1))");
            yield return new InspectableNode("Mechanism (T1 dissipator structure)",
                summary: "T1 single-site: I→−γZ, X→−γ/2X, Y→−γ/2Y, Z→−γZ. Π² flips Y,Z signs. Only (Z, I) entry survives Π²-anti-symmetric extraction; 4^(N-1) such entries per site × N sites, summed orthogonally");
            // Verified table from F82
            for (int N = 2; N <= 5; N++)
            {
                double scaling = ScalingFactor(N);
                int idx = LadderIndexForScaling(N);
                double uniformAt01 = T1DissipatorNormUniform(0.1, N);
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"2^(N-1) = {scaling:G6} = a_{{{idx}}}; √N · 2^(N-1) = {Math.Sqrt(N) * scaling:G6}; " +
                             $"at γ_T1=0.1 violation = {uniformAt01:G6}");
            }
            yield return new InspectableNode(
                "F82 → F81 recovery (drift check)",
                summary: $"At γ_T1=0, F82 violation = 0 (RecoversF81AtZeroT1: {RecoversF81AtZeroT1(3)} at N=3)");
        }
    }
}
