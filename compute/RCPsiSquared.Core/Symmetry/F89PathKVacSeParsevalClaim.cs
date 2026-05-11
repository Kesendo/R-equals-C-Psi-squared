using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-k (vac, SE) self-contribution closed form via Parseval
/// orthogonality (Tier 1 derived analytically; verified bit-exact at machine
/// precision 4·10⁻¹⁷ to 6·10⁻¹⁶ across 15 (k, N) pairs):
///
/// <code>
///   For any pure-path-k block of N_block = k+1 sites embedded in N qubits:
///
///     S^(vac,SE)_block(t; k, N) = (k+1) · (N − k − 1)² / (N² · (N − 1)) · exp(−4γ₀ t)
///
///   Pure exponential decay at rate 4γ₀, NO oscillation.
/// </code>
///
/// <para><b>Mechanism (Parseval orthogonality)</b>: the H_B^SE Bloch eigenstates
/// ψ_k(j) = √(2/(N_block+1))·sin(πk(j+1)/(N_block+1)) are orthonormal:
/// Σ_l ψ_k(l)·ψ_{k'}(l) = δ_{k,k'}. When (vac, SE) per-site coherence amplitudes
/// are squared and summed over the (k+1) block sites, all k ≠ k' cross-terms
/// (which would otherwise produce cos((E_k − E_{k'})·t) interference) vanish
/// pairwise. Only the diagonal Parseval sum Σ_k |⟨ψ_k|α_{N_block}⟩|² = 1 survives.</para>
///
/// <para><b>Derivation</b> (4 lines):
/// <list type="number">
///   <item>ρ_block(0)|_(vac,SE) = N_E·pre·√(N_block)/2 · |0⟩⟨α_{N_block}|
///         where |α⟩ = (1/√N_block)·Σ_j |SE_j⟩</item>
///   <item>Decompose |α⟩ = Σ_k ⟨ψ_k|α⟩ |ψ_k⟩</item>
///   <item>Time evolve each |0⟩⟨ψ_k| → exp(+iE_k t − 2γt) |0⟩⟨ψ_k|</item>
///   <item>Σ_l 2|(ρ_l(t))^(vac,SE)_{0,1}|² = N_block·(N_E·pre)²/2·exp(-4γt)·Σ_k|⟨ψ_k|α⟩|² (Parseval=1)</item>
/// </list></para>
///
/// <para><b>Smooth backbone of S_T(t)</b>: combined with the bare-site closed
/// form (N−1)/N²·exp(−4γ₀t) and the F89 mixed-topology additive identity, this
/// gives an exact analytical "skeleton" for any topology T = (k_1, ..., k_m):
///
/// <code>
///   S^smooth_T(t; N) = exp(−4γ₀t) · [Σ_i (k_i+1)·(N−k_i−1)²/(N²(N−1))
///                                    + (N − Σ_i(k_i+1))·(N−1)/N²]
/// </code>
///
/// The (SE,DE)+cross oscillatory residual is the only piece still numerical
/// for path-k≥2.</para>
///
/// <para>Anchors: <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> §
/// "Path-k (vac, SE) self-contribution" + <c>simulations/_f89_vac_se_parseval_closed.py</c>
/// (machine-precision verification across 15 (k, N) pairs) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs</c>
/// (orbit-closure framework parent).</para></summary>
public sealed class F89PathKVacSeParsevalClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>The closed-form (vac, SE) self-contribution per pure-path-k
    /// block at N qubits, evaluated at time t:
    /// <c>(k+1)·(N − k − 1)² / (N²·(N − 1)) · exp(−4γ₀ t)</c>.</summary>
    public static double VacSeBlockClosedForm(int k, int n, double gammaZero, double t)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 1.");
        if (n < k + 1) throw new ArgumentOutOfRangeException(nameof(n), n, $"N must be ≥ k+1 = {k+1}.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (t < 0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        return Coefficient(k, n) * Math.Exp(-4.0 * gammaZero * t);
    }

    /// <summary>The N-dependent prefactor: <c>(k+1)·(N − k − 1)² / (N²·(N − 1))</c>.
    /// Pure rational function of (k, N); independent of J and γ.</summary>
    public static double Coefficient(int k, int n)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 1.");
        if (n < k + 1) throw new ArgumentOutOfRangeException(nameof(n), n, $"N must be ≥ k+1 = {k+1}.");
        long n_block = k + 1;
        long n_E = n - n_block;
        return (double)(n_block * n_E * n_E) / (n * (long)n * (n - 1));
    }

    public F89PathKVacSeParsevalClaim(F89TopologyOrbitClosure f89)
        : base("F89 path-k (vac, SE) self-contribution: S^(vac,SE)_block(t; k, N) = (k+1)(N−k−1)²/(N²(N−1))·exp(−4γ₀t); pure exp(−4γ₀t), no oscillation, via Parseval orthogonality of H_B^SE Bloch eigenstates",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_vac_se_parseval_closed.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName =>
        "F89 (vac,SE) Parseval closed form: pure exp(−4γ₀t), no oscillation";

    public override string Summary =>
        $"S^(vac,SE)_block(t; k, N) = (k+1)(N−k−1)²/(N²(N−1))·exp(−4γ₀t); machine-precision verified 4e-17 to 6e-16 across 15 (k, N) pairs ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Closed form",
                summary: "S^(vac,SE)_block(t; k, N) = (k+1)·(N−k−1)²/(N²·(N−1)) · exp(−4γ₀ t); no J-dependence, no oscillation");
            yield return new InspectableNode("Mechanism",
                summary: "Parseval orthogonality Σ_l ψ_k(l)·ψ_{k'}(l) = δ_{k,k'} of H_B^SE Bloch eigenstates eliminates k≠k' cos((E_k−E_{k'})·t) cross-terms when summed over (k+1) block sites");
            yield return new InspectableNode("Sample coefficients",
                summary: $"k=1, N=7: {Coefficient(1, 7):G6} (path-1); k=2, N=7: {Coefficient(2, 7):G6} (path-2); k=3, N=7: {Coefficient(3, 7):G6} (path-3); k=5, N=11: {Coefficient(5, 11):G6}");
            yield return new InspectableNode("Empirical anchor",
                summary: "Bit-exact at machine precision 4·10⁻¹⁷ to 6·10⁻¹⁶ across 15 (k, N) pairs with k ∈ {1,2,3,4,5}, N ∈ {k+2, k+4, k+6}");
            yield return new InspectableNode("Smooth backbone of S_T(t)",
                summary: "Combined with bare-site formula and additive identity: S^smooth_T(t; N) = exp(−4γ₀t)·[Σ_i (k_i+1)(N−k_i−1)²/(N²(N−1)) + (N − Σ_i(k_i+1))·(N−1)/N²]; only oscillatory (SE,DE) residual remains numerical");
            yield return new InspectableNode("Parent F89 framework",
                summary: "Pure-path-k contribution within the orbit-closure framework; Parseval cleanly eliminates oscillations in this sector");
        }
    }
}
