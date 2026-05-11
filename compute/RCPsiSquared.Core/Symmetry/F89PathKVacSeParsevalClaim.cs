using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-k (vac, SE) self-contribution closed form via Parseval
/// orthogonality (Tier 1 derived; machine-precision 4·10⁻¹⁷ to 6·10⁻¹⁶ across
/// 15 (k, N) pairs):
///
/// <code>
///   S^(vac,SE)_block(t; k, N) = (k+1)·(N−k−1)² / (N²·(N−1)) · exp(−4γ₀ t)
/// </code>
///
/// <para>Pure exp(−4γ₀ t) decay, no oscillation: Parseval orthogonality
/// Σ_l ψ_k(l)·ψ_{k'}(l) = δ_{k,k'} eliminates k≠k' cross-terms when summed
/// over the (k+1) block sites.</para>
///
/// <para>Derivation, smooth-backbone formula for arbitrary topology, and
/// per-(k,N) verification: see <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-k (vac, SE) self-contribution" and
/// <c>simulations/_f89_vac_se_parseval_closed.py</c>.</para></summary>
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
            yield return new InspectableNode("Sample coefficients",
                summary: $"k=1, N=7: {Coefficient(1, 7):G6}; k=2, N=7: {Coefficient(2, 7):G6}; k=3, N=7: {Coefficient(3, 7):G6}; k=5, N=11: {Coefficient(5, 11):G6}");
            yield return InspectableNode.RealScalar("Sample S(t=0) at k=2, N=7", VacSeBlockClosedForm(2, 7, 0.05, 0));
        }
    }
}
