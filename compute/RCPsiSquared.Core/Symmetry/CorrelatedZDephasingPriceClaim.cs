using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The palindrome under a shared field (Tier 1 derived,
/// <c>experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md</c>): for Z dephasing with a real symmetric
/// Kossakowski matrix a over the Z_k (site rates on the diagonal, the fields' covariances off it)
/// and a real H commuting with X^⊗N, Π_Z·L·Π_Z⁻¹ = −L − {K, ·}, K = Σ_kl a_kl Z_l Z_k (= Σ_μ c_μ†c_μ for any
/// jumps realising a ⪰ 0), Π_Z(ρ) = ρᵀ·X^⊗N,
/// which is F1 when a is diagonal (K = σ·1, σ = tr a). The palindrome's price becomes an operator.
///
/// <para>What follows on the page: the dissipator residual M_D = −{K − σ, ·} with
/// ‖M_D‖²_F = 8·4^N·Σ_{k&lt;l} a_kl², adding in squares to the Hamiltonian's; for one covariance κ
/// common to every pair, K = σ − Nκ + 4κS_z², the generator on joint-popcount block (p, q) (ket and
/// bra popcounts, d = p − q, m = p + q − N) is the local one at the rates γ_k − κ minus 2κd², and
/// for a popcount-conserving X^⊗N-even H each block pairs with its image (q, N − p) about the pair
/// centre −σ(p, q), σ(p, q) = σ − Nκ + κ(d² + m²); the single
/// centre −σ fails for a positive semidefinite a whose covariances sum to a positive number (H
/// keeping the all-0 and all-1 states as eigenstates); for a general pattern and a popcount-conserving
/// X^⊗N-even H a block pairs wherever K takes one value on each of its two popcounts; the four edge blocks pay the local bill whatever the
/// covariances. Π_Z carries the spectator ladder W onto the sideways ladder S⁺, an operator identity
/// at every N; on the open XY chain W commutes for every real symmetric a, while S⁺ fails for every
/// correlated one, under one common covariance only by the block shift −8κ(d + 1).</para>
///
/// <para>Live witness: <c>inspect --root sharedfield</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/CorrelatedZDephasingPriceWitness.cs</c>), which
/// builds L from integer jumps through the Lindbladian builder and reads the residual in the Pauli
/// basis with the canonical Π, a different road from the page's integer verifier
/// <c>simulations/mirror_under_shared_field.py</c>.</para>
///
/// <para>Typed parents: <see cref="F1PalindromeIdentity"/> (the identity this generalises, K = σ·1
/// being its case) and <see cref="AbsorptionTheoremClaim"/> (its cell-rate reading at the two ends the
/// shared field interpolates, local dephasing and the collective jump √γ·Σ_k Z_k, proof §2). The scope
/// of <c>SidewaysSpinLadderClaim</c> (Tier 1 candidate, so not a typed parent here) is narrowed by the
/// ladder statement.</para></summary>
public sealed class CorrelatedZDephasingPriceClaim : Claim
{
    public F1PalindromeIdentity F1Palindrome { get; }
    public AbsorptionTheoremClaim AbsorptionTheorem { get; }

    /// <summary>The closed form ‖M_D‖²_F = 8·4^N·Σ_{k&lt;l} a_kl² of the dissipator residual.</summary>
    public static double ResidualNormSquared(int n, double sumSquaredCovariances)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (sumSquaredCovariances < 0) throw new ArgumentOutOfRangeException(nameof(sumSquaredCovariances), "a sum of squares is not negative");
        return 8.0 * Math.Pow(4, n) * sumSquaredCovariances;
    }

    /// <summary>σ(p, q) = σ − Nκ + κ(d² + m²), d = p − q, m = p + q − N, the rate-side value: for a
    /// popcount-conserving X^⊗N-even H, block (p, q) and its image pair about the eigenvalue centre
    /// −σ(p, q) under one common covariance κ, σ the summed site rate; σ(p, q) is also the
    /// mean (K_p + K_q)/2 of K on the block's two popcounts.</summary>
    public static double BlockSigma(int n, double sigma, double kappa, int p, int q)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (p < 0 || p > n) throw new ArgumentOutOfRangeException(nameof(p), p, "a popcount lies in [0, N]");
        if (q < 0 || q > n) throw new ArgumentOutOfRangeException(nameof(q), q, "a popcount lies in [0, N]");
        int d = p - q, m = p + q - n;
        return sigma - n * kappa + kappa * (d * d + m * m);
    }

    public CorrelatedZDephasingPriceClaim(F1PalindromeIdentity f1Palindrome, AbsorptionTheoremClaim absorptionTheorem)
        : base(
            "Correlated Z dephasing: Π_Z·L·Π_Z⁻¹ = −L − {K, ·}, K = Σ_kl a_kl Z_l Z_k, for a real symmetric Kossakowski " +
            "matrix a over the Z_k and a real H commuting with X^⊗N; for such an H that also conserves the popcount, under " +
            "one common covariance the spectrum pairs block by block, each joint-popcount block with its image about a centre of its own",
            Tier.Tier1Derived,
            "experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md")
    {
        F1Palindrome = f1Palindrome ?? throw new ArgumentNullException(nameof(f1Palindrome));
        AbsorptionTheorem = absorptionTheorem ?? throw new ArgumentNullException(nameof(absorptionTheorem));
    }

    public override string DisplayName =>
        "Correlated Z dephasing: the palindrome's price is the operator {Σc†c, ·}";

    public override string Summary =>
        "Π_Z·L·Π_Z⁻¹ = −L − {K, ·}, K = Σ_kl a_kl Z_l Z_k, for a real symmetric Kossakowski matrix a and a real X^⊗N-even H; " +
        "‖M_D‖² = 8·4^N·Σ_{k<l} a_kl²; one common covariance κ and such an H that also conserves the popcount: block (p, q) pairs with " +
        "(q, N − p) about −σ(p, q), σ(p, q) = σ − Nκ + κ(d² + m²), d = p − q, m = p + q − N; no single centre " +
        "serves the spectrum for a positive semidefinite a with a positive covariance sum (H keeping the all-0 and all-1 " +
        "states as eigenstates) " +
        "(Tier 1 derived; live witness: inspect --root sharedfield).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("price identity",
                summary: "for a real X^⊗N-even H, Π_Z L Π_Z⁻¹ = −L − {K, ·}, Π_Z(ρ) = ρᵀ·X^⊗N (F1's Π); F1 is the case K = σ·1 (a diagonal); the one-sided form with ρ ↦ ρ† gives the same blockwise spectral pairing for a complex, popcount-conserving X^⊗N-even H");
            yield return new InspectableNode("residual",
                summary: "M_D = −{K − σ, ·}, ‖M_D‖²_F = 8·4^N·Σ_{k<l} a_kl², H-independent, adding in squares to the Hamiltonian's residual");
            yield return new InspectableNode("shared field",
                summary: "a_kl = κ: K = σ − Nκ + 4κS_z²; on block (p, q) (ket and bra popcounts) the generator is the local one at γ_k − κ minus 2κ(p − q)²; each block pairs with (q, N − p) about −σ(p, q), σ(p, q) = σ − Nκ + κ(d² + m²), d = p − q, m = p + q − N, for a popcount-conserving X^⊗N-even H");
            yield return new InspectableNode("global break",
                summary: "for a positive semidefinite a whose covariances sum to a positive number and an H keeping the all-0 and all-1 states as eigenstates, no single centre serves the spectrum (the GHZ coherence's partner would grow)");
            yield return new InspectableNode("general pattern",
                summary: "for a popcount-conserving X^⊗N-even H a joint-popcount block pairs with its image about −(K_p + K_q)/2 wherever K takes one value on each of its two popcounts; the four edge blocks pay the local bill whatever the covariances");
            yield return new InspectableNode("ladders",
                summary: "Π_Z W Π_Z⁻¹ = S⁺ at every N, carried by the bra flip alone (a staggered Shiba map); on the open XY chain W commutes for every real symmetric a, while S⁺ fails for every correlated one, under one common covariance only by [L, S⁺] = −8κ(d + 1)·S⁺ on the source block (p, q), d = p − q");
        }
    }
}
