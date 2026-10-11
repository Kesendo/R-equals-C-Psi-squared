using System.Numerics;
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
/// keeping the all-0 and all-1 states as eigenstates), for a ⪰ 0 and a popcount-conserving X^⊗N-even H wherever K takes one
/// value above σ on some popcount (<see cref="FarEndBreakPopcounts"/>; under one common covariance every allowed κ ≠ 0), and, through
/// the third moment of the spectrum about −σ
/// (<see cref="ThirdMoment"/>), for a ⪰ 0 with a covariance and a traceless H commuting with X^⊗N whenever
/// Σ_{k≠l} a_kl·Tr(H²Z_kZ_l) ≤ 0, whatever the sign of the summed covariance (the pairing of the complex eigenvalues),
/// and, for a ⪰ 0 with a covariance and any H, at weak enough coupling, the rates included (an eigenvalue continued from
/// the cell of the configuration with the largest K lies below −2σ); for a general pattern and a popcount-conserving
/// X^⊗N-even H a block pairs wherever K takes one value on each of its two popcounts, and the popcounts that do
/// are settled for every pattern (<see cref="SingleValuePopcounts"/>: 0 and N always, 1 and N − 1 exactly when every
/// site carries the same summed covariance, N/2 exactly when the covariances are additive, a_kl = y_k + y_l, every
/// other popcount exactly when they are one number), K there taking the shared field's values at the mean
/// covariance; the four edge blocks pay the local bill whatever the
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

    /// <summary>The popcounts on which K = Σ_kl a_kl Z_l Z_k takes one value, read from the covariances alone (the
    /// off-diagonal part of a; its diagonal enters K only through σ and is ignored here): 0 and N always; 1 and N − 1
    /// exactly when every site carries the same summed covariance; N/2, for 2 ≤ N/2 ≤ N − 2, exactly when the
    /// covariances are additive, a_kl = y_k + y_l (k ≠ l); every other popcount exactly when they are one number. The
    /// page's proof writes the covariances as κ̄ + x_k + x_l + Ω_kl (Σx = 0, Ω with zero row sums) and cut(𝒜), the
    /// summed covariance across the excited set 𝒜, as κ̄·p(N − p) + (N − 2p)·x(𝒜) − 2·Ω(𝒜). The qualifying blocks are
    /// the pairs of these popcounts, 4, 16, 9 (even N) or (N + 1)², and on each (K_p + K_q)/2 is <see cref="BlockSigma"/>
    /// at the mean covariance κ = T/C(N, 2), T the summed covariance. Integer covariances (|a_kl| ≤ 10⁹), so every test
    /// is exact.</summary>
    public static IReadOnlyList<int> SingleValuePopcounts(long[,] covariance)
    {
        ArgumentNullException.ThrowIfNull(covariance);
        int n = covariance.GetLength(0);
        if (n < 1 || covariance.GetLength(1) != n)
            throw new ArgumentException("the covariances form a square N × N matrix, N ≥ 1", nameof(covariance));
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
            {
                if (covariance[k, l] != covariance[l, k])
                    throw new ArgumentException("a Kossakowski matrix is symmetric", nameof(covariance));
                if (k != l && (covariance[k, l] > 1_000_000_000L || covariance[k, l] < -1_000_000_000L))
                    throw new ArgumentOutOfRangeException(nameof(covariance), covariance[k, l], "|a_kl| ≤ 10⁹ (k ≠ l) keeps every sum exact");
            }
        if (n <= 2) return Enumerable.Range(0, n + 1).ToArray();

        bool constant = true, equalRows = true, additive = true;
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++)
                if (covariance[k, l] != covariance[0, 1]) constant = false;
        long RowSum(int k) { long r = 0; for (int l = 0; l < n; l++) if (l != k) r += covariance[k, l]; return r; }
        long r0 = RowSum(0);
        for (int k = 1; k < n; k++) if (RowSum(k) != r0) equalRows = false;
        var y2 = new long[n];                                   // 2·y_k, so the additive test stays in the integers
        y2[0] = covariance[0, 1] + covariance[0, 2] - covariance[1, 2];
        for (int k = 1; k < n; k++) y2[k] = 2 * covariance[0, k] - y2[0];
        for (int k = 0; k < n && additive; k++)
            for (int l = k + 1; l < n; l++)
                if (2 * covariance[k, l] != y2[k] + y2[l]) { additive = false; break; }

        var popcounts = new List<int>();
        for (int p = 0; p <= n; p++)
            if (p == 0 || p == n || constant
                || ((p == 1 || p == n - 1) && equalRows)
                || (2 * p == n && additive))
                popcounts.Add(p);
        return popcounts;
    }

    /// <summary>The popcounts whose far-end mode breaks the single mirror: for a popcount-conserving H commuting with
    /// X^⊗N and a popcount p on which K takes one value K_p (<see cref="SingleValuePopcounts"/>), P_p·X^⊗N (the
    /// bra-flipped populations of popcount p) is an eigenvector at −2K_p (H commutes with it, and each of its cells
    /// |i⟩⟨ī| pays 2K(i)), and for a ⪰ 0, wherever K_p > σ, its partner about −σ would grow, so not even the rates pair.
    /// On such a popcount K_p = σ − Nκ̄ + κ̄(N − 2p)², κ̄ the mean covariance, so K_p > σ exactly when
    /// T·((N − 2p)² − N) > 0, T the summed covariance: for T > 0 the qualifying popcounts with (N − 2p)² > N, the corners
    /// (the GHZ coherence) always among them; for T < 0 those with (N − 2p)² < N, the middle wherever it qualifies (one
    /// common covariance, an additive pattern at even N). Integer covariances (|a_kl| ≤ 10⁹), exact.</summary>
    public static IReadOnlyList<int> FarEndBreakPopcounts(long[,] covariance)
    {
        var qualifying = SingleValuePopcounts(covariance);
        int n = covariance.GetLength(0);
        long total = 0;
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++) total += covariance[k, l];
        return qualifying.Where(p => (BigInteger)total * ((long)(n - 2 * p) * (n - 2 * p) - n) > 0).ToArray();
    }

    /// <summary>The noise term of the third moment's bracket, 4·tr(C·(2a + Λ)·C) = 12·tr(ΛC²) + 8·tr(C³), Λ the diagonal of
    /// a and C its covariances (<see cref="ThirdMoment"/> is −4^N·[this − 3·tr(CΘ) − 6·ξᵀaξ]). For a ⪰ 0 with a covariance
    /// it is positive, 8·tr(CaC) + 4·tr(CΛC); an indefinite a can make it negative (zero diagonal and every covariance −1
    /// at N = 3 give −48). Integer entries, exact.</summary>
    public static BigInteger ThirdMomentNoiseTerm(long[,] a)
    {
        int n = SymmetricSize(a, nameof(a));
        BigInteger s = 0;
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
            {
                if (l == k) continue;
                for (int m = 0; m < n; m++)
                {
                    if (m == k) continue;
                    BigInteger w = 2 * (BigInteger)a[l, m] + (l == m ? a[l, l] : 0);
                    s += (BigInteger)a[k, l] * w * a[m, k];
                }
            }
        return 4 * s;
    }

    /// <summary>The third moment of the spectrum about F137's one candidate centre −σ, exact over the integers for every
    /// real symmetric a: Tr((L + σ)³) = −[4^N·W − 3·2^N·Σ_{k≠l} a_kl·T2_kl − 6·Σ_kl a_kl·T1_k·T1_l], W the
    /// <see cref="ThirdMomentNoiseTerm"/>, T2_kl = Tr(H²Z_kZ_l) and T1_k = Tr(H Z_k) for a TRACELESS Hermitian H (the
    /// page's Θ_kl = T2_kl/2^N and ξ_k = T1_k/2^N; a multiple of the identity in H changes nothing in L but would change
    /// T2). A palindrome about −σ makes it zero, so for a ⪰ 0 with a covariance and a traceless H commuting with X^⊗N
    /// (T1 = 0) a nonzero value, which every Σ_{k≠l} a_kl·T2_kl ≤ 0 forces, rules out the pairing of the complex
    /// eigenvalues; whether the rates alone pair it does not say.</summary>
    public static BigInteger ThirdMoment(long[,] a, long[,] twoSiteZZ, long[] siteZ)
    {
        int n = SymmetricSize(a, nameof(a));
        ArgumentNullException.ThrowIfNull(twoSiteZZ);
        ArgumentNullException.ThrowIfNull(siteZ);
        if (twoSiteZZ.GetLength(0) != n || twoSiteZZ.GetLength(1) != n)
            throw new ArgumentException($"Tr(H²Z_kZ_l) is an N × N table, N = {n}", nameof(twoSiteZZ));
        if (siteZ.Length != n)
            throw new ArgumentException($"Tr(H Z_k) has one entry per site, N = {n}", nameof(siteZ));
        BigInteger hTwo = 0, hOne = 0;
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
            {
                if (k != l) hTwo += (BigInteger)a[k, l] * twoSiteZZ[k, l];
                hOne += (BigInteger)a[k, l] * siteZ[k] * siteZ[l];
            }
        return -(BigInteger.Pow(4, n) * ThirdMomentNoiseTerm(a) - 3 * BigInteger.Pow(2, n) * hTwo - 6 * hOne);
    }

    private static int SymmetricSize(long[,] a, string name)
    {
        ArgumentNullException.ThrowIfNull(a, name);
        int n = a.GetLength(0);
        if (n < 1 || a.GetLength(1) != n)
            throw new ArgumentException("a Kossakowski matrix is square, N × N with N ≥ 1", name);
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++)
                if (a[k, l] != a[l, k])
                    throw new ArgumentException("a Kossakowski matrix is symmetric", name);
        return n;
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
        "states as eigenstates), nor, for a ⪰ 0 and a popcount-conserving X^⊗N-even H, wherever K takes one value above σ " +
        "on a popcount (a far-end mode decaying faster than 2σ; one allowed common covariance of either sign at every N), nor, for a ⪰ 0 " +
        "with a covariance and any H, at weak enough coupling, nor for its complex eigenvalues, for a ⪰ 0 with a covariance and a traceless " +
        "X^⊗N-even H, whenever Σ_{k≠l} a_kl·Tr(H²Z_kZ_l) ≤ 0, either sign of the sum (the third moment); for a general pattern popcount p sees one value of K exactly when p = 0, N, when p = 1, " +
        "N − 1 and every site carries the same summed covariance, when p = N/2 and the covariances are additive, and " +
        "otherwise when they are one number, K there the shared field's at the mean covariance " +
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
                summary: "for a positive semidefinite a whose covariances sum to a positive number and an H keeping the all-0 and all-1 states as eigenstates, no single centre serves the spectrum, the rates included (the GHZ coherence's partner would grow); for a ⪰ 0 with a covariance and a traceless H commuting with X^⊗N, none either whenever Σ_{k≠l} a_kl·Tr(H²Z_kZ_l) ≤ 0, whatever the sign of the summed covariance (the third moment, below; the complex eigenvalues only); for a ⪰ 0 with a covariance and any H, none either at weak enough coupling, the rates included (an eigenvalue continued from the cell of the configuration with the largest K lies below −2σ); where no rate route reaches, whether the rates pair is open; for a ⪰ 0 and a popcount-conserving X^⊗N-even H, none either wherever K takes one value K_p > σ on a popcount p (FarEndBreakPopcounts): P_p·X^⊗N is then an eigenvector at −2K_p, the rates included, which under one common covariance is every allowed κ ≠ 0 at every N (p = ⌊N/2⌋ for κ < 0, the GHZ case p = 0 for κ > 0); for one common covariance at uniform rates on the open XY and Heisenberg chains the third and fifth moments together reach the complex half of that by a second route at N = 2..7 (exact computer algebra, simulations/shared_field_odd_moments.py)");
            yield return new InspectableNode("the third moment",
                summary: "Tr((L + σ)³) = −4^N·[4·tr(C·(2a + Λ)·C) − 3·tr(CΘ) − 6·ξᵀaξ] for every real symmetric a = Λ + C and traceless Hermitian H, Θ_kl = Tr(H²Z_kZ_l)/2^N, ξ_k = Tr(H Z_k)/2^N (ThirdMoment, ThirdMomentNoiseTerm); C = 0 leaves F117's third-moment face 6·4^N·Σ_k γ_k ξ_k²; the noise term is 8·tr(CaC) + 4·tr(CΛC) > 0 for a ⪰ 0 with a covariance; an H commuting with X^⊗N has ξ = 0; on the open XY chain (H = J·Σ(XX + YY), the Pauli book) Θ is −2J² on the bonds, so the break holds whenever the bond covariances sum to zero or more");
            yield return new InspectableNode("general pattern",
                summary: "for a popcount-conserving X^⊗N-even H a joint-popcount block pairs with its image about −(K_p + K_q)/2 wherever K takes one value on each of its two popcounts; the four edge blocks pay the local bill whatever the covariances");
            yield return new InspectableNode("the four shapes",
                summary: "popcount p sees one value of K exactly when p = 0 or N; p = 1 or N − 1 and every site carries the same summed covariance; p = N/2 (2 ≤ p ≤ N − 2) and the covariances are additive, a_kl = y_k + y_l; or the covariances are one number (SingleValuePopcounts); so the qualifying blocks number 4, 16, 9 (even N) or (N + 1)², and on each (K_p + K_q)/2 = σ(p, q) at κ = T/C(N, 2), the mean covariance; whether any other block pairs is read, not proved");
            yield return new InspectableNode("ladders",
                summary: "Π_Z W Π_Z⁻¹ = S⁺ at every N, carried by the bra flip alone (a staggered Shiba map); on the open XY chain W commutes for every real symmetric a, while S⁺ fails for every correlated one, under one common covariance only by [L, S⁺] = −8κ(d + 1)·S⁺ on the source block (p, q), d = p − q");
        }
    }
}
