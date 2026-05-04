namespace RCPsiSquared.Core.Symmetry;

/// <summary>Closed-form predictor for the Π²-odd-fraction-within-memory of
/// popcount-coherence pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 (any popcount pair, any
/// HD), plus the Krawtchouk-form verifier.
///
/// <para>Π²-odd / memory = (1/2 − α · s) / (1 − s), with α = 0 at popcount-mirror,
/// α = Krawtchouk-rational at K-intermediate (one of n_p, n_q equals N/2 at even N),
/// α = 1/2 generic; Π²-odd / memory = 0 at the HD = N Π²-classical anchor (GHZ_N,
/// Bell-states, intra-complements). The HD = N anchor connects to F60: GHZ_N has
/// pair-CΨ = 0 (partial-trace) and Π²-EVEN-only content (projection): two
/// orthogonal readings of the same "classical" classification.</para>
///
/// <para>Canonical statement, structural reasons (Krawtchouk reflection
/// K_{N−n}(s; N) = (−1)^s K_n(s; N) at popcount-mirror; K_{N/2}(s; N) = 0 for odd s
/// at K-intermediate), and verified table: see
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> §Structural inheritance from F88.
/// Python verifier (213 cases at N = 2..7, max deviation 8.88e−16):
/// <c>simulations/_pi2_odd_general_closed_form.py</c>.</para>
/// </summary>
public static class PopcountCoherencePi2Odd
{
    /// <summary>Krawtchouk polynomial K_n(s; N) = Σ_k (−1)^k · C(s, k) · C(N − s, n − k).
    /// Integer-valued; valid k range is [max(0, n − (N − s)), min(s, n)] so both
    /// binomials have non-negative arguments.</summary>
    public static long Krawtchouk(int n, int s, int N)
    {
        int kLo = Math.Max(0, n - (N - s));
        int kHi = Math.Min(s, n);
        long total = 0;
        for (int k = kLo; k <= kHi; k++)
        {
            long term = Binomial(s, k) * Binomial(N - s, n - k);
            total += (k % 2 == 0) ? term : -term;
        }
        return total;
    }

    /// <summary>Binomial coefficient C(n, k); 0 if k out of [0, n].</summary>
    public static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n - k) k = n - k;
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    /// <summary>Static Frobenius² fraction of the kernel projection of |ψ⟩⟨ψ|
    /// with |ψ⟩ = (|p⟩ + |q⟩)/√2. Branches on inter-sector vs intra-sector
    /// because the kernel of L for Heisenberg + Z-dephasing is span{P_n}: when
    /// n_p = n_q, ρ_diag projects onto a single sector; when n_p ≠ n_q, onto two.
    /// HD/bit-position invariant.</summary>
    public static double StaticFraction(int N, int np, int nq)
    {
        if (np == nq) return 1.0 / Binomial(N, np);
        return 0.25 / Binomial(N, np) + 0.25 / Binomial(N, nq);
    }

    /// <summary>α = Π²-odd-Frobenius²-fraction of the kernel projection ρ_d0 of
    /// |ψ⟩⟨ψ|, computed via the explicit Krawtchouk sum
    /// Σ_{s odd} C(N, s)(A_s + B_s)² / Σ_s C(N, s)(A_s + B_s)² with
    /// A_s = K_{n_p}(s; N)/C(N, n_p), B_s = K_{n_q}(s; N)/C(N, n_q). Universal
    /// over all popcount pairs; HD/bit-position invariant.</summary>
    public static double AlphaKrawtchouk(int N, int np, int nq)
    {
        double cnp = Binomial(N, np);
        double cnq = Binomial(N, nq);
        double oddSum = 0, total = 0;
        for (int s = 0; s <= N; s++)
        {
            double aS = Krawtchouk(np, s, N) / cnp;
            double bS = Krawtchouk(nq, s, N) / cnq;
            double term = Binomial(N, s) * (aS + bS) * (aS + bS);
            total += term;
            if ((s & 1) == 1) oddSum += term;
        }
        return total > 0 ? oddSum / total : 0.0;
    }

    /// <summary>α via the three structural anchors. For α = 0 (popcount-mirror,
    /// covers intra-mirror at n_p = n_q = N/2) and α = 1/2 (generic) returns the
    /// closed-form value directly. For K-intermediate (one of n_p, n_q equals N/2
    /// at even N, no mirror), falls back to <see cref="AlphaKrawtchouk"/> since
    /// the value is a Krawtchouk-derived rational that varies with (n_p, n_q).
    /// Bit-exactly equivalent to <see cref="AlphaKrawtchouk"/> on all tested cases
    /// (N = 2..7, 213 popcount-coherence-pair configurations). The α = 1/2
    /// generic case is proven analytically only for the boundary pair (0, 1)
    /// (via Σ_s C(N, s) K_n(s; N) (−1)^s = 2^N · [n = N]); for other generic
    /// (n_p, n_q) it is verified bit-exact rather than formally derived.</summary>
    public static double AlphaAnchor(int N, int np, int nq)
    {
        if (np + nq == N) return 0.0;
        if (N % 2 == 0 && (np == N / 2 || nq == N / 2)) return AlphaKrawtchouk(N, np, nq);
        return 0.5;
    }

    /// <summary>Closed-form prediction Π²-odd-fraction-within-memory for
    /// |ψ⟩ = (|p⟩ + |q⟩)/√2 with popcount(p) = n_p, popcount(q) = n_q,
    /// HD(p, q) = hd. Returns 0 at the HD = N (Π²-classical) anchor;
    /// otherwise (1/2 − α·s) / (1 − s) via <see cref="AlphaKrawtchouk"/>
    /// (the universal verifier; equivalent to <see cref="AlphaAnchor"/>
    /// but unconditional, no anchor-classification branching).</summary>
    public static double Pi2OddInMemory(int N, int np, int nq, int hd)
    {
        if (hd == N) return 0.0;
        double s = StaticFraction(N, np, nq);
        double a = AlphaKrawtchouk(N, np, nq);
        return (0.5 - a * s) / (1.0 - s);
    }

    /// <summary>True iff n_p + n_q = N (popcount-mirror configuration; X-flip
    /// conjugation cancels all odd-|S| Pauli strings between sectors). Covers
    /// both inter-sector mirror (n_p ≠ n_q at odd N central pair) and intra-
    /// sector mirror (n_p = n_q = N/2 at even N).</summary>
    public static bool IsPopcountMirror(int N, int np, int nq) => np + nq == N;

    /// <summary>True iff even N and exactly one of {n_p, n_q} equals N/2 (no
    /// popcount-mirror). K_{N/2}(s; N) = 0 for odd s due to bit-flip symmetry
    /// of the half-popcount sector, producing a Krawtchouk-derived intermediate
    /// α (e.g. (N+2)/(4(N+1)) at adjacent near-mirror half, 5/13 at N=6
    /// popcount-(1, 3), etc.).</summary>
    public static bool IsKIntermediate(int N, int np, int nq) =>
        (N % 2 == 0) && (np + nq != N) && (np == N / 2 || nq == N / 2);

    /// <summary>True iff hd equals N (p and q are complementary bit patterns),
    /// forcing both static and memory Π²-odd content to zero (Y² = I cancellation
    /// in off-diagonal plus zero matching bits). Captures GHZ_N, Bell states
    /// at N = 2, intra-sector Singlet/Triplet at N = 2.</summary>
    public static bool IsHdComplement(int N, int hd) => hd == N;
}
