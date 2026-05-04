namespace RCPsiSquared.Core.Symmetry;

/// <summary>Closed-form predictor for the Π²-odd-fraction-within-memory of
/// popcount-coherence states |ψ⟩ = (|p⟩ + |q⟩)/√2 with popcount(p) = n_p,
/// popcount(q) = n_p + 1, plus the Krawtchouk-form verifier.
///
/// <para>Π²-odd / memory = (1 / 2 − α · s) / (1 − s), with HD-invariant
/// s = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q)) and α the Π²-odd-Frobenius²-fraction
/// of the kernel projection. α has a three-anchor closed form
/// (<see cref="AlphaThreeAnchor"/>) bit-exactly equivalent to the full
/// Krawtchouk sum (<see cref="AlphaKrawtchouk"/>) across N = 3..16.</para>
///
/// <para>See <c>docs/proofs/PROOF_F86_QPEAK.md</c> §Structural inheritance
/// from F88 for the canonical statement, structural reasons (Krawtchouk
/// reflection K_{N−n}(s; N) = (−1)^s K_n(s; N) at popcount-mirror;
/// K_{N/2}(s; N) = 0 for odd s at near-mirror near-half), and verified table.
/// Python counterpart: <c>simulations/_pi2_odd_landscape_sweep.py</c>.</para>
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

    /// <summary>Static Frobenius² fraction s = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q)) of
    /// the kernel projection of |ψ⟩⟨ψ| with |ψ⟩ = (|p⟩ + |q⟩)/√2,
    /// popcount(p) = n_p, popcount(q) = n_q. HD-invariant.</summary>
    public static double StaticFraction(int N, int np, int nq)
    {
        return 0.25 / Binomial(N, np) + 0.25 / Binomial(N, nq);
    }

    /// <summary>α = Π²-odd-Frobenius²-fraction of the kernel projection ρ_d0 of
    /// |ψ⟩⟨ψ|, computed via the explicit Krawtchouk sum
    /// Σ_{s odd} C(N, s)(A_s + B_s)² / Σ_s C(N, s)(A_s + B_s)² with
    /// A_s = K_{n_p}(s; N)/C(N, n_p), B_s = K_{n_q}(s; N)/C(N, n_q).
    /// HD/bit-position invariant.</summary>
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

    /// <summary>Three-anchor closed form for α(N, n_p, n_p + 1). Bit-exactly equivalent
    /// to <see cref="AlphaKrawtchouk"/> for all popcount-coherence pairs at N = 3..16.</summary>
    public static double AlphaThreeAnchor(int N, int np, int nq)
    {
        if (np + nq == N) return 0.0;
        if (IsNearMirrorHalf(N, np, nq)) return (double)(N + 2) / (4.0 * (N + 1));
        return 0.5;
    }

    /// <summary>Closed-form prediction Π²-odd-fraction-within-memory =
    /// (1/2 − α·s) / (1 − s) for |ψ⟩ = (|p⟩ + |q⟩)/√2 with
    /// popcount(p) = n_p, popcount(q) = n_q. Uses <see cref="AlphaKrawtchouk"/>.</summary>
    public static double Pi2OddInMemory(int N, int np, int nq)
    {
        double s = StaticFraction(N, np, nq);
        double a = AlphaKrawtchouk(N, np, nq);
        return (0.5 - a * s) / (1.0 - s);
    }

    /// <summary>True iff n_p + n_q = N (popcount-mirror configuration; X-flip
    /// conjugation cancels all odd-|S| Pauli strings between sectors).</summary>
    public static bool IsPopcountMirror(int N, int np, int nq) => np + nq == N;

    /// <summary>True iff even N and (n_p, n_q) is one of the two near-mirror near-half
    /// pairs {(N/2 − 1, N/2), (N/2, N/2 + 1)}. K_{N/2}(s; N) vanishes for odd s
    /// (bit-flip symmetry of the half-popcount sector), giving α = (N + 2)/(4(N + 1)).</summary>
    public static bool IsNearMirrorHalf(int N, int np, int nq) =>
        (N % 2 == 0) &&
        ((np == N / 2 - 1 && nq == N / 2) || (np == N / 2 && nq == N / 2 + 1));
}
