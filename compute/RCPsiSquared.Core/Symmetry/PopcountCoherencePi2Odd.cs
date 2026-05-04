namespace RCPsiSquared.Core.Symmetry;

/// <summary>Closed-form predictor for the Π²-odd-fraction-within-memory of
/// popcount-coherence pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 (any popcount pair, any
/// HD), plus the Krawtchouk-form verifier.
///
/// <para>Π²-odd / memory = (1/2 − α · s) / (1 − s), with all three α anchors in
/// closed form: α = 0 at popcount-mirror; α = C(N, N/2) / (2·(C(N, n_other) +
/// C(N, N/2))) at K-intermediate (one of n_p, n_q equals N/2 at even N, n_other
/// is the other); α = 1/2 generic. Π²-odd / memory = 0 at the HD = N
/// Π²-classical anchor (GHZ_N, Bell-states, intra-complements). The HD = N anchor
/// connects to F60: GHZ_N has pair-CΨ = 0 (partial-trace) and Π²-EVEN-only content
/// (projection): two orthogonal readings of the same "classical" classification.</para>
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

    /// <summary>α in closed form via the three structural anchors. All three
    /// branches are derived from the Krawtchouk identity
    /// Σ_s (−1)^s · C(N, s) · K_n(s; N) · K_m(s; N) = 2^N · C(N, n) · [n + m = N]
    /// (consequence of orthogonality + reflection K_n(s; N) = (−1)^s K_{N−n}(s; N)),
    /// which gives E − O := Σ_{s even} − Σ_{s odd} of C(N, s) · (A_s + B_s)² as
    /// a sum of three indicator-weighted terms; α = 0 (mirror), α = 1/2 (generic)
    /// follow by setting indicators to 0; α at K-intermediate is given by
    /// <see cref="AlphaKIntermediateClosed"/>. Bit-exactly equivalent to
    /// <see cref="AlphaKrawtchouk"/> at machine precision.</summary>
    public static double AlphaAnchor(int N, int np, int nq)
    {
        if (np + nq == N) return 0.0;
        if (N % 2 == 0 && (np == N / 2 || nq == N / 2)) return AlphaKIntermediateClosed(N, np, nq);
        return 0.5;
    }

    /// <summary>Closed form for α at K-intermediate configurations:
    /// α = C(N, N/2) / (2 · (C(N, n_other) + C(N, N/2))), where
    /// n_other ∈ {n_p, n_q} is the entry not equal to N/2. Derived from the
    /// Krawtchouk reflection-orthogonality lemma applied to E − O at the
    /// K-vanishing configuration. Recovers the adjacent near-mirror-half formula
    /// (N + 2)/(4·(N + 1)) when n_other = N/2 ± 1 (since C(N, N/2 ± 1) =
    /// C(N, N/2) · N/(N + 2)). Throws if (n_p, n_q) is not K-intermediate.</summary>
    public static double AlphaKIntermediateClosed(int N, int np, int nq)
    {
        if (!IsKIntermediate(N, np, nq))
            throw new ArgumentException($"({np}, {nq}) at N={N} is not K-intermediate");
        int nOther = (nq == N / 2) ? np : nq;
        double cHalf = Binomial(N, N / 2);
        double cOther = Binomial(N, nOther);
        return cHalf / (2.0 * (cOther + cHalf));
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

    // ──────────── Dicke-superposition extension (multi-state) ────────────
    //
    // Pair state (|p⟩+|q⟩)/√2 has total Π²-odd-of-ρ = 1/2 universally for HD < N.
    // Dicke superposition |ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2 (the canonical F86 K_CC_pr
    // probe at adjacent popcount-(n, n+1)) has a different total Π²-odd content
    // due to its bit-permutation symmetry. The static side is identical to the
    // pair-state formula (kernel only sees popcount weights, both have w_n = 1/2,
    // w_{n+1} = 1/2). The memory side differs: total Π²-odd-of-ρ has three anchors
    // in (N, n): 0 at popcount-mirror (X⊗N-symmetric → Π²-EVEN-only state),
    // 3/8 at K-intermediate (even N, n or n+1 = N/2), 1/2 generic.
    //
    // Verified bit-exact at N = 3..8 across all (n, n+1) pairs; analytical proof
    // of the 3/8 K-intermediate anchor remains open.

    /// <summary>True iff the Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 sits at the
    /// X⊗N-symmetric popcount-mirror configuration, where 2n + 1 = N (odd N only).
    /// At this anchor the state is X⊗N-symmetric, hence Π²-EVEN-only (X⊗N
    /// conjugation maps σ_α → (−1)^bit_b(α) σ_α, so eigenstates of X⊗N have zero
    /// expectation on Π²-odd Pauli strings). Π²-odd-fraction-within-memory = 0.
    ///
    /// <para>Inheritance: shares the X⊗N-symmetry root with the GHZ_N / Bell-state /
    /// HD = N pair-state anchor (<see cref="IsHdComplement"/>). Both classes are
    /// X⊗N-eigenstates and therefore F60-classical (pair-CΨ = 0) and F88-classical
    /// (Π²-odd = 0) simultaneously. Two structurally-distinct state classes
    /// (HD = N pair vs. Dicke-mirror multi-state) sharing the same Π²-classical
    /// origin via a single algebraic mechanism.</para></summary>
    public static bool IsDickeMirror(int N, int n) => 2 * n + 1 == N;

    /// <summary>True iff Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 sits at a
    /// K-intermediate configuration: even N and (n = N/2 − 1 or n = N/2).
    /// At this anchor total Π²-odd-of-ρ = 3/8 (verified bit-exact N = 4, 6, 8;
    /// analytical proof open).</summary>
    public static bool IsDickeKIntermediate(int N, int n) =>
        (N % 2 == 0) && (n == N / 2 - 1 || n == N / 2);

    /// <summary>Total Π²-odd-of-ρ for Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2.
    /// Three anchors:
    /// 0 at popcount-mirror (2n + 1 = N, odd N, state is X⊗N-symmetric);
    /// 3/8 at K-intermediate (even N, n ∈ {N/2 − 1, N/2});
    /// 1/2 generic. Differs from pair-state's universal 1/2 due to bit-permutation
    /// symmetry of Dicke states constraining off-diagonal Π²-odd content.</summary>
    public static double Pi2OddTotalDickeSuperposition(int N, int n)
    {
        if (IsDickeMirror(N, n)) return 0.0;
        if (IsDickeKIntermediate(N, n)) return 3.0 / 8.0;
        return 0.5;
    }

    /// <summary>Closed-form Π²-odd-fraction-within-memory for Dicke superposition
    /// |ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2. Static fraction s and α_static reuse the
    /// pair-state formulas (kernel projection only sees popcount weights;
    /// |D_n⟩'s diagonal projects to P_n / C(N, n) just like a single basis state).
    /// The memory side uses Pi2OddTotalDickeSuperposition's three-anchor structure:
    ///
    ///     Π²-odd / memory = (α_total_Dicke − α_static · s) / (1 − s)
    ///
    /// Verified bit-exact against numerical state-level reading at N = 3..8.</summary>
    public static double Pi2OddInMemoryDickeSuperposition(int N, int n)
    {
        double s = StaticFraction(N, n, n + 1);
        double aStatic = AlphaKrawtchouk(N, n, n + 1);
        double aTotal = Pi2OddTotalDickeSuperposition(N, n);
        return (aTotal - aStatic * s) / (1.0 - s);
    }
}
