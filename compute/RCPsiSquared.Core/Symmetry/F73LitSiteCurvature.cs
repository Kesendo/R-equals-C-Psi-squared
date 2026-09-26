using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The F73 lit-site corollary, computed exactly: one lit site l, the pure state
/// ψ = (|vac⟩ + |1_l⟩)/√2, and the Taylor coefficients at t = 0 of three readings of site l:
/// the F73 term p_l(t) = 2·|(ρ_l)_{0,1}(t)|², the excitation weight n_l(t) = (ρ_l)_{1,1}(t) and the
/// site purity P_l(t) = Tr(ρ_l(t)²).
///
/// <para>The corollary (docs/ANALYTICAL_FORMULAS.md F73, "Corollary: the lit site"): for any
/// Hermitian H with [H, N_total] = 0 and uniform Z-dephasing γ₀,
/// <code>
///   p_l(t) = ½·e^(−4γ₀t)·S(t), S γ₀-free       (every t: γ₀ is the F73 term's only rate)
///   p_l'(0)  = P_l'(0)  = −2γ₀                    (no H in it)
///   p_l''(0) = P_l''(0) = 8γ₀² − V_l,   V_l = Σ_{j≠l} |⟨1_j|H|1_l⟩|²
/// </code>
/// so at γ₀ = 0 the loss is (V_l/2)·t² + O(t⁴), and on the XXZ graph family
/// H = Σ_b J_b (X_a X_b + Y_a Y_b + Δ_b Z_a Z_b) + Σ_i h_i Z_i (Pauli form, hopping 2J_b)
/// V_l = 4·Σ_{b∋l} J_b², the loss 2·Σ_{b∋l} J_b²·t², which is 2·deg(l)·J²·t² at uniform J.
/// Δ_b and h_i are diagonal in the computational basis and a variance is blind to the diagonal.
/// The purity agrees with p_l through the third derivative and first departs at the fourth,
/// P_l'''' − p_l'''' = 12·V_l², since n_l'''(0) = 4γ₀·V_l: P_l is not ½e^(−4γ₀t) times a γ₀-free
/// function.</para>
///
/// <para>Everything here is exact over ℚ(i): ρ₀ has entries ½, H and γ₀ are rational, each
/// derivative is ρ⁽ᵏ⁾(0) = 𝓛ᵏ ρ₀ with 𝓛ρ = −i[H, ρ] + γ₀ Σ_i (Z_i ρ Z_i − ρ), computed exactly.
/// The all-t factorization of p_l is the derivation's; the tests check it through k = 6. Two ways out of F73's class are built in. The pairing term G(XX − YY) on a bond is
/// the mutation: on the gated chain it turns the lit bond's J_b² into (J_b − G)², so it breaks
/// the law for every G other than 0 and 2J_b. The transverse field h^x leaves the second
/// derivatives standing on the gated cases, the lit site included; so U(1) is sufficient for the
/// law, and the tests do not claim it is necessary.</para></summary>
public static class F73LitSiteCurvature
{
    /// <summary>One bond J·(X_a X_b + Y_a Y_b + Δ·Z_a Z_b) + G·(X_a X_b − Y_a Y_b), Pauli convention.
    /// G = 0 is the XXZ family inside F73's class; G ≠ 0 (J_x ≠ J_y, the pairing term) breaks
    /// [H, N_total] = 0.</summary>
    public sealed record Bond(int A, int B, BigRational J, BigRational Delta)
    {
        public BigRational G { get; init; } = BigRational.Zero;
    }

    /// <summary>Value, first and second derivative at t = 0 of the F73 term p_l and the purity P_l.</summary>
    public sealed record Taylor(
        BigRational P0Term, BigRational P1Term, BigRational P2Term,
        BigRational P0Purity, BigRational P1Purity, BigRational P2Purity);

    /// <summary>The full 2^N Hamiltonian, basis index bit i = site i, bit 1 = excited (Z = −1).</summary>
    public static BigRational[,] Hamiltonian(int n, IReadOnlyList<Bond> bonds,
        IReadOnlyList<BigRational>? hz = null, IReadOnlyList<BigRational>? hx = null)
    {
        int d = 1 << n;
        var h = new BigRational[d, d];
        for (int r = 0; r < d; r++) for (int c = 0; c < d; c++) h[r, c] = BigRational.Zero;
        for (int s = 0; s < d; s++)
        {
            foreach (var b in bonds)
            {
                int ba = (s >> b.A) & 1, bb = (s >> b.B) & 1;
                h[s, s] += b.J * b.Delta * (ba == bb ? 1 : -1);
                // XX + YY = 2(σ⁺σ⁻ + σ⁻σ⁺): hops a differing pair with amplitude 2J.
                if (ba != bb) h[s ^ (1 << b.A) ^ (1 << b.B), s] += 2 * b.J;
                // XX − YY = 2(σ⁺σ⁺ + σ⁻σ⁻): creates or annihilates an equal pair with amplitude 2G.
                else if (!b.G.IsZero) h[s ^ (1 << b.A) ^ (1 << b.B), s] += 2 * b.G;
            }
            for (int i = 0; i < n; i++)
            {
                if (hz is not null) h[s, s] += hz[i] * (((s >> i) & 1) == 0 ? 1 : -1);
                if (hx is not null) h[s ^ (1 << i), s] += hx[i];
            }
        }
        return h;
    }

    /// <summary>V_l = (H²)_{ll} − (H_{ll})² on the single-excitation block, read off H itself.</summary>
    public static BigRational SingleExcitationVariance(BigRational[,] h, int lit)
    {
        int d = h.GetLength(0), l = 1 << lit;
        var v = BigRational.Zero;
        for (int k = 0; k < d; k++) if (k != l && System.Numerics.BitOperations.PopCount((uint)k) == 1) v += h[k, l] * h[k, l];
        return v;
    }

    /// <summary>The closed form: V_l = 4·Σ_{b∋l} J_b² on the XXZ graph family (Pauli convention).</summary>
    public static BigRational PredictedVariance(IReadOnlyList<Bond> bonds, int lit)
    {
        var v = BigRational.Zero;
        foreach (var b in bonds) if (b.A == lit || b.B == lit) v += 4 * b.J * b.J;
        return v;
    }

    /// <summary>The corollary's prediction for the second derivative: 8γ₀² − V_l.</summary>
    public static BigRational PredictedSecondDerivative(BigRational gamma, BigRational variance) =>
        8 * gamma * gamma - variance;

    /// <summary>Exact Taylor data at t = 0 from 𝓛ρ₀ and 𝓛²ρ₀.</summary>
    public static Taylor TaylorAtZero(int n, BigRational[,] h, BigRational gamma, int lit)
    {
        var d = Derivatives(n, h, gamma, lit, 2);
        return new Taylor(d.Term[0], d.Term[1], d.Term[2], d.Purity[0], d.Purity[1], d.Purity[2]);
    }

    /// <summary>The k-th derivatives at t = 0, k = 0..order, of the three readings of site l.</summary>
    public sealed record DerivativeSeries(BigRational[] Term, BigRational[] Weight, BigRational[] Purity);

    /// <summary>Derivatives of p_l, n_l and P_l at t = 0 up to <paramref name="order"/>, by Leibniz
    /// on ρ_l⁽ᵏ⁾ = Tr_rest(𝓛ᵏρ₀).</summary>
    public static DerivativeSeries Derivatives(int n, BigRational[,] h, BigRational gamma, int lit, int order)
    {
        int d = 1 << n, l = 1 << lit;
        var half = new BigRational(1, 2);
        var rho = Zero(d);
        rho[0, 0] = half; rho[0, l] = half; rho[l, 0] = half; rho[l, l] = half;
        var red = new List<GaussianRational[,]> { Reduce(n, rho, lit) };
        for (int k = 1; k <= order; k++)
        {
            rho = Apply(n, h, gamma, rho);
            red.Add(Reduce(n, rho, lit));
        }
        var term = new BigRational[order + 1];
        var weight = new BigRational[order + 1];
        var purity = new BigRational[order + 1];
        for (int k = 0; k <= order; k++)
        {
            var pk = GaussianRational.Zero;
            var qk = GaussianRational.Zero;
            for (int i = 0; i <= k; i++)
            {
                GaussianRational binom = new BigRational(Binomial(k, i));
                // p = 2 c c̄ with c = (ρ_l)_{0,1}; P = Tr(ρ_l ρ_l).
                pk += binom * red[i][0, 1] * red[k - i][0, 1].Conjugate;
                qk += binom * TraceProduct(red[i], red[k - i]);
            }
            if (!pk.Im.IsZero || !qk.Im.IsZero || !red[k][1, 1].Im.IsZero)
                throw new InvalidOperationException("A Hermitian reduced state gave a non-real coefficient.");
            term[k] = 2 * pk.Re;
            purity[k] = qk.Re;
            weight[k] = red[k][1, 1].Re;
        }
        return new DerivativeSeries(term, weight, purity);
    }

    /// <summary>The k-th derivative at 0 of e^(−4γ₀t)·f(t), from f's derivatives (Leibniz).</summary>
    public static BigRational DampedDerivative(IReadOnlyList<BigRational> f, BigRational gamma, int k)
    {
        var sum = BigRational.Zero;
        for (int i = 0; i <= k; i++)
        {
            var power = BigRational.One;
            for (int m = 0; m < k - i; m++) power *= -4 * gamma;
            sum += new BigRational(Binomial(k, i)) * power * f[i];
        }
        return sum;
    }

    private static long Binomial(int k, int i)
    {
        long b = 1;
        for (int m = 1; m <= i; m++) b = b * (k - i + m) / m;
        return b;
    }

    private static GaussianRational[,] Zero(int d)
    {
        var m = new GaussianRational[d, d];
        for (int r = 0; r < d; r++) for (int c = 0; c < d; c++) m[r, c] = GaussianRational.Zero;
        return m;
    }

    /// <summary>𝓛ρ = −i(Hρ − ρH) + γ₀ Σ_i (Z_i ρ Z_i − ρ); the dissipator is −2γ₀·(disagreeing sites).</summary>
    private static GaussianRational[,] Apply(int n, BigRational[,] h, BigRational gamma, GaussianRational[,] rho)
    {
        int d = 1 << n;
        var outM = Zero(d);
        var minusI = new GaussianRational(BigRational.Zero, -1);
        for (int r = 0; r < d; r++)
            for (int c = 0; c < d; c++)
            {
                var comm = GaussianRational.Zero;
                for (int k = 0; k < d; k++)
                {
                    if (!h[r, k].IsZero) comm += (GaussianRational)h[r, k] * rho[k, c];
                    if (!h[k, c].IsZero) comm -= rho[r, k] * (GaussianRational)h[k, c];
                }
                int disagree = System.Numerics.BitOperations.PopCount((uint)(r ^ c));
                outM[r, c] = minusI * comm + (GaussianRational)(-2 * gamma * disagree) * rho[r, c];
            }
        return outM;
    }

    private static GaussianRational[,] Reduce(int n, GaussianRational[,] rho, int lit)
    {
        int d = 1 << n, bit = 1 << lit;
        var red = new GaussianRational[2, 2];
        for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) red[a, b] = GaussianRational.Zero;
        for (int s = 0; s < d; s++)
        {
            if ((s & bit) != 0) continue;
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                    red[a, b] += rho[s | (a * bit), s | (b * bit)];
        }
        return red;
    }

    private static GaussianRational TraceProduct(GaussianRational[,] a, GaussianRational[,] b)
    {
        var t = GaussianRational.Zero;
        for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++) t += a[i, k] * b[k, i];
        return t;
    }
}
