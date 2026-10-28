using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The cat's centred odd moments m_k = Tr(ρ·(L + σ)^k ρ) read exactly in the Pauli-string
/// basis, shared by the F168 witnesses (<see cref="CatTurnedFieldWitness"/>,
/// <see cref="CatTwistedCouplingWitness"/>).
///
/// <para>The state is held as its amounts a_P = Tr(P·ρ) on the 4^N strings, read from the masks of
/// <see cref="PauliMask"/> (P|b⟩ = i^{|x∧z|}·(−1)^{|z∧b|}·|b ⊕ x⟩); the generator acts word by word:
/// the commutator with a string T is the factor ±2 on the string T·P where the two anticommute
/// (<see cref="PauliMask.Multiply"/>, its phase power odd), the Z watching the charge −2γ_l at every
/// site where P carries X or Y. Every coefficient is rational, so every number is exact, and
/// m_k = 2^{−N}·Σ_P a_P·((L + σ)^k a)_P.</para></summary>
internal static class CatPauliMoments
{
    public static int Index(int n, PauliMask p) => (int)(p.X << n | p.Z);

    /// <summary>a_P = Tr(P·ρ) for the cat (|0…0⟩ + e^{iθ}|1…1⟩)/√2, (cos θ, sin θ) given.</summary>
    public static BigRational[] CatAmounts(int n, BigRational cosTheta, BigRational sinTheta)
    {
        int dim = 1 << (2 * n);
        ulong all = (1UL << n) - 1;
        var half = new BigRational(1, 2);
        // ⟨b|ρ|b'⟩ on the two corners, as (real, imaginary)
        (BigRational Re, BigRational Im) Rho(ulong b, ulong bPrime) =>
            b == bPrime ? (half, BigRational.Zero)
            : b == all ? (half * cosTheta, half * sinTheta)          // ⟨1…1|ρ|0…0⟩ = e^{iθ}/2
            : (half * cosTheta, -half * sinTheta);                   // ⟨0…0|ρ|1…1⟩ = e^{−iθ}/2
        var a = new BigRational[dim];
        for (int i = 0; i < dim; i++) a[i] = BigRational.Zero;
        for (ulong x = 0; x <= all; x++)
            for (ulong z = 0; z <= all; z++)
            {
                // Tr(Pρ) = Σ_b ⟨b ⊕ x|P|b⟩·⟨b|ρ|b ⊕ x⟩ over b, b ⊕ x ∈ {0…0, 1…1}
                int phase = BitOperations.PopCount(x & z) % 4;          // the power of i
                BigRational re = BigRational.Zero, im = BigRational.Zero;
                foreach (ulong b in new[] { 0UL, all })
                {
                    ulong target = b ^ x;
                    if (target != 0 && target != all) continue;
                    int sign = (BitOperations.PopCount(z & b) & 1) == 0 ? 1 : -1;
                    int r = phase == 0 ? 1 : phase == 2 ? -1 : 0, s = phase == 1 ? 1 : phase == 3 ? -1 : 0;
                    var (pr, pi) = Rho(b, target);
                    re += sign * (r * pr - s * pi);
                    im += sign * (r * pi + s * pr);
                }
                if (!im.IsZero) throw new InvalidOperationException("a Hermitian state has real Pauli amounts");
                a[Index(n, new PauliMask(x, z))] = re;
            }
        return a;
    }

    /// <summary>(L + σ)·v for H = Σ coefficient·T and Z watching at the given rates.</summary>
    public static BigRational[] Apply(int n, BigRational[] v, IReadOnlyList<(PauliMask T, BigRational Coefficient)> terms,
                                      IReadOnlyList<BigRational> rates, BigRational sigma)
    {
        int dim = v.Length;
        var w = new BigRational[dim];
        for (int i = 0; i < dim; i++) w[i] = BigRational.Zero;
        ulong all = (1UL << n) - 1;
        for (ulong x = 0; x <= all; x++)
            for (ulong z = 0; z <= all; z++)
            {
                var p = new PauliMask(x, z);
                var amount = v[Index(n, p)];
                if (amount.IsZero) continue;
                BigRational charge = sigma;
                for (int l = 0; l < n; l++)
                    if (((x >> l) & 1) == 1) charge -= 2 * rates[l];
                w[Index(n, p)] += charge * amount;
                foreach (var (t, coefficient) in terms)
                {
                    if (coefficient.IsZero || PauliMask.Commute(t, p)) continue;
                    var (q, k) = PauliMask.Multiply(t, p);                 // t·p = i^k·q, k odd
                    // −i[t, p] = −2i·i^k·q: +2 at k = 1, −2 at k = 3
                    BigRational factor = k == 1 ? new BigRational(2) : new BigRational(-2);
                    w[Index(n, q)] += factor * coefficient * amount;
                }
            }
        return w;
    }

    /// <summary>The odd centred moments m_1, m_3, …, up to kMax, of the state with amounts a.</summary>
    public static IReadOnlyDictionary<int, BigRational> OddMoments(int n, BigRational[] a,
        IReadOnlyList<(PauliMask T, BigRational Coefficient)> terms, IReadOnlyList<BigRational> rates, int kMax)
    {
        BigRational sigma = BigRational.Zero;
        for (int l = 0; l < n; l++) sigma += rates[l];
        var v = a;
        var scale = new BigRational(BigInteger.One, BigInteger.One << n);
        var moments = new SortedDictionary<int, BigRational>();
        for (int k = 1; k <= kMax; k++)
        {
            v = Apply(n, v, terms, rates, sigma);
            if (k % 2 == 0) continue;
            BigRational m = BigRational.Zero;
            for (int i = 0; i < a.Length; i++)
                if (!a[i].IsZero && !v[i].IsZero) m += a[i] * v[i];
            moments[k] = m * scale;
        }
        return moments;
    }
}
