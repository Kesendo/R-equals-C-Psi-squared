using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F168 (<c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>): the cat
/// (|0…0⟩ + |1…1⟩)/√2 under Z dephasing, a field in the equator plane at one angle and a coupling
/// commuting with X^⊗N and Σ Z_l; its centred odd moments m_k = Tr(ρ·(L + σ)^k ρ) read exactly.
///
/// <para>The page's gate builds density matrices in sympy and multiplies them. This witness takes
/// the other road: the state as its amounts a_P = Tr(P·ρ) on the 4^N Pauli strings, read from the
/// masks of <see cref="PauliMask"/> (P|b⟩ = i^{|x∧z|}·(−1)^{|z∧b|}·|b ⊕ x⟩), and the generator applied
/// word by word: the commutator with a string T is the factor ±2 on the string T·P where the two
/// anticommute (<see cref="PauliMask.Multiply"/>, its phase power odd), the watching the charge
/// −2γ_l at every site where P carries X or Y. Every coefficient is rational, so every number is
/// exact (<see cref="BigRational"/>), and m_k = 2^{−N}·Σ_P a_P·((L + σ)^k a)_P.</para>
///
/// <para>The setting is wider than the gate's: unequal rates, unequal field strengths, unequal
/// bond strengths J_b and anisotropies Δ_b, on the chain or the ring, and from N = 3 a three-body
/// term (X Y Z − Y X Z) on the first three sites, which keeps the excitation number and commutes
/// with X^⊗N while carrying an odd count of Y and of Z each, so it lies inside F168's class and
/// outside the narrower one of couplings that keep the whole mirror.</para>
///
/// <para>Children: the vanishing (odd m_k = 0 for k &lt; 2N + 1, the assert); the law (m_{2N+1}
/// against the claim's <see cref="CatTurnedFieldSelfEncounterClaim.LeadingOddMoment"/>, nonzero);
/// the strengths' control (the law with h^{2N} from one strength would differ); the zero point (the
/// field along Y at even N, sin(Nφ) = 0 off the trivial angles, where every computed odd moment
/// vanishes; at odd N no rational point of the circle other than ±X has sin(Nφ) = 0, and the field
/// along X is read instead).</para>
///
/// <para>Guard: N ∈ {1, …, 5} (4^N ≤ 1024 words). Typed home
/// <see cref="CatTurnedFieldSelfEncounterClaim"/>.</para></summary>
public sealed class CatTurnedFieldWitness : IInspectable
{
    public const int MaxN = 5;

    // no prefix of the squared strengths multiplies to 1, so the one-strength control can fail at every N
    private static readonly BigRational[] Strengths = { R(1), R(2, 3), R(5, 4), R(3, 2), R(5, 3) };
    private static readonly BigRational[] Rates = { R(1, 5), R(2, 7), R(3, 11), R(1, 2), R(2, 9) };
    private static readonly BigRational[] Couplings = { R(1), R(2, 3), R(5, 4), R(3, 7), R(1, 2) };
    private static readonly BigRational[] Anisotropies = { R(1, 2), R(3), R(-1), R(2, 5), R(7, 4) };
    private static readonly BigRational ThreeBody = R(1, 2);

    public int N { get; }
    public bool Ring { get; }
    public (BigRational C, BigRational S) Angle { get; } = (R(3, 5), R(4, 5));

    /// <summary>The odd centred moments m_1, m_3, …, m_{2N+3} at the setting's angle.</summary>
    public IReadOnlyDictionary<int, BigRational> OddMoments { get; }
    /// <summary>The claim's closed form at this setting.</summary>
    public BigRational LawValue { get; }
    /// <summary>The closed form with one strength h_0 for every bit, the control.</summary>
    public BigRational OneStrengthLawValue { get; }
    /// <summary>The angle of the zero-point reading, (0, 1) at even N, (1, 0) at odd N.</summary>
    public (BigRational C, BigRational S) ZeroAngle { get; }
    /// <summary>The odd moments m_1 … m_{2N+3} at the zero angle.</summary>
    public IReadOnlyDictionary<int, BigRational> ZeroAngleOddMoments { get; }

    public CatTurnedFieldWitness(int n = 3, bool ring = false)
    {
        if (n < 1 || n > MaxN) throw new ArgumentOutOfRangeException(nameof(n), n, $"N must lie in 1..{MaxN}");
        if (ring && n < 3) throw new ArgumentOutOfRangeException(nameof(ring), "a ring needs N ≥ 3");
        N = n;
        Ring = ring;
        int kMax = 2 * n + 3;
        OddMoments = Moments(Angle, kMax);
        BigRational sigma = BigRational.Zero, prod = BigRational.One;
        for (int l = 0; l < n; l++) { sigma += Rates[l]; prod *= Strengths[l] * Strengths[l]; }
        var sinN = CatTurnedFieldSelfEncounterClaim.SinOfMultiple(n, Angle.C, Angle.S);
        LawValue = CatTurnedFieldSelfEncounterClaim.LeadingOddMoment(n, sigma, prod, sinN * sinN);
        BigRational oneProd = BigRational.One;
        for (int l = 0; l < n; l++) oneProd *= Strengths[0] * Strengths[0];
        OneStrengthLawValue = CatTurnedFieldSelfEncounterClaim.LeadingOddMoment(n, sigma, oneProd, sinN * sinN);
        ZeroAngle = n % 2 == 0 ? (R(0), R(1)) : (R(1), R(0));
        ZeroAngleOddMoments = Moments(ZeroAngle, kMax);
    }

    private static BigRational R(long p, long q = 1) => new(new BigInteger(p), new BigInteger(q));

    private List<(PauliMask T, BigRational Coefficient)> Terms((BigRational C, BigRational S) angle)
    {
        var terms = new List<(PauliMask, BigRational)>();
        for (int l = 0; l < N; l++)
        {
            ulong bit = 1UL << l;
            if (!angle.C.IsZero) terms.Add((new PauliMask(bit, 0), Strengths[l] * angle.C));
            if (!angle.S.IsZero) terms.Add((new PauliMask(bit, bit), Strengths[l] * angle.S));
        }
        int bonds = Ring ? N : N - 1;
        for (int b = 0; b < bonds; b++)
        {
            ulong pair = (1UL << b) | (1UL << ((b + 1) % N));
            terms.Add((new PauliMask(pair, 0), Couplings[b]));
            terms.Add((new PauliMask(pair, pair), Couplings[b]));
            terms.Add((new PauliMask(0, pair), Anisotropies[b]));
        }
        if (N >= 3)
        {
            terms.Add((PauliMask.Parse("XYZ" + new string('I', N - 3)), ThreeBody));
            terms.Add((PauliMask.Parse("YXZ" + new string('I', N - 3)), -ThreeBody));
        }
        return terms;
    }

    private int Index(PauliMask p) => (int)(p.X << N | p.Z);

    /// <summary>a_P = Tr(P·ρ) for the cat (|0…0⟩ + |1…1⟩)/√2, from the masks.</summary>
    private BigRational[] CatAmounts()
    {
        int dim = 1 << (2 * N);
        ulong all = (1UL << N) - 1;
        var a = new BigRational[dim];
        for (int i = 0; i < dim; i++) a[i] = BigRational.Zero;
        for (ulong x = 0; x <= all; x++)
            for (ulong z = 0; z <= all; z++)
            {
                // ⟨b'|P|b⟩ over b, b' ∈ {0, all}: P|b⟩ = i^{|x∧z|}·(−1)^{|z∧b|}·|b ⊕ x⟩, summed and halved
                int phase = BitOperations.PopCount(x & z) % 4;      // the power of i
                int re = 0, im = 0;
                foreach (ulong b in new[] { 0UL, all })
                {
                    ulong target = b ^ x;
                    if (target != 0 && target != all) continue;
                    int sign = (BitOperations.PopCount(z & b) & 1) == 0 ? 1 : -1;
                    int p = phase;
                    int r = p == 0 ? 1 : p == 2 ? -1 : 0, s = p == 1 ? 1 : p == 3 ? -1 : 0;
                    re += sign * r;
                    im += sign * s;
                }
                if (im != 0) throw new InvalidOperationException("a Hermitian state has real Pauli amounts");
                a[Index(new PauliMask(x, z))] = new BigRational(re) / new BigRational(2);
            }
        return a;
    }

    private BigRational[] Apply(BigRational[] v, List<(PauliMask T, BigRational Coefficient)> terms, BigRational sigma)
    {
        int dim = v.Length;
        var w = new BigRational[dim];
        for (int i = 0; i < dim; i++) w[i] = BigRational.Zero;
        ulong all = (1UL << N) - 1;
        for (ulong x = 0; x <= all; x++)
            for (ulong z = 0; z <= all; z++)
            {
                var p = new PauliMask(x, z);
                var amount = v[Index(p)];
                if (amount.IsZero) continue;
                BigRational charge = sigma;
                for (int l = 0; l < N; l++)
                    if (((x >> l) & 1) == 1) charge -= 2 * Rates[l];
                w[Index(p)] += charge * amount;
                foreach (var (t, coefficient) in terms)
                {
                    if (PauliMask.Commute(t, p)) continue;
                    var (q, k) = PauliMask.Multiply(t, p);                 // t·p = i^k·q, k odd
                    // −i[t, p] = −2i·i^k·q: +2 at k = 1, −2 at k = 3
                    BigRational factor = k == 1 ? new BigRational(2) : new BigRational(-2);
                    w[Index(q)] += factor * coefficient * amount;
                }
            }
        return w;
    }

    private IReadOnlyDictionary<int, BigRational> Moments((BigRational C, BigRational S) angle, int kMax)
    {
        var terms = Terms(angle);
        BigRational sigma = BigRational.Zero;
        for (int l = 0; l < N; l++) sigma += Rates[l];
        var a = CatAmounts();
        var v = a;
        var scale = new BigRational(BigInteger.One, BigInteger.One << N);
        var moments = new SortedDictionary<int, BigRational>();
        for (int k = 1; k <= kMax; k++)
        {
            v = Apply(v, terms, sigma);
            if (k % 2 == 0) continue;
            BigRational m = BigRational.Zero;
            for (int i = 0; i < a.Length; i++)
                if (!a[i].IsZero && !v[i].IsZero) m += a[i] * v[i];
            moments[k] = m * scale;
        }
        return moments;
    }

    public bool VanishingHolds => OddMoments.Where(kv => kv.Key < 2 * N + 1).All(kv => kv.Value.IsZero);
    public bool LawHolds => OddMoments[2 * N + 1] == LawValue && !LawValue.IsZero;
    public bool StrengthControlDiffers => N >= 2 && OneStrengthLawValue != OddMoments[2 * N + 1];
    public bool ZeroAngleHolds => ZeroAngleOddMoments.Values.All(m => m.IsZero);

    public string DisplayName => $"CatTurnedFieldWitness (N={N}, {(Ring ? "ring" : "chain")}, cos φ = 3/5)";

    public string Summary =>
        $"N={N}: the cat under a field at cos φ = 3/5 with unequal rates, strengths, bonds and anisotropies" +
        $"{(N >= 3 ? " and a three-body term" : "")}: the odd moments below order {2 * N + 1} {(VanishingHolds ? "vanish" : "do NOT vanish")}, " +
        $"m_{2 * N + 1} = {OddMoments[2 * N + 1]} {(LawHolds ? "is" : "is NOT")} the claim's (σ/3)(2N+1)!·∏h_l²·sin²(Nφ), " +
        $"one strength for every bit {(StrengthControlDiffers ? "would give another number" : "gives the same number, control void")}, " +
        $"and at the zero angle every computed odd moment {(ZeroAngleHolds ? "vanishes" : "does NOT vanish")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the vanishing (the assert)",
                summary: string.Join(", ", OddMoments.Where(kv => kv.Key < 2 * N + 1).Select(kv => $"m_{kv.Key} = {kv.Value}")) +
                         $": {(VanishingHolds ? "every odd moment below 2N + 1 is exactly 0" : "a moment below 2N + 1 is NOT 0")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the law",
                summary: $"m_{2 * N + 1} = {OddMoments[2 * N + 1]} read word by word against (σ/3)(2N+1)!·∏h_l²·sin²(Nφ) = {LawValue}: " +
                         $"{(LawHolds ? "equal, and not zero" : "DIFFERENT or zero")}; the next, m_{2 * N + 3} = {OddMoments[2 * N + 3]}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the strengths' control",
                summary: $"with one strength h_0 for every bit the law would read {OneStrengthLawValue}: " +
                         $"{(StrengthControlDiffers ? "different from the moment, so the gate sees the product ∏h_l²" : "the same, the control is void")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the zero point",
                summary: $"at (cos φ, sin φ) = ({ZeroAngle.C}, {ZeroAngle.S}) " +
                         $"{(N % 2 == 0 ? "(the field along Y, Nφ a multiple of π)" : "(the field along X; at odd N no other rational point of the circle has sin(Nφ) = 0)")}: " +
                         string.Join(", ", ZeroAngleOddMoments.Select(kv => $"m_{kv.Key} = {kv.Value}")) +
                         $" ({(ZeroAngleHolds ? "all zero, the property holds" : "NOT all zero")}).",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
