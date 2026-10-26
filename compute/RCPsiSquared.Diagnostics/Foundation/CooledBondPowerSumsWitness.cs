using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for <see cref="CooledBondPowerSumsClaim"/>
/// (<c>experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md</c>): a chain, or a ring, of one class bond under loss
/// alone, its odd power sums p_k = Tr((2L + c)^k) read exactly.
///
/// <para>The page's gate builds 2L + c modulo primes, row-major, by Kronecker products, and interpolates
/// polynomials in the rates. This witness takes another road: the generator's image of every matrix unit
/// |i⟩⟨j| written down directly (the commutator with H row by row and column by column, loss at site l
/// clearing bit l of both indices), vectorized column-major, entries Gaussian integers held exactly in
/// checked 128-bit arithmetic, no prime and no interpolation. Site 0 is the most significant bit.</para>
///
/// <para>The bond is not one of the gate's: 2XX − XY + 3YX + YY − 2ZZ (Δ = −2, S = 3, A = −4, W = 5, no two of them
/// ±1 or equal, so a wrong power of any factor misses).
/// Children: the low sums on the chain (p₁, p₃, p₅ exactly zero at a profile with unequal end products);
/// the seventh against the claim's closed form (equal, and not zero), and the reversed profile computed from the
/// generator (−p₇, the page's p_k(Pγ) = −p_k(γ)); the equal end products (a profile
/// with γ₁γ₂ = γ_{N−1}γ_N, moved by the reversal from N = 4, where the closed form and the sum are both
/// zero); the control (a bond with A = 0, XX + 2XY + 2YX + ZZ, whose seventh is zero at the breaking
/// profile); the ring (p₁ … p₇ zero on the ring of N, p₉ not).</para>
///
/// <para>Guard: N ∈ {3, 4, 5} (4^N ≤ 1024 matrix units).</para></summary>
public sealed class CooledBondPowerSumsWitness : IInspectable
{
    public const int MaxN = 5;

    public static readonly CooledBondPowerSumsClaim.Bond MainBond = new(2, -1, 3, 1, -2);
    public static readonly CooledBondPowerSumsClaim.Bond ControlBond = new(1, 2, 2, 0, 1);
    private static readonly long[] BreakingRates = { 2, 1, 3, 5, 4 };
    private static readonly long[][] EqualEndRates = { new long[] { 2, 5, 2 }, new long[] { 1, 6, 2, 3 }, new long[] { 1, 6, 5, 2, 3 } };

    public int N { get; }
    public IReadOnlyList<long> Breaking { get; }
    public IReadOnlyList<long> EqualEnds { get; }

    /// <summary>p_1 … p_7 on the chain at the breaking profile.</summary>
    public IReadOnlyDictionary<int, BigInteger> ChainSums { get; }
    /// <summary>p_7 on the chain at the breaking profile reversed, from the generator.</summary>
    public BigInteger ReversedSeventh { get; }
    /// <summary>The claim's closed form at the breaking profile.</summary>
    public BigInteger Law { get; }
    /// <summary>p_7 on the chain at the equal-end profile, and the closed form there.</summary>
    public BigInteger EqualEndsSeventh { get; }
    public BigInteger EqualEndsLaw { get; }
    /// <summary>p_7 of the control bond (A = 0) at the breaking profile.</summary>
    public BigInteger ControlSeventh { get; }
    /// <summary>p_1 … p_9 on the ring at the breaking profile.</summary>
    public IReadOnlyDictionary<int, BigInteger> RingSums { get; }

    public CooledBondPowerSumsWitness(int n = 4)
    {
        if (n < 3 || n > MaxN) throw new ArgumentOutOfRangeException(nameof(n), n, $"N must lie in 3..{MaxN}");
        N = n;
        Breaking = BreakingRates.Take(n).ToArray();
        EqualEnds = EqualEndRates[n - 3];
        ChainSums = PowerSums(MainBond, Breaking, ring: false, kMax: 7);
        Law = CooledBondPowerSumsClaim.SeventhPowerSum(MainBond, Breaking);
        ReversedSeventh = PowerSums(MainBond, Breaking.Reverse().ToArray(), ring: false, kMax: 7)[7];
        EqualEndsSeventh = PowerSums(MainBond, EqualEnds, ring: false, kMax: 7)[7];
        EqualEndsLaw = CooledBondPowerSumsClaim.SeventhPowerSum(MainBond, EqualEnds);
        ControlSeventh = PowerSums(ControlBond, Breaking, ring: false, kMax: 7)[7];
        RingSums = PowerSums(MainBond, Breaking, ring: true, kMax: 9);
    }

    // ---- exact Gaussian integers ----

    private readonly record struct G(Int128 Re, Int128 Im)
    {
        public static readonly G Zero = new(0, 0);
        public bool IsZero => Re == 0 && Im == 0;
        public static G operator +(G a, G b) => checked(new G(a.Re + b.Re, a.Im + b.Im));
        public static G operator *(G a, G b) => checked(new G(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re));
    }

    /// <summary>The Hamiltonian, d × d, every edge of the chain (and the wrap edge on a ring) carrying the bond.</summary>
    private static G[,] Hamiltonian(int n, CooledBondPowerSumsClaim.Bond bond, bool ring)
    {
        int d = 1 << n;
        var h = new G[d, d];
        var edges = Enumerable.Range(0, n - 1).Select(l => (l, l + 1)).ToList();
        if (ring) edges.Add((n - 1, 0));
        var terms = new[] { ('X', 'X', bond.XX), ('X', 'Y', bond.XY), ('Y', 'X', bond.YX), ('Y', 'Y', bond.YY), ('Z', 'Z', bond.ZZ) };
        foreach (var (a, b) in edges)
            foreach (var (pa, pb, c) in terms)
            {
                if (c == 0) continue;
                for (int col = 0; col < d; col++)
                {
                    var (row1, f1) = Letter(pa, col, n, a);
                    var (row, f2) = Letter(pb, row1, n, b);
                    h[row, col] = h[row, col] + new G(c, 0) * f1 * f2;
                }
            }
        return h;
    }

    /// <summary>P on site l applied to the basis state x: the image and the factor, ⟨image|P|x⟩.</summary>
    private static (int Image, G Factor) Letter(char p, int x, int n, int l)
    {
        int mask = 1 << (n - 1 - l);
        bool one = (x & mask) != 0;
        return p switch
        {
            'X' => (x ^ mask, new G(1, 0)),
            'Y' => (x ^ mask, one ? new G(0, -1) : new G(0, 1)),     // Y|0⟩ = i|1⟩, Y|1⟩ = −i|0⟩
            _ => (x, one ? new G(-1, 0) : new G(1, 0)),
        };
    }

    /// <summary>The columns of 2L + c, column-major (|i⟩⟨j| at i + d·j), sparse.</summary>
    private static List<(int Row, G Value)>[] Generator(int n, CooledBondPowerSumsClaim.Bond bond, IReadOnlyList<long> rates, bool ring)
    {
        int d = 1 << n, dim = d * d;
        var h = Hamiltonian(n, bond, ring);
        var cols = new List<(int, G)>[dim];
        var acc = new Dictionary<int, G>();
        var minusTwoI = new G(0, -2);
        var plusTwoI = new G(0, 2);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
            {
                acc.Clear();
                void Add(int r, int c, G v) { int k = r + d * c; acc[k] = acc.TryGetValue(k, out var o) ? o + v : v; }
                // −2i·(H|i⟩⟨j| − |i⟩⟨j|H)
                for (int a = 0; a < d; a++) if (!h[a, i].IsZero) Add(a, j, minusTwoI * h[a, i]);
                for (int b = 0; b < d; b++) if (!h[j, b].IsZero) Add(i, b, plusTwoI * h[j, b]);
                // loss at site l: 2γ·σ⁻|i⟩⟨j|σ⁺ − γ·(n_l|i⟩⟨j| + |i⟩⟨j|n_l) + γ (the c)
                for (int l = 0; l < n; l++)
                {
                    long g = rates[l];
                    if (g == 0) continue;
                    int mask = 1 << (n - 1 - l);
                    bool bi = (i & mask) != 0, bj = (j & mask) != 0;
                    if (bi && bj) Add(i ^ mask, j ^ mask, new G(2 * g, 0));
                    Add(i, j, new G(g * (1 - (bi ? 1 : 0) - (bj ? 1 : 0)), 0));
                }
                cols[i + d * j] = acc.Where(kv => !kv.Value.IsZero).Select(kv => (kv.Key, kv.Value)).ToList();
            }
        return cols;
    }

    /// <summary>Tr((2L + c)^k) for k = 1 … kMax, exact; the imaginary part of every trace is checked to be zero.</summary>
    public static IReadOnlyDictionary<int, BigInteger> PowerSums(CooledBondPowerSumsClaim.Bond bond, IReadOnlyList<long> rates, bool ring, int kMax)
    {
        int n = rates.Count;
        var cols = Generator(n, bond, rates, ring);
        int dim = cols.Length;
        var p = new G[dim, dim];
        for (int c = 0; c < dim; c++) foreach (var (r, v) in cols[c]) p[r, c] = p[r, c] + v;
        var sums = new SortedDictionary<int, BigInteger>();
        for (int k = 1; k <= kMax; k++)
        {
            if (k > 1)
            {
                var next = new G[dim, dim];
                for (int c = 0; c < dim; c++)
                    foreach (var (m, v) in cols[c])
                        for (int r = 0; r < dim; r++)
                            if (!p[r, m].IsZero) next[r, c] = next[r, c] + p[r, m] * v;
                p = next;
            }
            var tr = G.Zero;
            for (int i = 0; i < dim; i++) tr = tr + p[i, i];
            if (tr.Im != 0) throw new InvalidOperationException($"Tr((2L + c)^{k}) has an imaginary part; the generator is not real");
            sums[k] = (BigInteger)tr.Re;
        }
        return sums;
    }

    public bool ChainLowVanish => ChainSums[1].IsZero && ChainSums[3].IsZero && ChainSums[5].IsZero;
    public bool SeventhIsLaw => ChainSums[7] == Law && !Law.IsZero;
    public bool ReversalFlipsSeventh => ReversedSeventh == -ChainSums[7] && !ReversedSeventh.IsZero;
    public bool EqualEndsSilent => EqualEndsSeventh.IsZero && EqualEndsLaw.IsZero;
    public bool EqualEndsMovedByReversal => !EqualEnds.SequenceEqual(EqualEnds.Reverse());
    public bool ControlSilent => ControlSeventh.IsZero && CooledBondPowerSumsClaim.SeventhPowerSum(ControlBond, Breaking).IsZero;
    public bool RingSilentThroughSeventh => RingSums[1].IsZero && RingSums[3].IsZero && RingSums[5].IsZero && RingSums[7].IsZero;
    public bool RingNinthBreaks => !RingSums[9].IsZero;

    public string DisplayName => $"CooledBondPowerSumsWitness (N={N}, {MainBond} under loss)";

    public string Summary =>
        $"N={N}, {MainBond} on every edge, loss at ({string.Join(", ", Breaking)}): on the chain p₁, p₃, p₅ " +
        $"{(ChainLowVanish ? "are exactly 0" : "are NOT all 0")}, p₇ = {ChainSums[7]} {(SeventhIsLaw ? "is" : "is NOT")} the claim's " +
        $"closed form; at ({string.Join(", ", EqualEnds)}), equal end products, p₇ = {EqualEndsSeventh} " +
        $"{(EqualEndsSilent ? "and the closed form are 0" : "or the closed form is NOT 0")}; the control with A = 0 " +
        $"{(ControlSilent ? "is silent" : "is NOT silent")}; on the ring p₁ … p₇ {(RingSilentThroughSeventh ? "vanish" : "do NOT all vanish")} " +
        $"and p₉ = {RingSums[9]} {(RingNinthBreaks ? "does not" : "DOES")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the low sums on the chain (the assert)",
                summary: string.Join(", ", new[] { 1, 3, 5 }.Select(k => $"p_{k} = {ChainSums[k]}")) +
                         $" at ({string.Join(", ", Breaking)}), where γ₁γ₂ ≠ γ_{{N−1}}γ_N: {(ChainLowVanish ? "every one exactly 0" : "NOT all 0")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the seventh against the closed form",
                summary: $"p₇ = {ChainSums[7]}, the claim's −7·2¹⁷·4^(N−3)·Δ·S·A·W·(γ₁γ₂ − γ_{{N−1}}γ_N) = {Law} " +
                         $"(Δ = {Signed(MainBond.ZZ)}, S = {Signed(MainBond.S)}, A = {Signed(MainBond.A)}, W = {Signed(MainBond.W)}): {(SeventhIsLaw ? "equal, and not zero" : "DIFFERENT or zero")}; " +
                         $"at the reversed profile the generator gives {ReversedSeventh}: {(ReversalFlipsSeventh ? "−p₇, the sign the reversal flips" : "NOT −p₇")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the equal end products",
                summary: $"at ({string.Join(", ", EqualEnds)}) " +
                         $"{(EqualEndsMovedByReversal ? "(the reversal moves the profile)" : "(the reversal keeps the profile)")}: " +
                         $"p₇ = {EqualEndsSeventh}, the closed form {EqualEndsLaw}: {(EqualEndsSilent ? "both 0, the seventh decides nothing here" : "NOT both 0")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the control",
                summary: $"{ControlBond} (A = {ControlBond.A}) at ({string.Join(", ", Breaking)}): p₇ = {ControlSeventh}: " +
                         $"{(ControlSilent ? "zero, as the factor A says" : "NOT zero")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the ring",
                summary: $"the ring of {N}, {MainBond} on every edge, loss at ({string.Join(", ", Breaking)}): " +
                         string.Join(", ", RingSums.Select(kv => $"p_{kv.Key} = {kv.Value}").Where((_, i) => i % 2 == 0)) +
                         $": {(RingSilentThroughSeventh ? "the odd sums through the seventh vanish" : "an odd sum through the seventh is NOT 0")}, " +
                         $"{(RingNinthBreaks ? "the ninth does not" : "the ninth is 0 too, the reading is void")}.",
                provenance: NodeProvenance.Live);
        }
    }

    private static string Signed(long v) => v < 0 ? $"−{-v}" : v.ToString();

    public InspectablePayload Payload => InspectablePayload.Empty;
}
