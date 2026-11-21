using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for F173 (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>, typed
/// <see cref="SignedMaskBalanceClaim"/>): the claim's readings, GF(2) on masks and signs and, where a flip set carries several
/// strings, Gaussian integers on its ratio and read set, held against
/// Theorem 2 of the complement connection walked on the bitstrings (<see cref="ComplementConnectionGraph"/>, every entry
/// a Gaussian integer, every edge cross-multiplied), Z dephasing on every site, templates placed with their identity
/// letters at every start of the open chain or at every start of the ring, wrapped; Theorem D also against the balance
/// (<see cref="SignedMaskBalanceClaim.Balanced"/>, Theorem A's elimination) on rings to N = 32, past the bitstrings' reach.
///
/// <para>Families, each a count of cells and readings and the readings where the two routes disagree (none expected):
/// two templates of width three with masks of different shapes, both signs of the second coefficient, at N* − 1, N*,
/// N* + 1 where c = 1 and at N = 3 to 6 where c = 0; two templates of unequal widths, where N* pairs each width with the
/// other template's quotient; the diagonal-cell masks at k = 3 and 4 against F115's own valuation
/// (<see cref="WindowedObstructionScan.IsHardPair"/>) and against Theorem 2 at N* − 1 and N*; one template on the ring
/// from N = 2w − 1 to the ring cap; XY + ZZ on the ring; Theorem D on two templates against Theorem 2 to the ring cap,
/// and on two templates up to five sites wide and three templates against the balance to N = 32; three templates at one
/// N against Theorem 2; Theorem E on two strings on one flip mask, chain and ring, Corollary E2 on the ring (also against
/// Theorem E's elimination past the bitstrings' reach) and on the products (aX + bY)⊗(cX + dY), a, b, c, d = ±1,
/// Corollary E3 on two flip shapes on the chain; Theorem B′ on three and four templates on the chain against the balance
/// and Theorem 2, and at two templates against Theorem B; Theorem E's ℤ/4 form, the refinement and its carrier against
/// Theorem 2's (★) and as a diagonal Clifford; Corollary E3′ against E3 at two templates and against Theorem 2 at three; the long chain against
/// Theorem D on the ring of 64; Theorem H at random rates, λ_min of the far edge's second-order generator
/// (<see cref="StrongDephasingEdges"/>, built from H's entries) against the weighted code distance
/// (<see cref="SignedMaskBalanceClaim.StrongDephasingBreak"/>), and the kernels of both second-order generators against
/// Theorem 2's two counts; Theorem H by the convolutional code's trellis against the enumeration, and Corollary H1's two
/// bounds, from below by disjoint translates of Theorem B's relation on the chain and from above by a cut on the ring
/// where the chain pairs; Corollary H2's cycle bound on Theorem E's cells and as the packing of XZ ± ZY's end squares;
/// Corollary H3's hop reading (<see cref="StrongDephasingEdges.ReadHops"/>) on Theorem E's cells, the whole far spectrum
/// coset by coset where the hops commute and the floor where every two anticommute; Corollary A″, ket and bra agreeing on
/// some sites, the balance with their Y and Z letters struck against Theorem 2 with the complement taken on the sites outside them
/// (<see cref="ComplementConnectionGraph.AgreeOn"/>), and its hypotheses on rows with diagonal strings and with two strings
/// on one flip set; the triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J, two of its premises exactly (the three bond strings anticommuting in pairs, H² a
/// multiple of 1) and its closed form (<see cref="SignedMaskBalanceClaim.TriangleFarEdgeAtLargeJ"/>) against the rates
/// compressed onto the eigenspaces of the commutator with H (degenerate perturbation in 1/J, the rule of
/// THE_SPREAD_IS_A_RESONANCE); the sums over disjoint bonds on
/// the ring of bonds X_l·Y_(l+1) (<see cref="SignedMaskBalanceClaim.DisjointBondSum"/>), each and its product with Z^⊗N commuting with H
/// exactly, with H² = Σ c_l² + 2·Q₂ from four sites on and the top terms' commutators at unit coefficients, and on the odd ring the operators with N − 1 X or Y letters commuting with H counted exactly, 2N + 2, between the kernel
/// modulo a prime and Corollary A″'s count plus those two; the worked examples of the page.
/// Twenty-eight controls show the
/// battery can fail: the widths swapped in
/// N*, the quotients swapped in c and the power 2^(v₂(N)) dropped from the ring law each disagree with Theorem 2 somewhere,
/// Theorem D with the parity of N for 2^(v₂(N)) or with i₀ at the greatest valuation disagrees with the balance, and
/// Theorem E with any one of its three conditions dropped, Corollary E2 without its read-set test or its ordering term and
/// Corollary E3 without its read-sum test or its ordering term disagree with Theorem 2, and Theorem B′ read off the
/// unreduced basis of the relations disagrees with the balance, the ℤ/4 carrier without its coset term fails (★) and moved
/// by a quarter on one component stops being a Clifford exactly where the flip sets' span has codimension two or more, and
/// Corollary E3′ with N₁* read away or with hol walked without its ordering term disagrees with Theorem 2, and Theorem H
/// with every rate read as one or with the rates in reverse site order misses λ_min of the far edge's generator, and
/// Corollary H1's packing bound with every offset of the relation counted (overlapping translates) exceeds the trellis,
/// and its cut bound read where the chain breaks falls short of it, and Corollary H2's bound without its double count
/// exceeds λ_min of the far edge's generator, and Theorem H read string by string in place of Corollary H3's reading
/// misses it, and Corollary A″ with no letter struck, with the next site's letters struck for one site, or with the
/// diagonal strings ignored where one is odd off the agreeing sites, disagrees with Theorem 2's count, and the triangle's
/// closed form handed the rate of a site on each bond in place of the site opposite it, or its odd term alone without
/// Corollary A″'s single sites, disagrees with that compression.</para></summary>
public sealed class SignedMaskBalanceWitness : IInspectable
{
    /// <summary>The largest ring read (the graph has 2^N vertices).</summary>
    public const int MaxRing = 12;

    public sealed record FamilyReading(string Name, int Cells, int Readings, int Mismatches, string? FirstMismatch, string Against = "Theorem 2")
    {
        public bool Agrees => Mismatches == 0;
    }

    public int RingCap { get; }
    public IReadOnlyList<FamilyReading> Families { get; }
    public IReadOnlyDictionary<int, int> F115FloorsAtK4 { get; }
    public int SwappedWidthsMismatches { get; }
    public int SwappedQuotientsMismatches { get; }
    public int DroppedRingPowerMismatches { get; }
    public int ParityOnlyMismatches { get; }
    public int MostValuedMismatches { get; }
    /// <summary>Theorem E's readings that fall outside its setting (an absent edge or a ratio that is no character).</summary>
    public int CharacterOutside { get; }
    /// <summary>Theorem E with condition (1), (2) or (3) dropped: readings against Theorem 2.</summary>
    public IReadOnlyList<int> DroppedConditionMismatches { get; }
    /// <summary>Corollary E2 without its read-set test: readings against Theorem 2.</summary>
    public int DroppedReadSetMismatches { get; }
    /// <summary>Corollary E2 with hol(g) walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedOrderingMismatches { get; }
    /// <summary>Corollary E3 without its read-sum test: readings against Theorem 2.</summary>
    public int DroppedReadSumMismatches { get; }
    /// <summary>Corollary E3 with hol(u₀) walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedE3OrderingMismatches { get; }
    /// <summary>Theorem B′ read off the unreduced basis of the relations: readings against the balance.</summary>
    public int UnreducedBasisMismatches { get; }
    /// <summary>Theorem E's ℤ/4 carrier without its coset term 2·⟨R(w), x₀⟩: paired readings where it fails (★).</summary>
    public int DroppedCosetTermMismatches { get; }
    /// <summary>Corollary E3′ with N₁* read away (hol on all of S): readings against Theorem 2.</summary>
    public int DroppedE3PrimeReadSumMismatches { get; }
    /// <summary>Corollary E3′ with hol walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedE3PrimeOrderingMismatches { get; }
    /// <summary>How often each branch of Corollary E3′ decided a three-shape chain read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.ChainShapeBranch, int> E3PrimeBranches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
    /// <summary>How many of the ℤ/4 readings pair, each carrier checked against (★).</summary>
    public int Z4Paired { get; private set; }
    /// <summary>How many readings of the four-string ℤ/4 family pair, each carrier checked against (★).</summary>
    public int Z4FourStringPaired { get; private set; }
    /// <summary>The ℤ/4 carriers of both families moved by one quarter on one component off the span W of the flip sets
    /// (a phase per component, so (★) still holds): paired readings where the moved quarters are no longer a ℤ/4 quadratic
    /// form on the bitstrings, which must be exactly those where W has codimension two or more.</summary>
    public int MovedCarrierNotQuadratic { get; private set; }
    /// <summary>How often each branch of Corollary E3 decided a chain read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.ChainShapeBranch, int> E3Branches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
    /// <summary>How often each branch of Corollary E2 decided a ring read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.RingShapeBranch, int> E2Branches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.RingShapeBranch, int>();
    /// <summary>How many of Theorem B′'s drawn cells break for large N.</summary>
    public int BPrimeBreaking { get; private set; }
    /// <summary>How many of those cells have an unreduced basis whose leading coefficient matrix falls short of rank k − 1
    /// (the rank test can fail).</summary>
    public int UnreducedBasisDeficient { get; private set; }
    /// <summary>Theorem H with the rates read away (each string weighted by c² alone): readings where the code distance
    /// then misses λ_min(𝓜) by more than the eigensolver's bound.</summary>
    public int RatesReadAwayMismatches { get; private set; }
    /// <summary>Theorem H with the rates handed over in reverse site order (γ_(N−1−l) on site l, the convention error a
    /// builder reading site l as bit N − 1 − l makes): readings where the code distance then misses λ_min(𝓜).</summary>
    public int RatesReversedMismatches { get; private set; }
    /// <summary>Theorem H: the largest |λ_min(𝓜) − code distance| over dim·eps·‖𝓜‖ (a backward-stable Hermitian solver).</summary>
    public double StrongDephasingWorstRatio { get; private set; }
    /// <summary>Theorem H's lemma: the largest modulus of an eigenvalue counted as a zero of 𝓛 or 𝓜, negative ones included,
    /// over ‖·‖, and the smallest one not counted, over ‖·‖ (the count's threshold, 10⁻⁹·‖·‖, lies between them); and the
    /// most negative eigenvalue of either, over ‖·‖, which positive semidefiniteness puts at rounding.</summary>
    public double KernelZeroMax { get; private set; }
    public double MostNegative { get; private set; }
    public double KernelGapMin { get; private set; } = double.PositiveInfinity;
    /// <summary>Theorem H: readings where the far edge misses −2σ at second order (λ_min(𝓜) &gt; 0).</summary>
    public int StrongDephasingBroken { get; private set; }
    /// <summary>Corollary H2: the frustrated relations read, the readings (a cell at one N) with at least one, the readings
    /// where the packing meets λ_min(𝓜), and the relations where k·λ_*, without the double count, exceeds λ_min(𝓜) (the
    /// control).</summary>
    public int CycleRelations { get; private set; }
    public int CyclePackingReadings { get; private set; }
    public int CyclePackingEqualities { get; private set; }
    public int CycleDoubleCountMismatches { get; private set; }
    /// <summary>Corollary H2: the largest packing over λ_min(𝓜) on Theorem E's cells; on the triangle with a potential, the
    /// readings with moduli that differ (k &gt; 0) and those where (3/2)·λ_* exceeds 2·μ_min, which no cycle without a
    /// potential can.</summary>
    public double CycleTightest { get; private set; }
    public int CyclePotentialReadings { get; private set; }
    public int CyclePotentialAboveTwoHops { get; private set; }
    /// <summary>The triangle at k = 0, moduli that agree: readings where (3/2)·λ_* exceeds 2·μ_min (none may).</summary>
    public int CycleEqualModuliAboveTwoHops { get; private set; }
    /// <summary>Corollary H3 on Theorem E's cells: the readings where the hops commute and the spectrum is read, those
    /// broken at second order, those whose relation signs differ between cosets, those where Theorem H read string by string
    /// (the control) misses λ_min; the readings where the hops do not commute, those where every two anticommute; the
    /// readings outside (an entry vanishing somewhere).</summary>
    public int HopReadings { get; private set; }
    public int HopBroken { get; private set; }
    public int HopSplit { get; private set; }
    public int HopControlMisses { get; private set; }
    public int HopNoncommuting { get; private set; }
    public int HopAllAnticommute { get; private set; }
    public int HopOutside { get; private set; }
    /// <summary>Theorem H by the trellis: the largest |trellis − enumeration| over 2N·eps times the larger value.</summary>
    public double TrellisWorstRatio { get; private set; }
    /// <summary>Corollary H1's packing bound: readings where it is met exactly (to 4N·eps).</summary>
    public int PackingEqualities { get; private set; }
    /// <summary>The packing bound with every offset of the relation counted, the translates overlapping: readings where it
    /// exceeds the trellis (the control).</summary>
    public int PackingOverlapMismatches { get; private set; }
    /// <summary>Corollary H1's cut bound: ring readings where it is met exactly (to 4N·eps), the readings where the ring
    /// breaks (λ &gt; 0) with the chain pairing, and the largest λ over twice the cut there.</summary>
    public int CutEqualities { get; private set; }
    public int CutBroken { get; private set; }
    public double CutTightest { get; private set; }
    /// <summary>The cut bound read on the rings where the chain of the same N breaks: readings where the trellis exceeds it
    /// (the control).</summary>
    public int CutWhereChainBreaksMismatches { get; private set; }
    /// <summary>Corollary A″ (ket and bra agreeing on a set Λ of sites): readings where the operators X on every site outside
    /// Λ times a diagonal D commuting with H exist, where they do not, and where they exist while Theorem A's far end is empty
    /// (a broken row with a mode at −2·(σ − γ(Λ)) at every scale of H); the readings outside the corollary's hypotheses (two
    /// flipping strings on one flip set, or a diagonal string odd in Z off Λ), where its balance must decline; and three
    /// controls against Theorem 2's count (<see cref="ComplementConnectionGraph.AgreeOn"/>), which must disagree: the signs read with no letter struck, a single site read with the
    /// letters of the next site struck (on the ring translation makes the two verdicts equal, so this one fires on the chain
    /// alone), and the balance read with the diagonal strings ignored where one is odd off Λ.</summary>
    public int AgreeOnNonzero { get; private set; }
    public int AgreeOnZero { get; private set; }
    public int AgreeOnWhereFarIsZero { get; private set; }
    public int AgreeOnOutside { get; private set; }
    public int AgreeOnUnstruckMismatches { get; private set; }
    public int AgreeOnWrongSiteMismatches { get; private set; }
    public int AgreeOnDiagonalIgnoredMismatches { get; private set; }

    /// <summary>The triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J (the page's proposition): the profiles where the closed form's odd term is
    /// the smaller and those where its single-site term is; the worst distance of the closed form
    /// from the rates compressed onto the eigenspaces of the commutator with H, over 64·eps·2σ; and two controls against that
    /// compression, which must disagree: the
    /// closed form handed the rate of a site on each bond in place of the site opposite it, and its odd term alone, without
    /// Corollary A″'s single sites.</summary>
    public int TriangleOddWins { get; private set; }
    public int TriangleCorollaryWins { get; private set; }
    public double TriangleWorstRatio { get; private set; }
    public int TriangleRotatedMismatches { get; private set; }
    public int TriangleOddOnlyMismatches { get; private set; }

    /// <summary>The sums over disjoint bonds (the page's paragraph on one rate): the top terms on the odd rings, each of which
    /// fails to commute with H alone and moves its one site of I or Z by one place with the letter exchanged, I for Z and Z for I;
    /// and the top terms on the even rings, each a colouring that commutes with H alone.</summary>
    public int OddTopTermsMoving { get; private set; }
    public int EvenTopTermsCommuting { get; private set; }
    public bool AllAgree => Families.All(f => f.Agrees);
    public bool ControlsBreak => SwappedWidthsMismatches > 0 && SwappedQuotientsMismatches > 0 && DroppedRingPowerMismatches > 0
                                 && ParityOnlyMismatches > 0 && MostValuedMismatches > 0
                                 && DroppedConditionMismatches.All(d => d > 0) && DroppedReadSetMismatches > 0 && DroppedOrderingMismatches > 0
                                 && DroppedReadSumMismatches > 0 && DroppedE3OrderingMismatches > 0
                                 && UnreducedBasisMismatches > 0 && DroppedCosetTermMismatches > 0
                                 && DroppedE3PrimeReadSumMismatches > 0 && DroppedE3PrimeOrderingMismatches > 0
                                 && RatesReadAwayMismatches > 0 && RatesReversedMismatches > 0
                                 && PackingOverlapMismatches > 0 && CutWhereChainBreaksMismatches > 0 && CycleDoubleCountMismatches > 0
                                 && HopControlMisses > 0 && MovedCarrierNotQuadratic > 0
                                 && AgreeOnUnstruckMismatches > 0 && AgreeOnWrongSiteMismatches > 0 && AgreeOnDiagonalIgnoredMismatches > 0
                                 && TriangleRotatedMismatches > 0 && TriangleOddOnlyMismatches > 0;

    // the GF(2) rank of a set of masks
    static int Rank(int[] vectors)
    {
        var basis = new List<int>();
        foreach (int v0 in vectors)
        {
            int v = v0;
            foreach (int q in basis) v = Math.Min(v, v ^ q);
            if (v != 0) { basis.Add(v); basis.Sort((x, y) => y.CompareTo(x)); }
        }
        return basis.Count;
    }

    static readonly string[] Windows = (from a in "XYZ" from m in "IXYZ" from b in "XYZ" select $"{a}{m}{b}").ToArray();

    /// <summary>Theorem 2 on the open chain of n sites, each template placed with its identity letters at every start
    /// 0 … n − w.</summary>
    public static bool PairsOnChain(int n, params (string T, long C)[] terms)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var (t, c) in terms)
            for (int p = 0; p + t.Length <= n; p++)
                placed.Add((PauliMask.Parse(new string('I', p) + t + new string('I', n - p - t.Length)), c));
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, placed));
        return near == far;
    }

    /// <summary>Theorem 2 on the ring of n sites, each template at every start, wrapped (n ≥ its width).</summary>
    public static bool PairsOnRing(int n, params (string T, long C)[] terms)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var (t, c) in terms)
        {
            if (t.Length > n) throw new ArgumentException($"{t} does not fit on the ring of {n} sites");
            for (int p = 0; p < n; p++)
            {
                var letters = Enumerable.Repeat('I', n).ToArray();
                for (int i = 0; i < t.Length; i++) letters[(p + i) % n] = t[i];
                placed.Add((PauliMask.Parse(new string(letters)), c));
            }
        }
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, placed));
        return near == far;
    }

    // Corollary E3's templates of the given widths: one string with a letter at each end, or two different strings on one
    // mask whose letters together reach both ends (coefficients 1 and ±1, each unordered pair once), kept where the flip
    // set's ratio is a constant times a character
    static List<((string T, long C)[] Strings, SignedMaskBalanceClaim.FlipRatio Ratio)> OneShapeTemplates(params int[] widths)
    {
        var result = new List<((string T, long C)[], SignedMaskBalanceClaim.FlipRatio)>();
        foreach (int w in widths)
        {
            for (ulong f = 1; f < 1UL << w; f++)
            {
                var onMask = new List<string> { "" };
                for (int j = 0; j < w; j++) onMask = onMask.SelectMany(s => (((f >> j) & 1) == 1 ? "XY" : "IZ").Select(d => s + d)).ToList();
                var candidates = new List<(string T, long C)[]>();
                foreach (var t in onMask.Where(t => t[0] != 'I' && t[^1] != 'I')) candidates.Add(new[] { (t, 1L) });
                for (int a = 0; a < onMask.Count; a++)
                    for (int b = a + 1; b < onMask.Count; b++)
                    {
                        string s = onMask[a], t = onMask[b];
                        if ((s[0] == 'I' && t[0] == 'I') || (s[^1] == 'I' && t[^1] == 'I')) continue;
                        foreach (long c in new[] { 1L, -1L }) candidates.Add(new[] { (s, 1L), (t, c) });
                    }
                foreach (var cand in candidates)
                {
                    var r = SignedMaskBalanceClaim.ReadRatio(f, SignedMaskBalanceClaim.PlaceStrings(cand, w, ring: false));
                    if (r.Kind == SignedMaskBalanceClaim.RatioKind.Character) result.Add((cand, r));
                }
            }
        }
        return result;
    }

    static bool SameShape(ulong a, ulong b) => a >> BitOperations.TrailingZeroCount(a) == b >> BitOperations.TrailingZeroCount(b);

    // the ring's bond l: X on site l, Y on site l + 1 mod n
    static PauliMask RingBond(int n, int l) => new((1UL << l) | (1UL << ((l + 1) % n)), 1UL << ((l + 1) % n));

    // z·i^k
    static GaussianInteger IPower(GaussianInteger z, int k)
    {
        for (int t = 0; t < (k & 3); t++) z = new GaussianInteger(-z.Im, z.Re);
        return z;
    }

    static void AddTerm(Dictionary<PauliMask, GaussianInteger> op, PauliMask s, GaussianInteger c)
    {
        var sum = (op.TryGetValue(s, out var a) ? a : GaussianInteger.Zero) + c;
        if (sum == GaussianInteger.Zero) op.Remove(s); else op[s] = sum;
    }

    // an operator, Pauli strings with Gaussian-integer coefficients, times one string on the right
    static Dictionary<PauliMask, GaussianInteger> TimesString(Dictionary<PauliMask, GaussianInteger> op, PauliMask m)
    {
        var result = new Dictionary<PauliMask, GaussianInteger>();
        foreach (var (s, c) in op)
        {
            var (prod, k) = PauliMask.Multiply(s, m);
            AddTerm(result, prod, IPower(c, k));
        }
        return result;
    }

    // [H, op] exactly, H = Σ c_b·b: [b, s] = 2·b·s = 2·i^k·(b s) where b and s anticommute, zero where they commute
    static Dictionary<PauliMask, GaussianInteger> Commutator(IReadOnlyList<(PauliMask S, long C)> h, Dictionary<PauliMask, GaussianInteger> op)
    {
        var result = new Dictionary<PauliMask, GaussianInteger>();
        foreach (var (b, cb) in h)
            foreach (var (s, c) in op)
            {
                if (PauliMask.Commute(b, s)) continue;
                var (prod, k) = PauliMask.Multiply(b, s);
                AddTerm(result, prod, IPower(c * new GaussianInteger(2 * cb, BigInteger.Zero), k));
            }
        return result;
    }

    // the rank modulo p of sparse integer columns (row → value), eliminating on the least row present
    static int SparseRankModP(IReadOnlyList<Dictionary<int, long>> columns, long p)
    {
        static long Inverse(long a, long p)
        {
            long result = 1, b = a % p, e = p - 2;
            for (; e > 0; e >>= 1, b = b * b % p) if ((e & 1) == 1) result = result * b % p;
            return result;
        }
        var pivots = new Dictionary<int, Dictionary<int, long>>();
        foreach (var column in columns)
        {
            var v = new Dictionary<int, long>();
            foreach (var (r, x) in column) { long y = ((x % p) + p) % p; if (y != 0) v[r] = y; }
            while (v.Count > 0)
            {
                int r = v.Keys.Min();
                if (!pivots.TryGetValue(r, out var b)) { pivots[r] = v; break; }
                long f = v[r] * Inverse(b[r], p) % p;
                foreach (var (rr, x) in b)
                {
                    long y = ((v.TryGetValue(rr, out var cur) ? cur : 0) - f * x % p + p) % p;
                    if (y == 0) v.Remove(rr); else v[rr] = y;
                }
            }
        }
        return pivots.Count;
    }

    public SignedMaskBalanceWitness(int ringCap = 10)
    {
        if (ringCap < 5 || ringCap > MaxRing) throw new ArgumentOutOfRangeException(nameof(ringCap), ringCap, $"the ring cap runs from 5 to {MaxRing}");
        RingCap = ringCap;
        var families = new List<FamilyReading>();
        int swappedWidths = 0, swappedQuotients = 0, droppedPower = 0, parityOnly = 0, mostValued = 0, characterOutside = 0, droppedReadSet = 0, droppedOrdering = 0, droppedReadSum = 0, droppedE3Ordering = 0, unreducedBasis = 0, droppedCosetTerm = 0, droppedE3PrimeReadSum = 0, droppedE3PrimeOrdering = 0;
        var droppedCondition = new int[3];

        // two templates of width three, masks of different shapes, both signs
        {
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            for (int i = 0; i < Windows.Length; i++)
                for (int j = i + 1; j < Windows.Length; j++)
                {
                    string a = Windows[i], b = Windows[j];
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                    foreach (long sign in new[] { 1L, -1L })
                    {
                        cells++;
                        var ns = breaks ? new[] { from!.Value - 1, from.Value, from.Value + 1 }.Where(n => n >= 3) : Enumerable.Range(3, 4);
                        foreach (int n in ns)
                        {
                            bool expected = !breaks || n < from;
                            bool pairs = PairsOnChain(n, (a, 1), (b, sign));
                            bool balanced = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { a, b }, n));
                            bool character = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { (a, 1L), (b, sign) }, n, ring: false))
                                             == SignedMaskBalanceClaim.CharacterVerdict.Pairs;
                            readings++;
                            if (pairs != expected || balanced != expected || character != expected) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} at N = {n}"; }
                            if (pairs != (SwappedQuotientC(a, b) == 0 || n < NStar(a, b))) swappedQuotients++;
                        }
                    }
                }
            families.Add(new FamilyReading("two templates of width three, masks of different shapes, both signs (the balance and Theorem E agreeing)", cells, readings, mismatches, first));
        }

        // unequal widths
        {
            var two = from x in "XYZ" from y in "XYZ" select $"{x}{y}";
            var longer = new[] { "XIXIX", "XXIX", "YZIX", "XYIZY" };
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            foreach (var a in two.Concat(Windows))
                foreach (var b in Windows.Concat(longer).Where(b => b.Length > a.Length))
                {
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                    cells++;
                    // read only where both templates are placed, and where c = 0 from just below the first relation to
                    // just past it, so that every reading but the lowest holds a cycle the law could close or break
                    int star = NStar(a, b);
                    var ns = Enumerable.Range(star - 1, 3).Where(n => n >= Math.Max(3, b.Length));
                    foreach (int n in ns)
                    {
                        bool pairs = PairsOnChain(n, (a, 1), (b, 1));
                        readings++;
                        if (pairs != (!breaks || n < from)) { mismatches++; first ??= $"{a} + {b} at N = {n}"; }
                        if (breaks && pairs != (n < SwappedWidthsFrom(a, b))) swappedWidths++;
                    }
                }
            families.Add(new FamilyReading("two templates of unequal widths, placed with their identity letters", cells, readings, mismatches, first));
        }

        // the F115 face: diagonal-cell masks at k = 3, 4
        {
            var floors = new Dictionary<int, int>();
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            foreach (int k in new[] { 3, 4 })
            {
                var templates = new List<string>();
                for (ulong m = 1; m < 1UL << k; m++)
                {
                    if (BitOperations.PopCount(m) % 2 != 0) continue;
                    var letters = Enumerable.Range(0, k).Select(j => ((m >> j) & 1) == 1 ? 'X' : 'Z').ToArray();
                    if (SignedMaskBalanceClaim.Sign(new string(letters)) == 0) letters[BitOperations.TrailingZeroCount(m)] = 'Y';
                    templates.Add(new string(letters));
                }
                foreach (var a in templates)
                    foreach (var b in templates)
                    {
                        ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                        if (string.CompareOrdinal(a, b) >= 0 || SameShape(ma, mb)) continue;
                        cells++;
                        var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                        if (breaks != WindowedObstructionScan.IsHardPair(ma, mb)) { mismatches++; first ??= $"{a} + {b}: valuation"; }
                        if (!breaks) continue;
                        if (k == 4) floors[from!.Value] = floors.GetValueOrDefault(from.Value) + 1;
                        foreach (int n in new[] { from!.Value - 1, from.Value })
                        {
                            readings++;
                            if (PairsOnChain(n, (a, 1), (b, 1)) != (n < from)) { mismatches++; first ??= $"{a} + {b} at N = {n}"; }
                        }
                    }
            }
            F115FloorsAtK4 = floors;
            families.Add(new FamilyReading("the F115 face: diagonal-cell masks at k = 3, 4, the valuation and the floor", cells, readings, mismatches, first));
        }

        // one template on the ring, and XY + ZZ
        {
            var templates = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows)
                .Concat(new[] { "XYYY", "XIIY", "XXIY", "YIIY", "XIIIY", "YIIIY", "XYIYX" })
                .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            int readings = 0, mismatches = 0;
            string? first = null;
            foreach (var t in templates)
                for (int n = Math.Max(3, 2 * t.Length - 1); n <= RingCap; n++)
                {
                    bool breaks = SignedMaskBalanceClaim.RingBreaks(t, n);
                    bool pairs = PairsOnRing(n, (t, 1));
                    readings++;
                    if (pairs == breaks) { mismatches++; first ??= $"{t} on the ring of {n}"; }
                    bool dropped = SignedMaskBalanceClaim.Sign(t) == 1 && SignedMaskBalanceClaim.Valuation(SignedMaskBalanceClaim.Mask(t)) >= 1;
                    if (pairs == dropped) droppedPower++;
                }
            families.Add(new FamilyReading("one template on the ring, N = 2w − 1 to the cap", templates.Count, readings, mismatches, first));

            int r2 = 0, m2 = 0;
            string? f2 = null;
            for (int n = 3; n <= RingCap; n++)
            {
                r2++;
                if (PairsOnRing(n, ("XY", 1), ("ZZ", 1)) != (n % 2 == 0)) { m2++; f2 ??= $"XY + ZZ on the ring of {n}"; }
            }
            families.Add(new FamilyReading("XY + ZZ on the ring, paired exactly at even N", 1, r2, m2, f2));
        }

        // Theorem D: two templates on the ring, against Theorem 2 to the ring cap and against the balance to N = 32
        {
            var twoSite = (from x in "XYZ" from y in "XYZ" select $"{x}{y}").Where(x => SignedMaskBalanceClaim.Mask(x) != 0);
            var pool = twoSite.Concat(Windows.Where(x => SignedMaskBalanceClaim.Mask(x) != 0)).ToList();
            var wide = pool.Concat(new[] { "XYYY", "XIIY", "XXIY", "YIIY", "XYIZ", "YZZX", "XIIIY", "YIIIY", "XYIYX" }).ToList();
            int cells = 0, readings = 0, mismatches = 0, rcells = 0, rreadings = 0, rmismatches = 0;
            string? first = null, rfirst = null;
            for (int i = 0; i < wide.Count; i++)
                for (int j = i + 1; j < wide.Count; j++)
                {
                    string a = wide[i], b = wide[j];
                    if (SameShape(SignedMaskBalanceClaim.Mask(a), SignedMaskBalanceClaim.Mask(b))) continue;
                    var ab = new[] { a, b };
                    int lo = Math.Max(3, 2 * Math.Max(a.Length, b.Length) - 1);
                    rcells++;
                    for (int n = lo; n <= 32; n++)
                    {
                        bool breaks = SignedMaskBalanceClaim.RingBreaks(ab, n);
                        rreadings++;
                        if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(ab, n)) == breaks) { rmismatches++; rfirst ??= $"{a} + {b} on the ring of {n}"; }
                        if (RingParityOnly(ab, n) != breaks) parityOnly++;
                        if (RingAtTheMostValued(ab, n) != breaks) mostValued++;
                    }
                    foreach (long sign in j < pool.Count ? new[] { 1L, -1L } : new[] { 1L })
                    {
                        cells++;
                        for (int n = lo; n <= RingCap; n++)
                        {
                            readings++;
                            if (PairsOnRing(n, (a, 1), (b, sign)) == SignedMaskBalanceClaim.RingBreaks(ab, n)) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                        }
                    }
                }
            families.Add(new FamilyReading("Theorem D: two templates on the ring, against Theorem 2 to the cap (both signs below four sites wide)", cells, readings, mismatches, first));
            families.Add(new FamilyReading("Theorem D: two templates up to five sites wide on the ring, against the balance to N = 32", rcells, rreadings, rmismatches, rfirst, "the balance"));

            // three templates, valuations 0 to 4, so that T > v₀ and T ≤ v₀ both meet ties and greater valuations
            var sixteen = new[] { "XY", "XX", "XZ", "XXX", "XIY", "XYZ", "YIX", "XZY", "ZXY", "XXIY", "XIIY", "XYYY", "YIIY", "XIIIY", "YIIIY", "XYIYX" };
            int tcells = 0, treadings = 0, tmismatches = 0;
            string? tfirst = null;
            for (int i = 0; i < sixteen.Length; i++)
                for (int j = i + 1; j < sixteen.Length; j++)
                    for (int l = j + 1; l < sixteen.Length; l++)
                    {
                        var abc = new[] { sixteen[i], sixteen[j], sixteen[l] };
                        if (abc.Select(x => SignedMaskBalanceClaim.Mask(x)).Select(m => m >> BitOperations.TrailingZeroCount(m)).Distinct().Count() < 3) continue;
                        tcells++;
                        for (int n = Math.Max(3, 2 * abc.Max(x => x.Length) - 1); n <= 32; n++)
                        {
                            bool breaks = SignedMaskBalanceClaim.RingBreaks(abc, n);
                            treadings++;
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(abc, n)) == breaks) { tmismatches++; tfirst ??= $"{string.Join(" + ", abc)} on the ring of {n}"; }
                            if (RingParityOnly(abc, n) != breaks) parityOnly++;
                            if (RingAtTheMostValued(abc, n) != breaks) mostValued++;
                        }
                    }
            families.Add(new FamilyReading("Theorem D: three templates of different shapes on the ring, against the balance to N = 32", tcells, treadings, tmismatches, tfirst, "the balance"));
        }

        // three templates of different shapes at one N, by the balance
        {
            var pool = new[] { "XIY", "XXX", "XIZ", "XXZ", "XYX", "YZX", "XZZ" };
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            for (int i = 0; i < pool.Length; i++)
                for (int j = i + 1; j < pool.Length; j++)
                    for (int l = j + 1; l < pool.Length; l++)
                    {
                        string a = pool[i], b = pool[j], c = pool[l];
                        var shapes = new[] { a, b, c }.Select(x => SignedMaskBalanceClaim.Mask(x)).Select(m => m >> BitOperations.TrailingZeroCount(m));
                        if (shapes.Distinct().Count() < 3) continue;
                        cells++;
                        for (int n = 4; n <= 6; n++)
                        {
                            readings++;
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { a, b, c }, n)) != PairsOnChain(n, (a, 1), (b, 1), (c, 1)))
                            { mismatches++; first ??= $"{a} + {b} + {c} on the chain of {n}"; }
                            if (n < 5) continue;
                            readings++;
                            bool ringPairs = PairsOnRing(n, (a, 1), (b, 1), (c, 1));
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(new[] { a, b, c }, n)) != ringPairs
                                || SignedMaskBalanceClaim.RingBreaks(new[] { a, b, c }, n) == ringPairs)
                            { mismatches++; first ??= $"{a} + {b} + {c} on the ring of {n}"; }
                        }
                    }
            families.Add(new FamilyReading("three templates of different shapes at one N, by the balance (on the ring also by Theorem D)", cells, readings, mismatches, first));
        }

        // Theorem B′: three and four templates of different shapes on the chain, drawn (seeded) from the templates two to
        // four sites wide with a letter at each end, against the balance at every N from 1 to max(Σw + max w, N* + 2)
        // and against Theorem 2 to N = 8, the reduced basis's leading coefficient matrix of full rank; the same reading off the
        // unreduced basis, the control; and two templates against Theorem B's closed form
        {
            var pool = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows)
                .Concat(from a in "XYZ" from m1 in "IXYZ" from m2 in "IXYZ" from b in "XYZ" select $"{a}{m1}{m2}{b}")
                .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            var rng = new Random(173);
            int cells = 0, readings = 0, mismatches = 0, graphReadings = 0, graphMismatches = 0, breaking = 0, unreducedDeficient = 0, longChainReadings = 0, longChainMismatches = 0;
            string? first = null, graphFirst = null, longChainFirst = null;
            foreach (int k in new[] { 3, 4 })
                for (int drawn = 0; drawn < (k == 3 ? 400 : 200);)
                {
                    var ts = Enumerable.Range(0, k).Select(_ => pool[rng.Next(pool.Count)]).ToArray();
                    if (ts.Select(SignedMaskBalanceClaim.Mask).Select(m => m >> BitOperations.TrailingZeroCount(m)).Distinct().Count() < k) continue;
                    drawn++;
                    cells++;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(ts);
                    longChainReadings++;
                    if (breaks != SignedMaskBalanceClaim.RingBreaks(ts, 64)) { longChainMismatches++; longChainFirst ??= string.Join(" + ", ts); }
                    var (_, unreducedFrom) = SignedMaskBalanceClaim.Chain(ts, reduced: false);
                    if (breaks) breaking++;
                    // reduced: the leading coefficient matrix (row j: the coefficient of x^(sdeg bⱼ − wᵢ) in bⱼᵢ) of full rank k − 1,
                    // the property the predictable degree rests on, read apart from the reduction's own pivots
                    if (LeadingRank(SignedMaskBalanceClaim.ChainRelationBasis(ts), ts) != k - 1)
                    { mismatches++; first ??= $"{string.Join(" + ", ts)}: the reduced basis is not reduced"; }
                    if (LeadingRank(SignedMaskBalanceClaim.ChainRelationBasis(ts, reduced: false), ts) != k - 1) unreducedDeficient++;
                    int widest = ts.Max(t => t.Length), top = Math.Max(ts.Sum(t => t.Length) + widest, (from ?? 0) + 2);
                    for (int n = 1; n <= top; n++)
                    {
                        bool law = from is not int f || n < f;
                        bool balanced = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(ts, n));
                        readings++;
                        if (balanced != law) { mismatches++; first ??= $"{string.Join(" + ", ts)} on the chain of {n}"; }
                        if (balanced != (unreducedFrom is not int g || n < g)) unreducedBasis++;
                        if (n > 8) continue;
                        graphReadings++;
                        if (PairsOnChain(n, ts.Select(t => (t, 1L)).ToArray()) != law) { graphMismatches++; graphFirst ??= $"{string.Join(" + ", ts)} on the chain of {n}"; }
                    }
                }
            families.Add(new FamilyReading("Theorem B′: three and four templates on the chain, against the balance at every N from 1 to past N*", cells, readings, mismatches, first, "the balance"));
            families.Add(new FamilyReading("Theorem B′: the same cells against Theorem 2 to N = 8", cells, graphReadings, graphMismatches, graphFirst));
            families.Add(new FamilyReading("Theorem B′: whether the long chain breaks, against Theorem D on the ring of 64 (T > v₀)", cells, longChainReadings, longChainMismatches, longChainFirst, "Theorem D"));
            BPrimeBreaking = breaking;
            UnreducedBasisDeficient = unreducedDeficient;

            var two = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows).Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            int pairs = 0, pairMismatches = 0;
            string? pairFirst = null;
            for (int i = 0; i < two.Count; i++)
                for (int j = i + 1; j < two.Count; j++)
                {
                    ulong ma = SignedMaskBalanceClaim.Mask(two[i]), mb = SignedMaskBalanceClaim.Mask(two[j]);
                    if (SameShape(ma, mb)) continue;
                    pairs++;
                    if (SignedMaskBalanceClaim.Chain(new[] { two[i], two[j] }) != SignedMaskBalanceClaim.Chain(two[i], two[j]))
                    { pairMismatches++; pairFirst ??= $"{two[i]} + {two[j]}"; }
                }
            families.Add(new FamilyReading("Theorem B′ at two templates against Theorem B's N*", pairs, pairs, pairMismatches, pairFirst, "Theorem B"));
        }

        // Theorem E: two strings on one flip mask in a window of two or three sites, against Theorem 2
        {
            var pairsE = new List<(string A, string B)>();
            foreach (int w in new[] { 2, 3 })
                for (ulong f = 1; f < 1UL << w; f++)
                {
                    var strs = new List<string> { "" };
                    for (int j = 0; j < w; j++)
                    {
                        var letters = ((f >> j) & 1) == 1 ? "XY" : "IZ";
                        strs = strs.SelectMany(s => letters.Select(c => s + c)).ToList();
                    }
                    for (int i = 0; i < strs.Count; i++)
                        for (int j = i + 1; j < strs.Count; j++) pairsE.Add((strs[i], strs[j]));
                }
            int cells = 0, readings = 0, mismatches = 0, e2Readings = 0, e2Mismatches = 0, e2Cells = 0, e2Far = 0, e2FarMismatches = 0;
            int z4Readings = 0, z4Paired = 0, z4Mismatches = 0, movedNotQuadratic = 0;
            string? first = null, e2First = null, e2FarFirst = null, z4First = null;
            var branches = new Dictionary<SignedMaskBalanceClaim.RingShapeBranch, int>();
            foreach (var (a, b) in pairsE)
                foreach (long sign in new[] { 1L, -1L })
                {
                    cells++;
                    var terms = new[] { (a, 1L), (b, sign) };
                    int w = a.Length;
                    foreach (bool ring in new[] { false, true })
                        for (int n = ring ? Math.Max(3, 2 * w - 1) : Math.Max(3, w); n <= (ring ? RingCap : 7); n++)
                        {
                            var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                            var v = SignedMaskBalanceClaim.DecideCharacter(placed);
                            if (v == SignedMaskBalanceClaim.CharacterVerdict.Outside) { characterOutside++; continue; }
                            bool pairs = ring ? PairsOnRing(n, (a, 1), (b, sign)) : PairsOnChain(n, (a, 1), (b, sign));
                            readings++;
                            if (pairs != (v == SignedMaskBalanceClaim.CharacterVerdict.Pairs)) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}"; }
                            // the ℤ/4 form: the refinement exists exactly where Theorem 2 pairs, is one (on W up to 128 elements), and
                            // its carrier i^(−Q₀ − 2⟨R(w), x₀⟩) satisfies Theorem 2's (★) entry by entry and is a diagonal Clifford (its
                            // quarters a ℤ/4 quadratic form on all bitstrings); without the coset term it must fail (★)
                            {
                                var form = SignedMaskBalanceClaim.QuadraticRefinement(placed);
                                z4Readings++;
                                string at = $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}";
                                bool bad = (form is not null) != pairs;
                                if (form is not null)
                                {
                                    z4Paired++;
                                    if (form.Count <= 128 && !IsRefinement(form)) bad = true;
                                    var carrier = SignedMaskBalanceClaim.CarrierQuarters(form, n);
                                    if (!StarHolds(placed, n, carrier) || !QuartersQuadratic(carrier, n)) bad = true;
                                    if (MovedCarrierBreaches(placed, n, carrier, form, ref movedNotQuadratic)) bad = true;
                                    if (!StarHolds(placed, n, SignedMaskBalanceClaim.CarrierQuarters(form, n, cosetTerm: false))) droppedCosetTerm++;
                                }
                                if (bad) { z4Mismatches++; z4First ??= at; }
                            }
                            foreach (var (flag, k) in new[] { ((false, true, true), 0), ((true, false, true), 1), ((true, true, false), 2) })
                                if (pairs != (SignedMaskBalanceClaim.DecideCharacter(placed, flag.Item1, flag.Item2, flag.Item3) == SignedMaskBalanceClaim.CharacterVerdict.Pairs))
                                    droppedCondition[k]++;
                            if (!ring) continue;
                            var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                            if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside) continue;
                            if (n == Math.Max(3, 2 * w - 1)) e2Cells++;
                            e2Readings++;
                            branches[branch] = branches.GetValueOrDefault(branch) + 1;
                            if (pairs == breaks) { e2Mismatches++; e2First ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                            if (pairs == SignedMaskBalanceClaim.RingOneShape(terms, n, readSetTest: false).Breaks) droppedReadSet++;
                        }
                    for (int n = RingCap + 1; n <= 40; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                        if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside) continue;
                        var v = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(terms, n, ring: true));
                        e2Far++;
                        if (breaks != (v == SignedMaskBalanceClaim.CharacterVerdict.Breaks)) { e2FarMismatches++; e2FarFirst ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                    }
                }
            // four strings on one flip set: the products (aX + bY)⊗(cX + dY), a, b, c, d = ±1, where hol(g) is the ordering term alone
            int pcells = 0, preadings = 0, pmismatches = 0;
            string? pfirst = null;
            foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
            {
                pcells++;
                var terms = new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) };
                for (int n = 3; n <= RingCap; n++)
                {
                    bool pairs = PairsOnRing(n, terms.Select(t => (t.Item1, t.Item2)).ToArray());
                    var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                    var v = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(terms, n, ring: true));
                    preadings++;
                    if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside || pairs == breaks || pairs != (v == SignedMaskBalanceClaim.CharacterVerdict.Pairs))
                    { pmismatches++; pfirst ??= $"({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y) on the ring of {n}"; }
                    if (pairs == SignedMaskBalanceClaim.RingOneShape(terms, n, orderingTerm: false).Breaks) droppedOrdering++;
                }
            }
            families.Add(new FamilyReading("Theorem E: two strings on one flip mask, windows of two and three sites, chain to 7 and ring to the cap", cells, readings, mismatches, first));
            families.Add(new FamilyReading("Corollary E2 and Theorem E on (aX + bY)⊗(cX + dY), four strings on one flip set, ring to the cap", pcells, preadings, pmismatches, pfirst));
            families.Add(new FamilyReading("Theorem E's ℤ/4 form: the refinement where Theorem 2 pairs, its carrier against (★) and as a diagonal Clifford", cells, z4Readings, z4Mismatches, z4First));
            Z4Paired = z4Paired;

            // the ℤ/4 form on four strings per flip set, where ω is real beside a nonempty read set: the sixteen products
            // (aX + bY)⊗(cX + dY) on the ring and X ⊗ 2·CZ on chain and ring, N = 3 to 8
            int fcells = 0, freadings = 0, fpaired = 0, fmismatches = 0;
            string? ffirst = null;
            var fours = new List<((string, long)[] Terms, string Name, bool[] Rings)>();
            foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
                fours.Add((new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) },
                           $"({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y)", new[] { true }));
            fours.Add((new[] { ("XII", 1L), ("XIZ", 1L), ("XZI", 1L), ("XZZ", -1L) }, "X ⊗ 2·CZ", new[] { false, true }));
            foreach (var (terms, name, rings) in fours)
            {
                fcells++;
                foreach (bool ring in rings)
                    for (int n = 3; n <= 8; n++)
                    {
                        var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                        bool pairs = ring ? PairsOnRing(n, terms) : PairsOnChain(n, terms);
                        var form = SignedMaskBalanceClaim.QuadraticRefinement(placed);
                        freadings++;
                        bool bad = (form is not null) != pairs;
                        if (form is not null)
                        {
                            fpaired++;
                            if (form.Count <= 128 && !IsRefinement(form)) bad = true;
                            var carrier = SignedMaskBalanceClaim.CarrierQuarters(form, n);
                            if (!StarHolds(placed, n, carrier) || !QuartersQuadratic(carrier, n)) bad = true;
                            if (MovedCarrierBreaches(placed, n, carrier, form, ref movedNotQuadratic)) bad = true;
                        }
                        if (bad) { fmismatches++; ffirst ??= $"{name} on the {(ring ? "ring" : "chain")} of {n}"; }
                    }
            }
            families.Add(new FamilyReading("Theorem E's ℤ/4 form on four strings per flip set, a real ω beside a nonempty read set", fcells, freadings, fmismatches, ffirst));
            Z4FourStringPaired = fpaired;
            MovedCarrierNotQuadratic = movedNotQuadratic;
            families.Add(new FamilyReading("Corollary E2: one flip shape on the ring, against Theorem 2 to the cap", e2Cells, e2Readings, e2Mismatches, e2First));
            families.Add(new FamilyReading("Corollary E2 against Theorem E's elimination, rings past the cap to N = 40", e2Cells, e2Far, e2FarMismatches, e2FarFirst, "Theorem E"));
            E2Branches = branches;
        }

        // Corollary E3: two templates of different flip shapes on the chain, each one string or two on one mask, with a
        // character ratio; (a) at least one with a read set or ω = ±i and (1) holding at every offset, N = 3 to 8; (b) one
        // template four sites wide, where D ≠ 0 and (1) holds at N*, read from N* − 1 (at least 4) to N* + 1
        {
            var small = OneShapeTemplates(1, 2, 3);
            var wide = OneShapeTemplates(4);
            int acells = 0, areadings = 0, amismatches = 0, bcells = 0, breadings = 0, bmismatches = 0, e3pReadings = 0, e3pMismatches = 0;
            string? afirst = null, bfirst = null, e3pFirst = null;
            var e3Branches = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
            string Name(IReadOnlyList<(string T, long C)> t) => string.Join(" ", t.Select((s, k) => k == 0 ? s.T : $"{s.C:+0;-0}{s.T}"));
            for (int i = 0; i < small.Count; i++)
                for (int j = i + 1; j < small.Count; j++)
                {
                    var (ta, ra) = small[i];
                    var (tb, rb) = small[j];
                    if (SameShape(ra.Mask, rb.Mask)) continue;
                    if (ra.Read == 0 && ra.Quarter % 2 == 0 && rb.Read == 0 && rb.Quarter % 2 == 0) continue;
                    int wa = ta[0].T.Length, wb = tb[0].T.Length;
                    int every = Math.Max(wa + wb - 1, 2 * Math.Max(wa, wb) - 1);
                    if (SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, every).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail) continue;
                    acells++;
                    for (int n = Math.Max(3, Math.Max(wa, wb)); n <= 8; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n);
                        bool pairs = PairsOnChain(n, ta.Concat(tb).ToArray());
                        areadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { amismatches++; afirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n, orderingTerm: false).Breaks) droppedE3Ordering++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new[] { ta, tb }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                    }
                }
            foreach (var (ta, ra) in small)
                foreach (var (tb, rb) in wide)
                {
                    if (SameShape(ra.Mask, rb.Mask)) continue;
                    var (nStar, readSum) = SignedMaskBalanceClaim.ChainTwoShapeData(ta, tb);
                    if (readSum == 0 || SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, nStar).Branch != SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero) continue;
                    bcells++;
                    for (int n = Math.Max(4, nStar - 1); n <= nStar + 1; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n);
                        bool pairs = PairsOnChain(n, ta.Concat(tb).ToArray());
                        breadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { bmismatches++; bfirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n, readSumTest: false).Breaks) droppedReadSum++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new[] { ta, tb }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                    }
                }
            // the field X ± Y beside the sixteen products (aX + bY)⊗(cX + dY): read sets on both, and an ordering term of +1 in
            // hol(u₀), where on the families above it is −1 wherever a read set is nonzero
            int pcells = 0, preadings = 0, pmismatches = 0;
            string? pfirst = null;
            foreach (long fy in new[] { 1L, -1L })
                foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
                {
                    pcells++;
                    var field = new[] { ("X", 1L), ("Y", fy) };
                    var product = new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) };
                    for (int n = 2; n <= 8; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(field, product, n);
                        bool pairs = PairsOnChain(n, field.Concat(product).ToArray());
                        preadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { pmismatches++; pfirst ??= $"X {fy:+0;-0}Y beside ({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y) at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(field, product, n, orderingTerm: false).Breaks) droppedE3Ordering++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new IReadOnlyList<(string, long)>[] { field, product }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"X {fy:+0;-0}Y beside the product at N = {n}"; }
                    }
                }
            families.Add(new FamilyReading("Corollary E3: two flip shapes on the chain, a read set or ω = ±i, (1) at every offset, N = 3 to 8", acells, areadings, amismatches, afirst));
            families.Add(new FamilyReading("Corollary E3 where the read sum D ≠ 0 brings the break forward to N*, one template four sites wide", bcells, breadings, bmismatches, bfirst));
            families.Add(new FamilyReading("Corollary E3 on the field X ± Y beside (aX + bY)⊗(cX + dY), an ordering term of +1, N = 2 to 8", pcells, preadings, pmismatches, pfirst));
            E3Branches = e3Branches;
            families.Add(new FamilyReading("Corollary E3′ at two templates against Corollary E3, branch for branch, on the three families above", acells + bcells + pcells, e3pReadings, e3pMismatches, e3pFirst, "Corollary E3"));

            // Corollary E3′ at three templates of different flip shapes, against Theorem 2: (a) drawn (seeded) from the
            // one-shape templates up to three sites wide, N from the widest to 8; (b) drawn from those and the four-site ones,
            // kept where (1) holds at every offset (read at N = 12), N from the widest to 9; the controls read N₁* away and
            // walk hol without its ordering term
            var rng = new Random(1173);
            var all = small.Concat(wide).ToList();
            var tBranches = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
            foreach (bool rich in new[] { false, true })
            {
                var from = rich ? all : small;
                int tcells = 0, treadings = 0, tmismatches = 0;
                string? tfirst = null;
                while (tcells < (rich ? 200 : 500))
                {
                    var ts = Enumerable.Range(0, 3).Select(_ => from[rng.Next(from.Count)]).ToArray();
                    if (ts.Select(t => t.Ratio.Mask >> BitOperations.TrailingZeroCount(t.Ratio.Mask)).Distinct().Count() < 3) continue;
                    var terms = ts.Select(t => (IReadOnlyList<(string T, long C)>)t.Strings).ToArray();
                    if (rich && SignedMaskBalanceClaim.ChainShapes(terms, 12).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail) continue;
                    tcells++;
                    int widest = ts.Max(t => t.Strings[0].T.Length);
                    for (int n = Math.Max(3, widest); n <= (rich ? 9 : 8); n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainShapes(terms, n);
                        bool pairs = PairsOnChain(n, ts.SelectMany(t => t.Strings).ToArray());
                        treadings++;
                        tBranches[branch] = tBranches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks || branch == SignedMaskBalanceClaim.ChainShapeBranch.Outside) { tmismatches++; tfirst ??= $"{string.Join(" with ", ts.Select(t => Name(t.Strings)))} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainShapes(terms, n, readSumTest: false).Breaks) droppedE3PrimeReadSum++;
                        if (pairs == SignedMaskBalanceClaim.ChainShapes(terms, n, orderingTerm: false).Breaks) droppedE3PrimeOrdering++;
                    }
                }
                families.Add(new FamilyReading(rich
                    ? "Corollary E3′: three flip shapes on the chain, (1) at every offset, templates up to four sites wide, N to 9"
                    : "Corollary E3′: three flip shapes on the chain, up to three sites wide, N to 8", tcells, treadings, tmismatches, tfirst));
            }
            E3PrimeBranches = tBranches;
        }

        // Theorem H at random rates: on the window pairs of the first family (one string per flip set), chain N = 3 to 7 and ring
        // N = 5 to 7, λ_min of 𝓜 built from H's entries against the weighted code distance; on those and on Theorem E's two
        // strings on one flip mask, chain to 7 and ring 2w − 1 to 7, the kernels of 𝓛 and 𝓜 against Theorem 2's two counts
        {
            var hRng = new Random(20261112);
            double eps = Math.Pow(2, -52);
            int codeCells = 0, codeReadings = 0, codeMismatches = 0, kernelCells = 0, kernelReadings = 0, kernelMismatches = 0;
            string? codeFirst = null, kernelFirst = null;
            double[] Eigs(MathNet.Numerics.LinearAlgebra.Matrix<Complex> m) =>
                m.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Select(z => z.Real).OrderBy(v => v).ToArray();
            // the zeros of a positive semidefinite matrix, counted below 10⁻⁹ of its largest eigenvalue; the separation recorded
            int Zeros(double[] ev)
            {
                double scale = Math.Max(ev[^1], double.Epsilon);
                int k = ev.Count(v => v < 1e-9 * scale);
                if (k > 0) KernelZeroMax = Math.Max(KernelZeroMax, Math.Max(Math.Abs(ev[0]), Math.Abs(ev[k - 1])) / scale);
                MostNegative = Math.Min(MostNegative, ev[0] / scale);
                if (k < ev.Length) KernelGapMin = Math.Min(KernelGapMin, ev[k] / scale);
                return k;
            }
            bool KernelsAgree(int n, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed, double[] gamma, out double[] farEigs)
            {
                var h = ComplementConnectionGraph.Columns(n, placed.Select(s => (new PauliMask(s.X, s.Z), (BigInteger)s.C)));
                var (nearCount, farCount, _) = ComplementConnectionGraph.EverySiteDephased(n, h);
                var (near, far) = StrongDephasingEdges.Build(n, h, gamma);
                farEigs = Eigs(far);
                return Zeros(Eigs(near)) == nearCount && Zeros(farEigs) == farCount;
            }
            double[] Rates(int n) => Enumerable.Range(0, n).Select(_ => 0.5 + hRng.NextDouble()).ToArray();

            for (int i = 0; i < Windows.Length; i++)
                for (int j = i + 1; j < Windows.Length; j++)
                {
                    string a = Windows[i], b = Windows[j];
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    foreach (long sign in new[] { 1L, -1L })
                    {
                        codeCells++;
                        kernelCells++;
                        var terms = new[] { (a, 1L), (b, sign) };
                        foreach (bool ring in new[] { false, true })
                            for (int n = ring ? 5 : 3; n <= 7; n++)
                            {
                                var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                                var gamma = Rates(n);
                                string at = $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}";
                                kernelReadings++;
                                if (!KernelsAgree(n, placed, gamma, out var farEigs)) { kernelMismatches++; kernelFirst ??= at; }
                                double lam = farEigs[0], bound = (1 << n) * eps * farEigs[^1];
                                double code = SignedMaskBalanceClaim.StrongDephasingBreak(placed, gamma);
                                codeReadings++;
                                StrongDephasingWorstRatio = Math.Max(StrongDephasingWorstRatio, Math.Abs(lam - code) / bound);
                                if (Math.Abs(lam - code) > bound) { codeMismatches++; codeFirst ??= at; }
                                if (lam > bound) StrongDephasingBroken++;
                                // the control: every rate read as one
                                if (Math.Abs(lam - SignedMaskBalanceClaim.StrongDephasingBreak(placed, Enumerable.Repeat(1.0, n).ToArray())) > bound) RatesReadAwayMismatches++;
                                // and the rates in reverse site order
                                if (Math.Abs(lam - SignedMaskBalanceClaim.StrongDephasingBreak(placed, gamma.Reverse().ToArray())) > bound) RatesReversedMismatches++;
                            }
                    }
                }
            foreach (int w in new[] { 2, 3 })
                for (ulong f = 1; f < 1UL << w; f++)
                {
                    var strs = new List<string> { "" };
                    for (int j = 0; j < w; j++) strs = strs.SelectMany(s => (((f >> j) & 1) == 1 ? "XY" : "IZ").Select(c => s + c)).ToList();
                    for (int i = 0; i < strs.Count; i++)
                        for (int j = i + 1; j < strs.Count; j++)
                            foreach (long sign in new[] { 1L, -1L })
                            {
                                kernelCells++;
                                var terms = new[] { (strs[i], 1L), (strs[j], sign) };
                                foreach (bool ring in new[] { false, true })
                                    for (int n = ring ? Math.Max(3, 2 * w - 1) : Math.Max(3, w); n <= 7; n++)
                                    {
                                        kernelReadings++;
                                        if (!KernelsAgree(n, SignedMaskBalanceClaim.PlaceStrings(terms, n, ring), Rates(n), out _))
                                        { kernelMismatches++; kernelFirst ??= $"{strs[i]} {sign:+0;-0}{strs[j]} on the {(ring ? "ring" : "chain")} of {n}"; }
                                    }
                            }
                }
            // XZ ± ZY on the chain of four to seven, one magnitude c for both, ten rate draws each: λ_min(𝓜) against the two end
            // squares q(c²/γ₀, 2c²/γ₁) + q(c²/γ_(N−1), 2c²/γ_(N−2)), q(u, v) = u + v − √(u² + v²)
            int endCells = 0, endMismatches = 0;
            string? endFirst = null;
            static double Q(double u, double v) => u + v - Math.Sqrt(u * u + v * v);
            for (int n = 4; n <= 7; n++)
                foreach (long sign in new[] { 1L, -1L })
                    for (int draw = 0; draw < 10; draw++)
                    {
                        long c = 1 + draw % 3;
                        var placed = SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", c), ("ZY", sign * c) }, n, ring: false);
                        var gamma = Rates(n);
                        var h = ComplementConnectionGraph.Columns(n, placed.Select(s => (new PauliMask(s.X, s.Z), (BigInteger)s.C)));
                        var farEigs = Eigs(StrongDephasingEdges.Build(n, h, gamma).Far);
                        double c2 = c * c, ends = Q(c2 / gamma[0], 2 * c2 / gamma[1]) + Q(c2 / gamma[n - 1], 2 * c2 / gamma[n - 2]);
                        endCells++;
                        // and Corollary H2's packing of the two squares, each the walk a, b, a, b (k = 4, every mask twice, s = 4)
                        double squares = 0;
                        foreach (var (sa, sb) in new[] { (1, 2), (1 << (n - 1), 1 << (n - 2)) })
                            squares += StrongDephasingEdges.CycleReading(n, h, gamma, new[] { sa, sb, sa, sb })!.Value.LambdaStar;
                        double endBound = (1 << n) * eps * farEigs[^1];
                        if (Math.Abs(farEigs[0] - ends) > endBound || Math.Abs(farEigs[0] - squares) > endBound) { endMismatches++; endFirst ??= $"XZ {sign:+0;-0}ZY, c = {c}, on the chain of {n}"; }
                    }
            families.Add(new FamilyReading("Theorem H: XZ ± ZY on the chain of four to seven, one magnitude for both strings, λ_min of the far edge's generator against the two end squares, q and Corollary H2's packing of the two squares, ten rate draws each", endCells, endCells, endMismatches, endFirst, "the end squares"));
            families.Add(new FamilyReading("Theorem H: λ_min of the far edge's second-order generator, from H's entries, against the weighted code distance, window pairs at random rates, chain N = 3 to 7 and ring N = 5 to 7", codeCells, codeReadings, codeMismatches, codeFirst, "the code distance"));
            families.Add(new FamilyReading("Theorem H's lemma: the kernels of the two second-order generators against Theorem 2's near and far counts, the window pairs and Theorem E's two strings on one flip mask", kernelCells, kernelReadings, kernelMismatches, kernelFirst));
        }

        // Corollary H2 on Theorem E's cells, two strings on one flip mask (2 or 3 sites wide, both signs), chain N = 3 to 7 and
        // ring N = 2w − 1 (from 3) to 7, rates of their own seed: every minimal relation among the placed flip sets whose walk
        // (in ascending order of the masks) is a simple cycle with holonomy ≠ 1 from every bitstring, its bound (k/2)·λ_* at most
        // λ_min(𝓜), λ_* at least the least hop times 2(1 − cos(θ_min/k)), and the relations sharing no flip set, taken greedily
        // by size of bound, summed; the control, k·λ_* without the double count, exceeds λ_min(𝓜) somewhere
        {
            var cRng = new Random(20261116);
            double eps = Math.Pow(2, -52);
            int cells = 0, mismatches = 0, hopMismatches = 0;
            string? first = null, hopFirst = null;
            foreach (int w in new[] { 2, 3 })
                for (ulong f = 1; f < 1UL << w; f++)
                {
                    var strs = new List<string> { "" };
                    for (int j = 0; j < w; j++) strs = strs.SelectMany(t => (((f >> j) & 1) == 1 ? "XY" : "IZ").Select(c => t + c)).ToList();
                    for (int i = 0; i < strs.Count; i++)
                        for (int j = i + 1; j < strs.Count; j++)
                            foreach (long sign in new[] { 1L, -1L })
                            {
                                cells++;
                                var terms = new[] { (strs[i], 1L), (strs[j], sign) };
                                foreach (bool ring in new[] { false, true })
                                    for (int n = ring ? Math.Max(3, 2 * w - 1) : Math.Max(3, w); n <= 7; n++)
                                    {
                                        var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                                        var gamma = Enumerable.Range(0, n).Select(_ => 0.5 + cRng.NextDouble()).ToArray();
                                        var h = ComplementConnectionGraph.Columns(n, placed.Select(p => (new PauliMask(p.X, p.Z), (BigInteger)p.C)));
                                        var far = StrongDephasingEdges.Build(n, h, gamma).Far.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian)
                                                                    .EigenValues.Select(z => z.Real).OrderBy(v => v).ToArray();
                                        double lam = far[0], tol = 4 * (1 << n) * eps * Math.Max(far[^1], 1);
                                        // Corollary H3 on the same reading: the whole spectrum where the hops commute, the
                                        // closed floor where every two anticommute; the control, Theorem H string by string
                                        var hops = StrongDephasingEdges.ReadHops(n, h, gamma);
                                        if (hops is null) HopOutside++;
                                        else if (hops.Spectrum is { } spec)
                                        {
                                            HopReadings++;
                                            bool off = spec.Length != far.Length || spec.Zip(far, (x, y) => Math.Abs(x - y)).Max() > tol;
                                            if (off) { hopMismatches++; hopFirst ??= $"{strs[i]} {sign:+0;-0}{strs[j]} on the {(ring ? "ring" : "chain")} of {n}"; }
                                            if (lam > tol) HopBroken++;
                                            if (hops.Split) HopSplit++;
                                            double naive = double.PositiveInfinity;
                                            for (int phi = 0; phi < 1 << n; phi++)
                                            {
                                                double tot = 0;
                                                foreach (var pl in placed)
                                                    if ((BitOperations.PopCount((ulong)phi & pl.X) & 1) != (BitOperations.PopCount(pl.Z) & 1))
                                                        tot += (double)pl.C * pl.C / Enumerable.Range(0, n).Where(l => ((pl.X >> l) & 1) == 1).Sum(l => gamma[l]);
                                                naive = Math.Min(naive, tot);
                                            }
                                            if (Math.Abs(2 * naive - lam) > tol) HopControlMisses++;
                                        }
                                        else
                                        {
                                            HopNoncommuting++;
                                            if (hops.LambdaMin is { } floor)
                                            {
                                                HopAllAnticommute++;
                                                if (Math.Abs(floor - lam) > tol) { hopMismatches++; hopFirst ??= $"{strs[i]} {sign:+0;-0}{strs[j]} on the {(ring ? "ring" : "chain")} of {n}, every two hops anticommuting"; }
                                            }
                                        }
                                        var masks = placed.Select(p => (int)p.X).Where(m => m != 0).Distinct().OrderBy(m => m).ToArray();
                                        var frustrated = new List<(int[] Masks, double Bound)>();
                                        bool bad = false;
                                        for (int sub = 1; sub < 1 << masks.Length; sub++)
                                        {
                                            var u = Enumerable.Range(0, masks.Length).Where(q => ((sub >> q) & 1) == 1).Select(q => masks[q]).ToArray();
                                            if (u.Length < 3 || u.Aggregate(0, (acc, m) => acc ^ m) != 0) continue;
                                            // minimal: rank k − 1, and k − 1 after dropping any one of them
                                            if (Rank(u) != u.Length - 1 || u.Any(e => Rank(u.Where(x => x != e).ToArray()) != u.Length - 1)) continue;
                                            var reading = StrongDephasingEdges.CycleReading(n, h, gamma, u);
                                            if (reading is not { } r || r.ThetaMin < 1e-9) continue;
                                            int k = u.Length;
                                            double bound = k / 2.0 * r.LambdaStar;
                                            CycleRelations++;
                                            if (bound > lam + tol) bad = true;
                                            if (r.LambdaStar < r.HopMin * 2 * (1 - Math.Cos(r.ThetaMin / k)) - tol) bad = true;
                                            if (k * r.LambdaStar > lam + tol) CycleDoubleCountMismatches++;   // the control
                                            frustrated.Add((u, bound));
                                        }
                                        if (frustrated.Count == 0) continue;
                                        double packed = 0;
                                        var used = new HashSet<int>();
                                        foreach (var (u, bound) in frustrated.OrderByDescending(x => x.Bound))
                                            if (u.All(m => !used.Contains(m))) { packed += bound; foreach (int m in u) used.Add(m); }
                                        if (packed > lam + tol) bad = true;
                                        if (Math.Abs(packed - lam) <= tol) CyclePackingEqualities++;
                                        CycleTightest = Math.Max(CycleTightest, packed / lam);
                                        CyclePackingReadings++;
                                        if (bad) { mismatches++; first ??= $"{strs[i]} {sign:+0;-0}{strs[j]} on the {(ring ? "ring" : "chain")} of {n}"; }
                                    }
                            }
                }
            // with a potential (moduli that differ on the cycle): the triangle 10·XXI + k·XXZ + 10·IXX + k·ZXX + 10·XIY on three
            // sites, k = 0 to 9, three rate draws each, the walk over the flip sets {0, 1}, {1, 2}, {0, 2}
            int potCells = 0, potReadings = 0, potBad = 0;
            string? potFirst = null;
            for (long kk = 0; kk <= 9; kk++)
            {
                potCells++;
                var placed = SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXI", 10L), ("XXZ", kk), ("IXX", 10L), ("ZXX", kk), ("XIY", 10L) }, 3, ring: false);
                var h = ComplementConnectionGraph.Columns(3, placed.Where(p => p.C != 0).Select(p => (new PauliMask(p.X, p.Z), (BigInteger)p.C)));
                for (int draw = 0; draw < 3; draw++)
                {
                    var gamma = Enumerable.Range(0, 3).Select(_ => 0.5 + cRng.NextDouble()).ToArray();
                    var far = StrongDephasingEdges.Build(3, h, gamma).Far.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian)
                                                .EigenValues.Select(z => z.Real).OrderBy(v => v).ToArray();
                    double lam = far[0], tol = 4 * 8 * eps * Math.Max(far[^1], 1);
                    var r = StrongDephasingEdges.CycleReading(3, h, gamma, new[] { 0b011, 0b110, 0b101 })!.Value;
                    potReadings++;
                    if (kk > 0) CyclePotentialReadings++;
                    if (1.5 * r.LambdaStar > lam + tol || r.LambdaStar < r.HopMin * 2 * (1 - Math.Cos(r.ThetaMin / 3)) - tol || r.ThetaMin < 1e-9)
                    { potBad++; potFirst ??= $"k = {kk}, draw {draw}"; }
                    if (1.5 * r.LambdaStar > 2 * r.HopMin + tol) { if (kk > 0) CyclePotentialAboveTwoHops++; else CycleEqualModuliAboveTwoHops++; }
                }
            }
            families.Add(new FamilyReading("Corollary H2 with a potential: the triangle 10·XXI + k·XXZ + 10·IXX + k·ZXX + 10·XIY on three sites, k = 0 to 9, three rate draws each, (3/2)·λ_* at most λ_min(𝓜) and λ_* above its floor", potCells, potReadings, potBad, potFirst, "the cycle bound"));
            families.Add(new FamilyReading("Corollary H3: Theorem E's two strings on one flip mask at random rates, chain N = 3 to 7 and ring N = 2w − 1 to 7, wherever no entry vanishes (every modulus there one number): the far edge's whole second-order spectrum against {2·Σ μ_m·y_m} on each coset where the hops commute, and λ_min against Σ μ − (Σ μ²)^(1/2) where every two anticommute", cells, HopReadings + HopAllAnticommute, hopMismatches, hopFirst, "the hop reading"));
            families.Add(new FamilyReading("Corollary H2: Theorem E's two strings on one flip mask at random rates, chain N = 3 to 7 and ring N = 2w − 1 to 7, every minimal relation whose walk is a simple cycle with holonomy ≠ 1 from every bitstring, (k/2)·λ_* and the packing of those sharing no flip set at most λ_min(𝓜); the readings, a cell at one N, those carrying such a relation", cells, CyclePackingReadings, mismatches, first, "the cycle bound"));
        }

        // Corollary H1 and the trellis, on the window pairs at random rates (a seed of their own): (a) the trellis against the
        // enumeration, chain N = 3 to 12 and ring N = 5 to 12; (b) on the chain where Theorem B breaks, the relation r = 1 at
        // N* and its translates by D = 1 + max deg q_i, pairwise disjoint and each unbalanced, against the trellis from N* to
        // 3N*: λ ≥ 2·Σ_t the least weight in translate t (the control: every offset 0 … N − N*, which overlap); (c) on the ring where the
        // chain of the same N pairs: λ ≤ 2·min over the N cuts of the weights of the strings whose window straddles the cut
        // (the control: the same bound on the rings where the chain breaks)
        {
            var dRng = new Random(20261115);
            double eps = Math.Pow(2, -52);
            double[] Rates(int n) => Enumerable.Range(0, n).Select(_ => 0.5 + dRng.NextDouble()).ToArray();
            int trellisCells = 0, trellisReadings = 0, trellisMismatches = 0;
            int packCells = 0, packReadings = 0, packMismatches = 0;
            int cutCells = 0, cutReadings = 0, cutMismatches = 0;
            string? trellisFirst = null, packFirst = null, cutFirst = null;
            static ulong PolyMod(ulong a, ulong b) { int db = 63 - BitOperations.LeadingZeroCount(b); while (a != 0 && 63 - BitOperations.LeadingZeroCount(a) >= db) a ^= b << (63 - BitOperations.LeadingZeroCount(a) - db); return a; }
            static ulong PolyGcd(ulong a, ulong b) { while (b != 0) (a, b) = (b, PolyMod(a, b)); return a; }
            static ulong PolyDiv(ulong a, ulong b)
            {
                ulong q = 0; int db = 63 - BitOperations.LeadingZeroCount(b);
                while (a != 0 && 63 - BitOperations.LeadingZeroCount(a) >= db) { int s = 63 - BitOperations.LeadingZeroCount(a) - db; q |= 1UL << s; a ^= b << s; }
                return q;
            }
            // the weight c²/γ(m) of the template t placed at start p (wrapped on n sites)
            static double Weight(string t, long c, int p, int n, double[] gamma)
            {
                double rate = 0;
                for (int i = 0; i < t.Length; i++) if (t[i] is 'X' or 'Y') rate += gamma[(p + i) % n];
                return (double)c * c / rate;
            }

            for (int i = 0; i < Windows.Length; i++)
                for (int j = i + 1; j < Windows.Length; j++)
                {
                    string a = Windows[i], b = Windows[j];
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                    foreach (long sign in new[] { 1L, -1L })
                    {
                        var terms = new[] { (a, 1L), (b, sign) };
                        string Name(bool ring, int n) => $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}";

                        // (a) the trellis against the enumeration
                        trellisCells++;
                        foreach (bool ring in new[] { false, true })
                            for (int n = ring ? 5 : 3; n <= 12; n++)
                            {
                                var gamma = Rates(n);
                                double enumerated = SignedMaskBalanceClaim.StrongDephasingBreak(SignedMaskBalanceClaim.PlaceStrings(terms, n, ring), gamma);
                                double trellis = SignedMaskBalanceClaim.StrongDephasingBreakTrellis(terms, gamma, ring);
                                // both sum at most 2N positive weights, each weight the same double in both (its rates summed in
                                // ascending site order): the summation order alone, a law of 2N·eps relative to the larger value
                                // (equal values, zeros among them, read 0; a value that is not finite reads NaN and counts as a miss)
                                double ratio = enumerated == trellis ? 0 : Math.Abs(enumerated - trellis) / (2 * n * eps * Math.Max(enumerated, trellis));
                                if (ratio > TrellisWorstRatio || double.IsNaN(ratio)) TrellisWorstRatio = ratio;
                                trellisReadings++;
                                if (!(ratio <= 1)) { trellisMismatches++; trellisFirst ??= Name(ring, n); }
                            }

                        // (b) the packing bound on the chain, where Theorem B breaks
                        if (breaks && from is int nStar)
                        {
                            packCells++;
                            ulong g = PolyGcd(ma, mb), qa = PolyDiv(ma, g), qb = PolyDiv(mb, g);
                            int step = 1 + Math.Max(63 - BitOperations.LeadingZeroCount(qa), 63 - BitOperations.LeadingZeroCount(qb));
                            // the relation r = 1: template a at the starts of q_b, template b at the starts of q_a, shifted by o
                            double Least(int o, double[] gamma)
                            {
                                double least = double.PositiveInfinity;
                                for (int k = 0; k < 64; k++)
                                {
                                    if (((qb >> k) & 1) == 1) least = Math.Min(least, Weight(a, 1, o + k, int.MaxValue, gamma));
                                    if (((qa >> k) & 1) == 1) least = Math.Min(least, Weight(b, sign, o + k, int.MaxValue, gamma));
                                }
                                return least;
                            }
                            for (int n = nStar; n <= 3 * nStar; n++)
                            {
                                var gamma = Rates(n);
                                double lam = SignedMaskBalanceClaim.StrongDephasingBreakTrellis(terms, gamma, ring: false);
                                double packed = 0, every = 0;
                                for (int o = 0; o + nStar <= n; o += step) packed += Least(o, gamma);
                                for (int o = 0; o + nStar <= n; o++) every += Least(o, gamma);
                                packReadings++;
                                if (lam < 2 * packed * (1 - 4 * n * eps)) { packMismatches++; packFirst ??= Name(false, n); }
                                if (packed > 0 && Math.Abs(lam - 2 * packed) <= 4 * n * eps * lam) PackingEqualities++;
                                if (lam < 2 * every * (1 - 4 * n * eps)) PackingOverlapMismatches++;   // the control
                            }
                        }

                        // (c) the cut bound on the ring, where the chain of the same N pairs
                        cutCells++;
                        for (int n = 5; n <= 12; n++)
                        {
                            bool chainPairs = !breaks || n < from;
                            var gamma = Rates(n);
                            double lam = SignedMaskBalanceClaim.StrongDephasingBreakTrellis(terms, gamma, ring: true);
                            double cut = double.PositiveInfinity;
                            for (int c = 0; c < n; c++)
                            {
                                // the strings whose window holds both site c − 1 and site c (mod n)
                                double straddle = 0;
                                foreach (var (t, coeff) in terms)
                                    for (int p = 0; p < n; p++)
                                    {
                                        bool left = false, right = false;
                                        for (int k = 0; k < t.Length; k++) { int s = (p + k) % n; left |= s == (c + n - 1) % n; right |= s == c; }
                                        if (left && right) straddle += Weight(t, coeff, p, n, gamma);
                                    }
                                cut = Math.Min(cut, straddle);
                            }
                            bool within = lam <= 2 * cut * (1 + 4 * n * eps);
                            if (chainPairs)
                            {
                                cutReadings++;
                                if (!within) { cutMismatches++; cutFirst ??= Name(true, n); }
                                if (Math.Abs(lam - 2 * cut) <= 4 * n * eps * Math.Max(lam, 2 * cut)) CutEqualities++;
                                if (lam > 0) { CutBroken++; CutTightest = Math.Max(CutTightest, lam / (2 * cut)); }
                            }
                            else if (!within) CutWhereChainBreaksMismatches++;   // the control
                        }
                    }
                }
            families.Add(new FamilyReading("Theorem H by the trellis: the code distance on the convolutional code's trellis against the enumeration over φ, window pairs at random rates, chain N = 3 to 12 and ring N = 5 to 12", trellisCells, trellisReadings, trellisMismatches, trellisFirst, "the enumeration"));
            families.Add(new FamilyReading("Corollary H1, the packing bound: where Theorem B breaks, the trellis at least twice the least weights of Theorem B's relation and its translates by 1 + max deg q_i, chain N = N* to 3N*", packCells, packReadings, packMismatches, packFirst, "the packing bound"));
            families.Add(new FamilyReading("Corollary H1, the cut bound: on the ring where the chain of the same N pairs, the trellis at most twice the weight over the lightest cut, ring N = 5 to 12", cutCells, cutReadings, cutMismatches, cutFirst, "the cut bound"));
        }

        // Corollary A″: ket and bra agreeing on a set Λ of sites, the balance with Λ's Y and Z letters struck against Theorem 2
        // with the complement taken on the sites outside Λ (ComplementConnectionGraph.AgreeOn), every single site and, to six
        // sites, every set of two
        {
            var two = (from x in "XYZ" from y in "XYZ" select $"{x}{y}").Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            var three = Windows.Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            var cellList = two.Concat(three).Select(t => new[] { (t, 1L) })
                .Concat(from a in two from b in three
                        where !SameShape(SignedMaskBalanceClaim.Mask(a), SignedMaskBalanceClaim.Mask(b))
                        select new[] { (a, 1L), (b, 1L) }).ToList();
            int ringTop = Math.Min(RingCap, 8);
            static IEnumerable<int> Sets(int n) => Enumerable.Range(0, n).Select(l => 1 << l)
                .Concat(n <= 6 ? from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (1 << a) | (1 << b) : Enumerable.Empty<int>());
            static int Components(int n, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed) =>
                1 << (n - Rank(placed.Where(p => p.X != 0).Select(p => (int)p.X).Distinct().ToArray()));
            static Dictionary<int, GaussianInteger>[] ColumnsOf(int n, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed) =>
                ComplementConnectionGraph.Columns(n, placed.Select(p => (new PauliMask(p.X, p.Z), (BigInteger)p.C)));
            static string Sites(int n, int set) => string.Join(",", Enumerable.Range(0, n).Where(l => ((set >> l) & 1) == 1));
            int cells = 0, readings = 0, mismatches = 0, nonzero = 0, zero = 0, whereFarIsZero = 0, unstruck = 0, wrongSite = 0;
            string? first = null;
            foreach (var cell in cellList)
            {
                cells++;
                int w = cell.Max(t => t.Item1.Length);
                foreach (bool ring in new[] { false, true })
                    for (int n = ring ? Math.Max(3, 2 * w - 1) : Math.Max(3, w); n <= (ring ? ringTop : 6); n++)
                    {
                        var placed = SignedMaskBalanceClaim.PlaceStrings(cell, n, ring);
                        var columns = ColumnsOf(n, placed);
                        int components = Components(n, placed);
                        int far = ComplementConnectionGraph.AgreeOn(n, columns, 0);
                        bool unstruckBalanced = SignedMaskBalanceClaim.Balanced(placed.Where(p => p.X != 0).Select(p => (p.X, BitOperations.PopCount(p.Z) & 1)));
                        foreach (int agree in Sets(n))
                        {
                            if (SignedMaskBalanceClaim.AgreeOnBalanced(placed, (ulong)agree) is not { } balanced) continue;
                            int walked = ComplementConnectionGraph.AgreeOn(n, columns, agree);
                            readings++;
                            if (walked != (balanced ? components : 0))
                            {
                                mismatches++;
                                first ??= $"{string.Join(" + ", cell.Select(t => t.Item1))} on the {(ring ? "ring" : "chain")} of {n}, agreeing on {Sites(n, agree)}";
                            }
                            if (walked > 0) nonzero++; else zero++;
                            if (walked > 0 && far == 0) whereFarIsZero++;
                            if (walked != (unstruckBalanced ? components : 0)) unstruck++;
                            if (BitOperations.PopCount((uint)agree) == 1
                                && SignedMaskBalanceClaim.AgreeOnBalanced(placed, 1UL << ((BitOperations.TrailingZeroCount(agree) + 1) % n)) is { } moved
                                && walked != (moved ? components : 0)) wrongSite++;
                        }
                    }
            }
            AgreeOnNonzero = nonzero;
            AgreeOnZero = zero;
            AgreeOnWhereFarIsZero = whereFarIsZero;
            AgreeOnUnstruckMismatches = unstruck;
            AgreeOnWrongSiteMismatches = wrongSite;
            families.Add(new FamilyReading("Corollary A″: ket and bra agreeing on some sites, the balance with their Y and Z letters struck against Theorem 2 with the complement taken on the sites outside them, one template of two or three sites and a two-site template beside a three-site one of a different shape, every single site (chain N = 3 to 6, ring N = 2w − 1 to the cap or 8, the smaller) and every set of two (to six sites)", cells, readings, mismatches, first));

            // its hypotheses: XY beside a field Z on site 0 alone (the field even off Λ exactly where Λ holds site 0) and XY
            // beside ZZ on every bond (odd off every Λ of one or two sites), against Theorem 2's count where the balance reads and
            // with the diagonal strings ignored where it declines; XX + XY and XZ + ZY, two strings on one flip set, where it must decline
            int hcells = 0, hreadings = 0, hmismatches = 0, outside = 0, diagonalIgnored = 0;
            string? hfirst = null;
            foreach (bool ring in new[] { false, true })
                for (int n = 3; n <= 6; n++)
                {
                    var rows = new List<(string Name, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> Placed, bool Shared)>
                    {
                        ("XY + Z on site 0", SignedMaskBalanceClaim.PlaceStrings(new[] { ("XY", 10L) }, n, ring).Append(new SignedMaskBalanceClaim.PlacedString(0, 1, 9)).ToList(), false),
                        ("XY + ZZ", SignedMaskBalanceClaim.PlaceStrings(new[] { ("XY", 10L), ("ZZ", 9L) }, n, ring), false),
                        ("XX + XY", SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 1L), ("XY", 1L) }, n, ring), true),
                        ("XZ + ZY", SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, n, ring), true),
                    };
                    foreach (var (name, placed, shared) in rows)
                    {
                        hcells++;
                        var columns = ColumnsOf(n, placed);
                        int components = Components(n, placed);
                        foreach (int agree in Sets(n))
                        {
                            var verdict = SignedMaskBalanceClaim.AgreeOnBalanced(placed, (ulong)agree);
                            int walked = ComplementConnectionGraph.AgreeOn(n, columns, agree);
                            string where = $"{name} on the {(ring ? "ring" : "chain")} of {n}, agreeing on {Sites(n, agree)}";
                            if (verdict is { } v)
                            {
                                hreadings++;
                                if (shared || walked != (v ? components : 0)) { hmismatches++; hfirst ??= where; }
                                continue;
                            }
                            outside++;
                            if (shared) continue;
                            bool flipsOnly = SignedMaskBalanceClaim.Balanced(placed.Where(p => p.X != 0).Select(p => (p.X, BitOperations.PopCount(p.Z & ~(ulong)agree) & 1)));
                            if (walked != (flipsOnly ? components : 0)) diagonalIgnored++;
                        }
                    }
                }
            AgreeOnOutside = outside;
            AgreeOnDiagonalIgnoredMismatches = diagonalIgnored;
            families.Add(new FamilyReading("Corollary A″'s hypotheses: XY beside a field Z on site 0 (read where site 0 is in the set) and beside ZZ on every bond (declined), and XX + XY and XZ + ZY, two strings on one flip set (declined); chain and ring N = 3 to 6, every single site and every set of two", hcells, hreadings, hmismatches, hfirst));
        }

        // the triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J (the page's proposition): two of its premises exactly, on integer bonds (the three bond strings
        // anticommute in pairs, so H² = |j|²·1 and the commutator with H obeys ad³ = 4|j|²·ad, its eigenvalues 0 and ±2|j|
        // alone); then the closed form against the strong-coupling rule computed directly (degenerate perturbation in 1/J, the
        // rule of THE_SPREAD_IS_A_RESONANCE), the rates compressed onto each eigenspace of the
        // commutator in the Pauli basis, the eigenprojections the polynomials 1 − ad²/(4|j|²) and (ad² ± 2|j|·ad)/(8|j|²), at
        // random signed bonds and rates spread over four decades
        {
            var rng = new Random(20261121);
            var bondStrings = new[] { PauliMask.Parse("XYI"), PauliMask.Parse("IXY"), PauliMask.Parse("YIX") };
            var strings = (from x in Enumerable.Range(0, 8) from z in Enumerable.Range(0, 8) select new PauliMask((ulong)x, (ulong)z)).ToArray();
            var index = Enumerable.Range(0, strings.Length).ToDictionary(i => strings[i], i => i);
            int premiseReadings = 3, premiseMismatches = 0, profiles = 0, mismatches = 0, oddWins = 0, corollaryWins = 0, rotated = 0, oddOnly = 0;
            double worst = 0;
            string? premiseFirst = null, first = null;
            for (int a = 0; a < 3; a++)
                if (PauliMask.Commute(bondStrings[a], bondStrings[(a + 1) % 3])) { premiseMismatches++; premiseFirst ??= $"the bonds {a} and {(a + 1) % 3} commute"; }
            for (int draw = 0; draw < 240; draw++)
            {
                var j = Enumerable.Range(0, 3).Select(_ => (long)rng.Next(1, 7) * (rng.Next(2) == 0 ? 1 : -1)).ToArray();   // every bond present: a ring
                long norm2 = j.Sum(v => v * v);
                var g = Enumerable.Range(0, 3).Select(_ => Math.Exp(Math.Log(0.01) + rng.NextDouble() * Math.Log(1e4))).ToArray();
                double sigma = g.Sum();

                // H² exactly, in Gaussian integers
                premiseReadings++;
                var square = new Dictionary<PauliMask, GaussianInteger>();
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                    {
                        var (prod, k) = PauliMask.Multiply(bondStrings[a], bondStrings[b]);
                        AddTerm(square, prod, IPower(new GaussianInteger(j[a] * j[b], BigInteger.Zero), k));
                    }
                if (square.Count != 1 || !square.TryGetValue(PauliMask.Identity, out var id) || id != new GaussianInteger(norm2, BigInteger.Zero))
                { premiseMismatches++; premiseFirst ??= $"H² at the bonds ({string.Join(", ", j)})"; }

                // the strong-coupling rule, computed directly
                profiles++;
                var ad = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.Dense(64, 64);
                for (int c = 0; c < 64; c++)
                    for (int a = 0; a < 3; a++)
                    {
                        if (j[a] == 0 || PauliMask.Commute(bondStrings[a], strings[c])) continue;
                        var (prod, k) = PauliMask.Multiply(bondStrings[a], strings[c]);   // [b, P] = 2·b·P where they anticommute
                        var turn = IPower(GaussianInteger.One, k);
                        ad[index[prod], c] += 2.0 * j[a] * new Complex((double)turn.Re, (double)turn.Im);
                    }
                var rates = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.DenseDiagonal(64, 64,
                    c => 2 * Enumerable.Range(0, 3).Where(l => ((strings[c].X >> l) & 1) == 1).Sum(l => g[l]));
                var ad2 = ad * ad;
                var one = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.DenseIdentity(64);
                double r1 = Math.Sqrt(norm2);
                var projectors = new[] { one - ad2 / (4.0 * norm2), (ad2 + ad * (2 * r1)) / (8.0 * norm2), (ad2 - ad * (2 * r1)) / (8.0 * norm2) };
                double top = projectors.Max(pk => (pk * rates * pk).Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Max(z => z.Real));
                double computed = 2 * sigma - top;
                double tol = 64 * 2.220446049250313e-16 * 2 * sigma;
                double closed = SignedMaskBalanceClaim.TriangleFarEdgeAtLargeJ(j[0], j[1], j[2], g[0], g[1], g[2]);
                double odd = SignedMaskBalanceClaim.TriangleKernelOddShortfall(j.Select(v => (double)v).ToArray(), new[] { g[2], g[0], g[1] });
                worst = Math.Max(worst, Math.Abs(computed - closed) / tol);
                if (Math.Abs(computed - closed) > tol) { mismatches++; first ??= $"bonds ({string.Join(", ", j)}), rates ({string.Join(", ", g.Select(v => v.ToString("G4")))})"; }
                if (2 * g.Min() < odd) corollaryWins++; else oddWins++;
                double onBond = Math.Min(2 * g.Min(), SignedMaskBalanceClaim.TriangleKernelOddShortfall(j.Select(v => (double)v).ToArray(), new[] { g[0], g[1], g[2] }));
                if (Math.Abs(computed - onBond) > tol) rotated++;
                if (Math.Abs(computed - odd) > tol) oddOnly++;
            }
            TriangleOddWins = oddWins;
            TriangleCorollaryWins = corollaryWins;
            TriangleWorstRatio = worst;
            TriangleRotatedMismatches = rotated;
            TriangleOddOnlyMismatches = oddOnly;
            families.Add(new FamilyReading("the triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J: the three bond strings anticommute in pairs and H² = |j|²·1, exactly, on integer bonds from ±1 to ±6", 1, premiseReadings, premiseMismatches, premiseFirst, "exact arithmetic"));
            families.Add(new FamilyReading("the triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J: the closed form against the rates compressed onto the three eigenspaces of the commutator with H (degenerate perturbation in 1/J), random signed bonds and rates over four decades", profiles, profiles, mismatches, first, "the compression"));
        }

        // sums over disjoint bonds (the page's paragraph on one rate): on the ring of bonds X_l·Y_(l+1), 3 to 11 sites, unit coefficients and two
        // random draws, every Q_k and Q_k·Z^⊗N commutes with H, in Gaussian integers, and H² = Σ c_l² + 2·Q₂ from four sites on;
        // on the odd rings each top term alone does not commute, its commutator moving its one site of I or Z by one place with
        // the letter exchanged, and at unit coefficients [H, S_l] = 2i·(T_(l−1) − T_(l+1)) and [H, T_l] = 2i·(S_(l−1) − S_(l+1))
        // exactly, S_l the product of the (N − 1)/2 disjoint bonds that cover every site but l and T_l = S_l·Z^⊗N; on the even
        // rings each top term commutes alone
        {
            var rng = new Random(11);
            int cells = 0, readings = 0, mismatches = 0, moving = 0, evenCommuting = 0;
            string? first = null;
            for (int n = 3; n <= 11; n++)
                for (int draw = 0; draw < 3; draw++)
                {
                    cells++;
                    var c = Enumerable.Range(0, n).Select(_ => draw == 0 ? 1L : rng.Next(1, 10) * (rng.Next(2) == 0 ? 1L : -1L)).ToArray();
                    var h = Enumerable.Range(0, n).Select(l => (S: RingBond(n, l), C: c[l])).ToArray();
                    var parity = new PauliMask(0, (1UL << n) - 1);
                    if (n >= 4)
                    {
                        // H² = Σ c_l² + 2·Q₂: two bonds anticommute where they meet and commute where they do not
                        readings++;
                        var square = new Dictionary<PauliMask, GaussianInteger>();
                        foreach (var (s1, c1) in h)
                            foreach (var (s2, c2) in h)
                            {
                                var (prod, k) = PauliMask.Multiply(s1, s2);
                                AddTerm(square, prod, IPower(new GaussianInteger(c1 * c2, BigInteger.Zero), k));
                            }
                        AddTerm(square, PauliMask.Identity, new GaussianInteger(-c.Sum(v => v * v), BigInteger.Zero));
                        foreach (var (s, q2) in SignedMaskBalanceClaim.DisjointBondSum(c, 2))
                            AddTerm(square, s, new GaussianInteger(-2 * q2, BigInteger.Zero));
                        if (square.Count != 0) { mismatches++; first ??= $"H² against Σ c² + 2·Q₂ on the ring of {n}"; }
                    }
                    if (draw == 0 && n % 2 == 1)
                    {
                        // the moves exactly, at unit coefficients: S[l] the term whose bonds cover every site but l, T[l] = S[l]·Z^⊗N with phase 1
                        var top = SignedMaskBalanceClaim.DisjointBondSum(c, (n - 1) / 2);
                        var sTerm = new PauliMask[n];
                        var tTerm = new PauliMask[n];
                        foreach (var (s, coefficient) in top)
                        {
                            int l = BitOperations.TrailingZeroCount(~s.X & ((1UL << n) - 1));
                            sTerm[l] = s;
                            var (prod, k) = PauliMask.Multiply(s, parity);
                            tTerm[l] = prod;
                            if (coefficient != BigInteger.One || k != 0) { mismatches++; first ??= $"a top term or its product with Z^⊗N off phase 1 on the ring of {n}"; }
                        }
                        var two = new GaussianInteger(BigInteger.Zero, 2);
                        for (int l = 0; l < n; l++)
                        {
                            int before = (l + n - 1) % n, after = (l + 1) % n;
                            foreach (var (term, target) in new[] { (sTerm, tTerm), (tTerm, sTerm) })
                            {
                                readings++;
                                var expected = new Dictionary<PauliMask, GaussianInteger>();
                                AddTerm(expected, target[before], two);
                                AddTerm(expected, target[after], -two);
                                var got = Commutator(h, new Dictionary<PauliMask, GaussianInteger> { [term[l]] = GaussianInteger.One });
                                if (got.Count != expected.Count || got.Any(kv => !expected.TryGetValue(kv.Key, out var e) || e != kv.Value))
                                { mismatches++; first ??= $"the move of the top term at {l} on the ring of {n}"; }
                            }
                        }
                    }
                    for (int k = 0; 2 * k <= n; k++)
                    {
                        var q = SignedMaskBalanceClaim.DisjointBondSum(c, k).ToDictionary(kv => kv.Key, kv => new GaussianInteger(kv.Value, BigInteger.Zero));
                        foreach (bool withParity in new[] { false, true })
                        {
                            var op = withParity ? TimesString(q, parity) : q;
                            readings++;
                            if (Commutator(h, op).Count != 0) { mismatches++; first ??= $"the ring of {n}, k = {k}{(withParity ? ", times Z^⊗N" : "")}"; }
                            if (2 * k == n - 1)
                                foreach (var (s, coefficient) in op)
                                {
                                    int l = BitOperations.TrailingZeroCount(~s.X & ((1UL << n) - 1));
                                    bool zAtDefect = ((s.Z >> l) & 1) == 1;
                                    var moved = Commutator(h, new Dictionary<PauliMask, GaussianInteger> { [s] = coefficient });
                                    bool ok = moved.Count > 0 && moved.Keys.All(o =>
                                    {
                                        ulong gap = ~o.X & ((1UL << n) - 1);
                                        if (BitOperations.PopCount(gap) != 1) return false;
                                        int m = BitOperations.TrailingZeroCount(gap);
                                        return (m == (l + 1) % n || m == (l + n - 1) % n) && (((o.Z >> m) & 1) == 1) != zAtDefect;
                                    });
                                    if (ok) moving++; else { mismatches++; first ??= $"a top term on the ring of {n} whose commutator does not move its agreeing site"; }
                                }
                            if (2 * k == n)
                                foreach (var (s, coefficient) in op)
                                {
                                    if (Commutator(h, new Dictionary<PauliMask, GaussianInteger> { [s] = coefficient }).Count == 0) evenCommuting++;
                                    else { mismatches++; first ??= $"a top term on the even ring of {n} that does not commute alone"; }
                                }
                        }
                    }
                }
            OddTopTermsMoving = moving;
            EvenTopTermsCommuting = evenCommuting;
            families.Add(new FamilyReading("sums over disjoint bonds on the ring of bonds X_l·Y_(l+1): every Q_k and Q_k·Z^⊗N commutes with H, N = 3 to 11, unit coefficients and two random draws, H² = Σ c² + 2·Q₂ from N = 4, on the odd rings each top term alone moves its agreeing site, and at unit coefficients [H, S_l] = 2i·(T_(l−1) − T_(l+1)) and [H, T_l] = 2i·(S_(l−1) − S_(l+1)) (exact)", cells, readings, mismatches, first, "exact arithmetic"));
        }

        // the operators with N − 1 X or Y letters at one rate: on the odd ring of bonds X_l·Y_(l+1) those that commute with H number at most
        // the kernel modulo p of i·[H, ·] on that plane (real integer entries), and at least Corollary A″'s count on each site
        // plus the two sums above, independent of the site-resolved ones at any nonzero coefficients, the commutator of each top
        // term holding two distinct strings, those of an S term with the letters of a T term and the reverse; where the two
        // bounds meet the count is exact
        {
            var rng = new Random(13);
            const long p = 1_000_000_007;
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            for (int n = 3; n <= 9; n += 2)
                for (int draw = 0; draw < 2; draw++)
                {
                    cells++;
                    readings++;
                    var c = Enumerable.Range(0, n).Select(_ => draw == 0 ? 1L : rng.Next(1, 10) * (rng.Next(2) == 0 ? 1L : -1L)).ToArray();
                    var columns = ComplementConnectionGraph.Columns(n, Enumerable.Range(0, n).Select(l => (RingBond(n, l), (BigInteger)c[l])));
                    int lower = Enumerable.Range(0, n).Sum(l => ComplementConnectionGraph.AgreeOn(n, columns, 1 << l)) + 2;
                    var rows = new Dictionary<PauliMask, int>();
                    var cols = new List<Dictionary<int, long>>();
                    ulong full = (1UL << n) - 1;
                    for (int l = 0; l < n; l++)
                    {
                        ulong x = full & ~(1UL << l);
                        for (ulong y = x; ; y = (y - 1) & x)   // the Y letters, every subset of the X/Y sites
                        {
                            foreach (ulong d in new[] { 0UL, 1UL << l })   // I or Z on site l
                            {
                                var s = new PauliMask(x, y | d);
                                var col = new Dictionary<int, long>();
                                for (int b = 0; b < n; b++)
                                {
                                    var bond = RingBond(n, b);
                                    if (PauliMask.Commute(bond, s)) continue;
                                    var (prod, k) = PauliMask.Multiply(bond, s);   // [bond, s] = 2·i^k·prod with k odd, so i·[bond, s] is real
                                    if (!rows.TryGetValue(prod, out int r)) rows[prod] = r = rows.Count;
                                    col[r] = (col.TryGetValue(r, out var v) ? v : 0) + 2 * c[b] * (k == 1 ? -1 : k == 3 ? 1 : 0);
                                }
                                cols.Add(col);
                            }
                            if (y == 0) break;
                        }
                    }
                    int upper = cols.Count - SparseRankModP(cols, p);
                    if (upper != lower || lower != 2 * n + 2) { mismatches++; first ??= $"the ring of {n}, draw {draw}: at most {upper}, at least {lower}"; }
                }
            families.Add(new FamilyReading("N − 1 X or Y letters at one rate: on the odd ring of bonds X_l·Y_(l+1) the operators with X or Y on N − 1 sites that commute with H, exactly 2N + 2 (Corollary A″'s 2N and the two sums over disjoint bonds), N = 3 to 9, unit coefficients and a random draw", cells, readings, mismatches, first, "both bounds"));
        }

        // the worked examples
        {
            var checks = new List<(string Name, bool Holds)>
            {
                ("IXZX + XIZX: N* = 6", SignedMaskBalanceClaim.Chain("IXZX", "XIZX") == (true, 6)),
                ("IXZX + XIZX pairs at N = 5", PairsOnChain(5, ("IXZX", 1), ("XIZX", 1))),
                ("IXZX + XIZX breaks at N = 6, 7", !PairsOnChain(6, ("IXZX", 1), ("XIZX", 1)) && !PairsOnChain(7, ("IXZX", 1), ("XIZX", 1))),
                ("XZ + XIXIX: N* = 6", SignedMaskBalanceClaim.Chain("XZ", "XIXIX") == (true, 6)),
                ("XY + XXIX: N* = 5", SignedMaskBalanceClaim.Chain("XY", "XXIX") == (true, 5)),
                ("+XIY + XXX: N* = 5", SignedMaskBalanceClaim.Chain("XIY", "XXX") == (true, 5)),
                ("+XIY + XIZ: N* = 5", SignedMaskBalanceClaim.Chain("XIY", "XIZ") == (true, 5)),
                ("+XIY + XXZ: N* = 4", SignedMaskBalanceClaim.Chain("XIY", "XXZ") == (true, 4)),
                ("+XIX + XXX never breaks", SignedMaskBalanceClaim.Chain("XIX", "XXX") == (false, null)),
                ("XX + XZ + XXX breaks from N = 3, where each of its three pairs pairs (Theorem B′, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XXX" }) == (true, 3)
                    && SignedMaskBalanceClaim.Chain("XX", "XZ") == (false, null) && SignedMaskBalanceClaim.Chain("XX", "XXX") == (false, null)
                    && SignedMaskBalanceClaim.Chain("XZ", "XXX") == (true, 4)
                    && !PairsOnChain(3, ("XX", 1), ("XZ", 1), ("XXX", 1)) && PairsOnChain(3, ("XZ", 1), ("XXX", 1))
                    && PairsOnChain(3, ("XX", 1), ("XZ", 1)) && PairsOnChain(3, ("XX", 1), ("XXX", 1))),
                ("XZ beside X·I·(X + Y) and X·X·(X + Y) breaks at N = 3 by the holonomy of (x, 1, 1), where each of its pairs pairs (Corollary E3′, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.ChainShapes(new IReadOnlyList<(string, long)>[] { new[] { ("XZ", 1L) }, new[] { ("XIX", 1L), ("XIY", 1L) }, new[] { ("XXX", 1L), ("XXY", 1L) } }, 3)
                        == (true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus)
                    && !PairsOnChain(3, ("XZ", 1), ("XIX", 1), ("XIY", 1), ("XXX", 1), ("XXY", 1))
                    && PairsOnChain(3, ("XZ", 1), ("XIX", 1), ("XIY", 1)) && PairsOnChain(3, ("XZ", 1), ("XXX", 1), ("XXY", 1))
                    && PairsOnChain(3, ("XIX", 1), ("XIY", 1), ("XXX", 1), ("XXY", 1))),
                ("XX + XZ + XIY breaks from N = 3, which its unreduced basis reads as 4 (Theorem B′ and its control, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XIY" }) == (true, 3)
                    && SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XIY" }, reduced: false) == (true, 4)
                    && !PairsOnChain(3, ("XX", 1), ("XZ", 1), ("XIY", 1))),
                ("XXI + XIX pairs at N = 4 to 6", Enumerable.Range(4, 3).All(n => PairsOnChain(n, ("XXI", 1), ("XIX", 1)))),
                ("XXI + XIY pairs at 3 and breaks at 4", PairsOnChain(3, ("XXI", 1), ("XIY", 1)) && !PairsOnChain(4, ("XXI", 1), ("XIY", 1))),
                ("XXI + XIX pairs on the ring of five", PairsOnRing(5, ("XXI", 1), ("XIX", 1))),
                ("+XXZ ± XZY and its seven relatives pair exactly at even N (Theorem D, N = 5 to 64)",
                    new[] { ("XXZ", "XZY"), ("XXZ", "YZX"), ("XZY", "YYZ"), ("XZY", "ZXX"), ("XZY", "ZYY"), ("YYZ", "YZX"), ("YZX", "ZXX"), ("YZX", "ZYY") }
                        .All(p => Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.RingBreaks(new[] { p.Item1, p.Item2 }, n) == (n % 2 == 1)))),
                ("XX + XIY breaks on every ring (Theorem D, N = 5 to 64)", Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIY" }, n))),
                ("XXX + XYX breaks on the ring exactly at 3 | N (Corollary E2, N = 5 to 63; Theorem 2 to 9)",
                    Enumerable.Range(5, 59).All(n => SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, n).Breaks == (n % 3 == 0))
                    && Enumerable.Range(5, 5).All(n => PairsOnRing(n, ("XXX", 1), ("XYX", 1)) == (n % 3 != 0))),
                ("XXX + XYX pairs on the chain, N = 3 to 7", Enumerable.Range(3, 5).All(n => PairsOnChain(n, ("XXX", 1), ("XYX", 1)))),
                ("XZ + ZY pairs on the ring and breaks on the chain, N = 3 to 8",
                    Enumerable.Range(3, 6).All(n => PairsOnRing(n, ("XZ", 1), ("ZY", 1)) && !PairsOnChain(n, ("XZ", 1), ("ZY", 1)))
                    && Enumerable.Range(3, 6).All(n => SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, n, true)) == SignedMaskBalanceClaim.CharacterVerdict.Pairs)),
                ("XXX + YYY pairs on every ring, relations at 3 | N (Corollary E2, N = 5 to 63)",
                    Enumerable.Range(5, 59).All(n => !SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("YYY", 1L) }, n).Breaks)),
                ("XIY + XZX pairs on the ring exactly at 4 | N (Corollary E2, N = 5 to 63; Theorem 2 to 9)",
                    Enumerable.Range(5, 59).All(n => SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, n).Breaks == (n % 4 != 0))
                    && Enumerable.Range(5, 5).All(n => PairsOnRing(n, ("XIY", 1), ("XZX", 1)) == (n % 4 == 0))),
                ("XXX + YYY beside the field X + Y on every site pairs on every chain, and XXX − YYY beside it breaks from N = 3 (Corollary E3 to 64, Theorem 2 to 8)",
                    Enumerable.Range(3, 62).All(n => !SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("X", 1L), ("Y", 1L) }, new[] { ("XXX", 1L), ("YYY", 1L) }, n).Breaks
                                                    && SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("X", 1L), ("Y", 1L) }, new[] { ("XXX", 1L), ("YYY", -1L) }, n).Breaks)
                    && Enumerable.Range(3, 6).All(n => PairsOnChain(n, ("X", 1), ("Y", 1), ("XXX", 1), ("YYY", 1)) && !PairsOnChain(n, ("X", 1), ("Y", 1), ("XXX", 1), ("YYY", -1)))),
                ("X·(X + Y)·X beside XIIX pairs at N = 3, where XIIX is not placed yet (Theorem 2), breaks at N* = 4 by the read sum, and from 5 by a square (Corollary E3 to 64, Theorem 2 to 8)",
                    SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XXX", 1L), ("XYX", 1L) }, new[] { ("XIIX", 1L) }, 4).Branch == SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero
                    && Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XXX", 1L), ("XYX", 1L) }, new[] { ("XIIX", 1L) }, n).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail)
                    && PairsOnChain(3, ("XXX", 1), ("XYX", 1)) && Enumerable.Range(4, 5).All(n => !PairsOnChain(n, ("XXX", 1), ("XYX", 1), ("XIIX", 1)))),
            };
            int bad = checks.Count(c => !c.Holds);
            families.Add(new FamilyReading("the worked examples of the page", checks.Count, checks.Count, bad, checks.FirstOrDefault(c => !c.Holds).Name, "the page"));
        }

        Families = families;
        SwappedWidthsMismatches = swappedWidths;
        SwappedQuotientsMismatches = swappedQuotients;
        DroppedRingPowerMismatches = droppedPower;
        ParityOnlyMismatches = parityOnly;
        MostValuedMismatches = mostValued;
        CharacterOutside = characterOutside;
        DroppedConditionMismatches = droppedCondition;
        DroppedReadSetMismatches = droppedReadSet;
        DroppedOrderingMismatches = droppedOrdering;
        DroppedReadSumMismatches = droppedReadSum;
        DroppedE3OrderingMismatches = droppedE3Ordering;
        UnreducedBasisMismatches = unreducedBasis;
        DroppedCosetTermMismatches = droppedCosetTerm;
        DroppedE3PrimeReadSumMismatches = droppedE3PrimeReadSum;
        DroppedE3PrimeOrderingMismatches = droppedE3PrimeOrdering;
    }

    // a ℤ/4 quadratic refinement of β(a, b) = ⟨R(a), b⟩ on all of W: Q(a ⊕ b) = Q(a) + Q(b) + 2β(a, b), β symmetric
    static bool IsRefinement(IReadOnlyDictionary<ulong, (int Q, ulong R)> form)
    {
        static int Beta(ulong r, ulong b) => BitOperations.PopCount(r & b) & 1;
        foreach (var (a, (qa, ra)) in form)
            foreach (var (b, (qb, rb)) in form)
                if (Beta(ra, b) != Beta(rb, a) || form[a ^ b].Q != (qa + qb + 2 * Beta(ra, b)) % 4) return false;
        return true;
    }

    // the carrier X^⊗N·diag(i^(d_x)) with d₀ = 0 is a Clifford (up to a global phase, which d₀ = 0 fixes) exactly when
    // Q = −d is a ℤ/4 quadratic refinement on all bitstrings: Q(0) = 0 and Q(x ⊕ e_l) − Q(x) − Q(e_l) = 2·B(x, e_l) at
    // every x and site l, B read off the basis vectors (symmetric by construction) and bilinear (Q is then fixed by its
    // values on the e_l and by B, as the refinement with those values is)
    static bool QuartersQuadratic(int[] d, int n)
    {
        int Q(int x) => ((-d[x]) % 4 + 4) % 4;
        int Step(int x, int l) => ((Q(x ^ (1 << l)) - Q(x) - Q(1 << l)) % 4 + 8) % 4;
        if (Q(0) != 0) return false;
        var col = new int[n];
        for (int l = 0; l < n; l++)
            for (int k = 0; k < n; k++)
            {
                int s = Step(1 << k, l);
                if ((s & 1) == 1) return false;
                if (s == 2) col[l] |= 1 << k;
            }
        for (int l = 0; l < n; l++)
            for (int x = 0; x < 1 << n; x++)
                if (Step(x, l) != 2 * (BitOperations.PopCount((uint)(x & col[l])) & 1)) return false;
        return true;
    }

    // the control: the carrier moved by one quarter on the component through the largest bitstring off W (a phase per
    // component, so (★) must still hold) stays quadratic exactly where W has codimension one (its parity then moves by a
    // linear form; at codimension two or more by the indicator of a smaller coset, of higher degree); true on a breach
    static bool MovedCarrierBreaches(IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed, int n, int[] carrier,
                                     IReadOnlyDictionary<ulong, (int Q, ulong R)> form, ref int notQuadratic)
    {
        int codim = n - BitOperations.Log2((uint)form.Count);
        if (codim == 0) return false;
        int x1 = Enumerable.Range(0, 1 << n).Reverse().First(x => !form.ContainsKey((ulong)x));
        var moved = (int[])carrier.Clone();
        foreach (var w in form.Keys) moved[x1 ^ (int)w] = (moved[x1 ^ (int)w] + 1) % 4;
        bool quadratic = QuartersQuadratic(moved, n);
        if (!quadratic) notQuadratic++;
        return !StarHolds(placed, n, moved) || quadratic != (codim == 1);
    }

    // Theorem 2's far condition (★) for X^⊗N·diag(i^(d_x)): H_x̄ȳ·i^(d_y) = i^(d_x)·H_xy at every entry, in Gaussian integers
    static bool StarHolds(IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed, int n, int[] d)
    {
        var h = ComplementConnectionGraph.Columns(n, placed.Select(s => (new PauliMask(s.X, s.Z), (BigInteger)s.C)));
        int full = (1 << n) - 1;
        static GaussianInteger Turn(GaussianInteger z, int k) { for (int j = 0; j < k; j++) z = z * GaussianInteger.I; return z; }
        for (int y = 0; y < h.Length; y++)
            foreach (var (x, hxy) in h[y])
            {
                if (Turn(ComplementConnectionGraph.Entry(h, x ^ full, y ^ full), d[y]) != Turn(hxy, d[x])) return false;
                if (Turn(hxy, d[y ^ full]) != Turn(ComplementConnectionGraph.Entry(h, x ^ full, y ^ full), d[x ^ full])) return false;
            }
        return true;
    }

    // the rank over GF(2) of a relation basis's leading coefficient matrix for the shift by the template widths
    static int LeadingRank(IReadOnlyList<ulong[]> basis, string[] ts)
    {
        var rows = new List<ulong>();
        foreach (var b in basis)
        {
            int top = Enumerable.Range(0, ts.Length).Where(i => b[i] != 0).Max(i => BitOperations.Log2(b[i]) + ts[i].Length);
            ulong row = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                int e = top - ts[i].Length;
                if (e >= 0 && ((b[i] >> e) & 1) == 1) row |= 1UL << i;
            }
            rows.Add(row);
        }
        int rank = 0;
        for (int bit = 0; bit < ts.Length; bit++)
        {
            int pivot = rows.FindIndex(rank, r => ((r >> bit) & 1) == 1);
            if (pivot < 0) continue;
            (rows[rank], rows[pivot]) = (rows[pivot], rows[rank]);
            for (int r = 0; r < rows.Count; r++) if (r != rank && ((rows[r] >> bit) & 1) == 1) rows[r] ^= rows[rank];
            rank++;
        }
        return rank;
    }

    // the controls of Theorem D, read the proof's way from one template i₀: 2^(v₂(N)) read off the parity of N alone,
    // and i₀ taken at the greatest valuation
    static bool RingParityOnly(string[] ts, int n) => RingRead(ts, n % 2 == 1 ? 1 : 64, least: true);
    static bool RingAtTheMostValued(string[] ts, int n) => RingRead(ts, 1 << BitOperations.TrailingZeroCount((uint)n), least: false);

    static bool RingRead(string[] ts, int big, bool least)
    {
        var read = ts.Select(t => (V: SignedMaskBalanceClaim.Valuation(SignedMaskBalanceClaim.Mask(t)), S: SignedMaskBalanceClaim.Sign(t))).ToList();
        var idx = Enumerable.Range(0, read.Count);
        int i0 = least ? idx.MinBy(i => read[i].V) : idx.MaxBy(i => read[i].V);
        var (v0, s0) = read[i0];
        if (v0 >= big) return read.Any(r => r.S == 1);
        return idx.Any(j => j != i0 && (read[j].S + s0 * (read[j].V == v0 ? 1 : 0)) % 2 == 1);
    }

    // the controls: Theorem B with its pieces swapped, and Theorem C without 2^(v₂(N))
    static int SwappedQuotientC(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return (SignedMaskBalanceClaim.Sign(a) * (BitOperations.PopCount(Div(m1, g)) & 1)
                + SignedMaskBalanceClaim.Sign(b) * (BitOperations.PopCount(Div(m2, g)) & 1)) % 2;
    }

    // N* = max(w₁ + deg q₂, w₂ + deg q₁) whatever c is
    static int NStar(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return Math.Max(a.Length + Deg(Div(m2, g)), b.Length + Deg(Div(m1, g)));
    }

    static int SwappedWidthsFrom(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return Math.Max(a.Length + Deg(Div(m1, g)), b.Length + Deg(Div(m2, g)));
    }

    static int Deg(ulong p) => 63 - BitOperations.LeadingZeroCount(p);

    static ulong Gcd(ulong a, ulong b)
    {
        while (b != 0)
        {
            ulong r = a;
            while (r != 0 && Deg(r) >= Deg(b)) r ^= b << (Deg(r) - Deg(b));
            (a, b) = (b, r);
        }
        return a;
    }

    static ulong Div(ulong a, ulong b)
    {
        ulong q = 0;
        while (a != 0 && Deg(a) >= Deg(b)) { int s = Deg(a) - Deg(b); q |= 1UL << s; a ^= b << s; }
        return q;
    }

    public string DisplayName => $"SignedMaskBalanceWitness (ring cap N = {RingCap})";

    public string Summary =>
        $"F173 on {Families.Count} families: {Families.Where(f => f.Against == "Theorem 2").Sum(f => f.Readings)} readings against Theorem 2 " +
        $"on the bitstrings, {Families.Where(f => f.Against == "the balance").Sum(f => f.Readings)} of Theorems D and B′ against the balance, " +
        $"{Families.Where(f => f.Against == "Theorem B").Sum(f => f.Readings)} of Theorem B′ against Theorem B, " +
        $"{Families.Where(f => f.Against == "Corollary E3").Sum(f => f.Readings)} of E3′ against E3, " +
        $"{Families.Where(f => f.Against == "Theorem D").Sum(f => f.Readings)} of the long chain against Theorem D, " +
        $"{Families.Where(f => f.Against == "the code distance").Sum(f => f.Readings)} of Theorem H's far edge against the code distance ({StrongDephasingBroken} broken at second order) and " +
        $"{Families.Where(f => f.Against == "the end squares").Sum(f => f.Readings)} against XZ ± ZY's two end squares, " +
        $"{Families.Where(f => f.Against == "the enumeration").Sum(f => f.Readings)} of the trellis against the enumeration, " +
        $"{Families.Where(f => f.Against == "the packing bound").Sum(f => f.Readings)} and {Families.Where(f => f.Against == "the cut bound").Sum(f => f.Readings)} against Corollary H1's two bounds, " +
        $"{Families.Where(f => f.Against == "the hop reading").Sum(f => f.Readings)} of the far edge against Corollary H3's hop reading, " +
        $"{Families.Where(f => f.Against == "the page").Sum(f => f.Readings)} worked checks, {Families.Sum(f => f.Mismatches)} disagreeing; the controls (widths swapped, quotients swapped, 2^(v₂(N)) dropped, " +
        $"the parity of N for 2^(v₂(N)), i₀ at the greatest valuation) disagree {SwappedWidthsMismatches}, {SwappedQuotientsMismatches}, " +
        $"{DroppedRingPowerMismatches}, {ParityOnlyMismatches}, {MostValuedMismatches} times; Theorem E with (1), (2), (3) dropped and E2 " +
        $"without its read-set test or its ordering term disagree {string.Join(", ", DroppedConditionMismatches)}, {DroppedReadSetMismatches}, " +
        $"{DroppedOrderingMismatches} times, E3 without its read-sum test or its ordering term {DroppedReadSumMismatches}, {DroppedE3OrderingMismatches} times, " +
        $"B′ off the unreduced basis {UnreducedBasisMismatches} times, the ℤ/4 carrier without its coset term {DroppedCosetTermMismatches} times, " +
        $"moved by a quarter on one component {MovedCarrierNotQuadratic} times no longer a Clifford, " +
        $"E3′ without N₁* or without its ordering term {DroppedE3PrimeReadSumMismatches}, {DroppedE3PrimeOrderingMismatches} times, " +
        $"Theorem H with the rates read away or reversed {RatesReadAwayMismatches}, {RatesReversedMismatches} times, " +
        $"Corollary H1's bounds with overlapping translates or where the chain breaks {PackingOverlapMismatches}, {CutWhereChainBreaksMismatches} times, " +
        $"Corollary H2 without its double count {CycleDoubleCountMismatches} times, Corollary H3's spectrum read as Theorem H string by string {HopControlMisses} times, " +
        $"Corollary A″ with no letter struck, with the next site's letters struck, or with the diagonal ignored {AgreeOnUnstruckMismatches}, {AgreeOnWrongSiteMismatches}, {AgreeOnDiagonalIgnoredMismatches} times; " +
        $"{Families.Where(f => f.Against == "exact arithmetic").Sum(f => f.Readings)} exact readings of two of the triangle's premises and of the sums over disjoint bonds, " +
        $"{Families.Where(f => f.Against == "the compression").Sum(f => f.Readings)} of the triangle at large J against the rates compressed onto the eigenspaces of the commutator with H, " +
        $"{Families.Where(f => f.Against == "both bounds").Sum(f => f.Readings)} of the operators with N − 1 X or Y letters at one rate; the triangle's closed form handed the rate of a site on each bond, " +
        $"or its odd term alone, disagree {TriangleRotatedMismatches}, {TriangleOddOnlyMismatches} times";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var f in Families)
                yield return new InspectableNode(
                    displayName: f.Name,
                    summary: $"{f.Cells} {(f.Cells == 1 ? "cell" : "cells")}, {f.Readings} readings, " +
                             (f.Agrees ? $"every one agreeing with {f.Against}" : $"{f.Mismatches} DISAGREEING, first {f.FirstMismatch}"),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the floors of F115's pairs at k = 4",
                summary: string.Join(", ", F115FloorsAtK4.OrderBy(kv => kv.Key).Select(kv => $"{kv.Value} pairs from N = {kv.Key}")) +
                         " (to be read beside THE_PALINDROME_CLASSIFIER's count of the pairs already hard at N = 5)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the controls, which must disagree",
                summary: $"N* with each width beside its own quotient: {SwappedWidthsMismatches} readings against Theorem 2; c with the " +
                         $"quotients swapped: {SwappedQuotientsMismatches}; the ring law without 2^(v₂(N)): {DroppedRingPowerMismatches}; " +
                         $"Theorem D with the parity of N in place of 2^(v₂(N)): {ParityOnlyMismatches}; Theorem D with i₀ at the greatest " +
                         $"valuation: {MostValuedMismatches} (the last two against the balance); Theorem E without the symmetric reading: " +
                         $"{DroppedConditionMismatches[0]}, without the read sets summed over the relations: {DroppedConditionMismatches[1]}, " +
                         $"without hol: {DroppedConditionMismatches[2]}; Corollary E2 without its read-set test: {DroppedReadSetMismatches}, " +
                         $"with hol(g) walked without its ordering term: {DroppedOrderingMismatches}; Corollary E3 without its read-sum test: " +
                         $"{DroppedReadSumMismatches}, with hol(u₀) walked without its ordering term: {DroppedE3OrderingMismatches}; Theorem B′ read " +
                         $"off the unreduced basis of the relations: {UnreducedBasisMismatches} (against the balance); Theorem E's ℤ/4 carrier " +
                         $"without its coset term: {DroppedCosetTermMismatches} (against (★)), moved by a quarter on one component: " +
                         $"{MovedCarrierNotQuadratic} no longer a Clifford (exactly where W has codimension two or more); Corollary E3′ without N₁*: " +
                         $"{DroppedE3PrimeReadSumMismatches}, with hol walked without its ordering term: {DroppedE3PrimeOrderingMismatches}; " +
                         $"Theorem H with every rate read as one: {RatesReadAwayMismatches}, with the rates in reverse site order: {RatesReversedMismatches} (both against λ_min of 𝓜); " +
                         $"Corollary H1's packing bound with every offset of the relation counted: {PackingOverlapMismatches}, its cut bound read where the chain breaks: {CutWhereChainBreaksMismatches} (both against the trellis); " +
                         $"Corollary H2 with k·λ_* in place of (k/2)·λ_*: {CycleDoubleCountMismatches} of {CycleRelations} relations (against λ_min of 𝓜); " +
                         $"Corollary H3's reading replaced by Theorem H read string by string: {HopControlMisses} of {HopReadings} readings (against λ_min of 𝓜); " +
                         $"Corollary A″ with no letter struck: {AgreeOnUnstruckMismatches}, with the next site's letters struck for one site: {AgreeOnWrongSiteMismatches} (on the chain; on the ring the two verdicts agree by translation), with the diagonal strings ignored where one is odd off Λ: {AgreeOnDiagonalIgnoredMismatches} (all against Theorem 2's count with the complement taken on the sites outside the agreeing ones); " +
                         $"the triangle's closed form handed the rate of a site on each bond in place of the site opposite it: {TriangleRotatedMismatches}, its odd term alone, without Corollary A″'s single sites: {TriangleOddOnlyMismatches} (both against the rates compressed onto the eigenspaces of the commutator with H)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "Theorem H's solver, measured",
                summary: $"λ_min(𝓜) against the code distance within {StrongDephasingWorstRatio:F3} of dim·eps·‖𝓜‖ at worst; the kernels counted " +
                         $"below 10⁻⁹·‖·‖, the largest zero at {KernelZeroMax:E1}·‖·‖, the smallest nonzero eigenvalue at {KernelGapMin:E1}·‖·‖, " +
                         $"the most negative at {MostNegative:E1}·‖·‖; the trellis against the enumeration within {TrellisWorstRatio:F3} of 2N·eps at worst; " +
                         $"the packing bound met exactly on {PackingEqualities} readings, the cut bound on {CutEqualities} (at random rates), " +
                         $"the ring broken with the chain pairing on {CutBroken} of the cut readings, λ at most {CutTightest:F3} of twice the cut",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "Corollary A″, ket and bra agreeing on some sites",
                summary: $"{AgreeOnNonzero} readings with operators X off the agreeing sites times a diagonal commuting with H (modes at −2·(σ − γ(Λ)) at every scale of H), " +
                         $"{AgreeOnZero} without; {AgreeOnWhereFarIsZero} of the first where Theorem A's far end is empty, a broken row with an exact mode at that slower rate; " +
                         $"{AgreeOnOutside} readings outside its hypotheses, where its balance declines",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the triangle X₀Y₁ + X₁Y₂ + X₂Y₀ at large J, and the sums over disjoint bonds",
                summary: $"the closed form against the rates compressed onto the eigenspaces of the commutator with H within {TriangleWorstRatio:F3} of 64·eps·2σ at worst; the closed form's odd term " +
                         $"the smaller on {TriangleOddWins} profiles, its single-site term on {TriangleCorollaryWins}; on the odd rings {OddTopTermsMoving} readings of top terms, " +
                         $"each failing to commute with H alone and moving its agreeing site one place with the letter exchanged; on the even rings {EvenTopTermsCommuting} readings of top terms, " +
                         $"each commuting with H alone",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "Theorem E's reach",
                summary: $"{CharacterOutside} readings outside its setting (an absent edge or a ratio that is no character); Corollary E2's " +
                         $"branches against Theorem 2: {string.Join(", ", E2Branches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}; " +
                         $"Corollary E3's: {string.Join(", ", E3Branches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}; " +
                         $"Corollary E3′'s at three shapes: {string.Join(", ", E3PrimeBranches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
