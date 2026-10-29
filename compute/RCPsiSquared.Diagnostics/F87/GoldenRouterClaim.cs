using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F116, the golden/metallic ceiling router (PROOF_CEILING_GOLDEN_ROUTER.md, reviewed bit-exact
/// 2026-06-22) as a standalone Tier1Derived Claim. This is the "C# witness first" inverse gap: the live
/// witness <see cref="Foundation.GoldenRouterWitness"/> (inspect --root router) and the reused helper
/// <see cref="KBodyPalindromeRouting"/> already carry the from-below truth; only the tiered Claim wrapper
/// was missing (open arc f116_golden_router_typed_claim).
///
/// <para><b>The result.</b> For the open-chain sliding-window Hamiltonian H = Σ_windows (XZX+XZY+YZX) (or
/// its X↔Y sibling YZY+XZY+YZX) at any N ≥ 3, under Z-dephasing at arbitrary site rates γ_l ≥ 0, there is
/// an invertible period-4 per-site product W = ⊗_l q_{l mod 4} with W L W⁻¹ = −L − 2σ, σ = Σ_l γ_l. The
/// per-site maps follow the [a, a, b, b] rhythm a = φX + Y, b = X − φY on the golden locus
/// α² − αβ − β² = 0, with h_l = q_l(Z) = (−1)^(l+1)·i·R(g_l) and q_l² = −(2+φ)·I (so cond(W) = 1). The
/// two-sided form W(ρ) = (2+φ)^(N/2)·P ρ Q reduces {W, [H,·]} = 0 to the chiral conditions
/// P H P⁻¹ = −H and Q H Q⁻¹ = −H, and the dephaser asks the two factors to differ in Z-parity at every
/// site: the golden chiral a^⊗N on one side and the Z pattern ⊗[I,I,Z,Z] on the other is the router. The
/// palindromizer is a per-site PRODUCT: the Z-middle ceiling cases are continuous-periodic-LOCAL, not
/// non-local. On a periodic ring the same router exists at every N divisible by four; for 4 ∤ N and c ≠ 0
/// no per-site class-swapping router exists.</para>
///
/// <para><b>Mechanism (the window lemma).</b> The three-site, template-summed anticommutator
/// {q_w⊗q_{w+1}⊗q_{w+2}, [XZX+XZY+YZX,·]₃} = 0 EXACTLY (in the ring ℤ[φ]+iℤ[φ]) at every window offset,
/// while no single template's anticommutator vanishes: the cancellation is CROSS-TEMPLATE inside one
/// window. Additivity over windows gives {W, A} = 0 at every N ≥ 3. (This is why the per-term certifier
/// could not see the router; the per-term <see cref="KBodyPalindromeRouting.Routes"/> declines these
/// cases.)</para>
///
/// <para><b>Exclusion (deductive).</b> {W,A}=0 on the identity column forces G = ⊗_l g_l to commute with
/// H, reducing per window to two bilinears K₁, K₂. From them: no uniform router (g = 0 only), no period-2,
/// and the golden locus α²−αβ−β²=0 is the only gate; period-3 is impossible for N ≥ 5 (an alternating
/// assignment cannot close an odd cycle; the N=4 period-3 solution is the now-explained 2026-06-07
/// small-N artifact); the discrete Klein candidates P1/P4/M2/M are all off-locus.</para>
///
/// <para><b>Metallic family.</b> On the soft line t₂ = t₃, the router transports verbatim with a = (r, 1),
/// b = (1, −r), r(c) = (c+√(c²+4))/2 the metallic mean of c = t₁/t₂ (golden c=1 → φ, silver c=2 → 1+√2,
/// bronze c=3, c=0 the 45° frame), q_l² = −(1+r²)·I, and W_c L_c W_c⁻¹ = −L_c − 2σ for every real c. The
/// window lemma is a polynomial identity in r (degree ≤ 5, verified exactly zero at 8 rational nodes), so
/// it holds for all real (indeed complex) c by derivation. The soft set of the weighted template on the open chain is
/// exact: every triple at N = 3, 4, and from N = 5 on exactly t₁ = 0 or t₂ = t₃ (F158's far end as two chain kernels over the
/// even and the odd sites, det K₃ = −t₁²(t₂−t₃)²(t₁²+(t₂+t₃)²)²; battery, PROOF §8 "the soft set"). On the
/// periodic ring it is exact as well: every triple at N = 4, t₁ = 0 or t₂ = t₃ at 4 | N ≥ 8, t₁ = 0 alone at
/// 4 ∤ N (the step w ↦ w + 2 closes into cycles; at N = 4 two windows share each end pair and cancel jointly;
/// battery, PROOF §8 "the soft set on rings").</para>
///
/// <para><b>Scope fences (the soft part, NOT Tier1Derived).</b> Rigidity (zero continuous physical moduli,
/// the router is unique modulo gauge) is PROVEN only at the golden and silver stations (Jacobian nullity =
/// gauge, N=5,6). The c=0 station's "8 physical moduli" is a finite-difference Jacobian count at N=5 ONLY
/// (spectral gap 3.5e9), NOT a closed-form/analytic result; 4 of 8 moduli are catalogued. Uniqueness is
/// scoped to INVERTIBLE W (the anticommutation equation alone carries abundant singular strata); within
/// the closed-form family at N=5 the router is essentially unique modulo per-site gauge, the four pattern
/// shifts, and an explicit order-32 sign group (128 gauge classes, every one a genuine palindromizer).</para>
///
/// <para>Tier1Derived: the existence + closed form + identity + exclusion + metallic-family carry no float
/// (exact ring / exact-Fraction interpolation); the end-to-end Lindbladian check is ~2e-15. Two typed
/// parents, both Tier1Derived (5 ≥ 5): <see cref="F1PalindromeIdentity"/> (the global palindrome
/// Π·L·Π⁻¹ = −L − 2Σγ·I that W realizes locally for the ceiling class) and
/// <see cref="WindowedConverseThresholdClaim"/> (the F87 two-reflection chiral spine 𝓕=F⊗F whose chiral
/// driving F H F = −H the two-sided form places on both sides, chirals of opposite Z-parity at every site). This router is the construction that
/// CLOSES <see cref="PalindromeSoftCertifierClaim"/>'s k=3 windowed soft-certifier ceiling at zero (the
/// 2 → 0 step); the certifier USES it as a helper, so the certifier is logically downstream, expressed
/// here as a see-cref, not a parent edge. Anchor: docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md +
/// ANALYTICAL_FORMULAS.md F116 + simulations/ceiling_golden_router.py + simulations/metallic_router_family.py
/// (both self-validating, exact ring arithmetic) + the live witness GoldenRouterWitness
/// (inspect --root router).</para></summary>
public sealed class GoldenRouterClaim : Claim
{
    /// <summary>The chain length used for the live window-summed routing re-check. The golden router is
    /// N-free; soft is established at N ≥ 4 for the Z-middle cases, so 4 suffices and is cheap.</summary>
    public const int CheckN = 4;

    /// <summary>Tolerance for the live-vs-closed-form metallic-mean MATCH (matches the witness/helper
    /// floor; the live ratio is a bisection root, the closed form the metallic quadratic).</summary>
    public const double MatchTolerance = 1e-9;

    /// <summary>The two typed parents, exposed for inspection.</summary>
    public F1PalindromeIdentity F1Palindrome { get; }
    public WindowedConverseThresholdClaim ChiralSpine { get; }
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>One self-check case tying the Claim to the from-below router machinery
    /// (<see cref="KBodyPalindromeRouting"/>), in the style of <see cref="WindowedConverseThresholdClaim"/>.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>One metallic-family reading: the weight c, the live router-derived frame ratio (root of
    /// ‖{W,S}‖_F = 0), the closed-form metallic mean, and whether they agree within
    /// <see cref="MatchTolerance"/>.</summary>
    public readonly record struct MetallicReading(double C, double LiveRatio, double ClosedForm)
    {
        public double Difference => Math.Abs(LiveRatio - ClosedForm);
        public bool Matches => Difference < MatchTolerance;
    }

    // --- live readings, recomputed once at construct time from the reused router helper ---

    /// <summary>The window-summed router certificate for the golden case XZX+XZY+YZX, or null. Expected
    /// <see cref="KBodyPalindromeRouting.GoldenDescription"/>.</summary>
    public string? GoldenRoutes { get; }

    /// <summary>True iff the per-term lens DECLINES the golden case (the documented cross-template coverage
    /// gap, NOT non-locality): the window-summed router certifies what the per-term lens cannot.</summary>
    public bool GoldenPerTermDeclines { get; }

    /// <summary>The window-summed router certificate for the X↔Y sibling YZY+XZY+YZX, or null. Expected
    /// <see cref="KBodyPalindromeRouting.GoldenMirrorDescription"/>.</summary>
    public string? SiblingRoutes { get; }

    /// <summary>The metallic family evaluated live at c ∈ {0, 1, 2, 3}: live root vs closed-form mean.</summary>
    public IReadOnlyList<MetallicReading> MetallicReadings { get; }

    public IReadOnlyList<BatteryCase> Battery => _battery.Value;
    private readonly Lazy<IReadOnlyList<BatteryCase>> _battery;
    public int PassCount => Battery.Count(c => c.Passes);
    public bool SelfCheckPasses => PassCount == Battery.Count;

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> Terms(params string[] labels) => labels.Select(T).ToList();

    public GoldenRouterClaim(F1PalindromeIdentity f1Palindrome, WindowedConverseThresholdClaim chiralSpine,
        PalindromeTwoEndCountClaim twoEndCount)
        : base(
            "F116 golden/metallic router: an invertible period-4 per-site product W = ⊗_l q_{l mod 4} " +
            "([a,a,b,b], a=φX+Y, b=X−φY on the golden locus α²−αβ−β²=0, q_l²=−(2+φ)I) palindromizes the " +
            "Z-middle sliding-window ceiling H=Σ(XZX+XZY+YZX) (and its X↔Y sibling) under Z-dephasing: " +
            "W L W⁻¹ = −L − 2σ at every N≥3 for arbitrary site rates γ_l≥0; the Z-middle ceiling cases are " +
            "continuous-periodic-LOCAL, not non-local. Generalises to the metallic family r(c)=(c+√(c²+4))/2 " +
            "as a polynomial identity in r (golden c=1, silver c=2, bronze c=3). Exclusion (no uniform/no " +
            "period-2/period-3-impossible N≥5/Klein off-locus) is deductive from the K₁,K₂ identity-column " +
            "functional",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md + docs/ANALYTICAL_FORMULAS.md F116 + " +
            "simulations/ceiling_golden_router.py + simulations/metallic_router_family.py")
    {
        F1Palindrome = f1Palindrome ?? throw new ArgumentNullException(nameof(f1Palindrome));
        ChiralSpine = chiralSpine ?? throw new ArgumentNullException(nameof(chiralSpine));
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));

        // Re-run the reused router machinery from below, once.
        var golden = Terms("XZX", "XZY", "YZX");
        var sibling = Terms("YZY", "XZY", "YZX");
        GoldenRoutes = KBodyPalindromeRouting.RoutesWindowSummed(golden, CheckN);
        GoldenPerTermDeclines = !KBodyPalindromeRouting.Routes(golden, CheckN);
        SiblingRoutes = KBodyPalindromeRouting.RoutesWindowSummed(sibling, CheckN);

        var readings = new List<MetallicReading>();
        foreach (double c in new[] { 0.0, 1.0, 2.0, 3.0 })
            readings.Add(new MetallicReading(c,
                KBodyPalindromeRouting.LiveMetallicRatio(c),
                KBodyPalindromeRouting.MetallicMean(c)));
        MetallicReadings = readings;

        LetterTableReadings = new[] { 0L, 1L, 2L, 3L }.Select(ReadLetterTable).ToList();
        ReflectedTableLeavesRemainder = ResidualEntries(new long[,] { { 1, 1 }, { -1, 0 } })
            .Any(e => !IsZero(DivideMonic(e, Locus(1)).Remainder));

        RhythmSolutions = Enumerable.Range(0, 16)
            .Select(bits => Enumerable.Range(0, 4).Select(l => (bits >> l) & 1).ToArray())
            .Where(z => Enumerable.Range(0, 4).All(w => z[w] + z[(w + 2) % 4] == 1))
            .Select(z => string.Concat(z.Select(b => b == 1 ? 'Z' : 'I')))
            .OrderBy(s => s, StringComparer.Ordinal)
            .ToList();

        RingRhythmCounts = Enumerable.Range(3, 7)
            .Select(n => Enumerable.Range(0, 1 << n)
                .Count(bits => Enumerable.Range(0, n).All(w => ((bits >> w) & 1) + ((bits >> ((w + 2) % n)) & 1) == 1)))
            .ToList();

        _splitRouterReadings = new Lazy<IReadOnlyList<SplitRouterReading>>(() => new[]
        {
            ReadSplitRouter("a^⊗N | S, open N=4", 4, ring: false, SplitFactor.Golden, SplitFactor.Rhythm),
            ReadSplitRouter("a^⊗N | S, open N=5", 5, ring: false, SplitFactor.Golden, SplitFactor.Rhythm),
            ReadSplitRouter("a^⊗N | S, periodic ring N=4", 4, ring: true, SplitFactor.Golden, SplitFactor.Rhythm),
            ReadSplitRouter("S | a^⊗N, open N=4", 4, ring: false, SplitFactor.Rhythm, SplitFactor.Golden),
            ReadSplitRouter("control a^⊗N | a^⊗N, open N=4", 4, ring: false, SplitFactor.Golden, SplitFactor.Golden),
            ReadSplitRouter("control a^⊗N·S | I, open N=4", 4, ring: false, SplitFactor.Both, SplitFactor.Identity),
            ReadSplitRouter("control a^⊗N | S, periodic ring N=5", 5, ring: true, SplitFactor.Golden, SplitFactor.Rhythm),
        });

        _chainDeterminantIdentityHolds = new Lazy<bool>(() =>
        {
            var grid = Enumerable.Range(-4, 9).Select(v => (long)v).ToArray();
            return grid.All(a => grid.All(b => grid.All(c =>
                Det(ChainConstraint(a, b, c)) == ChainDeterminantClosedForm(a, b, c))));
        });

        _reductionMismatches = new Lazy<int>(() =>
        {
            var grid = Enumerable.Range(-2, 5).Select(v => (long)v).ToArray();
            return (from n in new[] { 5, 6 } from a in grid from b in grid from c in grid
                    where FarKernelModP(n, a, b, c) != ChainKernelModP((n + 1) / 2, a, b, c) * ChainKernelModP(n / 2, a, b, c)
                    select n).Count();
        });

        (long, long) xPlusY = (1, 1), xMinus2Y = (1, -2), x = (1, 0);
        _softSetReadings = new Lazy<IReadOnlyList<SoftSetReading>>(() => new[]
        {
            new SoftSetReading("(3,4,5), N=3: (X+Y)⊗X⊗(X−2Y) commutes with H, invertible", "carries",
                ProductElementCarries(3, 4, 5, new[] { xPlusY, x, xMinus2Y }) ? "carries" : "fails"),
            new SoftSetReading("(3,4,5), N=4: (X+Y)⊗(X+Y)⊗(X−2Y)⊗(X−2Y) commutes with H, invertible", "carries",
                ProductElementCarries(3, 4, 5, new[] { xPlusY, xPlusY, xMinus2Y, xMinus2Y }) ? "carries" : "fails"),
            new SoftSetReading("control (3,4,5), N=3: (X+Y)⊗X⊗(X+Y), a root without its partner", "fails",
                ProductElementCarries(3, 4, 5, new[] { xPlusY, x, xPlusY }) ? "carries" : "fails"),
            new SoftSetReading("(3,4,5), N=5: far kernel mod p", "0", FarKernelModP(5, 3, 4, 5).ToString(Inv)),
            new SoftSetReading("(3,4,5), N=6: far kernel mod p", "0", FarKernelModP(6, 3, 4, 5).ToString(Inv)),
            new SoftSetReading("control on the line (2,1,1), N=5: far kernel mod p", "4", FarKernelModP(5, 2, 1, 1).ToString(Inv)),
            new SoftSetReading("control on the plane (0,2,7), N=5: far kernel mod p", "4", FarKernelModP(5, 0, 2, 7).ToString(Inv)),
        });

        _ringFourElementHolds = new Lazy<bool>(() =>
        {
            var grid = Enumerable.Range(-2, 5).Select(v => (long)v).ToArray();
            bool commutes = grid.All(a => grid.All(b => grid.All(c =>
                CommutesWithTemplate(4, a, b, c, ring: true, RingFourElement(a, b + c)))));
            bool determinant = Enumerable.Range(-16, 33).All(t =>
            {
                BigInteger closed = BigInteger.Pow(2, 16) * BigInteger.Pow(t, 8) * BigInteger.Pow((BigInteger)t * t + 1, 6);
                return GaussianDet(MatrixOf(RingFourElement(t, 1))) == new GaussianInteger(closed, 0);
            });
            return commutes && determinant;
        });

        _ringFourLineElementHolds = new Lazy<bool>(() =>
        {
            var grid = Enumerable.Range(-2, 5).Select(v => (long)v).ToArray();
            var g = RingFourLineElement();
            return grid.All(a => grid.All(b => CommutesWithTemplate(4, a, b, -b, ring: true, g)))
                && GaussianDet(MatrixOf(g)) == new GaussianInteger(BigInteger.Pow(2, 24), 0);
        });

        string Far(long t1, long t2, long t3, params int[] ns) =>
            string.Join(", ", ns.Select(n => FarKernelModP(n, t1, t2, t3, ring: true).ToString(Inv)));
        var xs = Enumerable.Repeat(x, 5).ToArray();
        (long, long) a2 = (2, 1), b2 = (1, -2), y = (0, 1);   // the metallic directions at c = 3/2 (r = 2)
        _ringSoftSetReadings = new Lazy<IReadOnlyList<SoftSetReading>>(() => new[]
        {
            new SoftSetReading("(1,2,3), off the line and the plane, ring N = 3, 5, 6, 7: far kernel mod p", "0, 0, 0, 0",
                Far(1, 2, 3, 3, 5, 6, 7)),
            new SoftSetReading("(2,1,1), on the line, ring N = 3, 5, 6, 7 (4 ∤ N): far kernel mod p", "0, 0, 0, 0",
                Far(2, 1, 1, 3, 5, 6, 7)),
            new SoftSetReading("control (2,1,1), on the line, ring N = 8 (4 | N): far kernel mod p", "4",
                Far(2, 1, 1, 8)),
            new SoftSetReading("(1,2,3), ring N = 4: far kernel mod p", "6", Far(1, 2, 3, 4)),
            new SoftSetReading("(3,2,2), on the line, ring N = 8: ⊗[a,a,b,b] with a = 2X+Y, b = X−2Y commutes with H, invertible",
                "carries", ProductElementCarries(3, 2, 2, new[] { a2, a2, b2, b2, a2, a2, b2, b2 }, ring: true) ? "carries" : "fails"),
            new SoftSetReading("(1,0,0), the line's end t₂ = t₃ = 0, ring N = 8: XXYYXXYY commutes with H, invertible",
                "carries", ProductElementCarries(1, 0, 0, new[] { x, x, y, y, x, x, y, y }, ring: true) ? "carries" : "fails"),
            new SoftSetReading("control (3,2,2), ring N = 6 (4 ∤ N): ⊗[a,a,b,b,a,a], the rhythm cannot close",
                "fails", ProductElementCarries(3, 2, 2, new[] { a2, a2, b2, b2, a2, a2 }, ring: true) ? "carries" : "fails"),
            new SoftSetReading("(0,2,7), on the plane, ring N = 5: X^⊗5 commutes with H, invertible", "carries",
                ProductElementCarries(0, 2, 7, xs, ring: true) ? "carries" : "fails"),
            new SoftSetReading("control (1,2,3), ring N = 5: X^⊗5", "fails",
                ProductElementCarries(1, 2, 3, xs, ring: true) ? "carries" : "fails"),
        });

        _battery = new Lazy<IReadOnlyList<BatteryCase>>(BuildBattery);
    }

    /// <summary>One exact reading of a two-sided map ρ ↦ P ρ Q against W∘L = (−L − 2σ)∘W (PROOF §2, "Why
    /// two-sided?"): the per-site parity rule needs P and Q to differ in Z-parity at every site, and the
    /// Hamiltonian leg needs both to be chiral; the golden chiral a^⊗N (lit at every site) and the Z pattern
    /// S = ⊗[I,I,Z,Z] (dead at every site) are such a pair.</summary>
    public readonly record struct SplitRouterReading(string Name, bool Routes);

    /// <summary>The split router read exactly at the rational metallic point c = 3/2, where r = 2 and
    /// a = 2X + Y (H = Σ 3·XZX + 2·XZY + 2·YZX, integer rates γ_l = l + 1), in Gaussian integers: the routing
    /// pairs, with the sides swapped, on the open chain and on the ring N = 4, and the three controls that must
    /// fail (a^⊗N on both sides, the identity column on one side, the ring at N = 5 where S is not chiral).</summary>
    public IReadOnlyList<SplitRouterReading> SplitRouterReadings => _splitRouterReadings.Value;
    private readonly Lazy<IReadOnlyList<SplitRouterReading>> _splitRouterReadings;

    /// <summary>On a periodic ring of N = 3..9 sites, the number of Z patterns with exactly one Z on every window's end
    /// pair (w, w+2 mod N). Since b ∝ a·Z, this is also the number of ways §4's alternation a ↔ b can close
    /// along the ring's parity cycles; it is nonzero exactly when 4 | N.</summary>
    public IReadOnlyList<int> RingRhythmCounts { get; }

    private enum SplitFactor { Identity, Golden, Rhythm, Both }

    private static SplitRouterReading ReadSplitRouter(string name, int n, bool ring, SplitFactor left, SplitFactor right)
    {
        int d = 1 << n;
        string rhythm = string.Concat(Enumerable.Range(0, n).Select(l => l % 4 < 2 ? 'I' : 'Z'));
        GMat Factor(SplitFactor f) => f switch
        {
            SplitFactor.Identity => GMat.Product(new string('I', n)),
            SplitFactor.Golden => GMat.Product(new string('a', n)),
            SplitFactor.Rhythm => GMat.Product(rhythm),
            _ => GMat.Mul(GMat.Product(new string('a', n)), GMat.Product(rhythm)),
        };
        var p = Factor(left);
        var q = Factor(right);

        // H = Σ_windows 3·XZX + 2·XZY + 2·YZX: c = 3/2 on the soft line, r = 2.
        var h = GMat.Zero(d);
        int windows = ring ? n : n - 2;
        foreach (var (template, weight) in new[] { ("XZX", 3L), ("XZY", 2L), ("YZX", 2L) })
            for (int w = 0; w < windows; w++)
            {
                var letters = Enumerable.Repeat('I', n).ToArray();
                for (int k = 0; k < 3; k++) letters[(w + k) % n] = template[k];
                h = GMat.Combine(h, 1, GMat.Product(new string(letters)), weight);
            }
        long[] gamma = Enumerable.Range(1, n).Select(l => (long)l).ToArray();
        long twoSigma = 2 * gamma.Sum();

        // L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ); Z_l is diagonal, site 0 the most significant bit.
        GMat L(GMat rho)
        {
            var comm = GMat.Combine(GMat.Mul(h, rho), 1, GMat.Mul(rho, h), -1);
            var outp = GMat.Zero(d);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                {
                    long rate = 0;
                    for (int l = 0; l < n; l++)
                        if ((((i ^ j) >> (n - 1 - l)) & 1) == 1) rate += 2 * gamma[l];
                    outp.Re[i, j] = comm.Im[i, j] - rate * rho.Re[i, j];                       // −i·comm, then −2γ_l
                    outp.Im[i, j] = -comm.Re[i, j] - rate * rho.Im[i, j];                      // per site where i, j differ
                }
            return outp;
        }

        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                var e = GMat.Zero(d);
                e.Re[i, j] = 1;
                var wle = GMat.Mul(GMat.Mul(p, L(e)), q);
                var we = GMat.Mul(GMat.Mul(p, e), q);
                var rhs = GMat.Combine(L(we), -1, we, -twoSigma);
                if (!GMat.Equal(wle, rhs)) return new SplitRouterReading(name, false);
            }
        return new SplitRouterReading(name, true);
    }

    /// <summary>The period-4 Z patterns (I/Z per site mod 4) that put exactly one Z on every window's end pair
    /// (w, w+2): the chiral S of the rhythm and its shifts. Periods 1 and 2 are included in the enumeration
    /// (they are period-4 patterns too) and fail, since their two end sites carry the same letter.</summary>
    public IReadOnlyList<string> RhythmSolutions { get; }

    /// <summary>One reading of the window's letter table at integer weight c (PROOF §2): the end letters of
    /// c·XZX + XZY + YZX carry M_c = [[c,1],[1,0]] over {X,Y}×{X,Y}. For the direction v = (r, 1), i.e.
    /// â = rX + Y, the eigen condition det[v, M_c v] is the locus polynomial r² − cr − 1; and conjugation by
    /// the uniform product â^⊗N flips the middle Z (−(r²+1)), reflects both end letters by
    /// R̃ = [[r²−1, 2r], [2r, 1−r²]] and scales every other site by r²+1, so â^⊗N H â^⊗N = −(r²+1)^N·H iff
    /// R̃M_cR̃ = (r²+1)²M_c. Computed in integer
    /// polynomials in r, that residual is exactly (r² − cr − 1)·K(r) with K = [[4r, 2(1−r²)], [2(1−r²), −4r]],
    /// whose entries 4r and 2(1−r²) never vanish together: the uniform product is chiral exactly on the
    /// metallic locus, where â is an eigenvector of M_c. Both the residual and (r² − cr − 1)·K are affine in
    /// c, so agreement at two values of c settles every c; the battery reads four.</summary>
    public readonly record struct LetterTableReading(long C, bool EigenEquationIsLocus, bool ChiralResidualIsLocusTimesK);

    /// <summary>The letter-table readings at c ∈ {0, 1, 2, 3} (exact, integer polynomial arithmetic).</summary>
    public IReadOnlyList<LetterTableReading> LetterTableReadings { get; }

    /// <summary>True iff the reflected table [[1,1],[−1,0]] (t₃ = −t₂, hard from N = 5) leaves a nonzero remainder modulo
    /// r² − r − 1: the control on which the letter-table check must fail.</summary>
    public bool ReflectedTableLeavesRemainder { get; }

    // Integer polynomials in r: index i holds the coefficient of r^i.
    private static long[] P(params long[] coefficients) => coefficients;
    private static long[] Locus(long c) => P(-1, -c, 1);

    private static long[] Add(long[] a, long[] b)
    {
        var s = new long[Math.Max(a.Length, b.Length)];
        for (int i = 0; i < a.Length; i++) s[i] += a[i];
        for (int i = 0; i < b.Length; i++) s[i] += b[i];
        return s;
    }

    private static long[] Scale(long[] a, long k) => a.Select(x => k * x).ToArray();

    private static long[] Mul(long[] a, long[] b)
    {
        var p = new long[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++)
                p[i + j] += a[i] * b[j];
        return p;
    }

    private static bool IsZero(long[] a) => a.All(x => x == 0);

    private static bool Equal(long[] a, long[] b) => IsZero(Add(a, Scale(b, -1)));

    /// <summary>Division by a monic polynomial over ℤ (exact, no rounding).</summary>
    private static (long[] Quotient, long[] Remainder) DivideMonic(long[] a, long[] monic)
    {
        int d = monic.Length - 1;
        var rem = (long[])a.Clone();
        var q = new long[Math.Max(1, a.Length - d)];
        for (int i = a.Length - 1; i >= d; i--)
        {
            long lead = rem[i];
            if (lead == 0) continue;
            q[i - d] = lead;
            for (int j = 0; j <= d; j++) rem[i - d + j] -= lead * monic[j];
        }
        return (q, rem);
    }

    /// <summary>R̃MR̃ − (r²+1)²M for an integer 2×2 letter table M, entry-wise as polynomials in r
    /// (row-major: entries 0..3 = (0,0), (0,1), (1,0), (1,1)).</summary>
    private static long[][] ResidualEntries(long[,] m)
    {
        var rt = new[,] { { P(-1, 0, 1), P(0, 2) }, { P(0, 2), P(1, 0, -1) } };   // R̃ = [[r²−1, 2r], [2r, 1−r²]]
        var n4 = Mul(P(1, 0, 1), P(1, 0, 1));                                      // (r²+1)²
        var entries = new long[4][];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                var s = P(0);
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++)
                        s = Add(s, Scale(Mul(rt[i, k], rt[l, j]), m[k, l]));
                entries[2 * i + j] = Add(s, Scale(n4, -m[i, j]));
            }
        return entries;
    }

    private static LetterTableReading ReadLetterTable(long c)
    {
        var locus = Locus(c);
        var table = new long[,] { { c, 1 }, { 1, 0 } };
        // det[v, M_c v] with v = (r, 1), M_c v read off the same table.
        var v = new[] { P(0, 1), P(1) };
        var mv = new[]
        {
            Add(Scale(v[0], table[0, 0]), Scale(v[1], table[0, 1])),
            Add(Scale(v[0], table[1, 0]), Scale(v[1], table[1, 1])),
        };
        var eigen = Add(Mul(v[0], mv[1]), Scale(Mul(mv[0], v[1]), -1));

        var k = new[] { P(0, 4), P(2, 0, -2), P(2, 0, -2), P(0, -4) };
        var residual = ResidualEntries(table);
        bool locusTimesK = residual.Select((e, i) => (e, i)).All(x =>
        {
            var (q, rem) = DivideMonic(x.e, locus);
            return IsZero(rem) && Equal(q, k[x.i]);
        });
        return new LetterTableReading(c, Equal(eigen, locus), locusTimesK);
    }

    // --- the soft set of the weighted template t₁·XZX + t₂·XZY + t₃·YZX (PROOF §8), exact and integer ---

    /// <summary>The window's two-site rule on the far end (PROOF §8, "the soft set"): a lit string carries X or Y
    /// on every site, the middle Z of a window term always turns, so the term's row depends only on whether the
    /// end pair (f_w, f_{w+2}) comes out II or ZZ. Over (XX, XY, YX, YY) that is one 2×4 rule
    /// P = [[t₁, t₂, t₃, 0], [0, t₃, t₂, −t₁]], whatever the middle and the rest carry, so the far kernel is the
    /// product of two chain kernels, one over the even sites and one over the odd ones. K₃ is a three-site
    /// chain's: P on sites (0, 1) and on (1, 2), an 8×8 integer matrix.</summary>
    private static long[,] ChainConstraint(long t1, long t2, long t3)
    {
        var p = new long[,] { { t1, t2, t3, 0 }, { 0, t3, t2, -t1 } };
        var a = new long[8, 8];
        for (int r = 0; r < 2; r++)
            for (int c = 0; c < 4; c++)
                for (int k = 0; k < 2; k++)
                {
                    a[2 * r + k, 2 * c + k] = p[r, c];       // P on (s0, s1), identity on s2: index 4·s0 + 2·s1 + s2
                    a[4 + 2 * k + r, 4 * k + c] = p[r, c];   // identity on s0, P on (s1, s2)
                }
        return a;
    }

    /// <summary>Exact determinant of an integer matrix (Bareiss, fraction-free).</summary>
    private static BigInteger Det(long[,] m)
    {
        int n = m.GetLength(0);
        var a = new BigInteger[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i, j] = m[i, j];
        BigInteger prev = 1, sign = 1;
        for (int k = 0; k < n - 1; k++)
        {
            if (a[k, k].IsZero)
            {
                int s = k + 1;
                while (s < n && a[s, k].IsZero) s++;
                if (s == n) return 0;
                for (int j = 0; j < n; j++) (a[k, j], a[s, j]) = (a[s, j], a[k, j]);
                sign = -sign;
            }
            for (int i = k + 1; i < n; i++)
                for (int j = k + 1; j < n; j++)
                    a[i, j] = (a[i, j] * a[k, k] - a[i, k] * a[k, j]) / prev;
            prev = a[k, k];
        }
        return sign * a[n - 1, n - 1];
    }

    /// <summary>The closed form of det K₃: −t₁²(t₂ − t₃)²(t₁² + (t₂ + t₃)²)², whose real zeros are exactly
    /// t₁ = 0 or t₂ = t₃.</summary>
    private static BigInteger ChainDeterminantClosedForm(long t1, long t2, long t3)
    {
        BigInteger a = t1, d = t2 - t3, s = (BigInteger)t1 * t1 + (BigInteger)(t2 + t3) * (t2 + t3);
        return -(a * a) * (d * d) * (s * s);
    }

    /// <summary>True iff det K₃ equals its closed form on the whole grid {−4..4}³. Both sides are polynomials of
    /// degree at most 8 in each variable, so agreement on nine values per variable is the identity.</summary>
    public bool ChainDeterminantIdentityHolds => _chainDeterminantIdentityHolds.Value;
    private readonly Lazy<bool> _chainDeterminantIdentityHolds;

    /// <summary>The reduction to two chain kernels read from below: the number of triples on the grid {−2..2}³, at
    /// N = 5 and 6, where the far kernel built from the Pauli strings of H differs from dim K_⌈N/2⌉ · dim K_⌊N/2⌋
    /// built from P alone (both modulo the same prime). Expected 0.</summary>
    public int ReductionMismatches => _reductionMismatches.Value;
    private readonly Lazy<int> _reductionMismatches;

    // Pauli letters 0 = I, 1 = X, 2 = Y, 3 = Z; a·b = i^Phase[a,b]·Letter[a,b].
    private static readonly int[,] LetterOf = { { 0, 1, 2, 3 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, { 3, 2, 1, 0 } };
    private static readonly int[,] PhaseOf = { { 0, 0, 0, 0 }, { 0, 0, 1, 3 }, { 0, 3, 0, 1 }, { 0, 1, 3, 0 } };

    private static (int Phase, int[] Letters) Product(int[] s, int[] t)
    {
        int e = 0;
        var u = new int[s.Length];
        for (int l = 0; l < s.Length; l++)
        {
            e += PhaseOf[s[l], t[l]];
            u[l] = LetterOf[s[l], t[l]];
        }
        return (e % 4, u);
    }

    private static List<(int[] Letters, long Coefficient)> Template(int n, long t1, long t2, long t3, bool ring = false)
    {
        var h = new List<(int[], long)>();
        int windows = ring ? n : n - 2;
        foreach (var (a, b, c) in new[] { (1, 1, t1), (1, 2, t2), (2, 1, t3) })
            for (int w = 0; w < windows && c != 0; w++)
            {
                var s = new int[n];
                s[w] = a; s[(w + 1) % n] = 3; s[(w + 2) % n] = b;
                h.Add((s, c));
            }
        return h;
    }

    private const long SoftPrime = 998244353;   // ≡ 1 mod 4, so i = √−1 lives in GF(p)

    private static long PowMod(long b, long e, long m)
    {
        long r = 1;
        b %= m;
        for (; e > 0; e >>= 1, b = b * b % m)
            if ((e & 1) == 1) r = r * b % m;
        return r;
    }

    private static long SqrtMinusOne(long p)
    {
        for (long g = 2; ; g++)
        {
            long x = PowMod(g, (p - 1) / 4, p);
            if (x * x % p == p - 1) return x;
        }
    }

    /// <summary>dim of the far kernel modulo p: the lit strings ({X, Y} on every site) whose combinations commute
    /// with H, read as the nullity of [H, ·] on their span over GF(p) with i = √−1. Reduction can only lower a
    /// rank, so this is an upper bound on the rational far kernel: a reading of 0 proves it empty.</summary>
    private static int FarKernelModP(int n, long t1, long t2, long t3, bool ring = false)
    {
        long p = SoftPrime, iu = SqrtMinusOne(p);
        long[] ipow = { 1, iu, p - 1, p - iu };
        var h = Template(n, t1, t2, t3, ring);
        int cols = 1 << n;
        var rows = new Dictionary<string, int>();
        var entries = new Dictionary<(int, int), long>();
        for (int col = 0; col < cols; col++)
        {
            var f = Enumerable.Range(0, n).Select(l => ((col >> l) & 1) == 0 ? 1 : 2).ToArray();
            foreach (var (term, c) in h)
            {
                int differ = Enumerable.Range(0, n).Count(l => term[l] != 0 && term[l] != f[l]);
                if (differ % 2 == 0) continue;                    // commutes: no contribution
                var (e, u) = Product(term, f);                    // [T, F] = 2·T·F when they anticommute
                string key = string.Concat(u);
                if (!rows.TryGetValue(key, out int r)) rows[key] = r = rows.Count;
                long v = ((2 * c) % p + p) % p * ipow[e] % p;
                entries[(r, col)] = ((entries.TryGetValue((r, col), out long old) ? old : 0) + v) % p;
            }
        }
        var m = new long[rows.Count, cols];
        foreach (var ((r, col), v) in entries) m[r, col] = v;
        return NullityModP(m);
    }

    /// <summary>dim K_m modulo p: the vectors on a chain of m sites (letter X or Y per site, bit l of the index)
    /// that P sends to zero on every neighbouring pair, built from P alone.</summary>
    private static int ChainKernelModP(int sites, long t1, long t2, long t3)
    {
        if (sites < 2) return 1 << sites;
        long p = SoftPrime;
        long[,] rule = { { t1, t2, t3, 0 }, { 0, t3, t2, -t1 } };
        int cols = 1 << sites, restCount = 1 << (sites - 2);
        var m = new long[(sites - 1) * 2 * restCount, cols];
        for (int j = 0; j + 1 < sites; j++)
            for (int col = 0; col < cols; col++)
            {
                int pair = 2 * ((col >> j) & 1) + ((col >> (j + 1)) & 1);
                int rest = (col & ((1 << j) - 1)) | ((col >> (j + 2)) << j);
                for (int r = 0; r < 2; r++)
                    m[(j * 2 + r) * restCount + rest, col] = ((rule[r, pair] % p) + p) % p;
            }
        return NullityModP(m);
    }

    private static int NullityModP(long[,] m)
    {
        long p = SoftPrime;
        int rowsN = m.GetLength(0), cols = m.GetLength(1), rank = 0;
        for (int col = 0; col < cols && rank < rowsN; col++)
        {
            int piv = rank;
            while (piv < rowsN && m[piv, col] == 0) piv++;
            if (piv == rowsN) continue;
            for (int j = 0; j < cols; j++) (m[rank, j], m[piv, j]) = (m[piv, j], m[rank, j]);
            long inv = PowMod(m[rank, col], p - 2, p);
            for (int i = 0; i < rowsN; i++)
            {
                if (i == rank || m[i, col] == 0) continue;
                long fct = m[i, col] * inv % p;
                for (int j = 0; j < cols; j++)
                    m[i, j] = ((m[i, j] - fct * m[rank, j]) % p + p) % p;
            }
            rank++;
        }
        return cols - rank;
    }

    /// <summary>True iff G = ⊗ g_l (g_l integer combinations of X and Y) commutes with the template exactly and
    /// squares to a nonzero multiple of the identity: an invertible element of F158's far kernel (G anticommutes
    /// with every Z_l, its letters lying in span{X, Y}), so the spectrum pairs at every positive rate.</summary>
    private static bool ProductElementCarries(long t1, long t2, long t3, (long X, long Y)[] g, bool ring = false)
    {
        int n = g.Length;
        var strings = new List<(int[] Letters, long Coefficient)> { (new int[0], 1) };
        foreach (var (x, y) in g)
            strings = strings.SelectMany(s => new[] { (s.Letters.Append(1).ToArray(), s.Coefficient * x),
                                                      (s.Letters.Append(2).ToArray(), s.Coefficient * y) })
                             .Where(s => s.Item2 != 0).ToList();
        // G² = Π (x² + y²)·I for a product of real XY-plane letters; read it off the factors
        bool invertible = g.All(f => f.X * f.X + f.Y * f.Y != 0);
        return CommutesWithTemplate(n, t1, t2, t3, ring, strings) && invertible;
    }

    /// <summary>True iff [H, G] = 0 exactly, G = Σ c_s·s a combination of Pauli strings with integer coefficients,
    /// H the weighted template on the open chain or on the ring; accumulated in Gaussian integers, keyed by the
    /// output string.</summary>
    private static bool CommutesWithTemplate(int n, long t1, long t2, long t3, bool ring,
                                             IReadOnlyList<(int[] Letters, long Coefficient)> strings)
    {
        var acc = new Dictionary<string, (long Re, long Im)>();
        void Add(string key, long c, int e)
        {
            var (re, im) = acc.TryGetValue(key, out var v) ? v : (0, 0);
            acc[key] = (e % 4) switch { 0 => (re + c, im), 1 => (re, im + c), 2 => (re - c, im), _ => (re, im - c) };
        }
        foreach (var (term, c) in Template(n, t1, t2, t3, ring))
            foreach (var (s, a) in strings)
            {
                var (e1, u) = Product(term, s);
                var (e2, _) = Product(s, term);
                string key = string.Concat(u);
                Add(key, c * a, e1);
                Add(key, -c * a, e2);
            }
        return acc.Values.All(v => v.Re == 0 && v.Im == 0);
    }

    // --- the soft set on the periodic ring (PROOF §8, "the soft set on rings") ---

    private static List<(int[] Letters, long Coefficient)> Strings(params (string Label, long Coefficient)[] terms) =>
        terms.Where(t => t.Coefficient != 0)
             .Select(t => (t.Label.Select(ch => "IXYZ".IndexOf(ch)).ToArray(), t.Coefficient)).ToList();

    /// <summary>The ring N = 4 element off t₁ = 0 and t₂ + t₃ = 0, with s = t₂ + t₃: windows w and
    /// w + 2 share the end pair {w, w + 2} with opposite middles, and off the line t₂ = t₃ F158's far end depends on
    /// t₁ and s only (on the line it grows).</summary>
    private static List<(int[] Letters, long Coefficient)> RingFourElement(long t1, long s) => Strings(
        ("XXYY", t1 * t1), ("XYYX", t1 * t1), ("XYXY", -s * s), ("YXYX", -s * s),
        ("YYYX", t1 * s), ("XYYY", t1 * s), ("XXYX", -t1 * s), ("XYXX", -t1 * s));

    /// <summary>The ring N = 4 element on the line t₃ = −t₂, the same for every t₁ and t₂.</summary>
    private static List<(int[] Letters, long Coefficient)> RingFourLineElement() => Strings(
        ("XXYY", 2), ("XYYX", 1), ("YXXY", 1), ("XXXX", 1), ("XYXY", 1), ("YXYX", 1), ("YYYY", 1));

    /// <summary>The 2^N × 2^N matrix of Σ c_s·s in Gaussian integers.</summary>
    private static GaussianInteger[,] MatrixOf(IReadOnlyList<(int[] Letters, long Coefficient)> strings)
    {
        int n = strings[0].Letters.Length, d = 1 << n;
        var m = new GaussianInteger[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = GaussianInteger.Zero;
        foreach (var (letters, c) in strings)
        {
            var g = GMat.Product(string.Concat(letters.Select(l => "IXYZ"[l])));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    m[i, j] += new GaussianInteger(c * g.Re[i, j], c * g.Im[i, j]);
        }
        return m;
    }

    /// <summary>Exact determinant over ℤ[i] (Bareiss, fraction-free; every division is exact and is checked).</summary>
    private static GaussianInteger GaussianDet(GaussianInteger[,] m)
    {
        int n = m.GetLength(0);
        var a = (GaussianInteger[,])m.Clone();
        GaussianInteger prev = GaussianInteger.One, sign = GaussianInteger.One;
        static bool IsZero(GaussianInteger z) => z.Re.IsZero && z.Im.IsZero;
        static GaussianInteger DivideExact(GaussianInteger x, GaussianInteger y)
        {
            BigInteger norm = y.Re * y.Re + y.Im * y.Im;
            var num = x * new GaussianInteger(y.Re, -y.Im);
            if (!(num.Re % norm).IsZero || !(num.Im % norm).IsZero)
                throw new InvalidOperationException("Bareiss division not exact in ℤ[i]");
            return new GaussianInteger(num.Re / norm, num.Im / norm);
        }
        for (int k = 0; k < n - 1; k++)
        {
            if (IsZero(a[k, k]))
            {
                int s = k + 1;
                while (s < n && IsZero(a[s, k])) s++;
                if (s == n) return GaussianInteger.Zero;
                for (int j = 0; j < n; j++) (a[k, j], a[s, j]) = (a[s, j], a[k, j]);
                sign = -sign;
            }
            for (int i = k + 1; i < n; i++)
                for (int j = k + 1; j < n; j++)
                    a[i, j] = DivideExact(a[i, j] * a[k, k] - a[i, k] * a[k, j], prev);
            prev = a[k, k];
        }
        return sign * a[n - 1, n - 1];
    }

    /// <summary>The ring N = 4 element read exactly. (i) It commutes with H at every triple of {−2..2}³: the entries
    /// of [H, G] are polynomials of degree at most three in each of t₁, t₂, t₃, so five values per variable make it
    /// an identity. (ii) det G = 2¹⁶·t₁⁸·s¹²·(t₁² + s²)⁶: G depends on t₁ and s only, each entry homogeneous of
    /// degree two, so det G is a binary form of degree 32, fixed by its values at s = 1 and the 33 integers
    /// t₁ = −16..16, Bareiss over ℤ[i]. Its real zeros are t₁ = 0 and s = 0. Its letters lie in span{X, Y}, so it
    /// anticommutes with every Z_l and is an invertible element of F158's far kernel.</summary>
    public bool RingFourElementHolds => _ringFourElementHolds.Value;
    private readonly Lazy<bool> _ringFourElementHolds;

    /// <summary>The ring N = 4 element on the line t₃ = −t₂: it commutes with H at every (t₁, t₂) of {−2..2}²
    /// (degree one in each), and its determinant is the constant 2²⁴.</summary>
    public bool RingFourLineElementHolds => _ringFourLineElementHolds.Value;
    private readonly Lazy<bool> _ringFourLineElementHolds;

    /// <summary>The soft set on the periodic ring read at integer triples: F158's far kernel modulo p (an upper
    /// bound, so a reading of 0 proves it empty against a near end holding the identity), and invertible product
    /// elements: the [a, a, b, b] rhythm and XXYY on the line at 4 | N, X^⊗N on the plane.</summary>
    public IReadOnlyList<SoftSetReading> RingSoftSetReadings => _ringSoftSetReadings.Value;
    private readonly Lazy<IReadOnlyList<SoftSetReading>> _ringSoftSetReadings;

    /// <summary>One reading of the soft set at an integer triple.</summary>
    public readonly record struct SoftSetReading(string Name, string Expected, string Actual);

    /// <summary>The soft set read at the Pythagorean triple (3, 4, 5), off the line and off the plane: there
    /// t₃² = t₁² + t₂², the root equation Q(r) = t₁t₃r² + (t₃² − t₁² − t₂²)r − t₁t₃ becomes t₁t₃(r² − 1), its roots
    /// r = ±1 are rational, and the element is in integers. Soft at N = 3, 4 (each parity chain at most two
    /// sites, so a root suffices), hard at N = 5, 6 (a chain of three sites, det K₃ ≠ 0, far kernel empty against
    /// a near end holding the identity). With the line and the plane as controls, where the far kernel reads four
    /// at N = 5.</summary>
    public IReadOnlyList<SoftSetReading> SoftSetReadings => _softSetReadings.Value;
    private readonly Lazy<IReadOnlyList<SoftSetReading>> _softSetReadings;

    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        cases.Add(new BatteryCase(
            Name: "golden case routes (window-summed)",
            Detail: "RoutesWindowSummed(XZX+XZY+YZX, N=4)",
            Expected: KBodyPalindromeRouting.GoldenDescription,
            Actual: GoldenRoutes ?? "null"));

        cases.Add(new BatteryCase(
            Name: "golden case declines per-term (cross-template gap)",
            Detail: "!Routes(XZX+XZY+YZX, N=4) (per-term lens cannot see cross-template cancellation)",
            Expected: bool.TrueString,
            Actual: GoldenPerTermDeclines.ToString()));

        cases.Add(new BatteryCase(
            Name: "X↔Y sibling routes (window-summed mirror)",
            Detail: "RoutesWindowSummed(YZY+XZY+YZX, N=4)",
            Expected: KBodyPalindromeRouting.GoldenMirrorDescription,
            Actual: SiblingRoutes ?? "null"));

        // The genuine metallic check: live router root = closed-form mean at the golden + silver stations.
        var golden = MetallicReadings.Single(r => r.C == 1.0);
        cases.Add(new BatteryCase(
            Name: "golden mean live = closed form (c=1 → φ)",
            Detail: $"|live {golden.LiveRatio.ToString("0.######", Inv)} − closed {golden.ClosedForm.ToString("0.######", Inv)}| < {MatchTolerance:0.0e+00}",
            Expected: "MATCH",
            Actual: golden.Matches ? "MATCH" : "MISMATCH"));

        var silver = MetallicReadings.Single(r => r.C == 2.0);
        cases.Add(new BatteryCase(
            Name: "silver mean live = closed form (c=2 → 1+√2)",
            Detail: $"|live {silver.LiveRatio.ToString("0.######", Inv)} − closed {silver.ClosedForm.ToString("0.######", Inv)}| < {MatchTolerance:0.0e+00}",
            Expected: "MATCH",
            Actual: silver.Matches ? "MATCH" : "MISMATCH"));

        // Where the direction comes from, exactly: the window's letter table and the uniform chiral.
        foreach (var t in LetterTableReadings)
            cases.Add(new BatteryCase(
                Name: $"c={t.C}: the locus is the letter table's eigen equation; â^⊗N is chiral exactly on it",
                Detail: $"det[v, M_c v] = r²−{t.C}r−1 for v = (r,1); R̃M_cR̃ − (r²+1)²M_c = (r²−{t.C}r−1)·K(r), " +
                        "K = [[4r, 2(1−r²)], [2(1−r²), −4r]] nowhere zero (integer polynomials in r)",
                Expected: "eigen, chiral iff locus",
                Actual: t.EigenEquationIsLocus && t.ChiralResidualIsLocusTimesK ? "eigen, chiral iff locus" : "FAIL"));

        // The control the divisibility test must fire on: the reflected line t₃ = −t₂ (M = [[1,1],[−1,0]],
        // hard from N = 5 by PROOF §8, the soft set) leaves a remainder modulo the golden locus.
        cases.Add(new BatteryCase(
            Name: "control: the reflected table [[1,1],[−1,0]] leaves a remainder mod r²−r−1",
            Detail: "R̃MR̃ − (r²+1)²M mod (r²−r−1) for the reflected line (hard from N = 5)",
            Expected: "remainder",
            Actual: ReflectedTableLeavesRemainder ? "remainder" : "FAIL"));

        // The rhythm as a count: S = ⊗[I,I,Z,Z] puts exactly one Z on every window's end pair, so it is
        // chiral; uniform and period-2 Z patterns put zero or two and commute (PROOF §2, "a second chiral").
        cases.Add(new BatteryCase(
            Name: "the rhythm is a chiral: one Z per window end pair has exactly the four period-4 solutions",
            Detail: "over the 16 Z patterns of period dividing 4, one Z on (w, w+2) at every offset w " +
                    "(period 3 cannot toggle an odd cycle)",
            Expected: "IIZZ, IZZI, ZIIZ, ZZII",
            Actual: string.Join(", ", RhythmSolutions)));

        // Why two-sided: the dephaser reads the product of the two factors' Z-parities at each site, so P and Q
        // must differ there; the golden chiral (lit everywhere) and the Z pattern (dead everywhere) do.
        string Verdicts(IEnumerable<SplitRouterReading> rs) => string.Join(", ", rs.Select(r => r.Routes ? "routes" : "fails"));
        cases.Add(new BatteryCase(
            Name: "why two-sided: ρ ↦ a^⊗N ρ S routes L exactly, either side, open N = 4, 5 and periodic ring N = 4",
            Detail: "W∘L = (−L − 2σ)∘W on every matrix unit, Gaussian integers, c = 3/2 (r = 2, a = 2X + Y), γ_l = l + 1: " +
                    string.Join("; ", SplitRouterReadings.Take(4).Select(r => r.Name)),
            Expected: "routes, routes, routes, routes",
            Actual: Verdicts(SplitRouterReadings.Take(4))));

        cases.Add(new BatteryCase(
            Name: "controls: equal Z-parity on both sides, a one-sided identity column, the periodic ring at N = 5",
            Detail: string.Join("; ", SplitRouterReadings.Skip(4).Select(r => r.Name)),
            Expected: "fails, fails, fails",
            Actual: Verdicts(SplitRouterReadings.Skip(4))));

        cases.Add(new BatteryCase(
            Name: "on a periodic ring the rhythm closes exactly at 4 | N",
            Detail: "Z patterns with one Z on every end pair (w, w+2 mod N), N = 3..9",
            Expected: "0, 4, 0, 0, 0, 4, 0",
            Actual: string.Join(", ", RingRhythmCounts)));

        // The soft set of the weighted template, exactly: N = 3, 4 every triple, N ≥ 5 the line and the plane.
        cases.Add(new BatteryCase(
            Name: "the soft set: det K₃ = −t₁²(t₂−t₃)²(t₁²+(t₂+t₃)²)², zero exactly on t₁ = 0 or t₂ = t₃",
            Detail: "the three-site chain constraint (P on two overlapping pairs), Bareiss over ℤ on the grid " +
                    "{−4..4}³, nine values per variable against degree ≤ 8",
            Expected: "identity",
            Actual: ChainDeterminantIdentityHolds ? "identity" : "FAIL"));
        cases.Add(new BatteryCase(
            Name: "the soft set: the far kernel is the product of two chain kernels (N = 5, 6, grid {−2..2}³)",
            Detail: "far kernel from the Pauli strings of H vs dim K_⌈N/2⌉ · dim K_⌊N/2⌋ from P alone, mismatches",
            Expected: "0",
            Actual: ReductionMismatches.ToString(Inv)));
        foreach (var r in SoftSetReadings)
            cases.Add(new BatteryCase(Name: "the soft set: " + r.Name, Detail: "PROOF §8, the soft set",
                                      Expected: r.Expected, Actual: r.Actual));

        // The soft set on the periodic ring: N = 4 every triple, 4 | N ≥ 8 the line and the plane, 4 ∤ N the plane.
        cases.Add(new BatteryCase(
            Name: "the ring soft set: N = 4 off t₁ = 0 and t₂ + t₃ = 0, an invertible element in the far kernel",
            Detail: "G = t₁²(XXYY+XYYX) − s²(XYXY+YXYX) + t₁s(YYYX+XYYY−XXYX−XYXX), s = t₂+t₃: [H, G] = 0 on " +
                    "{−2..2}³ (degree ≤ 3 per variable), det G = 2¹⁶t₁⁸s¹²(t₁²+s²)⁶ (binary form of degree 32, " +
                    "33 values, Bareiss over ℤ[i])",
            Expected: "identity",
            Actual: RingFourElementHolds ? "identity" : "FAIL"));
        cases.Add(new BatteryCase(
            Name: "the ring soft set: N = 4 on t₃ = −t₂, an invertible element in the far kernel",
            Detail: "2·XXYY + XYYX + YXXY + XXXX + XYXY + YXYX + YYYY: [H, G] = 0 on {−2..2}², det G = 2²⁴",
            Expected: "identity",
            Actual: RingFourLineElementHolds ? "identity" : "FAIL"));
        foreach (var r in RingSoftSetReadings)
            cases.Add(new BatteryCase(Name: "the ring soft set: " + r.Name, Detail: "PROOF §8, the soft set on rings",
                                      Expected: r.Expected, Actual: r.Actual));

        // The residual ‖{W,S}‖_F vanishes at the golden mean and is nonzero off it (the dip the witness draws).
        double phi = KBodyPalindromeRouting.MetallicMean(1.0);
        double atMean = KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, phi);
        double offMean = KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, phi + 0.25);
        cases.Add(new BatteryCase(
            Name: "residual ‖{W,S}‖_F vanishes at r = φ, nonzero off it",
            Detail: $"‖{{W,S}}‖_F(c=1, r=φ) = {atMean.ToString("0.##e+00", Inv)} (≈0), at r=φ+0.25 = {offMean.ToString("0.###", Inv)} (>0)",
            Expected: "vanishes",
            Actual: (atMean < MatchTolerance && offMean > MatchTolerance) ? "vanishes" : "FAIL"));

        return cases;
    }

    public override string DisplayName =>
        "F116: W L W⁻¹ = −L − 2σ (golden period-4 local palindromizer for the Z-middle ceiling, Tier1Derived)";

    public override string Summary
    {
        get
        {
            int matched = MetallicReadings.Count(r => r.Matches);
            return "an invertible period-4 per-site product W ([a,a,b,b] on the golden locus) palindromizes " +
                   "the Z-middle sliding-window ceiling H=Σ(XZX+XZY+YZX): W L W⁻¹ = −L − 2σ at every N≥3, " +
                   "arbitrary site γ; LOCAL, not non-local. Live: golden→" + (GoldenRoutes ?? "null") +
                   ", sibling→" + (SiblingRoutes ?? "null") + ", per-term declines=" + GoldenPerTermDeclines +
                   $", metallic r(c) live=closed at {matched}/{MetallicReadings.Count} c; {PassCount}/{Battery.Count} PASS ({Tier.Label()})";
        }
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the closed form (W, golden locus, two-sided form)",
                summary: "W = ⊗_l q_{l mod 4}, period-4 [a,a,b,b] with a = φX+Y, b = X−φY on the golden locus " +
                         "α²−αβ−β²=0 (slopes 1/φ and −φ, frame angle tan 2θ = 2), h_l = (−1)^(l+1)·i·R(g_l); " +
                         "q_l² = −(2+φ)·I so cond(W)=1. W L W⁻¹ = −L − 2σ, σ=Σ_l γ_l, at every N≥3 for arbitrary " +
                         "site-dependent γ_l≥0 (H = Σ_windows XZX+XZY+YZX and its X↔Y sibling). Two-sided form " +
                         "W(ρ) = (2+φ)^(N/2)·P ρ Q reduces {W,[H,·]}=0 to the chiral conditions P H P⁻¹ = −H, " +
                         "Q H Q⁻¹ = −H; G = PQ is a Hermitian product involution, an exact eigenmode at the floor " +
                         "L(G) = −2σ·G.");

            yield return new InspectableNode("why two-sided: the two chirals on the two sides",
                summary: "For ρ ↦ PρQ with product factors, the dephaser reads at each site the product of the two " +
                         "factors' Z-parities, and the class swap needs it to be −1 at every dephased site: there one " +
                         "factor carries a lit letter and the other a dead one. The Hamiltonian leg needs both " +
                         "factors chiral. So P = Q fails once any rate is nonzero, and no one-sided map reaches −L − 2σ (I is not chiral). The " +
                         "golden chiral a^⊗N is lit at every site and the Z pattern S = ⊗[I,I,Z,Z] dead at every " +
                         "site, so ρ ↦ a^⊗N ρ S is the router: per site §1's closed form with B signs [−,−,+,+] " +
                         "and h signs [+,−,+,−], a B-sign and an h-sign decoration of §5. On a periodic ring a^⊗N stays chiral at every N " +
                         "and S exists exactly at 4 | N; for 4 ∤ N and c ≠ 0, §4's alternation cannot " +
                         "close along the ring's parity cycles, so no per-site class-swapping router exists.");

            yield return new InspectableNode("the window lemma (cross-template, exact in the ring)",
                summary: "{q_w⊗q_{w+1}⊗q_{w+2}, [XZX+XZY+YZX,·]₃} = 0 exactly (ℤ[φ]+iℤ[φ]) at every window " +
                         "offset, while no single template's anticommutator vanishes: the cancellation is " +
                         "cross-template inside one window. Additivity over windows ⟹ {W,A}=0 at every N≥3 " +
                         "(single window included). This is why the per-term certifier cannot see the router; " +
                         "the per-term lens (KBodyPalindromeRouting.Routes) correctly returns false here.");

            yield return new InspectableNode("where the direction comes from: the window's letter table",
                summary: "The end letters of c·XZX+XZY+YZX carry the table M_c = [[c,1],[1,0]] over {X,Y}×{X,Y}; " +
                         "the middle Z is turned by every lit letter. A uniform product â^⊗N, â = αX+βY a unit " +
                         "direction, flips the middle Z and reflects both end letters, so it is chiral " +
                         "(â^⊗N H â^⊗N = −H) iff the " +
                         "reflection commutes with M_c, i.e. iff â is an EIGENVECTOR of M_c: α² − cαβ − β² = 0, " +
                         "the metallic locus, whose eigen equation r² = cr + 1 is the metallic mean (golden at " +
                         "c=1). Exact in integer polynomials in r: the residual is (r²−cr−1)·K(r) with K nowhere " +
                         "zero, live at c=0..3, with the reflected table (hard from N = 5) as the failing control. The uniform " +
                         "chiral alone gives spec(H) = −spec(H) but not the palindrome: conjugation by it fixes the " +
                         "dephaser; on one side, with the Z pattern on the other, it is the router (§2). The " +
                         "exclusion 'no uniform router' concerns the commuting G, not this chiral.");

            yield return new InspectableNode("the exclusion theorems (deductive from K₁, K₂)",
                summary: "{W,A}=0 on the identity column I^⊗N forces G = ⊗_l g_l to commute with H, reducing per " +
                         "window to two bilinears K₁ = α_w α_{w+2}+α_w β_{w+2}+β_w α_{w+2}, K₂ = α_w β_{w+2}+β_w " +
                         "α_{w+2}−β_w β_{w+2}. From them, deductively: no uniform router (g=0 only), no period-2; " +
                         "the golden locus α²−αβ−β²=0 is the only gate (kernel map exchanges the two locus " +
                         "directions, so each parity chain alternates a↔b); period-3 impossible for N≥5 (an " +
                         "alternating assignment cannot close an odd cycle; the N=4 period-3 solution is the " +
                         "now-explained 2026-06-07 small-N artifact); the discrete Klein candidates P1/P4/M2/M " +
                         "are all off-locus (the one-line reason the discrete-periodic search came back empty).");

            yield return new InspectableNode("the metallic family (polynomial identity in r)",
                summary: "on the soft line t₂=t₃, the router transports verbatim with a=(r,1), b=(1,−r), " +
                         "r(c)=(c+√(c²+4))/2 the metallic mean of c=t₁/t₂ (golden c=1→φ, silver c=2→1+√2, " +
                         "bronze c=3→(3+√13)/2, c=0→r=1 the 45° frame), q_l²=−(1+r²)·I, W_c L_c W_c⁻¹ = −L_c − 2σ " +
                         "for every real c. The window lemma is a polynomial identity in r (degree ≤ 5, verified " +
                         "exactly zero at 8 rational nodes via exact Fraction arithmetic), so it holds for all " +
                         "real, indeed all complex, c by derivation. The c≠0 identity-column determinant factors " +
                         "as c·(α²−cαβ−β²): the metallic locus is the gate.");

            yield return new InspectableNode("SCOPE FENCE: the c=0 station is soft (NOT Tier1Derived)",
                summary: "Rigidity (zero continuous physical moduli; the router unique modulo gauge) is PROVEN " +
                         "only at the golden and silver stations: finite-difference Jacobian nullity = gauge " +
                         "dimension at N=5,6. The c=0 station's '8 physical moduli' (nullity 16 = 8 gauge + 8 " +
                         "physical, spectral gap 3.5e9) is a FINITE-DIFFERENCE JACOBIAN COUNT AT N=5 ONLY, NOT a " +
                         "closed-form/analytic result; 4 of 8 moduli are catalogued (the [v,v,v̄,v̄] X-axis-mirror " +
                         "continuum + the period-2 Pauli-axis routers, the only place period-2 exists in the " +
                         "family). The first firing rung m* off the line and the plane is read at (1,2,1) only, " +
                         "not in closed form for general weights (the soft set itself is closed, see its node). " +
                         "The existence side carries no float anywhere. This soft sub-result is held below the " +
                         "Tier1Derived line.");

            yield return new InspectableNode("the soft set of t₁·XZX + t₂·XZY + t₃·YZX, exactly",
                summary: "N = 3, 4: the spectrum pairs for every real triple; N ≥ 5: iff t₁ = 0 or t₂ = t₃ (open chain, " +
                         "any positive site rates). On F158's far end a lit string carries X or Y on every site, the " +
                         "middle Z of a window always turns, so a window imposes one rule P = [[t₁,t₂,t₃,0],[0,t₃,t₂,−t₁]] " +
                         "on its end pair over (XX,XY,YX,YY), whatever the middle carries: the far kernel is the product " +
                         "of two chain kernels, over the even and the odd sites. A chain of at most two sites always holds " +
                         "a product of circle letters (a root of t₁t₃r² + (t₃²−t₁²−t₂²)r − t₁t₃, discriminant a sum of " +
                         "squares), an invertible element, so N = 3, 4 pair. A chain of three holds anything only where " +
                         "det K₃ = −t₁²(t₂−t₃)²(t₁²+(t₂+t₃)²)² vanishes; elsewhere the far end is empty while the near " +
                         "end holds the identity, so no N ≥ 5 pairs. Live: the determinant identity on {−4..4}³ and the " +
                         "Pythagorean witness (3,4,5), soft at N = 3, 4 in integers, far kernel 0 at N = 5, 6 (PROOF §8).");

            yield return new InspectableNode("the soft set on the periodic ring, exactly",
                summary: "N = 4: the spectrum pairs for every real triple; 4 | N ≥ 8: iff t₁ = 0 or t₂ = t₃; 4 ∤ N: iff " +
                         "t₁ = 0 (any positive site rates). Away from N = 4 each window still asks P of its own end pair, " +
                         "but the step w ↦ w + 2 mod N closes into cycles (one of length N at odd N, two of length N/2 at " +
                         "even N), each of at least three sites, so off the union det K₃ ≠ 0 empties the far end; on the " +
                         "line the words in a and b in which no letter follows itself close only along even cycles, so 4 | N. At " +
                         "N = 4 the windows w and w + 2 share the end pair {w, w + 2} with opposite middles and cancel " +
                         "jointly: off the line t₂ = t₃ the far end depends on t₁ and s = t₂ + t₃ only, and G = t₁²(XXYY+XYYX) − " +
                         "s²(XYXY+YXYX) + t₁s(YYYX+XYYY−XXYX−XYXX), det G = 2¹⁶t₁⁸s¹²(t₁²+s²)⁶, with a constant " +
                         "element on s = 0 (det 2²⁴) and X^⊗4 on t₁ = 0: the three carry the palindrome at every real triple. " +
                         "Live: both elements exact (commutation as a grid identity, determinants over ℤ[i]), far kernel " +
                         "mod p on rings N = 3..8 (PROOF §8, the soft set on rings).");

            yield return new InspectableNode("uniqueness, scoped to invertible W",
                summary: "Uniqueness statements are scoped to INVERTIBLE W: the anticommutation equation alone " +
                         "also carries abundant singular strata. Within the closed-form family at N=5, probed " +
                         "exhaustively (48 unbiased searches, every found zero matched), the invertible solution " +
                         "set is 4 cyclic shifts × an explicit order-32 sign-decoration group (128 gauge classes), " +
                         "every one a genuine palindromizer; the sign group is the diagonal ±1 (anti)commutant of " +
                         "A. Honest phrasing: essentially unique, modulo per-site gauge, the four pattern shifts, " +
                         "and the explicit structural sign group.");

            // Live battery (from-below re-run of the reused router machinery).
            foreach (var c in Battery)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));

            // The relationship to the soft-certifier (a see-cref, deliberately NOT a parent edge).
            yield return new InspectableNode("closes the soft-certifier ceiling (downstream, not a parent)",
                summary: "This router is the construction that closes PalindromeSoftCertifierClaim's k=3 windowed " +
                         "soft-certifier ceiling at zero (the 2→0 step of the 6→4→2→0 arc): the two Z-middle " +
                         "cases move from the ceiling to certified-positive witnesses (RoutingWindowSummed). The " +
                         "certifier USES this router as a helper, so it is logically DOWNSTREAM of this claim; it " +
                         "is referenced here, not wired as a parent (and the certifier is Tier1Candidate < " +
                         "Tier1Derived, so parenting on it would also violate strength inheritance). The live lab " +
                         "is GoldenRouterWitness (inspect --root router).");

            yield return new InspectableNode("parent: F1PalindromeIdentity",
                summary: "the global palindrome Π·L·Π⁻¹ = −L − 2Σγ·I; the golden router W realizes that same " +
                         "palindrome form LOCALLY (as a per-site product) for the Z-middle ceiling Hamiltonian " +
                         "class, where the standard Π argument does not directly give a per-site product.");

            yield return new InspectableNode("parent: PalindromeTwoEndCountClaim (F158)",
                summary: "the palindrome holds iff dim ker L = dim ker(L + 2σ), each end the commutator with H on " +
                         "a span of Pauli strings; the soft set of the weighted template is decided at its far end " +
                         "(the lit strings), which splits into two chain kernels, and any invertible element there " +
                         "carries the palindrome (the circle-letter products at N = 3, 4).");

            yield return new InspectableNode("parent: WindowedConverseThresholdClaim",
                summary: "the F87 two-reflection chiral spine (𝓕=F⊗F, R=I⊗F) and its chiral driving lemma " +
                         "F H F = −H (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2); the router's two-sided form " +
                         "puts a chiral on EACH side, of opposite Z-parity at every site (the golden a^⊗N and " +
                         "the Z pattern), the move that routes the Z-middle case where the single reflection " +
                         "F H F = −H fails on the cross-terms.");
        }
    }
}
