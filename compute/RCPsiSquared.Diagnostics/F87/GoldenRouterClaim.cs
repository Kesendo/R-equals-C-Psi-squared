using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

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
/// it holds for all real (indeed complex) c by derivation.</para>
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

    public GoldenRouterClaim(F1PalindromeIdentity f1Palindrome, WindowedConverseThresholdClaim chiralSpine)
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

    /// <summary>True iff the hard reflected table [[1,1],[−1,0]] (t₃ = −t₂) leaves a nonzero remainder modulo
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
        // hard by PROOF §8) leaves a remainder modulo the golden locus.
        cases.Add(new BatteryCase(
            Name: "control: the reflected table [[1,1],[−1,0]] leaves a remainder mod r²−r−1",
            Detail: "R̃MR̃ − (r²+1)²M mod (r²−r−1) for the hard reflected line",
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
                         "zero, live at c=0..3, with the hard reflected table as the failing control. The uniform " +
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
                         "family). The hard side off the line and the plane t1 = 0 (where X^N palindromizes) rests on the " +
                         "girth-ladder witness at (1,2,1) and " +
                         "the immediate-tilt scan, not on a closed-form p_{m*} for general weights. The existence " +
                         "side carries no float anywhere. This soft sub-result is held below the Tier1Derived line.");

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

            yield return new InspectableNode("parent: WindowedConverseThresholdClaim",
                summary: "the F87 two-reflection chiral spine (𝓕=F⊗F, R=I⊗F) and its chiral driving lemma " +
                         "F H F = −H (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2); the router's two-sided form " +
                         "puts a chiral on EACH side, of opposite Z-parity at every site (the golden a^⊗N and " +
                         "the Z pattern), the move that routes the Z-middle case where the single reflection " +
                         "F H F = −H fails on the cross-terms.");
        }
    }
}
