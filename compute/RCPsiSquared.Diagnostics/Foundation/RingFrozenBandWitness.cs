using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F170, the XY ring's frozen band at odd N (proof
/// <c>docs/proofs/PROOF_RING_FROZEN_BAND.md</c>, gate <c>simulations/ring_frozen_band.py</c>): on the
/// R₉₀ locus with γ̄ ≠ 0 every diagonal block (p, p), 1 ≤ p ≤ N − 1, carries the frozen root −4γ̄ exactly
/// ⌊N/2⌋ times for generic (J, δ), and N − 1 times at the uniform point for all but finitely many J, all of
/// it lowest weight on the bottom rung: nothing new is seeded above it. The reads below are taken at one
/// (J, γ̄) on one profile; an input on the exceptional set (a coupling where a higher rung gains modes or the
/// root turns defective) makes <see cref="Holds"/> false without contradicting the theorem.
///
/// <para><b>A different road from the gate, on purpose.</b> The Python gate builds every sector
/// combinatorially from the single-particle matrix. This witness never does for its counts: the spin ring
/// is a <see cref="PauliHamiltonian"/> of XX + YY bonds with the wrap bond (N − 1, 0) as one more Pauli
/// term, materialised as the full 2^N matrix, and every block is cut out of the Liouvillian by
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>, the same committed road the F140 witness takes.
/// The occupation pricing of the wrap bond is then not a modelling choice but what the Pauli algebra
/// does by itself, and the reduction to fixed-flux rings (Lemma 1) is read as a separate statement,
/// matrix against matrix (<see cref="SectorIdentityHolds"/>). The lowest-weight map Ψ(ρ) = Σ_l c_l ρ c_l†
/// is the transpose of F125's spectator W, from the committed builder
/// <see cref="SpectatorIntertwiner.BuildW"/>, so the typed parent's own object is what the lowest-weight
/// counts and the ladder check consume.</para>
///
/// <para><b>Counts and exact ranks, never a spectrum.</b> J and every γ_l sit on the dyadic grid
/// 1/<see cref="Grid"/>, so Grid·(L + 4γ̄) is an exact Gaussian-integer matrix, and every multiplicity is
/// its nullity over GF(p) at primes p ≡ 1 (mod 4), where i is a square root of −1. A rank modulo p can
/// only drop, so each count is an UPPER bound on the rational one. It is certified where it meets a proved
/// floor: ⌊N/2⌋ on the locus (F140's corner floor, carried up every rung by the ladder inside each
/// fixed-flux ring) and N − 1 at the uniform point (the traceless projectors on the plane waves, frozen at
/// every J); a lowest-weight count of 0 is certified outright. Where the first prime already meets the
/// floor the second is not asked, since it could only confirm. Reads above their floor run both primes and
/// stay upper bounds: the single-pair profile l = 1 at N = 9 (the others meet the floor), the N = 4
/// control, and the lowest-weight counts on rung 1 (where Ψ is the trace, so they equal the totals by trace
/// preservation and test Ψ not at all).</para>
///
/// <para><b>Algebraic multiplicity, where the block is small.</b> With A = Grid·(L + 4γ̄), the nullity of A²
/// modulo p bounds dim ker A² from above, and dim ker A² ≥ dim ker A; when it meets the certified kernel
/// dimension the two kernels coincide, so the root has no Jordan chain and its algebraic multiplicity is the
/// count. Read on every rung whose block has at most <see cref="AlgebraicCap"/> cells (every rung at N ≤ 5,
/// all but the two middle rungs at N = 7), for the locus and the uniform point.</para>
///
/// <para><b>Two-sided.</b> Beside the theorem's numbers sit zeros a broken builder or a mislabelled root
/// could not produce: the corner one grid step off the root, and the corner on a profile walked off the
/// locus at unchanged σ. N = 4 is read as a control: there two momenta can add up to π, so the resonance
/// lemma's argument is unavailable, and the middle rung reads 4 against ⌊N/2⌋ = 2 (an upper bound), the
/// extra two lowest weight on rung 2. The lemma's failure permits an excess and does not force one: on the
/// locus N = 6 reads ⌊N/2⌋ on every rung, and why N = 4 overfills is open. The resonance count itself is
/// the lemma written out over the momenta's definition, not a read of H, and cannot fail.</para>
///
/// <para><b>What pins Ψ.</b> The lowest-weight zeros above the bottom rung at odd N
/// cannot fail against an error in Ψ: a map that drops the Jordan-Wigner signs leaves them at 0, because
/// almost any linear map is injective on a frozen space of dimension ⌊N/2⌋ (N − 1 at the uniform point,
/// whose zeros are as blind), and only the N = 4 control
/// moves (its rung 2 reads 1 instead of 2). So Ψ is pinned directly, by the statement the band rests on:
/// Ψ·L₍p,p₎ = L'₍p−1,p−1₎·Ψ exactly in integers, L the spin ring and L' the spin ring with its wrap bond
/// reversed, since the rung p of the one and the rung p − 1 of the other are the same fixed-flux ring
/// (<see cref="LadderAcrossRings"/>); within one spin ring the residual is not zero
/// (<see cref="LadderWithinOneRing"/>). The sign-dropped Ψ fails it at 2J per entry. The same two-ring
/// statement in the raising direction, W from (p, q) to (p + 1, q + 1) with p ≡ q (mod 2), off-diagonal
/// blocks included, is gated in <c>SpectatorRingParityGateTests</c>; this check is its part-by-part
/// adjoint (the dissipator and the commutator with H are Hermitian superoperators) on the diagonal, where
/// the band lives. All three statements about the sign-dropped Ψ are gated in the tests.</para>
///
/// <para><b>The reversed ring is a consistency read at odd N.</b> With U = diag((−1)^l) the two fixed-flux
/// rings are conjugate up to J ↦ −J at odd N (the proof's §2), so <see cref="ReversedSpinRing"/> must repeat
/// <see cref="SpinRing"/>; it is kept because the ladder check needs both rings and because the sector
/// identity is read for both.</para>
///
/// <para><b>What stays outside.</b> The identity N·D̂ = ℓ² − Dmom on LW_ℓ ∩ V₀ and its floor
/// ℓ(N − ℓ)/N need the momentum operators n_k, defined over the cyclotomic field ℚ(ζ_N) rather than over
/// the Gaussian integers the ranks here use, so they live in the proof and in the gate's R3/R4 at rounding
/// precision; likewise the corner's large-J kernel span{P_k − P_{−k}} (R5). The defective corner pair at
/// N = 5 sits at a pair of couplings ±J* (read, not gated); A² shows the root non-defective at the first
/// coupling read here and is not asked at the second. The proof's sum formula and pencil count (PROOF_FROZEN_BAND_SO4 Corollary 4.2,
/// Propositions 5.1 and 5.3) are untyped. Even N stays a read.</para></summary>
public sealed class RingFrozenBandWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Smallest odd N with a balanced pair.</summary>
    public const int MinN = 3;

    /// <summary>Largest N whose every rung is read live: the middle rungs at N = 7 are 1225 × 1225.</summary>
    public const int MaxN = 7;

    /// <summary>The corner is swept over odd N up to here, both wraps (the 2^N Hamiltonian is 512 × 512).
    /// </summary>
    public const int CornerSweepMax = 9;

    /// <summary>The dyadic grid J and every γ_l are snapped onto.</summary>
    public const long Grid = 1024;

    /// <summary>Largest block, in cells, on which the algebraic multiplicity is read through A².</summary>
    public const int AlgebraicCap = 500;

    /// <summary>The composite N the single-pair profiles are read at.</summary>
    public const int ExceptionalN = 9;

    /// <summary>Two primes ≡ 1 (mod 4), different from the gate's.</summary>
    private static readonly long[] Primes = { 998244353L, 1004535809L };

    /// <summary>Half-profile offsets in grid units, mirrored antisymmetrically about γ̄; distinct
    /// magnitudes and mixed signs, so the profile is generic and not palindromic.</summary>
    private static readonly long[] HalfSeed = { 29, -13, 7, 19 };

    private const long OffLocusEps = 5;
    private const long SinglePairOffset = 37;

    /// <summary>One rung of one spin ring: its fixed-flux sign, the exact kernel reads, and the nullity of
    /// A² on the locus and at the uniform point (null above <see cref="AlgebraicCap"/>).</summary>
    public readonly record struct RungRead(
        int P, int SectorSign, int Locus, int LocusLowestWeight, int SecondCoupling, int Uniform,
        int UniformLowestWeight, int? LocusAlgebraic, int? UniformAlgebraic);

    /// <summary>The corner (1, 1) of one spin ring at one odd N.</summary>
    public readonly record struct CornerRead(int N, int Wrap, int Locus, int Uniform, int OffLocus, int WrongRoot);

    public int N { get; }
    public double J { get; }
    public double SecondJ { get; }
    public double GammaBar { get; }
    public IReadOnlyList<double> Gamma { get; }

    /// <summary>The theorem's generic count, ⌊N/2⌋.</summary>
    public int Generic => N / 2;

    /// <summary>The count at the uniform point, N − 1.</summary>
    public int AtUniform => N - 1;

    /// <summary>Every rung of the spin ring (wrap bond +J), whose rung p is the fixed-flux ring with
    /// s = (−1)^(p−1).</summary>
    public IReadOnlyList<RungRead> SpinRing { get; }

    /// <summary>Every rung of the spin ring with its wrap bond reversed (−J), whose rung p is the fixed-flux
    /// ring with the opposite sign; together the two cover both fixed-flux rings on every rung.</summary>
    public IReadOnlyList<RungRead> ReversedSpinRing { get; }

    /// <summary>The corner over odd N = 3..<see cref="CornerSweepMax"/>, both wraps.</summary>
    public IReadOnlyList<CornerRead> CornerSweep { get; }

    /// <summary>Lemma 1, entry for entry at this N, for both spin rings: every sector of the spin ring with
    /// wrap w equals the fixed-flux ring with s = w·(−1)^(p−1).</summary>
    public bool SectorIdentity { get; }

    /// <summary>The control for Lemma 1: at p = 2 the PERIODIC fermion ring differs from the spin ring.
    /// </summary>
    public bool PeriodicDiffersAtTwo { get; }

    /// <summary>Ordered pairs of allowed momenta adding up to π, both fixed-flux rings, at N and at N + 1.
    /// </summary>
    public (int Periodic, int PiFlux, int EvenPeriodic, int EvenPiFlux) Resonances { get; }

    /// <summary>The control N = 4, spin ring on a locus profile: total counts per rung and lowest-weight
    /// counts per rung. The read is [2, 4, 2] with [2, 2, 0].</summary>
    public (int[] Total, int[] LowestWeight) EvenControl { get; }

    /// <summary>The corner at N = 9 on a profile offset on the single pair (l, 8 − l), l = 0..3. Reads, not
    /// theorems: ⌊N/2⌋ + (gcd(2l + 1, N) − 1)/2, i.e. [4, 5, 4, 4].</summary>
    public IReadOnlyList<int> SinglePairCorners { get; }

    /// <summary>The ladder across the two rings: the largest entry of Ψ·L₍p,p₎ − L'₍p−1,p−1₎·Ψ over
    /// p = 2..N − 1, L the spin ring and L' the spin ring with its wrap reversed, and the other way round,
    /// in grid units. Exactly 0: the rung p of one and the rung p − 1 of the other are the same fixed-flux
    /// ring, where the η ladder commutes for any rates. This is what pins Ψ itself; the lowest-weight zeros
    /// at odd N cannot, since almost any map is injective on a frozen space that small (dropping Ψ's
    /// Jordan-Wigner signs leaves them at 0 and moves only the N = 4 control).</summary>
    public long LadderAcrossRings { get; }

    /// <summary>The same residual within ONE ring, the smallest over p = 2..N − 1: positive, because the
    /// rungs p and p − 1 of one spin ring see opposite fluxes. The ladder hands the band from ring to ring.
    /// </summary>
    public long LadderWithinOneRing { get; }

    public RingFrozenBandWitness(int n = 5, double gammaBar = 0.25, double j = 0.75)
    {
        if (n < MinN || n > MaxN || n % 2 == 0)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"N must be odd with {MinN} ≤ N ≤ {MaxN}: F170 is an odd-N theorem, and even N is read only " +
                "as the fixed N = 4 control inside this witness");
        if (!(gammaBar > 0.0) || gammaBar > 1.0)
            throw new ArgumentOutOfRangeException(nameof(gammaBar), gammaBar,
                "γ̄ must lie in (0, 1]. At γ̄ = 0 the statement is false (F140's zero-mean stratum, where the " +
                "corner carries N); a negative γ̄ (gain) is covered by the theorem but not read here");
        long gBarNum = Snap(gammaBar);
        long smallest = Math.Max(HalfSeed.Max(x => Math.Abs(x)), SinglePairOffset) + 1;
        if (gBarNum < smallest)
            throw new ArgumentOutOfRangeException(nameof(gammaBar), gammaBar,
                $"γ̄ snaps to {gBarNum}/{Grid}; the profiles here offset single rates by up to " +
                $"{smallest - 1}/{Grid}, so γ̄ must snap to at least {smallest}/{Grid} for every site rate to stay " +
                "positive, which is this witness's scope and not the theorem's");
        long jNum = Snap(j);
        if (jNum <= 0 || j > 2.0)
            throw new ArgumentOutOfRangeException(nameof(j), j, "J must lie in (0, 2] and not snap to 0");

        N = n;
        J = jNum / (double)Grid;
        SecondJ = 4.0 * J;
        GammaBar = gBarNum / (double)Grid;
        var locus = Locus(n, gBarNum);
        Gamma = locus.Select(g => g / (double)Grid).ToArray();

        SpinRing = Rungs(n, +1, jNum, gBarNum, locus);
        ReversedSpinRing = Rungs(n, -1, jNum, gBarNum, locus);

        var sweep = new List<CornerRead>();
        for (int m = 3; m <= CornerSweepMax; m += 2)
        {
            var g = Locus(m, gBarNum);
            var off = (long[])g.Clone();
            off[0] += OffLocusEps;
            off[1] -= OffLocusEps;
            var uniform = Enumerable.Repeat(gBarNum, m).ToArray();
            foreach (int wrap in new[] { +1, -1 })
            {
                var corner = Block(m, wrap, jNum, g, 1);
                sweep.Add(new CornerRead(m, wrap,
                    Nullity(corner, null, gBarNum, 0, m / 2),
                    Nullity(Block(m, wrap, jNum, uniform, 1), null, gBarNum, 0, m - 1),
                    Nullity(Block(m, wrap, jNum, off, 1), null, gBarNum, 0, 0),
                    Nullity(corner, null, gBarNum, 1, 0)));
            }
        }
        CornerSweep = sweep;

        SectorIdentity = SectorIdentityHolds(n, jNum, +1) && SectorIdentityHolds(n, jNum, -1);
        PeriodicDiffersAtTwo = !SectorMatches(n, 2, +1, +1, jNum);
        Resonances = (ResonantPairCount(n, +1), ResonantPairCount(n, -1),
                      ResonantPairCount(n + 1, +1), ResonantPairCount(n + 1, -1));

        var g4 = Locus(4, gBarNum);
        var tot4 = new int[3];
        var lw4 = new int[3];
        for (int p = 1; p <= 3; p++)
        {
            var b = Block(4, +1, jNum, g4, p);
            tot4[p - 1] = Nullity(b, null, gBarNum, 0, 2);
            lw4[p - 1] = Nullity(b, Psi(4, p), gBarNum, 0, 0);
        }
        EvenControl = (tot4, lw4);

        var single = new List<int>();
        for (int l = 0; l < ExceptionalN / 2; l++)
        {
            var g = Enumerable.Repeat(gBarNum, ExceptionalN).ToArray();
            g[l] += SinglePairOffset;
            g[ExceptionalN - 1 - l] -= SinglePairOffset;
            single.Add(Nullity(Block(ExceptionalN, +1, jNum, g, 1), null, gBarNum, 0, ExceptionalN / 2));
        }
        SinglePairCorners = single;

        long across = 0, within = long.MaxValue;
        for (int p = 2; p < n; p++)
        {
            across = Math.Max(across, LadderResidual(n, p, +1, -1, jNum, locus));
            across = Math.Max(across, LadderResidual(n, p, -1, +1, jNum, locus));
            within = Math.Min(within, LadderResidual(n, p, +1, +1, jNum, locus));
        }
        LadderAcrossRings = across;
        LadderWithinOneRing = within;
    }

    /// <summary>Every read meets the generic count, every control misses it where it must. False at an input
    /// on the theorem's exceptional set as well as under an error, so a false here is a question, not a
    /// refutation.</summary>
    public bool Holds =>
        SpinRing.Concat(ReversedSpinRing).All(r =>
            r.Locus == Generic && r.SecondCoupling == Generic && r.Uniform == AtUniform &&
            r.LocusLowestWeight == (r.P == 1 ? Generic : 0) &&
            r.UniformLowestWeight == (r.P == 1 ? AtUniform : 0) &&
            (r.LocusAlgebraic is null || r.LocusAlgebraic == Generic) &&
            (r.UniformAlgebraic is null || r.UniformAlgebraic == AtUniform)) &&
        CornerSweep.All(c => c.Locus == c.N / 2 && c.Uniform == c.N - 1 && c.OffLocus == 0 && c.WrongRoot == 0) &&
        SectorIdentity && PeriodicDiffersAtTwo &&
        LadderAcrossRings == 0 && LadderWithinOneRing > 0 &&
        Resonances.Periodic == 0 && Resonances.PiFlux == 0 &&
        Resonances.EvenPeriodic > 0 && Resonances.EvenPiFlux > 0 &&
        EvenControl.Total.SequenceEqual(new[] { 2, 4, 2 }) &&
        EvenControl.LowestWeight.SequenceEqual(new[] { 2, 2, 0 }) &&
        SinglePairCorners.Select((k, l) => k == SinglePairPrediction(ExceptionalN, l) &&
                                           k - ExceptionalN / 2 == BlindHeisenbergModes(ExceptionalN, l)).All(x => x);

    // ================= the statements, as static reads =================

    /// <summary>Ordered pairs (k, k') of allowed momenta of the fixed-flux ring with sign
    /// <paramref name="fermionWrap"/> with k + k' ≡ π (mod 2π). In units of π/N the momenta are 2m
    /// (s = +1) or 2m + 1 (s = −1), and π is N: zero at every odd N (Lemma 2), positive at every even N.
    /// </summary>
    public static int ResonantPairCount(int n, int fermionWrap)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "a ring needs N ≥ 2");
        if (fermionWrap != 1 && fermionWrap != -1)
            throw new ArgumentOutOfRangeException(nameof(fermionWrap), fermionWrap, "the wrap sign is ±1");
        int off = fermionWrap == 1 ? 0 : 1;
        int count = 0;
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
                if (((2 * a + off) + (2 * b + off) - n) % (2 * n) == 0) count++;
        return count;
    }

    /// <summary>F157's blind count of the seat <paramref name="j"/> on the uniform open Heisenberg chain, by
    /// enumerating the node condition (2j + 1)m ≡ N (mod 2N) over m = 1..N − 1, not by its gcd closed form,
    /// so that the match with the single-pair excess is a second route.</summary>
    public static int BlindHeisenbergModes(int n, int j)
    {
        int count = 0;
        for (int m = 1; m < n; m++)
            if (((long)m * (2 * j + 1) - n) % (2L * n) == 0) count++;
        return count;
    }

    /// <summary>The single-pair read at composite N: ⌊N/2⌋ + (gcd(2l + 1, N) − 1)/2.</summary>
    public static int SinglePairPrediction(int n, int l) =>
        n / 2 + ((int)BigInteger.GreatestCommonDivisor(2 * l + 1, n) - 1) / 2;

    /// <summary>Lemma 1 at one N: every popcount-p sector of the Hamiltonian of the spin ring with wrap
    /// <paramref name="spinWrap"/>, read off the full Pauli matrix, equals the fermion ring with wrap sign
    /// spinWrap·(−1)^(p−1) second-quantised with Slater signs. Compared exactly: the entries are 0 and ±J
    /// with J dyadic.</summary>
    public static bool SectorIdentityHolds(int n, long jNum, int spinWrap = +1)
    {
        for (int p = 1; p < n; p++)
            if (!SectorMatches(n, p, spinWrap, spinWrap * (p % 2 == 1 ? +1 : -1), jNum)) return false;
        return true;
    }

    private static bool SectorMatches(int n, int p, int spinWrap, int fermionWrap, long jNum)
    {
        var h = SpinRingHamiltonian(n, spinWrap, jNum / (double)Grid);
        var states = States(n, p);
        var fermion = FermionSector(n, states, fermionWrap, jNum / (double)Grid);
        for (int y = 0; y < states.Length; y++)
            for (int x = 0; x < states.Length; x++)
                if (h[states[y], states[x]] != fermion[y, x]) return false;
        return true;
    }

    /// <summary>The largest entry, in grid units, of Ψ·L₍p,p₎ − L'₍p−1,p−1₎·Ψ with L the spin ring of wrap
    /// <paramref name="sourceWrap"/> and L' that of wrap <paramref name="targetWrap"/>, both at the rates
    /// <paramref name="gammaNum"/>. Integer arithmetic throughout.</summary>
    public static long LadderResidual(int n, int p, int sourceWrap, int targetWrap, long jNum, long[] gammaNum,
        long[][]? psiOverride = null)
    {
        if (p < 2 || p > n - 1) throw new ArgumentOutOfRangeException(nameof(p), p, $"the ladder steps from rungs 2..{n - 1}");
        var src = Block(n, sourceWrap, jNum, gammaNum, p).Block;
        var dst = Block(n, targetWrap, jNum, gammaNum, p - 1).Block;
        var psi = psiOverride ?? Psi(n, p);
        int rows = psi.Length, cols = src.RowCount;
        long worst = 0;
        var vr = new long[rows];
        var vi = new long[rows];
        for (int j = 0; j < cols; j++)
        {
            Array.Clear(vr);
            Array.Clear(vi);
            for (int i = 0; i < cols; i++)
            {
                var z = src[i, j];
                if (z == Complex.Zero) continue;
                long zr = Exact(z.Real), zi = Exact(z.Imaginary);
                for (int k = 0; k < rows; k++)
                {
                    long s = psi[k][i];
                    if (s == 0) continue;
                    vr[k] += s * zr;
                    vi[k] += s * zi;
                }
            }
            for (int k = 0; k < rows; k++)
            {
                long s = psi[k][j];
                if (s == 0) continue;
                for (int r = 0; r < rows; r++)
                {
                    var z = dst[r, k];
                    if (z == Complex.Zero) continue;
                    vr[r] -= s * Exact(z.Real);
                    vi[r] -= s * Exact(z.Imaginary);
                }
            }
            for (int k = 0; k < rows; k++) worst = Math.Max(worst, Math.Max(Math.Abs(vr[k]), Math.Abs(vi[k])));
        }
        return worst;
    }

    /// <summary>The lowest-weight frozen count on the rung p of the spin ring with wrap
    /// <paramref name="wrap"/>: the nullity of [Grid·(L + 4γ̄); Ψ] modulo the primes, an upper bound, with
    /// Ψ replaceable so that a wrong Ψ can be shown to the reads.</summary>
    public static int LowestWeightCount(int n, int p, int wrap, long jNum, long[] gammaNum, long gBarNum,
        long[][]? psiOverride = null) =>
        Nullity(Block(n, wrap, jNum, gammaNum, p), psiOverride ?? Psi(n, p), gBarNum, 0, 0);

    // ================= the atoms =================

    private static long Snap(double v) => (long)Math.Round(v * Grid);

    private static long[] Locus(int n, long gBarNum)
    {
        var g = Enumerable.Repeat(gBarNum, n).ToArray();
        for (int i = 0; i < n / 2; i++)
        {
            long d = HalfSeed[i % HalfSeed.Length];
            g[i] = gBarNum + d;
            g[n - 1 - i] = gBarNum - d;
        }
        return g;
    }

    private static IReadOnlyList<RungRead> Rungs(int n, int wrap, long jNum, long gBarNum, long[] locus)
    {
        var uniform = Enumerable.Repeat(gBarNum, n).ToArray();
        var reads = new List<RungRead>();
        for (int p = 1; p < n; p++)
        {
            var psi = Psi(n, p);
            var atLocus = Block(n, wrap, jNum, locus, p);
            var atUniform = Block(n, wrap, jNum, uniform, p);
            int sign = wrap * (p % 2 == 1 ? 1 : -1);
            bool small = atLocus.Block.RowCount <= AlgebraicCap;
            reads.Add(new RungRead(p, sign,
                Nullity(atLocus, null, gBarNum, 0, n / 2),
                Nullity(atLocus, psi, gBarNum, 0, 0),
                Nullity(Block(n, wrap, 4 * jNum, locus, p), null, gBarNum, 0, n / 2),
                Nullity(atUniform, null, gBarNum, 0, n - 1),
                Nullity(atUniform, psi, gBarNum, 0, 0),
                small ? SquareNullity(atLocus, gBarNum, n / 2) : null,
                small ? SquareNullity(atUniform, gBarNum, n - 1) : null));
        }
        return reads;
    }

    /// <summary>The spin ring: (J/2)(XX + YY) on the N − 1 chain bonds and wrap·(J/2)(XX + YY) on the bond
    /// (N − 1, 0), as one Pauli Hamiltonian.</summary>
    private static ComplexMatrix SpinRingHamiltonian(int n, int wrap, double j)
    {
        var terms = new List<PauliTerm>(2 * n);
        for (int b = 0; b < n; b++)
        {
            int s1 = b, s2 = (b + 1) % n;
            Complex c = (b == n - 1 ? wrap : 1) * j / 2.0;
            terms.Add(PauliTerm.TwoSite(n, s1, PauliLetter.X, s2, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(n, s1, PauliLetter.Y, s2, PauliLetter.Y, c));
        }
        return new PauliHamiltonian(n, terms).ToMatrix();
    }

    private static (ComplexMatrix Block, int[] Rows, int[] Cols) Block(int n, int wrap, long jNum, long[] gammaNum, int p)
    {
        var h = SpinRingHamiltonian(n, wrap, jNum / (double)Grid);
        var states = States(n, p);
        int d = 1 << n;
        var flat = new int[states.Length * states.Length];
        var rows = new int[flat.Length];
        var cols = new int[flat.Length];
        int k = 0;
        foreach (int r in states)
            foreach (int c in states)
            {
                rows[k] = r;
                cols[k] = c;
                flat[k++] = r * d + c;
            }
        var gamma = gammaNum.Select(g => g / (double)Grid).ToArray();
        return (PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat), rows, cols);
    }

    /// <summary>Basis states of popcount p, ascending.</summary>
    private static int[] States(int n, int p) =>
        Enumerable.Range(0, 1 << n).Where(x => BitOperations.PopCount((uint)x) == p).ToArray();

    /// <summary>Site l is bit N − 1 − l (site 0 the leftmost Kronecker factor).</summary>
    private static int Bit(int n, int site) => 1 << (n - 1 - site);

    /// <summary>(−1) to the number of occupied sites before <paramref name="site"/>.</summary>
    private static int JordanWignerSign(int n, int state, int site)
    {
        int before = 0;
        for (int s = 0; s < site; s++) if ((state & Bit(n, s)) != 0) before++;
        return before % 2 == 0 ? 1 : -1;
    }

    private static double[,] FermionSector(int n, int[] states, int fermionWrap, double j)
    {
        var index = new Dictionary<int, int>();
        for (int i = 0; i < states.Length; i++) index[states[i]] = i;
        var hop = new double[n, n];
        for (int b = 0; b < n; b++)
        {
            int s1 = b, s2 = (b + 1) % n;
            double c = (b == n - 1 ? fermionWrap : 1) * j;
            hop[s1, s2] = hop[s2, s1] = c;
        }
        var f = new double[states.Length, states.Length];
        for (int x = 0; x < states.Length; x++)
        {
            int st = states[x];
            for (int s = 0; s < n; s++)
            {
                if ((st & Bit(n, s)) == 0) continue;
                int rest = st & ~Bit(n, s);
                for (int t = 0; t < n; t++)
                {
                    if (hop[t, s] == 0 || (rest & Bit(n, t)) != 0) continue;
                    int sign = JordanWignerSign(n, st, s) * JordanWignerSign(n, rest, t);
                    f[index[rest | Bit(n, t)], x] += sign * hop[t, s];
                }
            }
        }
        return f;
    }

    /// <summary>Ψ(ρ) = Σ_l c_l ρ c_l† from the block (p, p) to (p − 1, p − 1), as integer rows over the
    /// columns of the (p, p) block: the transpose of F125's spectator W from (p − 1, p − 1), taken from the
    /// committed builder <see cref="SpectatorIntertwiner.BuildW"/>, whose basis (popcount states ascending,
    /// ket outer, bra inner) is the one <see cref="Block"/> lays out. W's entries are ±1 and real, so the
    /// adjoint is the transpose.</summary>
    public static long[][] Psi(int n, int p)
    {
        var w = SpectatorIntertwiner.BuildW(n, p - 1, p - 1);
        var rows = new long[w.ColumnCount][];
        for (int i = 0; i < w.ColumnCount; i++)
        {
            rows[i] = new long[w.RowCount];
            for (int k = 0; k < w.RowCount; k++)
            {
                var z = w[k, i];
                if (z.Imaginary != 0 || z.Real != Math.Round(z.Real))
                    throw new InvalidOperationException($"W entry {z} is not an integer");
                rows[i][k] = (long)z.Real;
            }
        }
        return rows;
    }

    /// <summary>dim ker of [Grid·(L_block + 4γ̄) + shift; extra] over GF(p), the smallest over the primes
    /// asked; stops after the first prime when that already meets <paramref name="floor"/>, since a
    /// second prime can only confirm.</summary>
    private static int Nullity((ComplexMatrix Block, int[] Rows, int[] Cols) block, long[][]? extra,
        long gBarNum, long shift, int floor)
    {
        var b = block.Block;
        int d = b.RowCount;
        var re = new long[d][];
        var im = new long[d][];
        for (int i = 0; i < d; i++)
        {
            re[i] = new long[d];
            im[i] = new long[d];
            for (int k = 0; k < d; k++)
            {
                re[i][k] = Exact(b[i, k].Real) + (i == k ? 4 * gBarNum + shift : 0);
                im[i][k] = Exact(b[i, k].Imaginary);
            }
        }
        int best = int.MaxValue;
        foreach (long prime in Primes)
        {
            long r = SqrtMinusOne(prime);
            int extraRows = extra?.Length ?? 0;
            var a = new long[d + extraRows][];
            for (int i = 0; i < d; i++)
            {
                a[i] = new long[d];
                for (int k = 0; k < d; k++)
                    a[i][k] = Mod(Mod(re[i][k], prime) + Mod(r * Mod(im[i][k], prime), prime), prime);
            }
            for (int i = 0; i < extraRows; i++)
                a[d + i] = extra![i].Select(v => Mod(v, prime)).ToArray();
            best = Math.Min(best, d - RankModP(a, d, prime));
            if (best <= floor) break;
        }
        return best;
    }

    /// <summary>dim ker A² over GF(p), A = Grid·(L_block + 4γ̄), the smallest over the primes asked; an upper
    /// bound on the rational dim ker A², stopping once it meets <paramref name="floor"/>.</summary>
    private static int SquareNullity((ComplexMatrix Block, int[] Rows, int[] Cols) block, long gBarNum, int floor)
    {
        var b = block.Block;
        int d = b.RowCount;
        int best = int.MaxValue;
        foreach (long prime in Primes)
        {
            long r = SqrtMinusOne(prime);
            var a = new long[d][];
            for (int i = 0; i < d; i++)
            {
                a[i] = new long[d];
                for (int k = 0; k < d; k++)
                    a[i][k] = Mod(Mod(Exact(b[i, k].Real) + (i == k ? 4 * gBarNum : 0), prime)
                                  + Mod(r * Mod(Exact(b[i, k].Imaginary), prime), prime), prime);
            }
            var sq = new long[d][];
            for (int i = 0; i < d; i++)
            {
                var row = new long[d];
                var ai = a[i];
                for (int m = 0; m < d; m++)
                {
                    long x = ai[m];
                    if (x == 0) continue;
                    var am = a[m];
                    for (int k = 0; k < d; k++)
                        if (am[k] != 0) row[k] = (row[k] + x * am[k]) % prime;
                }
                sq[i] = row;
            }
            best = Math.Min(best, d - RankModP(sq, d, prime));
            if (best <= floor) break;
        }
        return best;
    }

    private static long Exact(double v)
    {
        double s = v * Grid;
        double r = Math.Round(s);
        if (s != r)
            throw new InvalidOperationException(
                $"block entry {v.ToString("R", Inv)} is not a multiple of 1/{Grid}: the dyadic grid was broken");
        return (long)r;
    }

    private static long Mod(long a, long p) { long m = a % p; return m < 0 ? m + p : m; }

    /// <summary>Forward elimination over GF(p), skipping the zeros of the pivot row; entries stay below
    /// p &lt; 2^30, so every product fits in a long.</summary>
    private static int RankModP(long[][] a, int cols, long p)
    {
        int rank = 0;
        var nz = new List<int>(cols);
        for (int col = 0; col < cols && rank < a.Length; col++)
        {
            int piv = -1;
            for (int row = rank; row < a.Length; row++) if (a[row][col] != 0) { piv = row; break; }
            if (piv < 0) continue;
            (a[rank], a[piv]) = (a[piv], a[rank]);
            var pr = a[rank];
            long inv = ModPow(pr[col], p - 2, p);
            nz.Clear();
            for (int k = col; k < cols; k++) if (pr[k] != 0) nz.Add(k);
            for (int row = rank + 1; row < a.Length; row++)
            {
                var ar = a[row];
                if (ar[col] == 0) continue;
                long f = ar[col] * inv % p;
                foreach (int k in nz) ar[k] = Mod(ar[k] - f * pr[k] % p, p);
            }
            rank++;
        }
        return rank;
    }

    private static long SqrtMinusOne(long p)
    {
        for (long t = 2; t < 200; t++)
        {
            long r = ModPow(t, (p - 1) / 4, p);
            if (r * r % p == p - 1) return r;
        }
        throw new InvalidOperationException($"no square root of −1 mod {p}");
    }

    private static long ModPow(long b, long e, long p)
    {
        long acc = 1;
        b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) acc = acc * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return acc;
    }

    // ================= the reading =================

    private static string Fmt(double v) => v.ToString("0.######", Inv);

    private static string Row(IEnumerable<int> xs) => "[" + string.Join(", ", xs) + "]";

    public string DisplayName =>
        $"RingFrozenBandWitness (F170: every rung of the XY ring carries the corner's ⌊N/2⌋, live; N={N})";

    public string Summary =>
        $"N={N}, J={Fmt(J)}, R₉₀ locus γ=[{string.Join(", ", Gamma.Select(Fmt))}] (γ̄={Fmt(GammaBar)}): the spin " +
        $"ring's rungs p = 1..{N - 1} carry −4γ̄ {Row(SpinRing.Select(r => r.Locus))} times on the locus, " +
        $"{Row(SpinRing.Select(r => r.Uniform))} at the uniform point, lowest weight " +
        $"{Row(SpinRing.Select(r => r.LocusLowestWeight))}, so nothing new is seeded above the corner; " +
        $"predicted ⌊N/2⌋ = {Generic} and N − 1 = {AtUniform}. The reversed wrap gives " +
        $"{Row(ReversedSpinRing.Select(r => r.Locus))}. Resonant momentum pairs: {Resonances.Periodic} and " +
        $"{Resonances.PiFlux} at N={N}, {Resonances.EvenPeriodic} and {Resonances.EvenPiFlux} at N={N + 1}. " +
        $"The algebraic multiplicity is read through A² on {SpinRing.Concat(ReversedSpinRing).Count(r => r.LocusAlgebraic is not null)} " +
        $"of {2 * (N - 1)} rungs. " +
        $"{(Holds ? "Every read meets the generic count." : "At least one read misses the generic count: an input on the exceptional set, or an error.")}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                $"the resonance lemma: {Resonances.Periodic + Resonances.PiFlux} pairs at N={N}, " +
                $"{Resonances.EvenPeriodic + Resonances.EvenPiFlux} at N={N + 1}",
                summary: "in units of π/N the momenta are even (periodic) or odd (π flux) integers, so two of " +
                         "them add up to an even number, and π is N: at odd N no pair reaches it. The lemma " +
                         "written out over the momenta's definition (ordered pairs with k + k' ≡ N mod 2N), not a " +
                         "read of H; it cannot fail",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the reduction: sectors equal entry for entry = {SectorIdentity}, control differs = " +
                $"{PeriodicDiffersAtTwo}",
                summary: "every popcount-p sector of the Pauli Hamiltonian of the spin ring with wrap w equals " +
                         "the fermion ring with wrap sign w·(−1)^(p−1), built with Slater signs, for w = ±1; the " +
                         "wrap hop passes the other p − 1 particles. Control: at p = 2 the periodic fermion ring " +
                         "is NOT the spin ring",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the ladder: Ψ·L − L'·Ψ = {LadderAcrossRings} across the two rings, ≥ {LadderWithinOneRing} within one",
                summary: "Ψ(ρ) = Σ_l c_l ρ c_l† intertwines the rung p of the spin ring with the rung p − 1 of the " +
                         "spin ring whose wrap bond is reversed, exactly, in integers, because the two are one " +
                         "fixed-flux ring; within one spin ring it does not, the two rungs seeing opposite fluxes. " +
                         "This pins Ψ itself, which the lowest-weight zeros below cannot do",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the rungs, spin ring",
                children: SpinRing.Select(RungNode).ToArray(), provenance: NodeProvenance.Live);
            yield return new InspectableNode("the rungs, wrap bond reversed",
                children: ReversedSpinRing.Select(RungNode).ToArray(), provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the corner over odd N = 3..{CornerSweepMax}, both wraps",
                summary: string.Join("; ", CornerSweep.Select(c =>
                    $"N={c.N} wrap {(c.Wrap > 0 ? "+" : "−")}: locus {c.Locus}, uniform {c.Uniform}, off locus " +
                    $"{c.OffLocus}, one grid step off the root {c.WrongRoot}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the even control N = 4: rungs {Row(EvenControl.Total)}, lowest weight {Row(EvenControl.LowestWeight)}",
                summary: "upper bounds modulo primes. At even N two momenta can add up to π, the resonance " +
                         "lemma's argument is unavailable, and here the middle rung reads 4 against ⌊N/2⌋ = 2, " +
                         "the extra two lowest weight on rung 2. The lemma's failure permits this and does not " +
                         "force it: N = 6 reads ⌊N/2⌋ on every rung of the locus, and why N = 4 overfills is open",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"the exceptional set at N = {ExceptionalN}: single-pair corners {Row(SinglePairCorners)}",
                summary: "an offset on the single pair (l, N − 1 − l) misses the second-order lift on the momenta " +
                         "whose order divides gcd(2l + 1, N), giving ⌊N/2⌋ + (gcd(2l + 1, N) − 1)/2; a read, " +
                         "certified where it meets the floor and an upper bound at l = 1. The excess equals F157's " +
                         $"blind Heisenberg count of the seat l, {Row(Enumerable.Range(0, ExceptionalN / 2).Select(l => BlindHeisenbergModes(ExceptionalN, l)))} " +
                         "by node enumeration: each missed momentum pair holds one blind mode, the same harmonic " +
                         "read at l + ½ (PROOF_RING_FROZEN_BAND §7)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("what stays outside",
                summary: "the identity N·D̂ = ℓ² − Dmom and its floor ℓ(N − ℓ)/N (irrational plane waves; gate " +
                         "R3/R4), the large-J corner span{P_k − P_{−k}} (R5), the defective corner pair at N = 5 " +
                         "(two couplings, not met here), the sum formula and pencil count of PROOF_FROZEN_BAND_SO4 " +
                         "(untyped), even N beyond the control");
        }
    }

    private IInspectable RungNode(RungRead r) =>
        new InspectableNode(
            $"rung {r.P} (fixed flux s = {(r.SectorSign > 0 ? "+1" : "−1")}): locus {r.Locus}, second coupling " +
            $"{r.SecondCoupling}, uniform {r.Uniform}; lowest weight {r.LocusLowestWeight} / {r.UniformLowestWeight}; " +
            $"dim ker A² {(r.LocusAlgebraic is null ? "not read (block too large)" : $"{r.LocusAlgebraic} / {r.UniformAlgebraic}")}",
            provenance: NodeProvenance.Live);

    public InspectablePayload Payload => InspectablePayload.Empty;
}
