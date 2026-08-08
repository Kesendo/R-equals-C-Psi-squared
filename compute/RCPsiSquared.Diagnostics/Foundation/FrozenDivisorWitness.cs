using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Which chain the frozen divisor is read on. The divisor BOUND holds on both; the
/// four-corner census does not (see <see cref="FrozenDivisorWitness.CornersAreTheOnlyCarriers"/>).
/// </summary>
public enum DivisorChain
{
    /// <summary>H = (J/4) Σ_b (XX + YY + ZZ): the chain the proof document's census was verified
    /// on.</summary>
    Heisenberg,

    /// <summary>H = (J/2) Σ_b (XX + YY): the same divisor, but the confinement to the corners is
    /// gone.</summary>
    Xy,
}

/// <summary>The live lab for the R90 frozen divisor (F140, proof
/// <c>docs/proofs/PROOF_R90_FROZEN_DIVISOR.md</c>, gate
/// <c>simulations/r90_frozen_divisor_gate.py</c>): on the mirror-balanced (anti-palindromic) locus of F91
/// (every reflection pair of site rates carrying the same total, γ_l + γ_{R(l)} = 2γ̄), the
/// single-excitation corner block of the Liouvillian carries
///
/// <code>
///     λ = −4γ̄   with multiplicity ≥ ⌊N/2⌋,   at EVERY coupling J
/// </code>
///
/// one frozen mode per balanced pair, and, for γ̄ ≠ 0, exactly ⌊N/2⌋ for all but finitely many J
/// (the proof's Sections 6 and 7). Two of the three exceptions are visible from here: J = 0, where
/// the multiplicity merely doubles to 2⌊N/2⌋ and stays semisimple, and the real exceptional
/// couplings of Section 9, where the root goes defective without its kernel dimension moving. The
/// third is a whole stratum rather than a point and this witness does not enter it: at γ̄ = 0 the
/// count is N for all but finitely many J (see the constructor, and MirrorWorld's
/// <c>Divisor</c>). The value hears only the MEAN rate: no J, no individual rate, nothing of
/// the Hamiltonian. Nothing symmetric is behind it; what pins it is a room
/// shortage of the cell mirror τQ : (a,b) ↦ (R(b), R(a)), which fixes exactly the 2⌊N/2⌋
/// anti-diagonal coherences, so an odd operator has more even rooms than odd ones and the surplus
/// must freeze. The diagonal cells' mirror pairs are EVEN, not odd, and tax away half:
/// ⌊N/2⌋ = 2⌊N/2⌋ − ⌊N/2⌋.
///
/// <para><b>Counts and exact ranks, never a spectrum.</b> The frozen root is an exact eigenvalue
/// at every coupling, and the other eigenvalues crowd it at spacing J^{2d}, so a floating-point
/// eigen- or SVD-count silently miscounts once the coupling is small and the chain long. Every
/// multiplicity here is instead <c>dim ker</c> by elimination over GF(p) at two primes p ≡ 1
/// (mod 4), where i is a genuine square root of −1: J and every γ_l are snapped onto the dyadic
/// grid 1/1024, so the block entries are exact in double and 4N·1024·(entries) are integers.</para>
///
/// <para><b>Independent of the MirrorWorld <c>Divisor</c>.</b> That object (run mode
/// <c>divisor N</c>) hand-builds the single-excitation handshake from the hopping matrix. This
/// witness goes the long way round through the repo's own committed primitives, and lands on the
/// same counts: <see cref="PauliHamiltonian"/> materialises the full 2^N Hamiltonian and
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> cuts the joint-popcount block out of the
/// Liouvillian. Two constructions, one number.</para>
///
/// <para><b>The census is a Heisenberg statement, and this witness is where that was found
/// (2026-07-25).</b> The proof document censuses every (p, q) block at N = 4, 5, 6 and finds only
/// the four corners carrying, but it runs that census on the Heisenberg chain alone. Drop the ZZ
/// term and the confinement goes: on the XY chain the two roots are carried, with the SAME
/// multiplicity ⌊N/2⌋, by many more blocks. Counting a block that carries either: nine of
/// twenty-five at N = 4 and twenty-four of thirty-six at N = 5, both of which this witness computes
/// itself, and twenty-one of forty-nine at N = 6, which is the gate's G2c rather than this
/// witness's (the full census stops at
/// <see cref="FullCensusMaxN"/>, so at N = 6 <see cref="CensusRan"/> is false). Never partially,
/// and every block carrying the UNFOLDED root −4γ̄ has p + q even (necessary, not sufficient:
/// plenty of even blocks carry nothing). The bound and the fold-parity split survive the change of
/// chain; "only the corners" does not. The invocation that shows it is
/// <c>inspect --root divisor --N 5 --chain xy</c>; read
/// <see cref="CornersAreTheOnlyCarriers"/> against <see cref="Chain"/>.</para>
///
/// <para><b>Two-sided by design.</b> Beside the ⌊N/2⌋ reads sit zeros that a broken rank routine
/// or a mislabelled root could not produce: the same block at a root one grid step away
/// (<see cref="WrongRootKernelDimension"/>) and the same profile walked OFF the locus at unchanged
/// σ and γ̄ (<see cref="OffLocusKernelDimension"/>, where partial balance yields nothing). The J = 0
/// corollary (<see cref="ZeroCouplingKernelDimension"/> = 2⌊N/2⌋, twice the multiplicity, from
/// which exactly ⌊N/2⌋ modes depart as the coupling turns on) is a third kind of read: not a zero,
/// not the theorem's number, and wrong under either mistake.</para>
///
/// <para><b>Scope.</b> This witness types the divisor BOUND (≥ ⌊N/2⌋ at every coupling, which
/// holds at every γ̄) and the fold-parity placement of the corners. Its READS are all taken at
/// γ̄ &gt; 0, the part of the taxed stratum γ̄ ≠ 0 where every site rate can also be positive.
/// Tightness is a theorem too (Sections 6 and 7, and it is the half that needs γ̄ ≠ 0), and the
/// ranks here meet it at every coupling sampled, but sampling is not what proves it. What no rank
/// read here can see: at the real exceptional couplings the root goes DEFECTIVE while its kernel
/// dimension does not move. What this witness deliberately does not REACH: the zero-mean stratum
/// γ̄ = 0, where the count is N rather than ⌊N/2⌋ and every ⌊N/2⌋ read below would be a mismatch.
/// The constructor closes that door rather than branching on it, because the whole read-set here
/// is built around one number; the stratum is pinned instead by the gate's G16 and by MirrorWorld's
/// <c>Divisor.Rooms()</c>, which carries both counts. The analytic half (the cofactor determinant,
/// the two boundary clocks, those defective couplings) is eigen-work and stays in the proof
/// document, as does the one item still open there, the upper half of the valuation
/// law.</para></summary>
public sealed class FrozenDivisorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Largest N the live build materialises: the 2^N Hamiltonian and the (N−1, N−1)
    /// corner keep the cost honest (N = 8 → a 256 × 256 H and 64 × 64 corner blocks).</summary>
    public const int MaxN = 8;

    /// <summary>Smallest N with a balanced pair AND a second site to unbalance for the off-locus
    /// control.</summary>
    public const int MinN = 3;

    /// <summary>Above this N the full (p, q) census is skipped: the exact rank is cubic and the
    /// middle blocks grow as C(N, N/2)².</summary>
    public const int FullCensusMaxN = 5;

    /// <summary>Above this block size the non-corner sample is skipped.</summary>
    public const int ControlBlockSizeCap = 500;

    /// <summary>The dyadic grid J and every γ_l are snapped onto. A power of two, so the snapped
    /// values are exact doubles and every block entry stays exact through the build.</summary>
    public const long Grid = 1024;

    /// <summary>Two primes ≡ 1 (mod 4), so −1 is a square and the Gaussian integers map into F_p.
    /// A rank mod p can only DROP at a bad prime, so the max over two pins the true rank.</summary>
    private static readonly long[] Primes = { 998244353L, 1004535809L };

    /// <summary>Half-profile numerators over <see cref="Grid"/>, mirrored antisymmetrically about
    /// the mean; deliberately irregular (no arithmetic progression, mixed signs, distinct
    /// magnitudes) so the locus is not accidentally uniform or palindromic.</summary>
    private static readonly long[] HalfSeed = { 29, -13, 7, 19 };

    /// <summary>A second half-profile with the same γ̄: the frozen root must not notice.</summary>
    private static readonly long[] SecondHalfSeed = { -5, 23, -31, 3 };

    /// <summary>The off-locus control's step, in grid units.</summary>
    private const long OffLocusEps = 5;

    public int N { get; }

    /// <summary>Which chain the reads below are taken on.</summary>
    public DivisorChain Chain { get; }

    /// <summary>The coupling, snapped onto the dyadic grid.</summary>
    public double J { get; }
    /// <summary>The second coupling the J-independence face is read at, snapped.</summary>
    public double SecondJ { get; }
    /// <summary>Mean dephasing rate γ̄ = σ/N, snapped; every site rate is γ̄ ± a grid step.</summary>
    public double GammaBar { get; }
    /// <summary>Total dephasing rate σ = Σ γ_l.</summary>
    public double Sigma { get; }
    /// <summary>The on-locus per-site rates, in SITE order (site 0 = leftmost).</summary>
    public IReadOnlyList<double> Gamma { get; }
    /// <summary>The off-locus control profile: same σ and same γ̄, two pair sums broken.</summary>
    public IReadOnlyList<double> OffLocusGamma { get; }

    /// <summary>How far the off-locus control walks the two pair sums against each other.
    /// </summary>
    public double OffLocusStep => OffLocusEps / (double)Grid;

    /// <summary>The frozen value, −4γ̄.</summary>
    public double Root => -4.0 * GammaBar;
    /// <summary>Its gamma-fold partner, r ↦ −r − 2σ, i.e. 4γ̄ − 2σ = (4 − 2N)γ̄.</summary>
    public double FoldRoot => 4.0 * GammaBar - 2.0 * Sigma;

    /// <summary>At N = 4 the fold sends the root to itself: (4 − 2N)γ̄ = −4γ̄ exactly when N = 4.
    /// There the four corners cannot be told apart by which root they carry, because there is only
    /// one root, and all four carry it.</summary>
    public bool RootsCoincide => N == 4;

    // ---- the room count (pure counting, no matrix) ----

    /// <summary>τQ's fixed cells among the off-diagonal coherences, the anti-diagonal (a, R(a)):
    /// 2⌊N/2⌋.</summary>
    public int FixedCells { get; }
    /// <summary>dim O₊ − dim O₋, which an involution forces to equal its fixed count.</summary>
    public int Surplus { get; }
    /// <summary>dim D₋, the diagonal cells' mirror pairs: EVEN under τQ, so they eat half.</summary>
    public int Tax { get; }
    /// <summary>Surplus − Tax = ⌊N/2⌋, the predicted frozen count.</summary>
    public int PredictedFrozen { get; }

    // ---- the exact reads ----

    /// <summary>dim ker(L₍₁,₁₎ + 4γ̄) at <see cref="J"/>, exact. The theorem on this stratum: ⌊N/2⌋.
    /// </summary>
    public int KernelDimension { get; }
    /// <summary>The same at <see cref="SecondJ"/>: the coupling cannot move the value.</summary>
    public int SecondCouplingKernelDimension { get; }
    /// <summary>The same at J = 0: the corollary, 2⌊N/2⌋ (twice, all semisimple; exactly ⌊N/2⌋
    /// modes depart as the coupling turns on, at order J^{2d_c} per pair).</summary>
    public int ZeroCouplingKernelDimension { get; }
    /// <summary>The same at a DIFFERENT on-locus profile with the same γ̄: the mode is blind to
    /// the individual rates, so this must again be ⌊N/2⌋.</summary>
    public int SecondProfileKernelDimension { get; }
    /// <summary>The same block read one grid step off the root: 0. The falsifier for the value.
    /// </summary>
    public int WrongRootKernelDimension { get; }
    /// <summary>The off-locus profile at the same root: 0. Partial balance yields nothing.
    /// </summary>
    public int OffLocusKernelDimension { get; }

    /// <summary>The four corner blocks p, q ∈ {1, N−1}: the fold count separating each from
    /// (1, 1), the kernel dimension at −4γ̄, and the kernel dimension at 4γ̄ − 2σ. Each corner
    /// carries ⌊N/2⌋ at exactly one of the two roots, chosen by fold parity (both chains).
    /// </summary>
    public IReadOnlyList<(int P, int Q, int Folds, int AtRoot, int AtFoldRoot)> Corners { get; }

    /// <summary>A non-corner sample block (2, 2) and its two reads; null when the block exceeds
    /// <see cref="ControlBlockSizeCap"/> or N is too small for (2, 2) to be off-corner. On the
    /// Heisenberg chain both reads are 0; on the XY chain they are not.</summary>
    public (int P, int Q, int Size, int AtRoot, int AtFoldRoot)? ControlBlock { get; }

    /// <summary>Every (p, q) block carrying either root, for N ≤ <see cref="FullCensusMaxN"/>;
    /// empty above it (see <see cref="CensusRan"/>).</summary>
    public IReadOnlyList<(int P, int Q, int AtRoot, int AtFoldRoot)> Census { get; }

    /// <summary>Whether the full (p, q) census was affordable at this N.</summary>
    public bool CensusRan { get; }

    /// <summary>Whether the census found the four corners and nothing else. TRUE on the Heisenberg
    /// chain, which is what the proof document verified; FALSE on the XY chain, where the quartic ZZ
    /// term that confines the divisor is missing (with it goes the ladder ρ ↦ Σ_l d†_l ρ d_l that
    /// carries the corner up the diagonal) and many more blocks carry the same root at the same
    /// multiplicity. Meaningless when <see cref="CensusRan"/> is false.</summary>
    public bool CornersAreTheOnlyCarriers { get; }

    // ---- the ladder (closed form, adopted from the proof's Section 8) ----

    /// <summary>Site distances d_c = N + 1 − 2c of the balanced pairs (c, R(c)), c = 1..⌊N/2⌋.
    /// </summary>
    public IReadOnlyList<int> PairDistances { get; }
    /// <summary>Departure orders J^{2d_c}, one per pair: the outermost pair holds longest.
    /// </summary>
    public IReadOnlyList<int> DepartureOrders { get; }
    /// <summary>Total valuation 2⌊N²/4⌋ = 2 Σ_c d_c.</summary>
    public int TotalValuation { get; }

    public FrozenDivisorWitness(
        int n = 6,
        double gammaBar = 0.09,
        double j = 0.75,
        DivisorChain chain = DivisorChain.Heisenberg)
    {
        if (n < MinN || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"N must satisfy {MinN} ≤ N ≤ {MaxN} (the live build materialises the 2^N Hamiltonian); got {n}");
        // γ̄ > 0 is this witness's scope, not the locus's. The locus itself is happy at γ̄ ≤ 0 (the
        // proof's Section 1 never uses positivity), and at γ̄ = 0 the corner kernel is N rather than
        // ⌊N/2⌋, and every read here is built around that one number, so it would mismatch. The
        // stratum belongs to the gate's G16 and to MirrorWorld's Divisor.Rooms(), which counts both.
        if (!(gammaBar > 0.0) || gammaBar > 1.0)
            throw new ArgumentOutOfRangeException(nameof(gammaBar), gammaBar,
                "γ̄ must lie in (0, 1]: it is the mean dephasing rate, and this witness reads the taxed stratum, " +
                "where the frozen count is ⌊N/2⌋. At γ̄ = 0 it is N instead; that stratum is pinned " +
                "by simulations/r90_frozen_divisor_gate.py (G16) and by MirrorWorld's Divisor, not here");
        if (!(j > 0.0) || j > 2.0)
            throw new ArgumentOutOfRangeException(nameof(j), j,
                "J must lie in (0, 2]: the J-independence face reads the same block at 4·J, which then stays " +
                "in a range where the scaled integers are comfortably exact");

        N = n;
        Chain = chain;
        int m = n / 2;

        // Snap onto the dyadic grid; everything downstream is exact.
        long gBarNum = SnapToGrid(gammaBar);
        long smallest = MaxHalfMagnitude(m);
        if (gBarNum <= smallest)
            throw new ArgumentOutOfRangeException(nameof(gammaBar), gammaBar,
                $"γ̄ snaps to {gBarNum}/{Grid}, at or below the largest half-profile step {smallest}/{Grid} this " +
                $"witness's profile uses at N = {n}; a site rate would come out non-positive, and positive " +
                $"rates are this witness's scope, not the locus's. Use " +
                $"γ̄ > {((smallest + 1) / (double)Grid).ToString("0.######", Inv)}.");
        long jNum = SnapToGrid(j);
        if (jNum <= 0)
            throw new ArgumentOutOfRangeException(nameof(j), j, $"J snaps to 0 on the 1/{Grid} grid");

        var gammaNum = Locus(n, gBarNum, HalfSeed);
        var secondGammaNum = Locus(n, gBarNum, SecondHalfSeed);

        // The off-locus control: +ε on site 0, −ε on site 1. σ and γ̄ are untouched, so the root
        // is the same number; only the pair balance breaks. That is the sharpest form of "partial
        // balance yields nothing".
        var offNum = (long[])gammaNum.Clone();
        offNum[0] += OffLocusEps;
        offNum[1] -= OffLocusEps;

        long sigmaNum = gammaNum.Sum();
        if (secondGammaNum.Sum() != sigmaNum || offNum.Sum() != sigmaNum)
            throw new InvalidOperationException(
                "the control profiles must share σ with the locus profile, or the roots are not comparable");

        J = jNum / (double)Grid;
        SecondJ = 4.0 * J;
        GammaBar = gBarNum / (double)Grid;
        Sigma = sigmaNum / (double)Grid;
        Gamma = gammaNum.Select(g => g / (double)Grid).ToArray();
        OffLocusGamma = offNum.Select(g => g / (double)Grid).ToArray();

        // The room count.
        FixedCells = 2 * m;
        Surplus = FixedCells;
        Tax = m;
        PredictedFrozen = Surplus - Tax;

        // The scale that turns every block entry into an integer. The Heisenberg chain's smallest
        // coefficient is J/4, so 4N·Grid·(entry) ∈ ℤ covers both chains; then
        // 4N·Grid·(−4γ̄) = −16·σNum and 4N·Grid·(4γ̄ − 2σ) = (16 − 8N)·σNum, and one grid step off
        // the root is 4N.
        long scale = 4L * n * Grid;
        long rootNum = -16 * sigmaNum;
        long foldRootNum = (16 - 8 * (long)n) * sigmaNum;
        long gridStep = 4L * n;

        KernelDimension = ExactKernelDim(CornerBlock(n, jNum, gammaNum, chain), scale, rootNum);
        SecondCouplingKernelDimension = ExactKernelDim(CornerBlock(n, 4 * jNum, gammaNum, chain), scale, rootNum);
        ZeroCouplingKernelDimension = ExactKernelDim(CornerBlock(n, 0, gammaNum, chain), scale, rootNum);
        SecondProfileKernelDimension = ExactKernelDim(CornerBlock(n, jNum, secondGammaNum, chain), scale, rootNum);
        WrongRootKernelDimension = ExactKernelDim(CornerBlock(n, jNum, gammaNum, chain), scale, rootNum + gridStep);
        OffLocusKernelDimension = ExactKernelDim(CornerBlock(n, jNum, offNum, chain), scale, rootNum);

        var h = Hamiltonian(n, J, chain);
        var gammaDouble = Gamma.ToArray();

        int Reads(int p, int q, long root) =>
            ExactKernelDim(PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaDouble, SectorFlat(n, p, q)), scale, root);

        // The four corners, each read at BOTH roots.
        var corners = new List<(int, int, int, int, int)>();
        foreach (int p in new[] { 1, n - 1 })
            foreach (int q in new[] { 1, n - 1 })
            {
                var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaDouble, SectorFlat(n, p, q));
                corners.Add((p, q,
                    (p == n - 1 ? 1 : 0) + (q == n - 1 ? 1 : 0),
                    ExactKernelDim(block, scale, rootNum),
                    ExactKernelDim(block, scale, foldRootNum)));
            }
        Corners = corners;

        // A non-corner sample: on the Heisenberg chain it carries nothing.
        if (n >= 4 && n - 1 != 2)
        {
            var controlFlat = SectorFlat(n, 2, 2);
            if (controlFlat.Length <= ControlBlockSizeCap)
            {
                var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaDouble, controlFlat);
                ControlBlock = (2, 2, controlFlat.Length,
                    ExactKernelDim(block, scale, rootNum),
                    ExactKernelDim(block, scale, foldRootNum));
            }
        }

        // The full (p, q) census, where it is affordable.
        var census = new List<(int, int, int, int)>();
        CensusRan = n <= FullCensusMaxN;
        if (CensusRan)
        {
            for (int p = 0; p <= n; p++)
                for (int q = 0; q <= n; q++)
                {
                    int atRoot = Reads(p, q, rootNum);
                    int atFold = Reads(p, q, foldRootNum);
                    if (atRoot > 0 || atFold > 0) census.Add((p, q, atRoot, atFold));
                }
        }
        Census = census;
        CornersAreTheOnlyCarriers =
            CensusRan && census.Count == 4 && census.All(c => IsCorner(c.Item1, c.Item2, n));

        // The ladder, closed form.
        var distances = new int[m];
        var orders = new int[m];
        for (int c = 1; c <= m; c++)
        {
            distances[c - 1] = n + 1 - 2 * c;
            orders[c - 1] = 2 * distances[c - 1];
        }
        PairDistances = distances;
        DepartureOrders = orders;
        TotalValuation = 2 * (n * n / 4);
    }

    /// <summary>True when every read meets what the proof claims FOR THIS CHAIN: the corner kernel
    /// is ⌊N/2⌋ at both nonzero couplings and at both profiles, the two falsifiers are zero, the
    /// J = 0 corollary doubles, and the four corners split by fold parity. The confinement to the
    /// corners is required only on the Heisenberg chain, which is where the proof verified it.
    /// </summary>
    public bool Holds =>
        KernelDimension == PredictedFrozen &&
        SecondCouplingKernelDimension == PredictedFrozen &&
        SecondProfileKernelDimension == PredictedFrozen &&
        ZeroCouplingKernelDimension == 2 * PredictedFrozen &&
        WrongRootKernelDimension == 0 &&
        OffLocusKernelDimension == 0 &&
        Corners.All(c => RootsCoincide
            ? c.AtRoot == PredictedFrozen && c.AtFoldRoot == PredictedFrozen
            : c.Folds % 2 == 0
                ? c.AtRoot == PredictedFrozen && c.AtFoldRoot == 0
                : c.AtRoot == 0 && c.AtFoldRoot == PredictedFrozen) &&
        (Chain != DivisorChain.Heisenberg ||
            ((ControlBlock is null || (ControlBlock.Value.AtRoot == 0 && ControlBlock.Value.AtFoldRoot == 0)) &&
             (!CensusRan || CornersAreTheOnlyCarriers)));

    // ================= the atoms =================

    private static bool IsCorner(int p, int q, int n) =>
        (p == 1 || p == n - 1) && (q == 1 || q == n - 1);

    /// <summary>γ_l = γ̄ + δ_l for l &lt; ⌊N/2⌋, γ_{N−1−l} = γ̄ − δ_l, middle = γ̄ at odd N: every
    /// reflection pair sums to 2γ̄ exactly, which is the R90 locus.</summary>
    private static long[] Locus(int n, long gBarNum, long[] halfSeed)
    {
        var g = new long[n];
        for (int i = 0; i < n; i++) g[i] = gBarNum;
        for (int i = 0; i < n / 2; i++)
        {
            long d = halfSeed[i % halfSeed.Length];
            g[i] = gBarNum + d;
            g[n - 1 - i] = gBarNum - d;
        }
        return g;
    }

    /// <summary>The largest half-profile step either seed asks for, so γ̄ can be required above it
    /// and every site rate stays positive.</summary>
    private static long MaxHalfMagnitude(int m)
    {
        long worst = 0;
        for (int i = 0; i < m; i++)
        {
            worst = Math.Max(worst, Math.Abs(HalfSeed[i % HalfSeed.Length]));
            worst = Math.Max(worst, Math.Abs(SecondHalfSeed[i % SecondHalfSeed.Length]));
        }
        return worst;
    }

    private static long SnapToGrid(double v) => (long)Math.Round(v * Grid);

    private static ComplexMatrix Hamiltonian(int n, double j, DivisorChain chain) =>
        (chain == DivisorChain.Heisenberg
            ? PauliHamiltonian.HeisenbergChain(n, j)
            : PauliHamiltonian.XYChain(n, j)).ToMatrix();

    /// <summary>The (1, 1) corner block of the Liouvillian at coupling jNum/Grid and the given
    /// integer rate profile, through the repo's committed builders.</summary>
    private static ComplexMatrix CornerBlock(int n, long jNum, long[] gammaNum, DivisorChain chain)
    {
        var h = Hamiltonian(n, jNum / (double)Grid, chain);
        var gamma = gammaNum.Select(g => g / (double)Grid).ToArray();
        return PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, SectorFlat(n, 1, 1));
    }

    /// <summary>Every Liouville-space flat index r·d + c of the joint-popcount sector
    /// (pRow, pCol).</summary>
    private static int[] SectorFlat(int n, int pRow, int pCol)
    {
        int d = 1 << n;
        var flat = new List<int>();
        for (int r = 0; r < d; r++)
        {
            if (BitOperations.PopCount((uint)r) != pRow) continue;
            for (int c = 0; c < d; c++)
                if (BitOperations.PopCount((uint)c) == pCol)
                    flat.Add(r * d + c);
        }
        return flat.ToArray();
    }

    /// <summary>dim ker(block − root·I), exact: scale·(block) − rootNum·I has Gaussian-integer
    /// entries, and the rank is taken over F_p at two primes with i a genuine square root of −1.
    /// </summary>
    private static int ExactKernelDim(ComplexMatrix block, long scale, long rootNum)
    {
        int d = block.RowCount;
        if (d == 0) return 0;
        var re = new long[d, d];
        var im = new long[d, d];
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++)
            {
                re[i, k] = ExactScale(block[i, k].Real, scale) - (i == k ? rootNum : 0L);
                im[i, k] = ExactScale(block[i, k].Imaginary, scale);
            }

        int rank = 0;
        foreach (long p in Primes) rank = Math.Max(rank, RankModP(re, im, d, p));
        return d - rank;
    }

    /// <summary>scale · v as an exact integer. Every input is dyadic with denominator dividing
    /// 4·<see cref="Grid"/>, so the product is integral; anything else means a value slipped off
    /// the grid and the exact path would silently round.</summary>
    private static long ExactScale(double v, long scale)
    {
        double s = v * scale;
        double r = Math.Round(s);
        if (s != r)
            throw new InvalidOperationException(
                $"block entry {v.ToString("R", Inv)} is not an exact multiple of 1/{scale}: the dyadic " +
                "grid was broken and the exact rank would round silently");
        if (Math.Abs(r) > 9.0e15)
            throw new InvalidOperationException($"scaled entry {r} exceeds the exact-integer range of double");
        return (long)r;
    }

    /// <summary>Rank over F_p with i realised as a square root of −1 (p ≡ 1 mod 4): a + b·i ↦
    /// a + b·r.</summary>
    private static int RankModP(long[,] re, long[,] im, int d, long p)
    {
        long r = SqrtMinusOne(p);
        var a = new long[d, d];
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++)
                a[i, k] = ((re[i, k] % p + r % p * (im[i, k] % p)) % p + p) % p;

        int rank = 0;
        for (int col = 0; col < d && rank < d; col++)
        {
            int piv = -1;
            for (int row = rank; row < d; row++)
                if (a[row, col] != 0) { piv = row; break; }
            if (piv < 0) continue;
            for (int k = 0; k < d; k++) (a[rank, k], a[piv, k]) = (a[piv, k], a[rank, k]);
            long inv = ModPow(a[rank, col], p - 2, p);
            for (int k = col; k < d; k++) a[rank, k] = a[rank, k] * inv % p;
            for (int row = rank + 1; row < d; row++)
            {
                if (a[row, col] == 0) continue;
                long f = a[row, col];
                for (int k = col; k < d; k++)
                    a[row, k] = ((a[row, k] - f * a[rank, k]) % p + p) % p;
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
        throw new InvalidOperationException($"no square root of −1 mod {p} (is p ≡ 1 mod 4?)");
    }

    private static long ModPow(long b, long e, long p)
    {
        long acc = 1;
        b = ((b % p) + p) % p;
        while (e > 0)
        {
            if ((e & 1) == 1) acc = acc * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return acc;
    }

    // ================= the reading =================

    private string Fmt(double v) => v.ToString("0.######", Inv);

    private string GammaList(IReadOnlyList<double> g) => "[" + string.Join(", ", g.Select(Fmt)) + "]";

    private string ChainName => Chain == DivisorChain.Heisenberg ? "Heisenberg" : "XY";

    public string DisplayName =>
        $"FrozenDivisorWitness (F140: λ = −4γ̄ frozen at least ⌊N/2⌋ times at every J, live; N={N}, "
      + $"{ChainName})";

    public string Summary =>
        $"N={N} on the {ChainName} chain, R90 locus γ={GammaList(Gamma)} (γ̄={Fmt(GammaBar)}, σ={Fmt(Sigma)}): " +
        $"the corner block holds λ = {Fmt(Root)} = −4γ̄ with exact kernel dimension {KernelDimension} at " +
        $"J={Fmt(J)} and {SecondCouplingKernelDimension} at J={Fmt(SecondJ)}, predicted ⌊N/2⌋ = " +
        $"{PredictedFrozen} from the room count 2⌊N/2⌋ − ⌊N/2⌋ = {Surplus} − {Tax}" +
        $"{(Holds ? "" : ", MISMATCH")}; a second locus profile with the same γ̄ gives " +
        $"{SecondProfileKernelDimension} (blind to the individual rates), J = 0 gives " +
        $"{ZeroCouplingKernelDimension} = 2⌊N/2⌋ (the corollary), and the falsifiers are " +
        $"{WrongRootKernelDimension}/{OffLocusKernelDimension} (root one grid step off / off the locus at " +
        $"unchanged σ). " +
        (CensusRan
            ? $"The full (p,q) census finds {Census.Count} carrying block(s), " +
              (CornersAreTheOnlyCarriers
                  ? "exactly the four corners."
                  : "MORE than the four corners: the confinement is a Heisenberg statement, and this chain " +
                    "does not obey it.")
            : "The full census is skipped at this N.") +
        $" {(Holds ? "Every read meets the theorem." : "At least one read does NOT meet the theorem.")}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: $"the room shortage: {Surplus} − {Tax} = {PredictedFrozen} = ⌊N/2⌋",
                summary: $"the cell mirror τQ : (a,b) ↦ (R(b), R(a)) fixes exactly the {FixedCells} = 2⌊N/2⌋ " +
                         "anti-diagonal coherences (a, R(a)), one per site, and an involution's even rooms " +
                         $"outnumber its odd rooms by exactly its fixed count, so dim O₊ − dim O₋ = {Surplus}. " +
                         $"The diagonal cells' mirror pairs (dim D₋ = {Tax}) are EVEN under τQ, not odd, and " +
                         "eat half. An operator odd under that mirror must send the bigger half into the " +
                         "smaller; what will not fit freezes. Pure counting: no matrix, no coupling, no γ.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"the exact multiplicity: dim ker(L₍₁,₁₎ + 4γ̄) = {KernelDimension} at both couplings",
                summary: $"J = {Fmt(J)} → {KernelDimension}, J = {Fmt(SecondJ)} → {SecondCouplingKernelDimension}, " +
                         $"and a DIFFERENT on-locus profile at the same γ̄ → {SecondProfileKernelDimension}: the " +
                         "value hears the mean dephasing rate and nothing else, not the coupling and not the " +
                         "individual rates. Read by elimination over GF(p) at two primes ≡ 1 (mod 4), never by " +
                         "an eigensolver: the departing modes sit at spacing J^{2d} from the root, so a " +
                         "floating-point count would miscount exactly where the theorem is sharpest.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"the falsifiers: {WrongRootKernelDimension} one grid step off the root, " +
                             $"{OffLocusKernelDimension} off the locus",
                summary: $"λ = {Fmt(Root)} + 1/{Grid} carries {WrongRootKernelDimension}: the value is the value, " +
                         "not a neighbourhood. Off the locus at UNCHANGED σ and γ̄ (two pair sums moved by " +
                         $"±{Fmt(OffLocusStep)} against each other, so the root is literally the same number) → " +
                         $"{OffLocusKernelDimension}: partial balance yields nothing, so the mode does not fade " +
                         "in, it snaps in. Together these two say the reads above are about this locus and this " +
                         "value, not about the rank routine.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the four corners, and which root each one carries",
                summary: string.Join("; ", Corners.Select(c =>
                             $"({c.P},{c.Q}) folds={c.Folds}: {c.AtRoot} at −4γ̄, {c.AtFoldRoot} at 4γ̄−2σ")) +
                         $". Each one-sided gamma fold sends a rate r ↦ −r − 2σ and flips one popcount index " +
                         $"to N − it, so an even number of folds returns to {Fmt(Root)} and an odd number lands " +
                         $"on {Fmt(FoldRoot)}" +
                         (RootsCoincide
                             ? ", except at N = 4, where (4 − 2N)γ̄ = −4γ̄ and the fold sends the root to " +
                               "ITSELF, so all four corners carry the one root and the parity is invisible here."
                             : ".") +
                         " The parity law survives the change of chain; what does not is the confinement below.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: CensusRan
                    ? $"the full (p,q) census on the {ChainName} chain: {Census.Count} carrying block(s)"
                    : $"the non-corner sample (the full census is skipped above N = {FullCensusMaxN})",
                summary: (CensusRan
                             ? string.Join("; ", Census.Select(c =>
                                   $"({c.P},{c.Q}){(IsCorner(c.P, c.Q, N) ? "" : "*")}: {c.AtRoot}/{c.AtFoldRoot}")) +
                               " (starred = NOT a corner). "
                             : "") +
                         (ControlBlock is { } cb
                             ? $"The non-corner block ({cb.P},{cb.Q}), size {cb.Size}, carries {cb.AtRoot} at " +
                               $"−4γ̄ and {cb.AtFoldRoot} at 4γ̄ − 2σ. "
                             : "") +
                         (Chain == DivisorChain.Heisenberg
                             ? "On the Heisenberg chain only the four corners carry, which is the proof " +
                               "document's census: what confines the divisor is ZZ as a QUARTIC term, not " +
                               "the diagonal it puts on h."
                             : "On the XY chain the confinement is GONE. The same root is carried at the same " +
                               "multiplicity ⌊N/2⌋ by blocks that are not corners" +
                               (CensusRan
                                   ? $" ({Census.Count} of the {(N + 1) * (N + 1)} blocks here, against the " +
                                     "four on Heisenberg), never partially, and the ones carrying the " +
                                     "UNFOLDED root −4γ̄ all have p + q even"
                                   : "") +
                               ". The bound and the fold parity survive the change of chain; \"only the " +
                               "corners\" is a Heisenberg statement, and this witness is where that was " +
                               "found."),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"the J = 0 corollary: {ZeroCouplingKernelDimension} = 2⌊N/2⌋, and the ladder that " +
                             "empties it",
                summary: $"at zero coupling the root carries {ZeroCouplingKernelDimension} = twice the " +
                         $"multiplicity, all semisimple; exactly {PredictedFrozen} modes depart as the coupling " +
                         "turns on, one per pair, and a mode cannot move until the coupling has walked the " +
                         "excitation from one site of its pair to the other: " +
                         string.Join(", ", PairDistances.Select((d, i) =>
                             $"pair {i + 1} at distance {d} departs at order J^{DepartureOrders[i]}")) +
                         $", total valuation 2⌊N²/4⌋ = {TotalValuation}. Distance buys immunity: the two ends " +
                         "of the chain are the last to let go. The orders are the proof's Section 8, adopted " +
                         "as closed form; what is recomputed here is the count they start from.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the remaining open ink (scope of this witness)",
                summary: "this witness types the divisor BOUND (≥ ⌊N/2⌋ at every coupling) and the " +
                         "fold-parity placement of the corners; the bound holds at every γ̄, " +
                         "and every read here is taken at γ̄ > 0. Tightness " +
                         "(exactly ⌊N/2⌋ for all but finitely many J, and only for γ̄ ≠ 0) " +
                         "is a theorem too, by the proof's Sections 6 and 7; the ranks " +
                         "above meet it at every coupling sampled, but sampling is not what proves it. " +
                         "What this witness does not REACH: γ̄ = 0, where the even defect is absent, the " +
                         "whole block is τQ-odd and the count is N. The constructor closes that door; " +
                         "MirrorWorld's Divisor.Rooms() and the gate's G16 carry the stratum. " +
                         "What no rank read here can see: at the real exceptional couplings the root goes " +
                         "DEFECTIVE (one 2×2 Jordan block) while its kernel dimension does not move, so " +
                         "this witness cannot tell that failure from a healthy coupling. Still open in the " +
                         "proof, and about the ladder rather than the multiplicity: the upper half of the " +
                         "valuation law, the J → 0 order being exactly 2⌊N²/4⌋ rather than at least it.",
                provenance: NodeProvenance.Stored);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
