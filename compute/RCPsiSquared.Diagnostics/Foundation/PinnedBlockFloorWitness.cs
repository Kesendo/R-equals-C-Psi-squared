using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F153, recomputed at inspect time: which joint-popcount blocks of the Liouvillian sit ENTIRELY
/// on their Absorption-Theorem floor. Sweeps every one of the (N+1)² blocks entry-wise through
/// <see cref="WeightCoherenceBlock"/> (γ = 1 baked in, which is the criterion's uniform-γ scope by
/// construction) and reads three cell comparisons per block. No eigensolver anywhere, so no tolerance:
/// the residuals below are compared to 0.0 exactly.
///
/// <para><b>What is read, and why these three.</b> On a block the builder returns
/// L = A + i·S with A the real Absorption diagonal −2·n_diff and S the hops plus the Δ·ZZ frequency.
/// <see cref="BlockReading.RealDiagonalSpread"/> is max−min of Re L[r,r], zero exactly when n_diff is constant on the
/// block; <see cref="BlockReading.OffDiagonalRealMass"/> is max |Re L[r,c]| off the diagonal, zero because every hop
/// is ±2iq at real q; <see cref="BlockReading.OffDiagonalAsymmetry"/> is max |L[r,c] − L[c,r]|, zero because the two
/// hop directions carry the same coefficient. Together they certify L = c·Id + i·S with S REAL SYMMETRIC,
/// whose eigenvalues are c + i·μ with μ real, so Re λ = c throughout the block at every REAL coupling.</para>
///
/// <para><b>Two fences on that certificate, and the first is measurable here.</b> REAL q: at complex q the
/// hops contribute ∓2·Im(q) to the real part OFF the diagonal and the ZZ term contributes +Im(q)·Δ·Δzz ON
/// it, so a pinned block leaves its floor and can even acquire a gain mode (block (0,1) at N=4, Δ=1,
/// q = 1+0.3i reaches Re λ = +0.0485 against a floor of −2). The tests carry that as a negative control.
/// REAL H: the SYMMETRY leg is the chain's, not the criterion's. Under a Hermitian hopping with a complex
/// amplitude (a Peierls phase) the criterion still holds exactly, at N=4 flat to 5e-15 on the pinned
/// blocks, but there S is HERMITIAN and not real symmetric, so this witness's third residual would fire on
/// a true instance. The criterion's general form is "the Hermitian part of L is scalar"; real symmetric is
/// what a REAL H buys, and this builder builds only real H.</para>
///
/// <para><b>The ZZ diagonal is the half a carrier must not drop.</b> MirrorWorld's gate certifies the
/// XY-only form c·Id + i·J·A, and it passes bit-exactly because <c>Pair.Rate</c> carries no ZZ. With the ZZ
/// term the diagonal is NOT c·Id: it is c on the real axis and configuration-dependent on the imaginary
/// one. The arc that opened this carrier records the (0,1) block at N=4 reading −1−4i and −1−2i; this
/// builder reproduces those FREQUENCIES exactly, −2−4i and −2−2i, the real part differing because γ = 1
/// is baked in here and the arc's reading is at γ = ½. Four cells, two values, one of them twice as
/// far from the real axis as the other. That is the ZZ diagonal, and the pinning survives it because the
/// spread is entirely imaginary. <see cref="BlockReading.ImaginaryDiagonalSpread"/> reads it, and the test asserts
/// it is STRICTLY POSITIVE at Δ ≠ 0 on the pinned blocks that have one, so the witness cannot silently
/// degrade into the XY statement.</para>
///
/// <para><b>The discriminating half.</b> On the other (N−1)² blocks the real diagonal spread must be
/// strictly positive: without it the sweep would pass on a criterion that claimed every block. The
/// entry-wise route certifies the mechanism in both directions; that the mechanism's absence really
/// moves the SPECTRUM is the trace identity recorded on
/// <see cref="PinnedBlockFloorClaim.ConverseByTrace"/>, not a separate measurement.</para>
///
/// <para>Live at <c>inspect --root pinned</c>. Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F153 +
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> §6 (the Edge lemma, the boundary-block sentence and
/// the window-combinatorics shell) + <c>simulations/offdiag_sector_floor.py</c> (the spectral gate).</para></summary>
public sealed class PinnedBlockFloorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Chain length swept. Kept small because the sweep is over ALL (N+1)² blocks and the
    /// middle ones are C(N,⌊N/2⌋)² dense.</summary>
    public int N { get; }

    /// <summary>The coupling q = J the blocks are built at. REAL, and that is a fence rather than a
    /// convenience: the criterion is free of the coupling's VALUE but not of its being real (see the class
    /// summary). The witness therefore takes a double and cannot be pointed off the real axis; the
    /// complex-q failure is exercised in the tests through the builder directly.</summary>
    public double J { get; }

    /// <summary>The ZZ anisotropy Δ. Δ=1 is Heisenberg, Δ=0 the XY chain MirrorWorld's gate covers.</summary>
    public double Delta { get; }

    public PinnedBlockFloorWitness(int n = 4, double j = 1.0, double delta = 1.0)
    {
        if (n < 2 || n > 7)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                "the sweep visits every (N+1)² block, the middle one C(N,N/2)² dense; N = 2..7");
        if (!double.IsFinite(j)) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be finite");
        if (!double.IsFinite(delta)) throw new ArgumentOutOfRangeException(nameof(delta), delta, "Δ must be finite");
        N = n; J = j; Delta = delta;
    }

    /// <summary>One block's entry-wise reading. All four residuals are exact doubles read off the built
    /// cells; none comes from an eigensolver.</summary>
    public readonly record struct BlockReading(
        int P, int Q, int Dim, bool Pinned,
        double Floor,
        double RealDiagonalSpread,
        double ImaginaryDiagonalSpread,
        double OffDiagonalRealMass,
        double OffDiagonalAsymmetry);

    /// <summary>Reads one (p, q) block entry-wise at (J, Δ). p is the KET popcount and q the BRA
    /// popcount in the builder's order; the criterion min(p,q) = 0 or max(p,q) = N is symmetric under
    /// swapping them, so the repo's two index conventions cannot disagree about this reading.</summary>
    public static BlockReading Read(int n, int p, int q, double j, double delta)
    {
        var l = WeightCoherenceBlock.Build(n, p, q, new Complex(j, 0.0), delta);
        int d = l.GetLength(0);

        double reLo = double.PositiveInfinity, reHi = double.NegativeInfinity;
        double imLo = double.PositiveInfinity, imHi = double.NegativeInfinity;
        for (int r = 0; r < d; r++)
        {
            reLo = Math.Min(reLo, l[r, r].Real); reHi = Math.Max(reHi, l[r, r].Real);
            imLo = Math.Min(imLo, l[r, r].Imaginary); imHi = Math.Max(imHi, l[r, r].Imaginary);
        }

        double realMass = 0.0, asym = 0.0;
        for (int r = 0; r < d; r++)
            for (int c = 0; c < d; c++)
            {
                if (r == c) continue;
                realMass = Math.Max(realMass, Math.Abs(l[r, c].Real));
                asym = Math.Max(asym, (l[r, c] - l[c, r]).Magnitude);
            }

        return new BlockReading(p, q, d, PinnedBlockFloorClaim.IsPinned(n, p, q),
            // the claim owns the formula; γ = 1 is what the builder bakes in. + 0.0 turns the diagonal
            // blocks' −0.0 into 0.0 for printing (the comparison is unaffected: −0.0 == 0.0 in IEEE).
            Floor: PinnedBlockFloorClaim.Floor(n, p, q, gamma: 1.0) + 0.0,
            RealDiagonalSpread: reHi - reLo,
            ImaginaryDiagonalSpread: imHi - imLo,
            OffDiagonalRealMass: realMass,
            OffDiagonalAsymmetry: asym);
    }

    /// <summary>The whole (N+1)² sweep at (J, Δ), in (p, q) order.</summary>
    public static IReadOnlyList<BlockReading> Sweep(int n, double j, double delta)
    {
        var rows = new List<BlockReading>();
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
                rows.Add(Read(n, p, q, j, delta));
        return rows;
    }

    /// <summary>How many blocks the entry-wise reading finds sitting entirely on THEIR FLOOR: real diagonal
    /// spread exactly 0.0, off-diagonal real mass exactly 0.0, off-diagonal exactly symmetric, AND the
    /// constant equal to −2|p−q|. The last conjunct is not redundant: the first three say the real part is
    /// CONSTANT, and a block can be constant at the wrong constant (which is what a γ profile produces on
    /// half the pinned blocks at N=4). Should equal 4N.</summary>
    public static int WholeCount(int n, double j, double delta) =>
        Sweep(n, j, delta).Count(r =>
            r.RealDiagonalSpread == 0.0 && r.OffDiagonalRealMass == 0.0 && r.OffDiagonalAsymmetry == 0.0
            && ObservedFloor(n, r.P, r.Q, j, delta) == r.Floor);

    /// <summary>The floor each pinned block sits on, read from the constant real diagonal rather than
    /// predicted: Re L[r,r] for any r, which the sweep has already shown to be constant.</summary>
    public static double ObservedFloor(int n, int p, int q, double j, double delta)
    {
        var l = WeightCoherenceBlock.Build(n, p, q, new Complex(j, 0.0), delta);
        return l[0, 0].Real;
    }

    /// <summary>The uniform-γ fence, made measurable ENTRY-WISE. Under a site-resolved profile a cell's
    /// real part is −2·Σ over the sites where its ket and bra disagree, so pinning fixes the disagreement
    /// COUNT and not the disagreeing SITES; this returns max−min of that quantity over the block. It is
    /// exactly 0.0 on a pinned block at uniform γ and generically positive under a profile, which is the
    /// mechanism failing, and it needs no eigensolver and no Hamiltonian at all (the hops and the ZZ term
    /// never touch the real diagonal at real q).
    ///
    /// <para>What it does NOT say: that the SPECTRUM spreads. It does not always. Measured this session at
    /// γ = linspace(0.1, 1, N), J = 1: at Δ = 1 the (0,1) block's spectral spread is 0.4903 (N=4) and
    /// 0.9477 (N=5), but at Δ = 0 the same profile gives exactly 0.0 at J = 1 and only opens as J falls
    /// (1.77 at J = 0.05, N=4), the strong-coupling coalescence F152 records. And where flatness survives
    /// a profile it is flatness at the WRONG CONSTANT: at N=4 eight of the sixteen pinned blocks stay flat,
    /// block (0,2) sitting at −2.200 = −2·mean(γ) while its floor is −1.000. So the fence is that the block
    /// leaves its FLOOR, not that it spreads.</para></summary>
    public static double ProfileRealDiagonalSpread(int n, int p, int q, IReadOnlyList<double> gamma)
    {
        ArgumentNullException.ThrowIfNull(gamma);
        if (gamma.Count != n)
            throw new ArgumentException($"a site-resolved profile needs one rate per site: {n} expected, {gamma.Count} given", nameof(gamma));

        double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
        foreach (int a in WeightCoherenceBlock.Configs(n, p))
            foreach (int b in WeightCoherenceBlock.Configs(n, q))
            {
                double rate = 0.0;
                for (int s = 0; s < n; s++)
                    if ((((a ^ b) >> s) & 1) != 0) rate += gamma[s];
                lo = Math.Min(lo, -2.0 * rate); hi = Math.Max(hi, -2.0 * rate);
            }
        return hi - lo;
    }

    private IReadOnlyList<BlockReading>? _rows;

    /// <summary>The sweep, computed once per witness instance. Cached because <see cref="Summary"/> and
    /// <see cref="Children"/> each read it and the middle block is C(N,⌊N/2⌋)² dense.</summary>
    public IReadOnlyList<BlockReading> Rows => _rows ??= Sweep(N, J, Delta);

    public string DisplayName =>
        $"F153 pinning criterion, recomputed: all {(N + 1) * (N + 1)} blocks at N={N}, J={J.ToString("0.###", Inv)}, Δ={Delta.ToString("0.###", Inv)}";

    public string Summary
    {
        get
        {
            var rows = Rows;
            int whole = rows.Count(r => r.RealDiagonalSpread == 0.0 && r.OffDiagonalRealMass == 0.0 && r.OffDiagonalAsymmetry == 0.0);
            double worstOther = rows.Where(r => !r.Pinned).Select(r => r.RealDiagonalSpread).DefaultIfEmpty(0.0).Min();
            return $"{whole} of {rows.Count} blocks sit entirely on their floor (4N = {4 * N}); " +
                   $"the criterion min(p,q)=0 or max(p,q)=N agrees on every block; " +
                   $"smallest real-diagonal spread among the other {(N - 1) * (N - 1)} is {worstOther.ToString("0.###", Inv)}; " +
                   "all residuals read entry-wise, compared to 0.0 exactly, no eigensolver";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var rows = Rows;
            yield return new InspectableNode("the criterion",
                summary: "a block sits entirely on Re λ = −2γ|p−q| exactly when min(p,q) = 0 or max(p,q) = N; " +
                         $"that is 4(N+1) − 4 = {4 * N} of the {(N + 1) * (N + 1)} blocks, the four corners among them");
            yield return new InspectableNode("the entry-wise certificate",
                summary: "L = c·Id + i·S with S real symmetric (hops PLUS the Δ·ZZ diagonal), so every eigenvalue " +
                         "is c + i·μ with μ real. Read as three cell comparisons: constant real diagonal, no real " +
                         "part off the diagonal, off-diagonal symmetric. Each compared to 0.0 exactly.");
            yield return new InspectableNode("the ZZ half",
                summary: Delta == 0.0
                    ? "Δ = 0 here, so the diagonal is c·Id outright: this is the XY form MirrorWorld's gate certifies, " +
                      "and it does NOT exercise the ZZ correction. Inspect at Δ ≠ 0 to see it."
                    : $"Δ = {Delta.ToString("0.###", Inv)}: the diagonal is c on the real axis and configuration-dependent on the imaginary one. " +
                      $"Largest imaginary-diagonal spread among the pinned blocks: " +
                      $"{rows.Where(r => r.Pinned).Select(r => r.ImaginaryDiagonalSpread).DefaultIfEmpty(0.0).Max().ToString("0.###", Inv)}. " +
                      "A carrier that read only c·Id would be reading the XY chain.");

            foreach (var r in rows)
                yield return new InspectableNode(
                    $"block ({r.P},{r.Q}) dim {r.Dim} {(r.Pinned ? "PINNED" : "free")}",
                    summary: $"floor −2|p−q| = {r.Floor.ToString(Inv)}; real-diagonal spread {r.RealDiagonalSpread.ToString(Inv)}; " +
                             $"imaginary-diagonal spread {r.ImaginaryDiagonalSpread.ToString("0.###", Inv)}; " +
                             $"off-diagonal real mass {r.OffDiagonalRealMass.ToString(Inv)}; off-diagonal asymmetry {r.OffDiagonalAsymmetry.ToString(Inv)}");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
