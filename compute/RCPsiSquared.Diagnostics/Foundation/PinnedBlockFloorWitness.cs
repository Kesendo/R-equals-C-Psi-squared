using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F153, recomputed at inspect time: which joint-popcount blocks of the Liouvillian sit ENTIRELY
/// on their Absorption-Theorem floor. Sweeps every one of the (N+1)² blocks entry-wise through
/// <see cref="WeightCoherenceBlock"/> (γ = 1 baked in, which is the criterion's uniform-γ scope by
/// construction) and reads three cell comparisons per block. No eigensolver in THAT sweep, so no tolerance
/// there: its residuals are compared to 0.0 exactly. One member added later, <see cref="ProfileReSpan"/>,
/// does call one: the FENCE's verdict never depends on it, but the flat-set law around the fence is an
/// eigensolver reading throughout. See the γ-fence paragraph below.
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
/// <para><b>The γ fence, and it is exact.</b> The members from <see cref="CellRate"/> down carry the
/// criterion under a site-resolved γ, where the block leaves its floor. That half was long thought to
/// need a spectrum; it does not. Under Z-dephasing, Herm(L_block) = −2·diag(rate) for any Hermitian H,
/// any graph, any Δ and any profile, and it is DIAGONAL, so the real parts follow from it and the trace
/// gives the converse; the SPECTRAL reading needs the number-conserving H the criterion already assumes,
/// since only then is the index set an invariant block. The blocks there come from
/// <see cref="PerBlockLiouvillianBuilder"/> rather than from <see cref="WeightCoherenceBlock"/>,
/// deliberately, though only half of that buys independence: off the diagonal the comparison tests that
/// the Hamiltonian cancels in the Hermitian part, while on the diagonal both sides implement the same
/// sum in opposite orders and the check is a parity one. Live under the node "the γ fence, spectral
/// half".</para>
///
/// <para>Live at <c>inspect --root pinned</c>. Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F153 +
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> §6 (the Edge lemma, the boundary-block sentence and
/// the window-combinatorics shell) + <c>simulations/offdiag_sector_floor.py</c> (the spectral gate at
/// UNIFORM γ; its γ is a scalar, which is why the profile side is gated here instead).</para></summary>
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
    /// four of the sixteen pinned blocks at N=4, those with free index N/2). Should equal 4N.</summary>
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
    /// <para>What it does NOT say: that the SPECTRUM spreads. It does not always. In the PAULI book at
    /// γ = linspace(0.1, 1, N), J = 1: at Δ = 1 the (0,1) block's spectral spread is 0.4903 (N=4) and
    /// 0.9477 (N=5), but at Δ = 0 the same profile gives MACHINE ZERO at J = 1 and only opens as J falls
    /// (1.77 at J = 0.05, N=4), which F152 records as a strong-coupling coalescence, though the ramp is on the
    /// R₉₀ locus and an OFF-locus profile at the same Δ does not coalesce at any coupling measured. And where flatness survives
    /// a profile it is not always flatness at the WRONG CONSTANT: at N=4, J=1 eight pinned blocks stay
    /// flat, but only the four with free index N/2 are at the wrong constant, block (0,2) sitting at
    /// −2.200 while its floor is −1.000; the four one-cell blocks are flat AT their floor. On a PINNED block that constant is the
    /// TRACE's, −2·|p−q|·γ̄, and not −2·mean(γ) = −1.100; the two coincide only at |p−q| = 1. Which
    /// blocks those are obeys a law that is NOT a |p−q| parity, spelled out on
    /// <see cref="ProfileReSpan"/>. So the fence is that the block leaves its FLOOR, not that it
    /// spreads.</para>
    ///
    /// <para>THE SPECTRAL CONSEQUENCE IS NOW EXACT TOO, and this member is no longer the only entry-wise
    /// route to it: <see cref="ProfileHermitianPartResidual"/> shows the Hermitian part of ANY block is
    /// −2·diag(rate) and DIAGONAL, which decides every real part without an eigensolver, the trace
    /// supplying the converse. What that turns this member into is the fence's exact predicate rather
    /// than a proxy for it.</para></summary>
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

    // ---- the SPECTRAL half of the γ fence, and it turns out to need almost no eigensolver ----

    /// <summary>γ INDEX CONVENTION for everything below, and it is NOT the one
    /// <see cref="ProfileRealDiagonalSpread"/> uses: <paramref name="gamma"/>[l] is SITE l, site 0 the
    /// leftmost Kronecker factor, hence bit n−1−l. That is the repo's canonical convention, the one
    /// <see cref="PerBlockLiouvillianBuilder"/> documents and gates. The older member above reads bit l
    /// for site l instead, and the disagreement is UNOBSERVABLE there rather than a bug: it returns a
    /// max−min over the block's cells, and the cell set of a (p, q) block is closed under every
    /// permutation of the bits, so relabelling the sites permutes the rates without moving their spread.
    /// It becomes observable the moment a reading is entry-wise, which is what these members are, so
    /// they take the canonical one and this paragraph is the fence.</summary>
    public static double CellRate(int n, int ket, int bra, IReadOnlyList<double> gamma)
    {
        ArgumentNullException.ThrowIfNull(gamma);
        if (gamma.Count != n)
            throw new ArgumentException($"a site-resolved profile needs one rate per site: {n} expected, {gamma.Count} given", nameof(gamma));
        int diff = ket ^ bra;
        double rate = 0.0;
        for (int l = 0; l < n; l++)
            if (((diff >> (n - 1 - l)) & 1) != 0) rate += gamma[l];
        return rate;
    }

    /// <summary>The block's cells as (ket, bra) config pairs, kets outer and bras inner, both in
    /// ascending-mask order: <see cref="WeightCoherenceBlock"/>'s basis order, so a reading here can be
    /// laid beside that builder's cell for cell.</summary>
    public static IReadOnlyList<(int Ket, int Bra)> Cells(int n, int p, int q)
    {
        var cells = new List<(int, int)>();
        foreach (int ket in WeightCoherenceBlock.Configs(n, p))
            foreach (int bra in WeightCoherenceBlock.Configs(n, q))
                cells.Add((ket, bra));
        return cells;
    }

    /// <summary>The same cells as Liouville-space flat indices <c>ket·2^n + bra</c>, which is what
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> takes.</summary>
    public static IReadOnlyList<int> FlatIndices(int n, int p, int q)
    {
        int d = 1 << n;
        return Cells(n, p, q).Select(c => c.Ket * d + c.Bra).ToArray();
    }

    /// <summary>H = J·Σ_bond (X·X + Y·Y + Δ·Z·Z) on the open chain: the PAULI book, which is the book
    /// <see cref="WeightCoherenceBlock"/> builds in and therefore the book F153's measured spectral
    /// numbers are quoted in. It is NOT the spin book of <c>BlockSpectrumWitness</c>, whose J is four
    /// times this one; the entry's numbers do not reproduce in that book, which is how the collision was
    /// found.</summary>
    public static ComplexMatrix ChainHamiltonian(int n, double j, double delta)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < n - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.X, b + 1, PauliLetter.X, j));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Y, b + 1, PauliLetter.Y, j));
            terms.Add(PauliTerm.TwoSite(n, b, PauliLetter.Z, b + 1, PauliLetter.Z, j * delta));
        }
        return new PauliHamiltonian(n, terms).ToMatrix();
    }

    /// <summary>The (p, q) block under a site-resolved γ, built by the GENERAL Liouvillian block builder
    /// from a Hilbert-space H rather than by the combinatorial block builder this witness otherwise uses.
    /// The second route is the point: every residual below compares this construction against a
    /// prediction assembled from the cells alone, so no gate here compares a code path with itself.</summary>
    public static ComplexMatrix ProfileBlock(int n, int p, int q, double j, double delta, IReadOnlyList<double> gamma) =>
        ProfileBlock(n, p, q, ChainHamiltonian(n, j, delta), gamma);

    /// <summary>The same block from an ALREADY BUILT H. Sweeping every block otherwise rebuilds the full
    /// 2^N × 2^N Hamiltonian once per block, which is why this overload exists; the block builds themselves
    /// dominate at larger N, and the fence node hoists those separately.</summary>
    public static ComplexMatrix ProfileBlock(int n, int p, int q, ComplexMatrix h, IReadOnlyList<double> gamma) =>
        PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, FlatIndices(n, p, q));

    /// <summary>THE MASTER INVARIANT, and it is what makes the spectral half exact: for ANY Hermitian H,
    /// any graph, any Δ, any γ profile and any block at all,
    /// <para>Herm(L_block) = (L + L†)/2 = −2·diag(rate), rate(a,b) = Σ over the sites where a and b disagree.</para>
    /// The dissipator is real and diagonal in the coherence basis, and −i·ad_H is anti-Hermitian as a
    /// superoperator for Hermitian H, so it contributes nothing. This returns the worst entry-wise
    /// deviation, and it is 0.0 EXACTLY, not to a tolerance.
    ///
    /// <para>NOT NEW, and the entry says which parts are whose: PROOF_CODIM1_BY_ADDITIVITY §6 owns the
    /// rate window Re λ = v†Av/v†v, "exact for every eigenvalue, defective or not", with the Edge lemma
    /// as its zero-width case; PROOF_ABSORPTION_THEOREM owns the anti-Hermiticity for every Hermitian H;
    /// and <c>BlockSpectrumWitness.HermitianPartResidual</c> already gates this identity under a profile
    /// for the (0,1) block. What is added here is the step to a GENERAL block and from a disagreement
    /// COUNT to a site-resolved sum, plus the trace converse that turns the window's zero-width case
    /// into an equivalence.</para>
    ///
    /// <para>Four steps, the last of which is what makes it survive restriction: the inner product is
    /// Hilbert-Schmidt; ad_H is HS-self-adjoint because H is Hermitian, so −i·ad_H is anti-Hermitian;
    /// {|a⟩⟨b|} is HS-orthonormal, so the matrix conjugate transpose IS the HS adjoint; and the block is
    /// a PRINCIPAL submatrix (one index list for rows and columns), so its Hermitian part is the
    /// Hermitian part's submatrix. Herm being DIAGONAL is what removes the eigensolver: its own spectrum
    /// is read off. If rate is constant then Herm = c·Id, so L − c·Id is anti-Hermitian and L is NORMAL,
    /// no Jordan block available and every Re λ exactly c. Conversely the trace gives Σ Re λ = −2·Σ rate
    /// (<see cref="ProfileTraceResidual"/>), so if every Re λ sits on the floor −2·min(rate) then
    /// mean = min, and rate ≥ min pointwise forces rate ≡ min. Together: a block sits ENTIRELY on its
    /// floor exactly when its cell rates are constant, an entry-wise statement and never a measurement.</para>
    ///
    /// <para>TWO PREMISES THAT ARE NOT FREE. The dissipator must be Z-DEPHASING, the jumps diagonal in
    /// the computational basis; under amplitude damping the Hermitian part is neither diagonal nor this.
    /// And the SPECTRAL half needs a NUMBER-CONSERVING H, because only then is the index set an invariant
    /// block whose submatrix eigenvalues are the Liouvillian's: with a transverse field the entry-wise
    /// identity still holds bit-exactly while the floor read off that submatrix is a number the dynamics
    /// does not have. The test file gates exactly that pair.</para>
    ///
    /// <para>ONE FLOATING-POINT FENCE, and it is the reason the gates use DYADIC profiles. The criterion
    /// keys on rate CONSTANCY, not on γ uniformity, and sums of different site subsets that are equal in
    /// exact arithmetic need not be bit-equal in doubles. So a comparison to 0.0 here can read the
    /// SUMMATION ORDER rather than the physics, exactly as the generator residual next door does. The
    /// answer is dyadic γ, where every subset sum is exact, and never a widened threshold.</para></summary>
    public static double ProfileHermitianPartResidual(int n, int p, int q, double j, double delta, IReadOnlyList<double> gamma) =>
        ProfileHermitianPartResidual(n, p, q, ChainHamiltonian(n, j, delta), gamma);

    /// <summary>The same residual from an ALREADY BUILT H, for sweeps.</summary>
    public static double ProfileHermitianPartResidual(int n, int p, int q, ComplexMatrix h, IReadOnlyList<double> gamma) =>
        HermitianPartResidual(n, ProfileBlock(n, p, q, h, gamma), Cells(n, p, q), gamma);

    /// <summary>The residual off an ALREADY BUILT block and its cells, so a sweep can read both residuals
    /// from one build.</summary>
    private static double HermitianPartResidual(int n, ComplexMatrix l, IReadOnlyList<(int Ket, int Bra)> cells, IReadOnlyList<double> gamma)
    {
        int d = cells.Count;
        double worst = 0.0;
        for (int r = 0; r < d; r++)
            for (int c = 0; c < d; c++)
            {
                var herm = 0.5 * (l[r, c] + Complex.Conjugate(l[c, r]));
                if (r == c) herm += 2.0 * CellRate(n, cells[r].Ket, cells[r].Bra, gamma);
                worst = Math.Max(worst, herm.Magnitude);
            }
        return worst;
    }

    /// <summary>Σ Re L[r,r] + 2·Σ_cells rate, read off the diagonal without an eigensolver. Exactly 0.0,
    /// and it is the exact route to Σ Re λ because the trace of a matrix is the sum of its eigenvalues
    /// whether or not the matrix is diagonalizable, defective blocks included.</summary>
    public static double ProfileTraceResidual(int n, int p, int q, double j, double delta, IReadOnlyList<double> gamma) =>
        ProfileTraceResidual(n, p, q, ChainHamiltonian(n, j, delta), gamma);

    /// <summary>The same residual from an ALREADY BUILT H, for sweeps.</summary>
    public static double ProfileTraceResidual(int n, int p, int q, ComplexMatrix h, IReadOnlyList<double> gamma) =>
        TraceResidual(n, ProfileBlock(n, p, q, h, gamma), Cells(n, p, q), gamma);

    /// <summary>The trace residual off an ALREADY BUILT block and its cells.</summary>
    private static double TraceResidual(int n, ComplexMatrix l, IReadOnlyList<(int Ket, int Bra)> cells, IReadOnlyList<double> gamma)
    {
        double traceRe = 0.0, rateSum = 0.0;
        for (int r = 0; r < cells.Count; r++)
        {
            traceRe += l[r, r].Real;
            rateSum += CellRate(n, cells[r].Ket, cells[r].Bra, gamma);
        }
        return traceRe + 2.0 * rateSum;
    }

    /// <summary>max − min of the cell rates in the CANONICAL site convention: the exact predicate the
    /// fence turns on. Zero exactly ⟺ the block sits entirely on its floor, at any coupling and any Δ.
    /// Sibling of <see cref="ProfileRealDiagonalSpread"/>, which is the same quantity scaled by −2 and
    /// read in the other site convention (see <see cref="CellRate"/> for why that is invisible there).</summary>
    public static double CellRateSpread(int n, int p, int q, IReadOnlyList<double> gamma)
    {
        double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
        foreach (var (ket, bra) in Cells(n, p, q))
        {
            double rate = CellRate(n, ket, bra, gamma);
            lo = Math.Min(lo, rate); hi = Math.Max(hi, rate);
        }
        return hi - lo;
    }

    /// <summary>The Re-span of the block's SPECTRUM under a profile. THE ONLY EIGENSOLVER IN THIS WITNESS,
    /// and it is here to reproduce the entry's measured numbers rather than to decide anything: what the
    /// fence asserts is settled exactly by <see cref="CellRateSpread"/>. Read it as a number with an
    /// error model, never as the criterion.
    ///
    /// <para>Note what the span is NOT. A block can be spectrally FLAT under a profile while sitting at
    /// the WRONG constant, which is a coalescence and not the criterion: on a PINNED block that constant
    /// is the trace's, −2·|p−q|·γ̄, strictly BELOW the floor −2·min(rate), so a flattened block is never a
    /// block meeting the criterion. FLAT and ON THE FLOOR are two predicates and merging them makes the
    /// law unreadable. (Off the pinned set the trace constant is −2γ̄·(p + q − 2pq/N), which collapses to
    /// −2γ̄·|p−q| exactly when min(p,q) = 0 or max(p,q) = N; the short form is a CONSEQUENCE of pinning,
    /// not a general fact, and this member is general.)</para>
    ///
    /// <para><b>Which blocks flatten, and it is NOT a |p−q| parity.</b> Two earlier readings of this were
    /// wrong, the second one worse than the first, so the law is stated with its mechanism and its
    /// measured range. Flattening needs an antilinear involution of the block about the trace constant,
    /// i.e. a pairing of cells sending rate ↦ 2·mean(rate) − rate, and there are two independent sources.
    /// (1) THE FREE INDEX IS N/2. On a pinned block one index is 0 or N and the other is free; when the
    /// free one is N/2 the bitwise COMPLEMENT of the free config stays inside the block and sends
    /// rate ↦ Σγ − rate = 2·mean(rate), for ANY profile. That is four blocks at even N and NONE at odd N.
    /// At N=4 they carry |p−q| = 2 and at N=6 |p−q| = 3, which is where the parity reading died: "even"
    /// was an N=4 accident. (2) γ ANTI-PALINDROMIC, γ_l + γ_{N−1−l} constant, the R₉₀ locus of the
    /// parameter Klein group (F91): the chain reflection does the same for EVERY pinned block, so all 4N
    /// flatten. Source (2) being the REFLECTION, it needs H reflection-symmetric too (palindromic bond
    /// couplings); source (1) is the global complement X^N, which commutes with each XX+YY+ZZ bond and is
    /// therefore blind to the bond profile. Gated as a pair. Measured N = 3..6 at Δ = 1; both need J above a threshold that is profile- and
    /// N-dependent and is NOT gated. What Δ decides is whether there is an onset AT ALL, and it decides that for
    /// the COLLAPSE and not for the involution. That reading is OFF-LOCUS and the clause was missing until
    /// 2026-08-15: OFF the locus at Δ = 0 the palindrome still holds while the span SATURATES, so no coupling
    /// collapses it, and at Δ = 0.5 the same happens at N = 6 but not at N = 4; ON the locus every pinned block
    /// of dimension above one collapses, past its onset, at every Δ measured. The gate next door reads the off-locus case, which
    /// is why it pairs Δ = 0 with N = 4 and Δ = 0.5 with N = 6. The ramp profile
    /// and <see cref="DyadicRampProfile"/> both lie on the R₉₀ locus, which is the source that CAN reach the
    /// whole 4N; whether it does depends on the onset, and the linspace ramp at J = 1 is below it.</para>
    ///
    /// <para><b>The involution buys a PALINDROME; flatness is its unbroken case.</b> What follows at once
    /// from rate ↦ 2·mean(rate) − rate is that the Re-spectrum is SYMMETRIC about the trace constant, and
    /// that holds at EVERY coupling: the symmetry defect is machine zero at J = 0.05 exactly as at J = 10
    /// while the span moves from order one to zero. Flatness is the case where every eigenvalue sits ON
    /// the axis, and that is what the J threshold is the onset of. Below it the block is symmetric and
    /// wide, above it symmetric and collapsed.</para>
    ///
    /// <para>NOT F140, though it is next door and on the same locus. F140 pins an EIGENVALUE, −4γ̄, with
    /// multiplicity at least ⌊N/2⌋, exactly ⌊N/2⌋ for all but finitely many J at γ̄ ≠ 0, INSIDE the
    /// N²-dimensional (1,1)-type corner block on the Heisenberg
    /// chain; those are not pinned blocks in F153's sense, and its p + q even rule is an XY-chain
    /// observation, necessary and not sufficient, not a counterpart of anything above. Its proof states
    /// that THERE the spectrum is not palindromic about the root and no multiset symmetry argument
    /// applies, so the involutions above must not be carried across. What the two share is the locus, and
    /// that sharing is F91's.</para></summary>
    public static (double MinRe, double MaxRe) ProfileReSpan(int n, int p, int q, double j, double delta, IReadOnlyList<double> gamma)
    {
        var l = ProfileBlock(n, p, q, j, delta, gamma);
        double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
        foreach (var z in l.Evd().EigenValues)
        {
            if (z.Real < lo) lo = z.Real;
            if (z.Real > hi) hi = z.Real;
        }
        return (lo, hi);
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

    /// <summary>A dyadic ramp, one rate per site: the profile the fence node runs on. Dyadic on purpose,
    /// so every subset sum of it is exact and a comparison to 0.0 reads the physics rather than the
    /// summation order (see <see cref="ProfileHermitianPartResidual"/>'s floating-point fence).</summary>
    public static double[] DyadicRampProfile(int n) =>
        Enumerable.Range(0, n).Select(k => 0.25 + 0.25 * k).ToArray();

    private InspectableNode? _fence;

    /// <summary>The γ fence, live, and the point of the node is that its verdict costs no eigensolver.
    /// The master invariant Herm(L) = −2·diag(rate) is recomputed entry-wise over every block, and the
    /// fence's predicate is then a rate spread rather than a spectrum.</summary>
    private InspectableNode TheGammaFenceNode()
    {
        var gamma = DyadicRampProfile(N);
        var h = ChainHamiltonian(N, J, Delta);
        double worstHerm = 0.0, worstTrace = 0.0;
        int oneCellFlat = 0, pinnedLeaving = 0;
        for (int p = 0; p <= N; p++)
            for (int q = 0; q <= N; q++)
            {
                // ONE block build per (p, q), read twice. The two residuals each build their own otherwise,
                // and at N = 7 the middle block is 1225 × 1225, which dwarfs the Hamiltonian this hoists.
                var l = ProfileBlock(N, p, q, h, gamma);
                var cells = Cells(N, p, q);
                worstHerm = Math.Max(worstHerm, HermitianPartResidual(N, l, cells, gamma));
                worstTrace = Math.Max(worstTrace, Math.Abs(TraceResidual(N, l, cells, gamma)));
                if (!PinnedBlockFloorClaim.IsPinned(N, p, q)) continue;
                if (CellRateSpread(N, p, q, gamma) == 0.0) oneCellFlat++; else pinnedLeaving++;
            }

        var kids = new List<IInspectable>
        {
            new InspectableNode("the master invariant, entry-wise",
                summary: $"Herm(L_block) + 2·diag(rate) over all {(N + 1) * (N + 1)} blocks: {worstHerm.ToString(Inv)}. " +
                         "The dissipator is real and diagonal in the coherence basis and −i·ad_H is anti-Hermitian " +
                         "as a superoperator, so the Hermitian part is −2·diag(rate) for ANY Hermitian H, any graph, " +
                         "any Δ and any profile. It is DIAGONAL, which is what removes the eigensolver: its spectrum " +
                         "is read off, and Re λ is its Rayleigh quotient on the right eigenvector.", provenance: NodeProvenance.Live),
            new InspectableNode("the trace, the converse's exact route",
                summary: $"Re tr L + 2·Σ rate: {worstTrace.ToString(Inv)}. Σ Re λ = Re tr L holds for defective " +
                         "blocks too, the trace being the characteristic polynomial's, so 'entirely on the floor ⇒ " +
                         "rates constant' needs no eigensolver: mean = min plus rate ≥ min pointwise forces rate ≡ min.", provenance: NodeProvenance.Live),
            new InspectableNode("the fence, and its only exception",
                summary: $"{pinnedLeaving} of the {4 * N} pinned blocks leave their floor under this profile; the " +
                         $"other {oneCellFlat} are the four ONE-CELL blocks, which have no " +
                         "rate to spread. That is the sharp form: the criterion turns on the cell rates being " +
                         "CONSTANT, which on a pinned block of dimension above one happens exactly when γ is uniform.", provenance: NodeProvenance.Live),
            new InspectableNode("what the fence does NOT say",
                summary: "not that the block SPREADS. A block can stay spectrally flat under a profile while sitting " +
                         "at the WRONG constant: on a pinned block the trace's, −2·|p−q|·γ̄, and NOT −2·mean(γ), two " +
                         "expressions that coincide only at |p−q| = 1. Which blocks flatten is NOT a |p−q| parity. " +
                         "Flattening needs an antilinear involution about that constant, and there are two sources: " +
                         "the FREE INDEX being N/2, where the complement of the free config stays in the block and " +
                         "sends rate ↦ Σγ − rate under ANY profile (four blocks at even N, none at odd N); and γ " +
                         "ANTI-PALINDROMIC, the R₉₀ locus of the parameter Klein group, where the chain reflection " +
                         "does the same for every pinned block and all 4N flatten. What the involution buys at once " +
                         "is a PALINDROME of the Re-spectrum about that constant, holding at EVERY coupling; " +
                         "flatness is its UNBROKEN case, and the J threshold is that onset. Measured N = 3..6 at " +
                         "Δ = 1. What Δ decides is whether the collapse has an onset AT ALL, not the involution, " +
                         "and that reading is OFF-LOCUS (clause added 2026-08-15): off the locus at Δ = 0 the " +
                         "palindrome holds while the span saturates, and at Δ = 0.5 it saturates at N = 6 but not " +
                         "at N = 4; ON the locus every pinned block of dimension above one collapses, past its " +
                         "onset, at every Δ measured. The profile printed above " +
                         "lies on the R₉₀ locus, the source that CAN reach the whole 4N; whether it does at this " +
                         "coupling is the onset's business, and this node does not read it, its verdict being a " +
                         "rate spread and not a spectrum.", provenance: NodeProvenance.Live),
        };
        return new InspectableNode("the γ fence, spectral half (a site-resolved profile)",
            summary: $"γ = [{string.Join(", ", gamma.Select(x => x.ToString("0.##", Inv)))}] by site, dyadic so every " +
                     "subset sum is exact. The spectral half of the fence reduces to an ENTRY-WISE statement and is " +
                     $"decided at 0.0 exactly: worst Hermitian-part residual {worstHerm.ToString(Inv)}, worst trace " +
                     $"residual {worstTrace.ToString(Inv)}. An eigensolver appears nowhere in this verdict.",
            children: kids,
            provenance: NodeProvenance.Live);
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var rows = Rows;
            yield return new InspectableNode("the criterion",
                summary: "a block sits entirely on Re λ = −2γ|p−q| exactly when min(p,q) = 0 or max(p,q) = N; " +
                         $"that is 4(N+1) − 4 = {4 * N} of the {(N + 1) * (N + 1)} blocks, the four one-cell blocks among them");
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

            yield return _fence ??= TheGammaFenceNode();

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
