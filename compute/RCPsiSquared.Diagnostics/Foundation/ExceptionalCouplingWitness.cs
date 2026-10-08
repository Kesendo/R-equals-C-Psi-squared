using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The exceptional couplings of F50, live (<c>inspect --root exceptional</c>). The count
/// d_real(Re = −2γ) = 2N of <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> holds at every ratio γ/J
/// outside a finite set E(N, G): by the Absorption Theorem a real eigenvalue −2γ needs ⟨n_XY⟩ = 1, which
/// only the diagonal joint-popcount blocks (p, p) can reach without pure weight 1, and in them
/// P_p(γ) = det(B_pp(γ) + 2γ·I) is a polynomial with leading coefficient Π(2 − 2·Hamming) ≠ 0. E(N, G)
/// is the set of positive roots of the P_p ("The count at exceptional couplings").
///
/// <para><b>What is recomputed, and on which road.</b> The Python gate
/// (<c>simulations/f50_exceptional_couplings.py</c>) builds the blocks by hand and takes sympy's
/// determinant over ℤ[i][γ]. This witness deliberately walks the other road: the Hamiltonian is the
/// repo's <see cref="PauliHamiltonian.Bilinear"/> in the PAULI book (H = J·Σ(XX+YY+ZZ), J = 1, so γ/J is
/// the gate's unit; the spin book's J is four times this), the blocks are
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> on <see cref="SectorBlock.SectorFlatIndices"/>, the
/// polynomial is recovered EXACTLY from dim + 1 determinants over ℤ[i] (Berkowitz,
/// <see cref="GaussianMatrixCharpoly"/>) by integer Lagrange interpolation, and the positive roots are
/// COUNTED by Sturm's theorem (<see cref="CentreLineExactCount"/>'s chain) and ISOLATED by Sturm counts on
/// dyadic intervals, every sign read in exact integer arithmetic. The checks on this route that can fail: the
/// builder's dissipator diagonal equals −2·Hamming of the cell exactly (the leading coefficient Π(2 − 2·Hamming)
/// of Re P then follows as det D, an identity rather than a check), the isolated roots are as many as the Sturm
/// count says, and at the rational points the realified nullity is even. At a rational point the measured count
/// and nullity are taken at the rational itself, since a dyadic root is read at 2 − 2^−41 and each extra eigenvalue of
/// an EP2 there moves 4√δ ≈ 2.7·10⁻⁶ from −2γ, outside the window: a residual that is a function of the reading offset
/// alone, and a count below the nullity is now refused as that signature. The root polynomial is
/// Q_p = gcd(Re P, Im P), whose real roots are P's with the same multiplicity. For N ≤ 4 on the graphs the gate
/// pins, the pinned minimal polynomials are certified by EXACT DIVISION: q^m divides Q_p and q^(m+1) does not,
/// and the positive roots of the q's, with their multiplicities, add up to every positive root of Q_p, so E is
/// pinned without a float (that the q's are irreducible and pairwise distinct is the table's word, checked once
/// with sympy, not recomputed here).</para>
///
/// <para><b>What is measured.</b> The count at −2γ at a generic ratio and at every point of E, read off
/// the eigensolver on every (p, q) block (window 10⁻⁶ on |λ + 2γ|), and for the chain the spectral gap
/// on both sides of the smallest point: 2γ just below it, a real mode below 2γ just above it, which is
/// the two-route reading min E(chain) = 1/Q*_gap(N) the proof records at N = 2 to 5. Those are
/// readings, not the theorem.</para>
///
/// <para>Guard: N ≤ <see cref="MaxN"/> (the (2,2) block at N = 4 is 36-dimensional and takes 37
/// Berkowitz charpolys; N = 5's 100-dimensional blocks are left to the gate's <c>--n5</c>).</para></summary>
public sealed class ExceptionalCouplingWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest N the live build admits.</summary>
    public const int MaxN = 4;

    /// <summary>The generic ratio γ/J at which the count 2N (+ δ_G) is read beside the exceptional points.</summary>
    public const double GenericGammaOverJ = 1.3;

    /// <summary>The eigensolver window on |λ + 2γ| for the measured counts.</summary>
    public const double CountWindow = 1e-6;

    /// <summary>One diagonal or even off-diagonal block: its polynomial P(γ) = det(B(γ) + 2γ·I) over ℤ[i], ascending
    /// in γ, as a real and an imaginary part (the leading coefficient Π(2 − 2·Hamming) is real, the lower
    /// coefficients need not be), the root polynomial Q = gcd(Re P, Im P) over ℤ (or the nonzero part when the other
    /// vanishes), whose real roots are the real roots of P, and Q's positive roots (exact counts, isolated values).
    /// Multiplicities are read in Q, as the gate reads them.</summary>
    public sealed record BlockPolynomial(
        int PKet, int QBra, int Dim,
        BigInteger[] RealPart, BigInteger[] ImagPart,
        BigInteger[] Coefficients,
        BigInteger PredictedLeading,
        int DistinctPositiveRoots,
        int PositiveRootsWithMultiplicity,
        IReadOnlyList<double> Roots)
    {
        /// <summary>Re P has degree Dim with the predicted leading coefficient, and Im P has lower degree.</summary>
        public bool LeadingMatches => RealPart.Length == Dim + 1 && RealPart[Dim] == PredictedLeading && ImagPart.Length <= Dim;
        public bool RootsAllIsolated => Roots.Count == DistinctPositiveRoots;
    }

    /// <summary>A pinned minimal polynomial q (ascending in x = γ/J) of a point of E in block (p, p): q^m | Q_p
    /// exactly, q^(m+1) ∤ Q_p, with q's own positive-root count (the number of points of E it carries).</summary>
    public sealed record PinnedCertificate(
        int P, BigInteger[] MinimalPolynomial, int Multiplicity, bool PowerDivides, bool NextPowerFails, int PositiveRoots)
    {
        public bool Holds => PowerDivides && NextPowerFails;
    }

    public int N { get; }
    public string Topology { get; }
    public IReadOnlyList<(int A, int B)> Edges { get; }

    public IReadOnlyList<BlockPolynomial> DiagonalBlocks { get; }
    public IReadOnlyList<BlockPolynomial> OffDiagonalEvenBlocks { get; }

    /// <summary>E(N, G): the distinct positive roots over all diagonal blocks, ascending, as γ/J.</summary>
    public IReadOnlyList<double> ExceptionalSet { get; }

    /// <summary>The measured count at −2γ at <see cref="GenericGammaOverJ"/>.</summary>
    public int GenericCount { get; }

    /// <summary>The measured count at −2γ at every point of E, in the order of <see cref="ExceptionalSet"/>; at a rational
    /// point the eigensolver runs at the rational itself, elsewhere at the 2^−41-accurate reading.</summary>
    public IReadOnlyList<int> CountsAtPoints { get; }

    /// <summary>The measured nullity of L + 2γ at every point of E (singular values below 10⁻⁸·σ_max over all blocks),
    /// the gate's G5 companion, evaluated where the count is. Where it equals the count the eigenvalue is semisimple to
    /// the eigensolver, and the window count is then a law: a semisimple eigenvalue moves by O(ε‖L‖) under the solver's
    /// backward error, and by O(δ) from the 2^−41 reading offset at an irrational point, both far inside 10⁻⁶. Where the
    /// count exceeds it the point is defective, and the window is a reading rather than a law (an EP of order k splits
    /// as ε^(1/k), about 10⁻⁷ for an EP2 at machine ε); the rational such points are read exactly in
    /// <see cref="ExactAtRationalPoints"/>. A count BELOW the nullity is the third case and the signature of a reading
    /// offset at a defective point (each extra eigenvalue of an EP2 moves 4√δ from −2γ, 2.7·10⁻⁶ at δ = 2^−41): the
    /// constructor throws on it, since it is the error this witness once carried.</summary>
    public IReadOnlyList<int> NullitiesAtPoints { get; }

    /// <summary>The pinned certificates (empty when the gate pins nothing for this (N, graph)).</summary>
    public IReadOnlyList<PinnedCertificate> Pinned { get; }

    /// <summary>Whether a pinned table exists for this (N, graph).</summary>
    public bool HasPinnedTable { get; }

    /// <summary>With a pinned table: every certificate holds and, block by block, the positive roots of the pinned
    /// q's are exactly as many as the Sturm count of P_p, so E is pinned without a float.</summary>
    public bool PinnedCoversAllRoots { get; }

    /// <summary>At a RATIONAL point of E the multiplicity of −2γ on the diagonal blocks, exactly: algebraic (the
    /// generalized kernel) and geometric (the kernel), both over ℚ(i) by realification to ℚ and
    /// <see cref="BigRationalLinearAlgebra.Nullspace"/>, summed over 1 ≤ p ≤ N−1. These ARE the extra modes, since off
    /// the diagonal blocks −2γ carries only the commutant. The window count is no law at a defective point (an EP
    /// splits as √ε), which is where this exact route replaces it.</summary>
    public sealed record ExactMultiplicity(double Point, BigInteger Numerator, BigInteger Denominator, int Algebraic, int Geometric);

    /// <summary>One entry per rational point of E (the N = 2 chain's and the N = 4 ring's γ/J = 2, the defective ones).</summary>
    public IReadOnlyList<ExactMultiplicity> ExactAtRationalPoints { get; }

    /// <summary>Theorem D of the plane-crossing section (the inertia identity n(γ) = #{r_j(γ) &lt; 1}) makes Theorem C an
    /// equality: in block (p, p), #E_p = C(N, p) − c_p − m_p with c_p the components of the exclusion graph and m_p the
    /// Krein debt, the non-stationary modes of the block inside the half-plane Re λ &gt; −2γ at the Hamiltonian end.
    /// <see cref="Debt"/> is read EXACTLY from the Sturm count with multiplicity, no eigensolver; <see cref="DebtRead"/>
    /// is the second route, the eigensolver's count of modes with Re λ + 2γ above 10⁻⁹ (Pauli J = 1) at J/γ = 30 and 60
    /// (the same at both, the law), minus c_p. A limit-1 branch's height rises to 1 as (γ/J)², so Re λ + 2γ falls as
    /// γ³/J²: at least 5.8·10⁻⁷ at J/γ = 60, 2.7 decades above the window, which is 5 decades above rounding
    /// eps·‖B‖ ≈ 10⁻¹⁴. <see cref="StationaryRead"/> is the eigensolver's second route to c_p, the modes with |λ| below
    /// 10⁻⁹ (the stationary modes, one per component). The routes meeting is the handshake; the 2γ regime of D06 exists
    /// iff every debt is zero.</summary>
    public sealed record KreinDebt(int P, int Binomial, int Components, int ExceptionalCountWithMultiplicity, int Debt, int DebtRead, int StationaryRead)
    {
        public bool RoutesMeet => Debt == DebtRead && Components == StationaryRead;
    }

    /// <summary>The Krein debt m_p of every diagonal block 1 ≤ p ≤ N − 1, exact and read.</summary>
    public IReadOnlyList<KreinDebt> Debts { get; }

    /// <summary>The 2γ regime exists (every mode slower than 2γ for γ below min E) iff m_p = 0 in every block, which by
    /// Theorem D is #E_p = C(N, p) − c_p in every block: decided here by the exact count.</summary>
    public bool GapRegimeExists => Debts.All(d => d.Debt == 0);

    /// <summary>The chain's gap reading on both sides of min E: (gap − 2γ just below, 2γ − gap just above,
    /// |Im| of the slow mode just above). Null off the chain or when E is empty.</summary>
    public (double BelowOffset, double AboveOffset, double AboveImag)? GapSides { get; }

    /// <summary>J/γ at the smallest point of E (the handover ratio); null when E is empty.</summary>
    public double? HandoverQ => ExceptionalSet.Count == 0 ? null : 1.0 / ExceptionalSet[0];

    /// <summary>The chain's gap reading holds: gap = 2γ below min E (within 10⁻⁹) and a real mode below 2γ above it.</summary>
    public bool GapSidesHold => GapSides is { } g && Math.Abs(g.BelowOffset) < 1e-9 && g.AboveOffset > 1e-8 && g.AboveImag < 1e-9;

    public ExceptionalCouplingWitness(int n = 3, string? topology = "chain")
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"N must lie in 2..{MaxN}.");
        N = n;
        Topology = (topology ?? "chain").ToLowerInvariant();
        if (Topology == "ring" && n < 3)
            throw new ArgumentException("a ring needs N ≥ 3 (at N = 2 the wrap bond doubles the chain's bond).", nameof(topology));
        Edges = Topology switch
        {
            "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray(),
            "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray(),
            "star" => Enumerable.Range(1, n - 1).Select(i => (0, i)).ToArray(),
            "complete" => (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
            _ => throw new ArgumentException($"--topology takes chain, ring, star or complete; got \"{topology}\".", nameof(topology)),
        };

        var H = BuildHamiltonian(n, Edges);
        int d = 1 << n;

        var diag = new List<BlockPolynomial>();
        var off = new List<BlockPolynomial>();
        for (int p = 0; p <= n; p++)
            for (int q = p; q <= n; q++)
            {
                if (((q - p) & 1) == 1) continue;
                var bp = BlockPolynomialOf(H, n, d, p, q);
                if (p == q) diag.Add(bp); else off.Add(bp);
            }
        DiagonalBlocks = diag;
        OffDiagonalEvenBlocks = off;

        var debts = new List<KreinDebt>();
        foreach (var b in diag.Where(b => b.PKet >= 1 && b.PKet <= n - 1))
        {
            int components = ExclusionGraphComponents(n, b.PKet, Edges);
            int binomial = (int)Binomial(n, b.PKet);
            int debt = binomial - components - b.PositiveRootsWithMultiplicity;
            var (inside30, stationary30) = ModesInsideTheHalfPlane(H, n, b.PKet, 1.0 / 30);
            var (inside60, stationary60) = ModesInsideTheHalfPlane(H, n, b.PKet, 1.0 / 60);
            if (inside30 != inside60 || stationary30 != stationary60)
                throw new InvalidOperationException(
                    $"block ({b.PKet},{b.PKet}): the eigensolver reads {inside30}/{stationary30} modes inside the half-plane/stationary at J/γ = 30 and " +
                    $"{inside60}/{stationary60} at 60; both counts are constant below min E (Theorem D), so a change is a window failure, not a result.");
            debts.Add(new KreinDebt(b.PKet, binomial, components, b.PositiveRootsWithMultiplicity, debt, inside30 - components, stationary30));
        }
        Debts = debts;

        var points = new SortedSet<double>();
        foreach (var b in diag) foreach (var r in b.Roots) points.Add(r);
        // every root is read at the midpoint of an aligned 2^−40 cell, so two blocks sharing a root (p and N−p are
        // mirror images) read it identically and the set merges them exactly; a residual difference below 1e-9
        // would be a construction finding, so it is flagged
        var merged = new List<double>();
        foreach (var r in points)
            if (merged.Count == 0 || r - merged[^1] > 1e-9) merged.Add(r);
            else throw new InvalidOperationException($"two isolated roots differ by less than 1e-9 but are not identical: {merged[^1]:R} vs {r:R}");
        ExceptionalSet = merged;

        if (ExceptionalSet.Any(r => Math.Abs(r - GenericGammaOverJ) < 1e-3))
            throw new InvalidOperationException($"the generic ratio {GenericGammaOverJ} lies within 1e-3 of a point of E({n}, {Topology}); pick another.");

        // the rational points first: there the measured readings are taken AT the rational, not at the dyadic reading
        // (a dyadic root such as 2 is read at 2 − 2^−41, and each extra eigenvalue of an EP2 there moves 4√δ ≈ 2.7·10⁻⁶
        // from −2γ, outside the window: a residual that is a function of the reading offset and of nothing else, so the
        // exact point is used; exact in double for a dyadic rational, which 2 is)
        var exact = new List<ExactMultiplicity>();
        var evalAt = ExceptionalSet.ToArray();
        for (int i = 0; i < evalAt.Length; i++)
            if (RationalPoint(evalAt[i], diag) is { } ab)
            {
                evalAt[i] = (double)ab.Num / (double)ab.Den;
                int alg = 0, geo = 0;
                for (int p = 1; p < n; p++)
                {
                    var (a, g) = ExactMultiplicityOnBlock(H, n, d, p, ab.Num, ab.Den);
                    alg += a; geo += g;
                }
                exact.Add(new ExactMultiplicity(ExceptionalSet[i], ab.Num, ab.Den, alg, geo));
            }
        ExactAtRationalPoints = exact;

        GenericCount = CountAtMinusTwoGamma(H, n, GenericGammaOverJ);
        CountsAtPoints = evalAt.Select(g => CountAtMinusTwoGamma(H, n, g)).ToArray();
        NullitiesAtPoints = evalAt.Select(g => NullityAtMinusTwoGamma(H, n, g)).ToArray();
        for (int i = 0; i < evalAt.Length; i++)
            if (CountsAtPoints[i] < NullitiesAtPoints[i])
                throw new InvalidOperationException(
                    $"at γ/J = {evalAt[i]:R} the window count {CountsAtPoints[i]} is below the nullity {NullitiesAtPoints[i]}: a defective " +
                    "point evaluated off its exact position (each extra eigenvalue of an EP2 moves 4√δ from −2γ), which is a reading error, not a result.");

        // at N = 3 the ring is K₃ and the star is the chain (isomorphic graphs share every block polynomial)
        string pinnedKey = (n, Topology) switch { (2, _) => "chain", (3, "ring") => "complete", (3, "star") => "chain", _ => Topology };
        HasPinnedTable = PinnedTable.TryGetValue((n, pinnedKey), out var table);
        var certs = new List<PinnedCertificate>();
        bool covers = HasPinnedTable;
        if (HasPinnedTable)
        {
            foreach (var (p, entries) in table!)
            {
                var block = diag.Single(b => b.PKet == p);
                int rootsPinned = 0;
                foreach (var grp in entries.GroupBy(e => (Key: string.Join(",", e.Poly), e.Multiplicity)))
                {
                    var q = grp.First().Poly;
                    int m = grp.Key.Multiplicity;
                    var pow = Power(q, m);
                    bool divides = ExactDivide(block.Coefficients, pow) is not null;
                    bool nextFails = ExactDivide(block.Coefficients, Multiply(pow, q)) is null;
                    int posRoots = DistinctPositiveRoots(q);
                    var cert = new PinnedCertificate(p, q, m, divides, nextFails, posRoots);
                    certs.Add(cert);
                    covers &= cert.Holds && posRoots == grp.Count();
                    rootsPinned += posRoots;
                }
                covers &= rootsPinned == block.DistinctPositiveRoots;
                covers &= entries.Sum(e => e.Multiplicity) == block.PositiveRootsWithMultiplicity;
            }
            // a block the table does not list must have no positive root
            foreach (var b in diag)
                if (!table!.ContainsKey(b.PKet)) covers &= b.DistinctPositiveRoots == 0;
        }
        Pinned = certs;
        PinnedCoversAllRoots = covers;

        if (Topology == "chain" && ExceptionalSet.Count > 0)
            GapSides = ReadGapSides(H, n, ExceptionalSet[0]);
    }

    // ---------------------------------------------------------------- the Hamiltonian, Pauli book

    private static ComplexMatrix BuildHamiltonian(int n, IReadOnlyList<(int A, int B)> edges)
    {
        var bonds = edges.Select(e => new Bond(e.A, e.B, 1.0)).ToArray();
        var terms = new (PauliLetter, PauliLetter, Complex)[]
        {
            (PauliLetter.X, PauliLetter.X, Complex.One),
            (PauliLetter.Y, PauliLetter.Y, Complex.One),
            (PauliLetter.Z, PauliLetter.Z, Complex.One),
        };
        return PauliHamiltonian.Bilinear(n, bonds, terms).ToMatrix();
    }

    /// <summary>The flat indices of the gate's block (p, q) = (popcount ket, popcount bra); the sector builder
    /// labels (PCol = bra, PRow = ket).</summary>
    private static int[] Flat(int n, int pKet, int qBra) => SectorBlock.SectorFlatIndices(n, qBra, pKet);

    // ---------------------------------------------------------------- the exact polynomial

    private static BlockPolynomial BlockPolynomialOf(ComplexMatrix H, int n, int d, int pKet, int qBra)
    {
        // the γ-free part exactly in ℤ[i] (every entry of −i(H⊗I − I⊗Hᵀ) is an integer multiple of i here) and the
        // coefficient of γ on the diagonal, 2 − 2·Hamming, read from the cell rather than from the builder
        var (A0, delta, _) = ExactBlock(H, n, d, pKet, qBra);
        int m = delta.Length;
        BigInteger predicted = BigInteger.One;
        foreach (var dk in delta) predicted *= dk;

        // P(g) = det(A0 + g·diag(delta)) ∈ ℤ[i][g], sampled at g = 0..m, real and imaginary parts interpolated over ℤ
        var re = new BigInteger[m + 1];
        var im = new BigInteger[m + 1];
        for (int g = 0; g <= m; g++)
        {
            var A = (GaussianInteger[,])A0.Clone();
            for (int k = 0; k < m; k++) A[k, k] = A[k, k] + new GaussianInteger(g * delta[k], BigInteger.Zero);
            var cp = GaussianMatrixCharpoly.Characteristic(A);   // det(λI − A), lowest first
            var det = (m & 1) == 0 ? cp[0] : -cp[0];             // det A = (−1)^m · cp[0]
            re[g] = det.Re;
            im[g] = det.Im;
        }
        var realPart = SeatBlindnessDeltaLocusWitness.InterpolateAtZeroToN(re);
        var imagPart = SeatBlindnessDeltaLocusWitness.InterpolateAtZeroToN(im);
        // the real roots of P are the common real roots of its two parts
        BigInteger[] rootPoly = imagPart.Length == 0 ? CentreLineExactCount.Primitive(realPart)
                              : realPart.Length == 0 ? CentreLineExactCount.Primitive(imagPart)
                              : CentreLineExactCount.Gcd(realPart, imagPart);

        int distinct = DistinctPositiveRoots(rootPoly);
        int withMult = PositiveRootsWithMultiplicity(rootPoly);
        var roots = IsolatePositiveRoots(rootPoly);
        if (roots.Count != distinct)
            throw new InvalidOperationException(
                $"block ({pKet},{qBra}): Sturm counts {distinct} distinct positive roots, isolation found {roots.Count}.");
        return new BlockPolynomial(pKet, qBra, m, realPart, imagPart, rootPoly, predicted, distinct, withMult, roots);
    }

    /// <summary>The exact ℤ[i] block B(0) = −i(H⊗I − I⊗Hᵀ) restricted, and the γ-coefficient 2 − 2·Hamming per cell.
    /// The Hamming count is NOT computed by hand alone: the builder's dissipator diagonal at γ = 1 minus the γ = 0 block
    /// must equal −2·Hamming(x, y) of the cell, entry for entry and exactly, with every off-diagonal entry unchanged.
    /// That is the check on this route that can fail (a wrong rate convention or a wrong block), since the leading
    /// coefficient of det(A0 + γ·D) being det D is an arithmetic identity whatever D holds.</summary>
    private static (GaussianInteger[,] A0, BigInteger[] Delta, int[] Flat) ExactBlock(ComplexMatrix H, int n, int d, int pKet, int qBra)
    {
        var flat = Flat(n, pKet, qBra);
        int m = flat.Length;
        var M0 = PerBlockLiouvillianBuilder.BuildBlockZ(H, new double[n], flat);
        var M1 = PerBlockLiouvillianBuilder.BuildBlockZ(H, Enumerable.Repeat(1.0, n).ToArray(), flat);
        var A0 = new GaussianInteger[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                A0[i, j] = new GaussianInteger(ExactInteger(M0[i, j].Real, i, j), ExactInteger(M0[i, j].Imaginary, i, j));
                if (i != j && M1[i, j] != M0[i, j])
                    throw new InvalidOperationException($"the dissipator changed off-diagonal entry ({i},{j}) of block ({pKet},{qBra}); Z-dephasing is diagonal in this basis.");
            }
        var delta = new BigInteger[m];
        for (int k = 0; k < m; k++)
        {
            int x = flat[k] / d, y = flat[k] % d;
            var diss = M1[k, k] - M0[k, k];
            if (diss.Imaginary != 0.0 || diss.Real != -2.0 * BitOperations.PopCount((uint)(x ^ y)))
                throw new InvalidOperationException(
                    $"the builder's dissipator on cell |{x}⟩⟨{y}| of block ({pKet},{qBra}) is {diss.Real:R} + {diss.Imaginary:R}i at γ = 1, " +
                    $"not −2·Hamming = {-2 * BitOperations.PopCount((uint)(x ^ y))}; the rate convention or the block is wrong.");
            delta[k] = 2 + new BigInteger(diss.Real);
        }
        return (A0, delta, flat);
    }

    /// <summary>Is the isolated root a rational a/b? By the rational-root theorem b divides the leading coefficient of
    /// the carrying block's root polynomial, so its divisors are walked in ascending order; a candidate within 2^−30 of
    /// the reading is decided EXACTLY, as a root of every diagonal block's root polynomial that carries the reading. A
    /// rational with a denominator beyond about 2^40 could not be proposed from a 2^−41 reading and would be missed; at
    /// N ≤ 4 the pinned tables certify that x − 2 is the only linear factor.</summary>
    private static (BigInteger Num, BigInteger Den)? RationalPoint(double r, IReadOnlyList<BlockPolynomial> diag)
    {
        var carrying = diag.Where(blk => blk.Roots.Any(x => Math.Abs(x - r) < 1e-9)).ToList();
        if (carrying.Count == 0) return null;
        foreach (var b in Divisors(BigInteger.Abs(carrying[0].Coefficients[^1])))
        {
            double bd = (double)b;
            double a = Math.Round(r * bd);
            if (Math.Abs(a / bd - r) > Math.Pow(2.0, -30)) continue;
            var num = new BigInteger(a);
            if (carrying.All(blk => EvaluateScaled(blk.Coefficients, num, b).IsZero)) return (num, b);
        }
        return null;
    }

    /// <summary>The positive divisors of n, ascending (trial division; the leading coefficients here are products of
    /// small even numbers).</summary>
    private static IEnumerable<BigInteger> Divisors(BigInteger n)
    {
        var large = new List<BigInteger>();
        for (BigInteger k = 1; k * k <= n; k++)
            if ((n % k).IsZero) { yield return k; if (k * k != n) large.Add(n / k); }
        for (int i = large.Count - 1; i >= 0; i--) yield return large[i];
    }

    /// <summary>q(a/b)·b^deg, an integer.</summary>
    private static BigInteger EvaluateScaled(BigInteger[] q, BigInteger num, BigInteger den)
    {
        int deg = q.Length - 1;
        BigInteger v = BigInteger.Zero, numPow = BigInteger.One;
        var denPow = new BigInteger[deg + 1];
        denPow[0] = BigInteger.One;
        for (int j = 1; j <= deg; j++) denPow[j] = denPow[j - 1] * den;
        for (int j = 0; j <= deg; j++) { v += q[j] * numPow * denPow[deg - j]; numPow *= num; }
        return v;
    }

    /// <summary>(algebraic, geometric) multiplicity of the eigenvalue −2γ of block (p,p) at γ/J = a/b, exactly over
    /// ℚ(i): A = B + 2γ·I = X + iY is realified to R = [[X, −Y], [Y, X]] over ℚ, whose nullity is twice A's, and the
    /// generalized kernel is read off R, R², … until the nullity stops growing.</summary>
    private static (int Algebraic, int Geometric) ExactMultiplicityOnBlock(ComplexMatrix H, int n, int d, int p, BigInteger num, BigInteger den)
    {
        var (A0, delta, _) = ExactBlock(H, n, d, p, p);
        int m = delta.Length;
        var g = new BigRational(num, den);
        var R = new BigRational[2 * m, 2 * m];
        for (int i = 0; i < 2 * m; i++) for (int j = 0; j < 2 * m; j++) R[i, j] = BigRational.Zero;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                var x = new BigRational(A0[i, j].Re);
                var y = new BigRational(A0[i, j].Im);
                if (i == j) x = x + g * new BigRational(delta[i]);
                R[i, j] = x; R[i + m, j + m] = x;
                R[i, j + m] = -y; R[i + m, j] = y;
            }
        int Nullity(BigRational[,] mat) => BigRationalLinearAlgebra.Nullspace(mat).Count;
        int geo2 = Nullity(R);
        if ((geo2 & 1) == 1) throw new InvalidOperationException("the realified nullity is odd; the realification is wrong.");
        var P = R; int prev = geo2;
        for (int k = 2; k <= 2 * m; k++)
        {
            P = BigRationalLinearAlgebra.Multiply(P, R);
            int nk = Nullity(P);
            if (nk == prev) break;
            prev = nk;
        }
        return (prev / 2, geo2 / 2);
    }

    private static BigInteger ExactInteger(double v, int i, int j)
    {
        double r = Math.Round(v);
        if (r != v)
            throw new InvalidOperationException($"block entry ({i},{j}) = {v:R} is not an integer; the Pauli-book H with J = 1 must give integer entries.");
        return new BigInteger(r);
    }

    // ---------------------------------------------------------------- exact root counting and isolation

    /// <summary>P(y²): its real roots are 0 (iff P(0) = 0) and ±√r for every positive root r of P.</summary>
    private static BigInteger[] EvenLift(BigInteger[] p)
    {
        var t = CentreLineExactCount.Trim(p);
        var q = new BigInteger[2 * t.Length - 1];
        for (int k = 0; k < t.Length; k++) q[2 * k] = t[k];
        return q;
    }

    /// <summary>Distinct positive roots of an integer polynomial, by Sturm's theorem on P(y²).</summary>
    public static int DistinctPositiveRoots(BigInteger[] p)
    {
        var t = CentreLineExactCount.Trim(p);
        if (t.Length <= 1) return 0;
        int all = CentreLineExactCount.DistinctRealRoots(EvenLift(t));
        return (all - (t[0].IsZero ? 1 : 0)) / 2;
    }

    /// <summary>Positive roots with multiplicity: the real roots of P(y²) with multiplicity, minus twice the
    /// order of P at 0, halved.</summary>
    public static int PositiveRootsWithMultiplicity(BigInteger[] p)
    {
        var t = CentreLineExactCount.Trim(p);
        if (t.Length <= 1) return 0;
        int ord = 0;
        while (t[ord].IsZero) ord++;
        return (CentreLineExactCount.RealRootsWithMultiplicity(EvenLift(t)) - 2 * ord) / 2;
    }

    /// <summary>The distinct positive roots, ascending, each isolated by Sturm counts on dyadic intervals of
    /// the squarefree part and read at the midpoint of an interval of width 2^−40 (every sign exact).</summary>
    public static IReadOnlyList<double> IsolatePositiveRoots(BigInteger[] p)
    {
        var t = CentreLineExactCount.Trim(p);
        if (t.Length <= 1) return Array.Empty<double>();
        var g = CentreLineExactCount.Gcd(t, CentreLineExactCount.Derivative(t));
        var s = ExactDivide(t, g) ?? throw new InvalidOperationException("P / gcd(P, P′) is not exact; the gcd is wrong.");
        s = CentreLineExactCount.Primitive(s);
        if (s.Length <= 1) return Array.Empty<double>();

        var chain = new List<BigInteger[]> { s, CentreLineExactCount.Primitive(CentreLineExactCount.Derivative(s)) };
        while (chain[^1].Length > 1)
        {
            var r = CentreLineExactCount.PositivePseudoRemainder(chain[^2], chain[^1]);
            if (r.Length == 0) break;
            chain.Add(CentreLineExactCount.Primitive(CentreLineExactCount.Negate(r)));
        }

        // Cauchy bound on |root|: 1 + max|c_k| / |c_n|
        BigInteger maxAbs = BigInteger.Zero;
        for (int k = 0; k < s.Length - 1; k++) maxAbs = BigInteger.Max(maxAbs, BigInteger.Abs(s[k]));
        BigInteger bound = 2 + maxAbs / BigInteger.Abs(s[^1]);

        // the interval (0, 2^b] with 2^b ≥ bound is halved k times to width 2^(b−k): every root is read at the midpoint
        // of an ALIGNED cell of width 2^−IsolationDepth, so two polynomials sharing a root read it identically
        int b = (int)bound.GetBitLength();
        int depth = IsolationDepth + b;
        var roots = new List<double>();
        Isolate(chain, BigInteger.Zero, BigInteger.One << b, 0, depth, roots);
        roots.Sort();
        return roots;
    }

    private const int IsolationDepth = 40;

    /// <summary>Roots of the squarefree part in (lo, hi] with lo, hi in units of 2^−k.</summary>
    private static void Isolate(List<BigInteger[]> chain, BigInteger lo, BigInteger hi, int k, int depth, List<double> roots)
    {
        int count = SignChanges(chain, lo, k) - SignChanges(chain, hi, k);
        if (count == 0) return;
        if (k > depth + IsolationDepth)
            throw new InvalidOperationException($"{count} roots of the squarefree part within 2^−{k - depth + IsolationDepth}; isolation gives up.");
        if (count == 1 && k >= depth)
        {
            roots.Add((double)(lo + hi) / Math.Pow(2.0, k + 1));
            return;
        }
        var mid = lo + hi;                         // (lo + hi) / 2 in units of 2^−(k+1)
        Isolate(chain, 2 * lo, mid, k + 1, depth, roots);
        Isolate(chain, mid, 2 * hi, k + 1, depth, roots);
    }

    /// <summary>Sign changes of the Sturm chain at t = num / 2^k, zeros skipped.</summary>
    private static int SignChanges(List<BigInteger[]> chain, BigInteger num, int k)
    {
        int changes = 0, last = 0;
        BigInteger den = BigInteger.One << k;
        foreach (var q in chain)
        {
            if (q.Length == 0) continue;
            // q(num/den) · den^(deg) = Σ c_j · num^j · den^(deg − j)
            BigInteger v = BigInteger.Zero, numPow = BigInteger.One;
            int deg = q.Length - 1;
            var denPows = new BigInteger[deg + 1];
            denPows[0] = BigInteger.One;
            for (int j = 1; j <= deg; j++) denPows[j] = denPows[j - 1] * den;
            for (int j = 0; j <= deg; j++) { v += q[j] * numPow * denPows[deg - j]; numPow *= num; }
            int sgn = v.Sign;
            if (sgn == 0) continue;
            if (last != 0 && sgn != last) changes++;
            last = sgn;
        }
        return changes;
    }

    // ---------------------------------------------------------------- exact polynomial arithmetic

    /// <summary>a / b over ℤ[x], or null when the division is not exact (a nonzero remainder or a non-integral
    /// quotient coefficient).</summary>
    public static BigInteger[]? ExactDivide(BigInteger[] a, BigInteger[] b)
    {
        a = CentreLineExactCount.Trim(a);
        b = CentreLineExactCount.Trim(b);
        if (b.Length == 0) throw new DivideByZeroException("division by the zero polynomial");
        if (a.Length == 0) return Array.Empty<BigInteger>();
        if (a.Length < b.Length) return null;
        var r = (BigInteger[])a.Clone();
        var q = new BigInteger[a.Length - b.Length + 1];
        for (int k = q.Length - 1; k >= 0; k--)
        {
            var top = r[k + b.Length - 1];
            if (top.IsZero) continue;
            var qq = BigInteger.DivRem(top, b[^1], out var rem);
            if (!rem.IsZero) return null;
            q[k] = qq;
            for (int j = 0; j < b.Length; j++) r[k + j] -= qq * b[j];
        }
        for (int j = 0; j < b.Length - 1; j++) if (!r[j].IsZero) return null;
        return q;
    }

    public static BigInteger[] Multiply(BigInteger[] a, BigInteger[] b)
    {
        var r = new BigInteger[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++) r[i + j] += a[i] * b[j];
        return r;
    }

    private static BigInteger[] Power(BigInteger[] a, int m)
    {
        var r = new BigInteger[] { BigInteger.One };
        for (int i = 0; i < m; i++) r = Multiply(r, a);
        return r;
    }

    // ---------------------------------------------------------------- the measured side

    /// <summary>The number of eigenvalues of L within <see cref="CountWindow"/> of −2γ (complex distance), summed over
    /// every (p, q) block, at uniform γ = gammaOverJ in the Pauli book with J = 1.</summary>
    public static int CountAtMinusTwoGamma(ComplexMatrix H, int n, double gammaOverJ)
    {
        int count = 0;
        foreach (var ev in AllEigenvalues(H, n, gammaOverJ))
            if ((ev + 2 * gammaOverJ).Magnitude < CountWindow) count++;
        return count;
    }

    /// <summary>The nullity of L + 2γ·I read block by block: singular values below 10⁻⁸ times the largest singular
    /// value over all blocks (the gate's G5 threshold on the full L).</summary>
    public static int NullityAtMinusTwoGamma(ComplexMatrix H, int n, double gammaOverJ)
    {
        var gammas = Enumerable.Repeat(gammaOverJ, n).ToArray();
        var spectra = new List<double[]>();
        double sigmaMax = 0;
        foreach (var r in JointPopcountSectorBuilder.Build(n).SectorRanges)
        {
            var flat = SectorBlock.SectorFlatIndices(n, r.PCol, r.PRow);
            var B = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammas, flat);
            for (int i = 0; i < flat.Length; i++) B[i, i] += 2 * gammaOverJ;
            var s = B.Svd(computeVectors: false).S.Select(c => c.Real).ToArray();
            spectra.Add(s);
            sigmaMax = Math.Max(sigmaMax, s.Length == 0 ? 0 : s[0]);
        }
        return spectra.Sum(s => s.Count(v => v < 1e-8 * sigmaMax));
    }

    /// <summary>c_p: the components of the exclusion graph on the popcount-p configurations (two configurations
    /// adjacent when one bond moves one excitation between them), by union-find. One on every connected graph.</summary>
    public static int ExclusionGraphComponents(int n, int p, IReadOnlyList<(int A, int B)> edges)
    {
        var confs = Enumerable.Range(0, 1 << n).Where(c => BitOperations.PopCount((uint)c) == p).ToArray();
        var parent = confs.ToDictionary(c => c, c => c);
        int Find(int c)
        {
            while (parent[c] != c) { parent[c] = parent[parent[c]]; c = parent[c]; }
            return c;
        }
        foreach (var c in confs)
            foreach (var (a, b) in edges)
                if (((c >> a) & 1) != ((c >> b) & 1))
                    parent[Find(c)] = Find(c ^ (1 << a) ^ (1 << b));
        return confs.Select(Find).Distinct().Count();
    }

    /// <summary>n(γ) of Theorem C read by the eigensolver on block (p, p), the eigenvalues with Re λ + 2γ &gt; 10⁻⁹
    /// (Pauli J = 1), and beside it the stationary count, the eigenvalues with |λ| &lt; 10⁻⁹ (c_p of them). A complex
    /// pair running along the line (Re λ = −2γ exactly, Theorem A's boundary case) is not counted; a branch with
    /// r_j(0⁺) = 1 sits inside with Re λ + 2γ of order γ³/J², about 4.6·10⁻⁶ at J/γ = 30.</summary>
    public static (int Inside, int Stationary) ModesInsideTheHalfPlane(ComplexMatrix H, int n, int p, double gammaOverJ)
    {
        var gammas = Enumerable.Repeat(gammaOverJ, n).ToArray();
        var flat = SectorBlock.SectorFlatIndices(n, p, p);
        var B = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammas, flat);
        var ev = B.Evd().EigenValues;
        return (ev.Count(e => e.Real + 2 * gammaOverJ > 1e-9), ev.Count(e => e.Magnitude < 1e-9));
    }

    private static IEnumerable<Complex> AllEigenvalues(ComplexMatrix H, int n, double gamma)
    {
        var gammas = Enumerable.Repeat(gamma, n).ToArray();
        foreach (var r in JointPopcountSectorBuilder.Build(n).SectorRanges)
        {
            var flat = SectorBlock.SectorFlatIndices(n, r.PCol, r.PRow);
            var B = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammas, flat);
            foreach (var ev in B.Evd().EigenValues) yield return ev;
        }
    }

    private static (double, double, double) ReadGapSides(ComplexMatrix H, int n, double gStar)
    {
        const double eps = 1e-6;
        var below = SlowestNonzero(H, n, gStar * (1 - eps));
        var above = SlowestNonzero(H, n, gStar * (1 + eps));
        return (-below.Real - 2 * gStar * (1 - eps), 2 * gStar * (1 + eps) + above.Real, Math.Abs(above.Imaginary));
    }

    private static Complex SlowestNonzero(ComplexMatrix H, int n, double gamma)
    {
        Complex best = default; double bestRe = double.PositiveInfinity;
        foreach (var ev in AllEigenvalues(H, n, gamma))
        {
            if (ev.Magnitude <= 1e-9) continue;
            if (Math.Abs(ev.Real) < bestRe) { bestRe = Math.Abs(ev.Real); best = ev; }
        }
        return best;
    }

    // ---------------------------------------------------------------- the pinned table (the gate's, x = γ/J)

    private static BigInteger[] Poly(params long[] ascending) => ascending.Select(c => new BigInteger(c)).ToArray();

    private static readonly BigInteger[] Deg12 = Poly(256, 0, 1280, 0, -2240, 0, -1696, 0, 68, 0, 132, 0, 9);
    private static readonly BigInteger[] Deg8 = Poly(384, 0, -656, 0, -84, 0, 36, 0, 3);
    private static readonly BigInteger[] Deg6 = Poly(-16, 0, -20, 0, 4, 0, 3);
    private static readonly BigInteger[] Q4a = Poly(-4, 0, 4, 0, 1);
    private static readonly BigInteger[] Q4b = Poly(-16, 0, 4, 0, 1);
    private static readonly BigInteger[] Q4c = Poly(-4, 0, -4, 0, 1);

    /// <summary>(N, graph) → p → the pinned (minimal polynomial, multiplicity in γ) entries, one per point of E in
    /// that block (a polynomial with two positive roots that are both points appears twice).</summary>
    private static readonly Dictionary<(int, string), Dictionary<int, (BigInteger[] Poly, int Multiplicity)[]>> PinnedTable = new()
    {
        [(2, "chain")] = new() { [1] = new[] { (Poly(-2, 1), 1) } },
        [(3, "chain")] = new()
        {
            [1] = new[] { (Poly(-4, 0, 1, 0, 1), 1), (Poly(-3, 0, 1), 1) },
            [2] = new[] { (Poly(-4, 0, 1, 0, 1), 1), (Poly(-3, 0, 1), 1) },
        },
        [(3, "complete")] = new() { [1] = new[] { (Poly(-3, 0, 1), 2) }, [2] = new[] { (Poly(-3, 0, 1), 2) } },
        [(4, "chain")] = new()
        {
            [1] = new[] { (Q4a, 1), (Q4b, 1), (Q4c, 1) },
            [2] = new[] { (Deg12, 1), (Deg8, 1), (Deg6, 1), (Deg12, 1), (Deg8, 1) },
            [3] = new[] { (Q4a, 1), (Q4b, 1), (Q4c, 1) },
        },
        [(4, "ring")] = new()
        {
            [1] = new[] { (Poly(-2, 1), 2) },
            [2] = new[] { (Poly(-64, 0, 28, 0, 3), 1), (Poly(-16, 0, 3), 1) },
            [3] = new[] { (Poly(-2, 1), 2) },
        },
        [(4, "star")] = new()
        {
            [1] = new[] { (Poly(-3, 0, 6, 0, 1), 2) },
            [2] = new[] { (Poly(-4, 0, 3), 2), (Poly(-3, 0, 1), 2), (Poly(-16, 0, 3), 1) },
            [3] = new[] { (Poly(-3, 0, 6, 0, 1), 2) },
        },
        [(4, "complete")] = new() { [2] = new[] { (Poly(-16, 0, 3), 3) } },
    };

    // ---------------------------------------------------------------- IInspectable

    public string DisplayName => $"F50 exceptional couplings, N={N} {Topology}";

    public string Summary
    {
        get
        {
            string pts = ExceptionalSet.Count == 0 ? "∅" : "{" + string.Join(", ", ExceptionalSet.Select(r => r.ToString("0.000000", Inv))) + "}";
            string pinned = HasPinnedTable ? (PinnedCoversAllRoots ? "pinned exactly" : "PINNED TABLE NOT CERTIFIED") : "unpinned";
            return $"E({N}, {Topology}) = {pts} as γ/J (Pauli book), {pinned}; count at −2γ: {GenericCount} at γ/J = {GenericGammaOverJ.ToString(Inv)}, " +
                   $"[{string.Join(", ", CountsAtPoints)}] at E" +
                   (HandoverQ is { } q && Topology == "chain" ? $"; J/γ at min E = {q.ToString("0.000000", Inv)}" : "");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            bool leadsOk = DiagonalBlocks.Concat(OffDiagonalEvenBlocks).All(b => b.LeadingMatches);
            yield return new InspectableNode("the polynomials, exactly",
                summary: $"P_p(γ) = det(B_pp(γ) + 2γ·I) ∈ ℤ[i][γ] recovered from dim+1 Berkowitz determinants over ℤ[i] by integer " +
                         $"interpolation on the repo's own blocks (PauliHamiltonian.Bilinear, PerBlockLiouvillianBuilder); " +
                         $"real roots = real roots of gcd(Re P, Im P). The construction check that can fail: the builder's " +
                         $"dissipator diagonal equals −2·Hamming of every cell exactly, on all {DiagonalBlocks.Count + OffDiagonalEvenBlocks.Count} " +
                         $"even blocks (it did, or this object would not exist); the leading coefficient Π(2 − 2·Hamming) of Re P then " +
                         $"follows as det D, read here as " + (leadsOk ? "consistent" : "INCONSISTENT, a Berkowitz or interpolation fault"),
                provenance: NodeProvenance.Live);

            foreach (var b in DiagonalBlocks)
                yield return new InspectableNode($"block ({b.PKet},{b.QBra}), dim {b.Dim}",
                    summary: $"leading {b.RealPart[^1]} = Π(2 − 2·Hamming) {(b.LeadingMatches ? "exactly" : "MISMATCH")}; " +
                             $"root polynomial Q = {(b.ImagPart.Length == 0 ? "Re P (Im P ≡ 0)" : $"gcd(Re P, Im P), degree {b.Coefficients.Length - 1}")}; " +
                             $"positive roots: {b.DistinctPositiveRoots} distinct, {b.PositiveRootsWithMultiplicity} with multiplicity (Sturm); " +
                             (b.Roots.Count == 0 ? "none" : "at γ/J = " + string.Join(", ", b.Roots.Select(r => r.ToString("0.0000000000", Inv)))) +
                             (b.PKet > 0 && b.PKet < N ? $"; C(N,p) − c_p = {Binomial(N, b.PKet) - Debts.Single(d => d.P == b.PKet).Components} (the Schur bound), the shortfall is the Krein debt m_p = {Debts.Single(d => d.P == b.PKet).Debt} (Theorem D)" : ""),
                    provenance: NodeProvenance.Live);

            yield return new InspectableNode("the Krein debt, exactly (Theorem D: #E_p = C(N,p) − c_p − m_p)",
                summary: string.Join("; ", Debts.Select(d =>
                             $"p={d.P}: C(N,p) − c_p = {d.Binomial} − {d.Components}, #E_p = {d.ExceptionalCountWithMultiplicity} (Sturm, with multiplicity), " +
                             $"m_p = {d.Debt} exact, {d.DebtRead} read by the eigensolver at J/γ = 30 and 60, c_p = {d.StationaryRead} read as the stationary modes{(d.RoutesMeet ? "" : " ROUTES DISAGREE")}")) +
                         $"; the 2γ regime {(GapRegimeExists ? "exists: every mode slower than 2γ for γ below min E, and the threshold is 1/min E" : "does NOT exist: a block holds a mode inside the half-plane at every small γ")}",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("the even off-diagonal blocks",
                summary: OffDiagonalEvenBlocks.Count == 0 ? "none at this N" :
                         "no positive root on any of them (distance ≥ 2 cannot reach ⟨n_XY⟩ = 1): " +
                         string.Join("; ", OffDiagonalEvenBlocks.Select(b => $"({b.PKet},{b.QBra}) {b.DistinctPositiveRoots}")),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("E(N, G)",
                summary: ExceptionalSet.Count == 0 ? "empty" :
                         $"{ExceptionalSet.Count} points, γ/J ∈ {{{string.Join(", ", ExceptionalSet.Select(r => r.ToString("0.0000000000", Inv)))}}}",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("the pinned minimal polynomials",
                summary: !HasPinnedTable ? "no pinned table for this (N, graph); the counts above are the live reading" :
                         (PinnedCoversAllRoots ? $"every pinned q^m (table key {(N, Topology) switch { (2, _) => "chain", (3, "ring") => "K₃", (3, "star") => "chain", _ => Topology }}) divides the root polynomial Q_p exactly and q^(m+1) does not, and the pinned positive roots " +
                                                 "exhaust the Sturm count of every block: E is pinned without a float"
                                               : "A PINNED ENTRY FAILS: " + string.Join("; ", Pinned.Where(c => !c.Holds).Select(c => $"p={c.P} q=[{string.Join(",", c.MinimalPolynomial)}]^{c.Multiplicity}"))),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("the multiplicity at the rational points, exactly",
                summary: ExactAtRationalPoints.Count == 0 ? "no rational point in E at this (N, graph); the window count below is the reading" :
                         string.Join("; ", ExactAtRationalPoints.Select(e =>
                             $"γ/J = {e.Numerator}/{e.Denominator}: extra algebraic {e.Algebraic}, geometric {e.Geometric} over ℚ(i) by realification " +
                             $"(generalized kernel of B_pp + 2γ summed over the diagonal blocks){(e.Algebraic > e.Geometric ? ", DEFECTIVE, where the window count is no law" : ", semisimple")}")),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("the count at −2γ (measured)",
                summary: $"eigensolver on every (p,q) block, |λ + 2γ| < {CountWindow.ToString(Inv)}: {GenericCount} at the generic γ/J = {GenericGammaOverJ.ToString(Inv)} " +
                         $"(2N = {2 * N}{(N == 3 && Topology is "complete" or "ring" ? ", K₃ generic 2N+2 = 8" : "")}); at E: [{string.Join(", ", CountsAtPoints)}], " +
                         $"nullity of L + 2γ beside it: [{string.Join(", ", NullitiesAtPoints)}], both evaluated AT the rational where the point is one; " +
                         $"count = nullity reads semisimple and makes the window a law (O(ε‖L‖) solver motion plus O(δ) from the 2^−41 reading), " +
                         $"count > nullity reads defective and the window is a reading there; count < nullity would be a defective point evaluated off its position and is refused",
                provenance: NodeProvenance.Live);

            if (Topology == "chain")
                yield return new InspectableNode("the handover: min E(chain) = 1/Q*_gap",
                    summary: GapSides is { } g
                        ? $"J/γ at the smallest point = {HandoverQ!.Value.ToString("0.0000000", Inv)}; the spectrum on both sides: gap − 2γ = {g.BelowOffset.ToString("0.0e0", Inv)} just below, " +
                          $"2γ − gap = {g.AboveOffset.ToString("0.0e0", Inv)} just above with |Im| = {g.AboveImag.ToString("0.0e0", Inv)}: " +
                          (GapSidesHold ? "gap 2γ below, a real mode below 2γ above (the two-route reading beside the theorem: a gap below 2γ is a real mode of a diagonal block)" : "THE GAP READING FAILS")
                        : "E is empty",
                    provenance: NodeProvenance.Live);
        }
    }

    private static long Binomial(int n, int k)
    {
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }
}
