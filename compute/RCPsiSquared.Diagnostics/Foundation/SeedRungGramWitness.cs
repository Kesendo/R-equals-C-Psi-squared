using System.Globalization;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F143, recomputed at inspect time: the seed-rung Gram matrix G = WᵀW, its closed form
/// (1/M)·(𝟏𝟏ᵀ + (I + R)/2), and its whole spectrum with multiplicities. Two halves that meet, and they are
/// kept apart on purpose because only one of them can be exact.
///
/// <para><b>The exact half, integer, compared to 0.</b> Ĝ := 2M·G has integer entries
/// (<see cref="SeedRungGramClaim.ScaledGram"/>), so the three eigenrelations are integer identities and
/// <see cref="ExactResidual"/> evaluates them in <c>long</c> arithmetic: Ĝ·(e_k − e_{M−k}) = 0 on the ⌊N/2⌋
/// chiral differences, Ĝ·(o_i − o_j) = 2·(o_i − o_j) on the ⌈N/2⌉ − 1 differences of R-orbit vectors, and
/// Ĝ·𝟏 = 2M·𝟏. Those N vectors exhaust the space and <see cref="CertifiedRank"/> shows them independent by
/// an exact GF(p) rank, so the spectrum is certified outright: no eigensolver, and nothing to tolerate.</para>
///
/// <para><b>The measured half, and its law rather than its threshold.</b> That the DST-I modes
/// v_k(l) = √(2/M)·sin(πkl/M) produce that matrix is where floating point enters and cannot be argued away.
/// <see cref="ClosedFormResidual"/> rebuilds W from the sines, forms WᵀW, and reads
/// max |2M·(WᵀW)_{ac} − Ĝ_{ac}|; <see cref="ResidualErrorRatio"/> divides that by ε·M. The law is the
/// FLATNESS of that ratio across N, not any one value of it: it reaches 2.09 over N = 2..40 and does not
/// grow, while a fixed tolerance that happens to pass says nothing about whether the identity is the reason
/// it passed. Two mechanisms feed the residual, the reduction of sine arguments running up to ≈πN and the
/// summation of N terms; which dominates is not settled here, and the gate does not depend on it.</para>
///
/// <para><b>The mode ORDER cannot bias either half, between the two orders anyone uses.</b> Reversing the
/// mode order is exactly R, and R Ĝ R = Ĝ, so ascending and descending in energy give the same matrix and
/// the eigensolver's convention in the fence below is not a hidden knob. That is a statement about the
/// REVERSAL and not about relabellings in general: a generic permutation moves the 3s off the diagonal and
/// the anti-diagonal and does change the matrix.</para>
///
/// <para><b>The uniform-hopping fence, and it is narrower than its usual phrasing.</b>
/// <see cref="ReadHopping"/> puts bond disorder on the hopping matrix, diagonalises it, and reads SEVERAL
/// numbers rather than one, because they do not fail together. On a profile with NO ZERO BOND the chiral
/// symmetry of W and G's annihilation of the R-odd sector stay at the rounding level, the first non-kernel
/// eigenvalue stays away from zero, and only the middle band's common value 1/M splits, taking the
/// entry-wise closed form with it. λ_max sits at 1 throughout and is reported rather than used as a
/// control, for the reason given at <see cref="ReadHopping"/>. The F143 entry's prose says the closed form
/// wants uniform hopping, which is true; this reading is what shows how much of the statement that leaves
/// untouched.</para>
///
/// <para>Live at <c>inspect --root seedrung</c>. Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F143 +
/// <c>docs/proofs/PROOF_FROZEN_BAND_SO4.md</c> §6 + <c>docs/proofs/PROOF_R90_FROZEN_DIVISOR.md</c> Lemma 5
/// + <c>simulations/eta_ceiling_reduction.py</c> block V5.</summary>
public sealed class SeedRungGramWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The chain length the detail rows are read at. The sweep behind the summary runs
    /// N = 2..<see cref="SweepMax"/> regardless, because the point of the residual is its behaviour ACROSS
    /// N and a single N cannot show it.</summary>
    public int N { get; }

    /// <summary>The largest N the residual sweep visits. 40 matches the Python gate's <c>--deep</c> reach
    /// in block V5, and the whole sweep is N² work at each step, so it costs nothing.</summary>
    public const int SweepMax = 40;

    public SeedRungGramWitness(int n = 9)
    {
        if (n < 2 || n > 200)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                "N ≥ 2 for a kernel to exist and a gap to sit above it; the upper bound is only the " +
                "detail listing's, the closed form has no N ceiling");
        N = n;
    }

    /// <summary>Maximum absolute residual of the three eigenrelations, in exact integer arithmetic.
    /// Zero, and it is a genuine zero rather than a rounded one.</summary>
    public static long ExactResidual(int n)
    {
        var g = SeedRungGramClaim.ScaledGram(n);
        int m = SeedRungGramClaim.TransformLength(n);
        long worst = 0;

        // (a) the kernel: one chiral difference per 2-cycle of R.
        foreach (var (k, kbar) in ChiralPairs(n))
        {
            var x = new long[n];
            x[k - 1] = 1; x[kbar - 1] = -1;
            worst = Math.Max(worst, MaxAbs(Apply(g, x)));
        }

        // (b) the middle band: differences of R-orbit vectors, each R-even and summing to zero.
        var orbits = OrbitVectors(n);
        for (int i = 1; i < orbits.Count; i++)
        {
            var d = new long[n];
            for (int a = 0; a < n; a++) d[a] = orbits[i][a] - orbits[0][a];
            var image = Apply(g, d);
            for (int a = 0; a < n; a++) image[a] -= 2 * d[a];
            worst = Math.Max(worst, MaxAbs(image));
        }

        // (c) the top: the all-ones vector at 2M.
        var ones = new long[n];
        Array.Fill(ones, 1L);
        var top = Apply(g, ones);
        for (int a = 0; a < n; a++) top[a] -= 2L * m;
        worst = Math.Max(worst, MaxAbs(top));

        return worst;
    }

    /// <summary>The eigenvectors the exact half uses, counted: ⌊N/2⌋ kernel differences, ⌈N/2⌉ − 1 orbit
    /// differences, one all-ones vector. They must sum to N or the certification has a hole, and this
    /// returns the sum so the caller can say so.</summary>
    public static int CertifiedDimension(int n)
    {
        var (frozen, middle, top) = SeedRungGramClaim.Multiplicities(n);
        return frozen + middle + top;
    }

    /// <summary>The N certified eigenvectors, in the order the exact half builds them: the chiral
    /// differences, then the orbit differences, then 𝟏. Integer, so nothing here rounds.</summary>
    public static IReadOnlyList<long[]> CertifiedVectors(int n)
    {
        var vectors = new List<long[]>();
        foreach (var (k, kbar) in ChiralPairs(n))
        {
            var x = new long[n];
            x[k - 1] = 1; x[kbar - 1] = -1;
            vectors.Add(x);
        }
        var orbits = OrbitVectors(n);
        for (int i = 1; i < orbits.Count; i++)
        {
            var d = new long[n];
            for (int a = 0; a < n; a++) d[a] = orbits[i][a] - orbits[0][a];
            vectors.Add(d);
        }
        var ones = new long[n];
        Array.Fill(ones, 1L);
        vectors.Add(ones);
        return vectors;
    }

    /// <summary>The rank of <see cref="CertifiedVectors"/>, computed EXACTLY over GF(p). Counting N
    /// relations is not certifying a spectrum: if two of them were dependent, the multiplicities would be
    /// asserted rather than shown. Integer inputs of size ≤ 2, so a single prime is decisive and a float
    /// rank would be the wrong tool here for the usual reason — it would need a threshold.</summary>
    public static int CertifiedRank(int n)
    {
        const long p = 1_000_003L;
        var rows = CertifiedVectors(n).Select(v => v.Select(x => ((x % p) + p) % p).ToArray()).ToList();

        int rank = 0;
        for (int col = 0; col < n && rank < rows.Count; col++)
        {
            int pivot = -1;
            for (int r = rank; r < rows.Count; r++)
                if (rows[r][col] != 0) { pivot = r; break; }
            if (pivot < 0) continue;

            (rows[rank], rows[pivot]) = (rows[pivot], rows[rank]);
            long inv = ModInverse(rows[rank][col], p);
            for (int c = col; c < n; c++) rows[rank][c] = rows[rank][c] * inv % p;
            for (int r = 0; r < rows.Count; r++)
            {
                if (r == rank || rows[r][col] == 0) continue;
                long factor = rows[r][col];
                for (int c = col; c < n; c++)
                    rows[r][c] = ((rows[r][c] - factor * rows[rank][c]) % p + p) % p;
            }
            rank++;
        }
        return rank;
    }

    private static long ModInverse(long a, long p)
    {
        long result = 1, e = p - 2, b = a % p;
        while (e > 0)
        {
            if ((e & 1) == 1) result = result * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return result;
    }

    /// <summary>max |2M·(WᵀW)_{ac} − Ĝ_{ac}| with W rebuilt from the DST-I sines: the one statement of this
    /// entry that floating point has to carry.</summary>
    public static double ClosedFormResidual(int n) => ResidualAgainstClosedForm(ModeAmplitudeSquares(n), n);

    /// <summary>The residual divided by ε·M, the error model. Reported instead of gated on a number: what
    /// a gate should assert is that this ratio stays bounded across N, not that some residual is below a
    /// constant that happens to hold at one N.</summary>
    public static double ResidualErrorRatio(int n) =>
        ClosedFormResidual(n) / (Eps * SeedRungGramClaim.TransformLength(n));

    /// <summary>How far W is from doubly stochastic, max_l |Σ_k W_{lk} − 1|. Exactly 1 in the algebra
    /// (Σ_k (2/M)sin²(πkl/M) = (N + 1)/M), so this is a second ε-level read on the same rebuild and it
    /// isolates the sines from the squaring.</summary>
    public static double DoublyStochasticResidual(int n)
    {
        var w = ModeAmplitudeSquares(n);
        double worst = 0.0;
        for (int i = 0; i < n; i++)
        {
            double rowSum = 0.0, columnSum = 0.0;
            for (int j = 0; j < n; j++) { rowSum += w[i, j]; columnSum += w[j, i]; }
            worst = Math.Max(worst, Math.Max(Math.Abs(rowSum - 1.0), Math.Abs(columnSum - 1.0)));
        }
        return worst;
    }

    /// <summary>What a bond profile does to the separate statements of this entry, read off the modes an
    /// eigensolver returns rather than off the sine formula. Reported as several numbers because they do
    /// NOT stand or fall together, which is the finding this reading exists to carry.
    ///
    /// <para><b>The hypothesis, which is not "any profile".</b> K h K = −h pairs the eigenVALUES at ±ε;
    /// concluding v_{M−k} ∝ K v_k, and with it W R = W, needs the eigenVECTOR to be fixed by its
    /// eigenvalue, i.e. a SIMPLE spectrum. On an open nearest-neighbour chain with every bond nonzero h is
    /// a Jacobi matrix and its spectrum is simple, so the hypothesis is exactly "no zero bond". It is not
    /// decoration: at N = 6 with bonds [1,1,0,1,1] the chain falls into two identical trimers and EVERY
    /// eigenvalue of h doubles. There the question is not merely false but ILL-POSED, because W stops being
    /// a function of h: MathNet returns trimer-localised vectors and reads the residuals as O(1), while the
    /// equally legal basis (u^A ± u^B)/√2 of the same eigenspaces reads them as exactly 0. Both are correct
    /// answers to a question that has none, so the test asserts the DEGENERACY, which is a function of the
    /// eigenvalues alone, and then DEMONSTRATES the basis-dependence by building the second basis rather
    /// than pinning either number. Two
    /// further hypotheses ride along and are the chain's, not this reading's: nearest-neighbour BIPARTITE
    /// hopping with NO on-site terms (a longitudinal field breaks K h K = −h outright) and an OPEN
    /// boundary, an even ring being bipartite but carrying doubly degenerate modes, which is the same
    /// failure as a zero bond.</para>
    ///
    /// <para><see cref="HoppingReading.ChiralSymmetryResidual"/> is max |W_{l,k} − W_{l,M−k}|, the
    /// involution as a MEASUREMENT rather than as an index map.
    /// <see cref="HoppingReading.ROddAnnihilationResidual"/> is max ‖G x‖ over the ⌊N/2⌋ chiral differences
    /// x = e_k − e_{M−k}, which reads the kernel's IDENTITY and not merely its size: W R = W forces
    /// G R = G, so every R-odd vector is annihilated. <see cref="HoppingReading.FirstNonKernelEigenvalue"/>
    /// closes the other side, since an annihilation test alone certifies nullity ≥ ⌊N/2⌋ and never ≤; it is
    /// the (⌊N/2⌋+1)-th eigenvalue of G and must stay bounded AWAY from zero.
    /// <see cref="HoppingReading.MiddleBandSpread"/> is the one that MOVES: the ⌈N/2⌉ − 1 middle
    /// eigenvalues share the value 1/M only on the uniform chain, and disorder splits them.
    /// <see cref="HoppingReading.TopEigenvalue"/> is reported and deliberately NOT used as a control: W is
    /// doubly stochastic for the squared amplitudes of ANY orthonormal basis, so λ_max = 1 cannot fail on
    /// physics and would only ever catch the eigensolver. What double stochasticity does NOT give is
    /// SIMPLICITY of that 1, which is an irreducibility statement and fails on the disconnected profile
    /// above; this reading does not certify it and the prose here does not claim it.
    /// <see cref="HoppingReading.HoppingSpectralGap"/> is the smallest spacing of h's own eigenvalues,
    /// carried because it is the error model: an eigenvector's backward error is O(ε‖h‖/gap), not O(ε‖h‖),
    /// and at bonds [1,1,1e-10,1,1] (all nonzero, so inside the scope) the chiral residual reaches 6.7e-6.
    /// A constant times ε·N would be a number wearing an error model's clothes.</para>
    ///
    /// <para>So the fence is narrower than "the closed form needs uniform hopping" makes it sound. On a
    /// chain with no zero bond the frozen COUNT and the kernel's identity as the R-odd sector are the
    /// involution's and survive; what uniform hopping buys is the middle band's common VALUE, and with it
    /// the entry-wise closed form. That matters downstream, because the reading F144 and F146 consume is
    /// the count. One clause against reading that as more than it is: the count survives, but its MOAT is
    /// no longer GUARANTEED. Under disorder the first non-kernel eigenvalue is not bounded below by 1/M any
    /// more; under J_l = 1 + 0.3·sin(l+1) it reads 0.074 against a uniform 0.100 at N = 9 and 0.171 against
    /// 0.167 at N = 5. So the BOUND goes, not the margin, and a NUMERIC rank read of that count is
    /// protected only on the uniform chain.</para>
    /// </summary>
    public static HoppingReading ReadHopping(int n, IReadOnlyList<double> bonds)
    {
        ArgumentNullException.ThrowIfNull(bonds);
        if (bonds.Count != n - 1)
            throw new ArgumentException(
                $"an open chain of {n} sites has {n - 1} bonds; {bonds.Count} given", nameof(bonds));

        var h = Matrix<double>.Build.Dense(n, n);
        for (int l = 0; l < n - 1; l++)
        {
            h[l, l + 1] = bonds[l];
            h[l + 1, l] = bonds[l];
        }

        var (v, hSpectrum) = Diagonalise(h, n);

        double hoppingGap = double.PositiveInfinity;
        for (int i = 1; i < n; i++) hoppingGap = Math.Min(hoppingGap, hSpectrum[i] - hSpectrum[i - 1]);

        var w = SquaredAmplitudes(n, v);
        double chiral = ChiralSymmetryOfModes(n, w);

        var g = new double[n, n];
        for (int a = 0; a < n; a++)
            for (int c = 0; c < n; c++)
            {
                double entry = 0.0;
                for (int l = 0; l < n; l++) entry += w[l, a] * w[l, c];
                g[a, c] = entry;
            }

        // The kernel's IDENTITY, read directly: G applied to each chiral difference, not the size of the
        // lowest eigenvalues. A small eigenvalue says a kernel is there; only this says it is the R-odd one.
        double annihilation = 0.0;
        foreach (var (k, kbar) in ChiralPairs(n))
            for (int a = 0; a < n; a++)
                annihilation = Math.Max(annihilation, Math.Abs(g[a, k - 1] - g[a, kbar - 1]));

        var gMatrix = Matrix<double>.Build.DenseOfArray(g);
        var spectrum = gMatrix.Evd().EigenValues.Select(z => z.Real).OrderBy(x => x).ToArray();
        int frozen = SeedRungGramClaim.FrozenNullity(n);
        int middle = SeedRungGramClaim.MiddleBandMultiplicity(n);

        double firstNonKernel = frozen < n ? spectrum[frozen] : double.NaN;
        double middleSpread = middle == 0
            ? 0.0
            : spectrum.Skip(frozen).Take(middle).Max() - spectrum.Skip(frozen).Take(middle).Min();

        return new HoppingReading(chiral, annihilation, firstNonKernel, middleSpread, spectrum[^1],
            ResidualAgainstClosedForm(w, n), hoppingGap);
    }

    /// <summary>The chain's single-excitation matrix, diagonalised, with the eigenvectors put in ASCENDING
    /// eigenvalue order explicitly rather than trusting the solver's convention. The chiral pairing
    /// k ↔ M − k only means anything in a monotone order, and the ambiguity that remains, ascending versus
    /// descending, is the reversal R, which the closed form commutes with.</summary>
    public static (double[,] Modes, double[] Spectrum) DiagonaliseChain(int n, IReadOnlyList<double> bonds)
    {
        ArgumentNullException.ThrowIfNull(bonds);
        if (bonds.Count != n - 1)
            throw new ArgumentException(
                $"an open chain of {n} sites has {n - 1} bonds; {bonds.Count} given", nameof(bonds));

        var h = Matrix<double>.Build.Dense(n, n);
        for (int l = 0; l < n - 1; l++)
        {
            h[l, l + 1] = bonds[l];
            h[l + 1, l] = bonds[l];
        }
        return Diagonalise(h, n);
    }

    private static (double[,] Modes, double[] Spectrum) Diagonalise(Matrix<double> h, int n)
    {
        var evd = h.Evd();
        var order = Enumerable.Range(0, n)
            .OrderBy(k => evd.EigenValues[k].Real)
            .ToArray();

        var modes = new double[n, n];
        var spectrum = new double[n];
        for (int k = 0; k < n; k++)
        {
            spectrum[k] = evd.EigenValues[order[k]].Real;
            for (int l = 0; l < n; l++) modes[l, k] = evd.EigenVectors[l, order[k]];
        }
        return (modes, spectrum);
    }

    /// <summary>W_{lk} = v_k(l)² from a mode matrix, so a caller can supply an eigenbasis of its own. That
    /// is not a convenience: on a DEGENERATE h the eigenbasis is not determined, so W is not a function of
    /// h, and the only way to say that rather than assert it is to build two valid bases and compare.
    /// </summary>
    public static double[,] SquaredAmplitudes(int n, double[,] modes)
    {
        ArgumentNullException.ThrowIfNull(modes);
        var w = new double[n, n];
        for (int l = 0; l < n; l++)
            for (int k = 0; k < n; k++)
                w[l, k] = modes[l, k] * modes[l, k];
        return w;
    }

    /// <summary>max |W_{l,k} − W_{l,M−k}|, the chiral symmetry of a given W.</summary>
    public static double ChiralSymmetryOfModes(int n, double[,] w)
    {
        ArgumentNullException.ThrowIfNull(w);
        double worst = 0.0;
        for (int l = 0; l < n; l++)
            for (int k = 1; k <= n; k++)
                worst = Math.Max(worst, Math.Abs(w[l, k - 1] - w[l, SeedRungGramClaim.ChiralPartner(n, k) - 1]));
        return worst;
    }

    /// <summary>The separate statements under a given bond profile. See <see cref="ReadHopping"/> for what
    /// each one is, which of them can fail on physics, and which hypothesis they all sit under.</summary>
    public readonly record struct HoppingReading(
        double ChiralSymmetryResidual,
        double ROddAnnihilationResidual,
        double FirstNonKernelEigenvalue,
        double MiddleBandSpread,
        double TopEigenvalue,
        double ClosedFormResidual,
        double HoppingSpectralGap);

    /// <summary>The entry-wise deviation from the closed form under a bond profile, i.e.
    /// <see cref="HoppingReading.ClosedFormResidual"/>, on the 2M scaling. Kept as its own name because it
    /// is the number the F143 entry's prose refers to.</summary>
    public static double DisorderedResidual(int n, IReadOnlyList<double> bonds) =>
        ReadHopping(n, bonds).ClosedFormResidual;

    /// <summary>A fixed, deterministic bond profile for the disorder fence: J_l = 1 + 0.3·sin(l + 1). No
    /// randomness, so the number in the inspect view is the number a test asserts.</summary>
    public static double[] DefaultDisorder(int n)
    {
        var bonds = new double[n - 1];
        for (int l = 0; l < n - 1; l++) bonds[l] = 1.0 + 0.3 * Math.Sin(l + 1);
        return bonds;
    }

    private const double Eps = 2.220446049250313e-16;

    /// <summary>W_{lk} = v_k(l)² = (2/M)·sin²(πkl/M), sites and modes both one-based, stored zero-based.
    /// </summary>
    private static double[,] ModeAmplitudeSquares(int n)
    {
        int m = SeedRungGramClaim.TransformLength(n);
        var w = new double[n, n];
        for (int l = 1; l <= n; l++)
            for (int k = 1; k <= n; k++)
            {
                double amp = Math.Sqrt(2.0 / m) * Math.Sin(Math.PI * k * l / m);
                w[l - 1, k - 1] = amp * amp;
            }
        return w;
    }

    private static double ResidualAgainstClosedForm(double[,] w, int n)
    {
        int m = SeedRungGramClaim.TransformLength(n);
        var g = SeedRungGramClaim.ScaledGram(n);
        double worst = 0.0;
        for (int a = 0; a < n; a++)
            for (int c = 0; c < n; c++)
            {
                double entry = 0.0;
                for (int l = 0; l < n; l++) entry += w[l, a] * w[l, c];
                worst = Math.Max(worst, Math.Abs(2.0 * m * entry - g[a, c]));
            }
        return worst;
    }

    /// <summary>The 2-cycles of R on the modes, one-based, each listed once with k &lt; M − k. There are
    /// ⌊N/2⌋ of them, and the self-paired middle mode of an odd chain is not among them.</summary>
    public static IReadOnlyList<(int Mode, int Partner)> ChiralPairs(int n)
    {
        var pairs = new List<(int, int)>();
        for (int k = 1; k <= n; k++)
        {
            int partner = SeedRungGramClaim.ChiralPartner(n, k);
            if (k < partner) pairs.Add((k, partner));
        }
        return pairs;
    }

    /// <summary>One integer vector per R-orbit, each with entry sum 2 so that differences of them are
    /// R-even AND orthogonal to 𝟏: a chiral pair gives e_k + e_{M−k}, the self-paired middle mode gives
    /// 2·e_{M/2}. There are ⌈N/2⌉ of them and their differences span the middle band.</summary>
    public static IReadOnlyList<long[]> OrbitVectors(int n)
    {
        var orbits = new List<long[]>();
        foreach (var (k, kbar) in ChiralPairs(n))
        {
            var o = new long[n];
            o[k - 1] = 1; o[kbar - 1] = 1;
            orbits.Add(o);
        }
        if (n % 2 == 1)
        {
            int middle = SeedRungGramClaim.TransformLength(n) / 2;
            var o = new long[n];
            o[middle - 1] = 2;
            orbits.Add(o);
        }
        return orbits;
    }

    private static long[] Apply(long[,] g, long[] x)
    {
        int n = x.Length;
        var y = new long[n];
        for (int a = 0; a < n; a++)
        {
            long sum = 0;
            for (int c = 0; c < n; c++) sum += g[a, c] * x[c];
            y[a] = sum;
        }
        return y;
    }

    private static long MaxAbs(long[] x)
    {
        long worst = 0;
        foreach (var v in x) worst = Math.Max(worst, Math.Abs(v));
        return worst;
    }

    /// <summary>One N of the residual sweep.</summary>
    public readonly record struct SweepRow(int N, long Exact, double Rebuild, double Ratio);

    /// <summary>The residual across N = 2..<see cref="SweepMax"/>. Computed once per instance.</summary>
    public IReadOnlyList<SweepRow> Sweep => _sweep ??= BuildSweep();

    private IReadOnlyList<SweepRow>? _sweep;

    private static IReadOnlyList<SweepRow> BuildSweep()
    {
        var rows = new List<SweepRow>();
        for (int n = 2; n <= SweepMax; n++)
            rows.Add(new SweepRow(n, ExactResidual(n), ClosedFormResidual(n), ResidualErrorRatio(n)));
        return rows;
    }

    public string DisplayName =>
        $"F143 seed-rung Gram matrix, recomputed: spec(G) certified in integers at N={N}, rebuild residual " +
        $"swept N = 2..{SweepMax}";

    public string Summary
    {
        get
        {
            var rows = Sweep;
            long worstExact = rows.Max(r => r.Exact);
            double worstRatio = rows.Max(r => r.Ratio);
            var (frozen, middle, _) = SeedRungGramClaim.Multiplicities(N);
            string gapReading = middle > 0
                ? $"gap {SeedRungGramClaim.Gap(N).ToString("0.######", Inv)} = 1/M"
                : $"gap {SeedRungGramClaim.Gap(N).ToString("0.######", Inv)}, the middle band being EMPTY " +
                  "here, so it runs straight to the simple 1 and is not 1/M";
            return $"at N={N}: kernel {frozen} = ⌊N/2⌋, middle band {middle} = ⌈N/2⌉−1, top 1, " +
                   $"total {CertifiedDimension(N)} = N; {gapReading}. " +
                   $"All three eigenrelations exact in integers across N = 2..{SweepMax}, worst residual {worstExact}. " +
                   $"The one measured statement, WᵀW against the closed form, stays at " +
                   $"{worstRatio.ToString("0.##", Inv)}·ε·M at worst over the same range";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var rows = Sweep;

            yield return new InspectableNode("the closed form",
                summary: $"2M·G = 2·𝟏𝟏ᵀ + I + R at M = {SeedRungGramClaim.TransformLength(N)}; entries 2 " +
                         "generically, 3 on the diagonal and the anti-diagonal, 4 where they meet at the " +
                         "self-paired middle mode of an odd chain");

            yield return new InspectableNode("the kernel, listed",
                summary: string.Join(", ", ChiralPairs(N).Select(p => $"e_{p.Mode} − e_{p.Partner}")) +
                         (N % 2 == 1
                             ? $"; mode {SeedRungGramClaim.TransformLength(N) / 2} is self-paired and contributes nothing"
                             : "; N even, so every mode is in a 2-cycle"));

            yield return InspectableNode.RealScalar("gap above the frozen eigenvalue", SeedRungGramClaim.Gap(N));

            yield return new InspectableNode("the exact half",
                summary: $"Ĝ·(e_k − e_{{M−k}}) = 0, Ĝ·(o_i − o_j) = 2(o_i − o_j), Ĝ·𝟏 = 2M·𝟏, all in long " +
                         $"arithmetic; worst residual over N = 2..{SweepMax} is {rows.Max(r => r.Exact)}, " +
                         "compared to 0 with no tolerance because there is none to give. The N relations " +
                         $"are independent, GF(p) rank {CertifiedRank(N)} of {N} at N={N}, so the " +
                         "multiplicities are shown and not merely counted");

            yield return new InspectableNode("the measured half",
                summary: $"W rebuilt from the sines, WᵀW read against the closed form: " +
                         $"{ClosedFormResidual(N).ToString("0.###e+0", Inv)} at N={N}, " +
                         $"{ResidualErrorRatio(N).ToString("0.##", Inv)}·ε·M; the doubly-stochastic row and " +
                         $"column sums deviate by {DoublyStochasticResidual(N).ToString("0.###e+0", Inv)}");

            var disordered = ReadHopping(N, DefaultDisorder(N));
            yield return new InspectableNode("the uniform-hopping fence, split",
                summary: $"with bond disorder J_l = 1 + 0.3·sin(l+1) at N={N}, every bond nonzero so the " +
                         $"h-spectrum stays simple (its smallest spacing is " +
                         $"{disordered.HoppingSpectralGap.ToString("0.####", Inv)}, which is the error model " +
                         "these residuals live on, an eigenvector's backward error being O(ε‖h‖/gap)). " +
                         $"WHAT SURVIVES: chiral symmetry of W {disordered.ChiralSymmetryResidual.ToString("0.###e+0", Inv)}, " +
                         $"G annihilating every chiral difference {disordered.ROddAnnihilationResidual.ToString("0.###e+0", Inv)} " +
                         "(the kernel's IDENTITY, not just its size), and the first NON-kernel eigenvalue " +
                         $"still away from zero at {disordered.FirstNonKernelEigenvalue.ToString("0.####", Inv)}, " +
                         "which closes the other side of the count. " +
                         $"WHAT BREAKS: the middle-band spread is {disordered.MiddleBandSpread.ToString("0.####", Inv)} " +
                         $"against a uniform-chain gap of {SeedRungGramClaim.Gap(N).ToString("0.####", Inv)}, and " +
                         $"the closed-form residual {disordered.ClosedFormResidual.ToString("0.####", Inv)} on 2M·G, i.e. " +
                         $"{(disordered.ClosedFormResidual / (2.0 * SeedRungGramClaim.TransformLength(N))).ToString("0.####", Inv)} " +
                         $"on G (the uniform chain reads {ClosedFormResidual(N).ToString("0.###e+0", Inv)} there). " +
                         $"NOT A CONTROL: λ_max = {disordered.TopEigenvalue.ToString("0.######", Inv)} is reported " +
                         "and cannot fail on physics, W being doubly stochastic for the squared amplitudes of " +
                         "ANY orthonormal basis; and it is the VALUE, not the multiplicity, so simplicity of " +
                         "the 1 is not certified here" +
                         (SeedRungGramClaim.MiddleBandMultiplicity(N) == 0
                             ? ". At this N the fence is DEGENERATE and shows nothing: there is no middle " +
                               "band to split, and a 2-site chain's modes do not depend on its one bond at all"
                             : string.Empty));

            yield return new InspectableNode("the Heisenberg fence",
                summary: "mode reflection an involution on the open XY chain's modes 1..N at modulus N+1: " +
                         $"{SeedRungGramClaim.ModeReflectionIsAnInvolution(1, N, N + 1)}; on the Heisenberg " +
                         $"chain's modes 0..N−1 at modulus N: " +
                         $"{SeedRungGramClaim.ModeReflectionIsAnInvolution(0, N - 1, N)}, because k = 0 is " +
                         "sent to N, outside the range. THREE things still carry there: the ⌊N/2⌋ count, " +
                         "the gap at M = N, and the chirality in the SITE index, where ker(B Bᵀ) is exactly " +
                         "the R-odd site space on both chains");

            foreach (var r in rows)
                yield return new InspectableNode(
                    $"N = {r.N}",
                    summary: $"exact residual {r.Exact}; rebuild residual {r.Rebuild.ToString("0.###e+0", Inv)} " +
                             $"= {r.Ratio.ToString("0.##", Inv)}·ε·M; multiplicities " +
                             $"{SeedRungGramClaim.Multiplicities(r.N)}");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
