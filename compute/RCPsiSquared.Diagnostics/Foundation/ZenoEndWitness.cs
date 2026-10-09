using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F50's Zeno end, live (<c>inspect --root zeno</c>): Theorem E of
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c> ("The Zeno end and the dispersion", "Which filling lasts
/// longest"), recomputed at inspect time on the repo's own Hamiltonian (<see cref="PauliHamiltonian"/>, Pauli book
/// H = Σ_b (XX + YY + Δ·ZZ), J = 1) and its own Liouvillian blocks (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>
/// on <see cref="SectorBlock.SectorFlatIndices"/>), uniform rates γ = 1, x = J/γ.
///
/// <para><b>What is checked exactly, every block 1 ≤ p ≤ N − 1, each a comparison that can fail.</b>
/// (a) The exclusion Laplacian L_p, built from the bond list, is the ferromagnet: 2·L_p = n_B − (Σ_b XX + YY + ZZ)_p,
/// the right side read off the repo's Pauli matrices, over the integers. (b) The Zeno generator does not see Δ:
/// (A²)_PP, from H's entries as 2·δ_zx·Σ_y |H_yx|² − 2|H_zx|² (the double commutator's diagonal, y ≠ x), equals 8·L_p
/// at the given Δ over the rationals; that the diagonal of H cancels is structural, so this row checks the
/// construction, and the Liouvillian reading below is where the blindness is tested. (c) The lift intertwines:
/// L_p·F_p = F_p·L_site with F_p[S][i] = 1 if i ∈ S, over the integers. (d) The token gap is λ₁: the number of
/// eigenvalues of L_p below a rational just under λ₁ is 1 and below one just over it is 1 + mult(λ₁), read by
/// Jacobi's rule on the leading principal minors of the integer matrix b·L_p − a·I (Bareiss), the two rationals
/// verified the same way on the site Laplacian; the gap is Theorem E (b) on the chain and the star, trees, and on the
/// ring the outside route the proof cites (Caputo, Liggett and Richthammer), read here; the multiplicity at λ₁ is
/// Theorem E (c) on the chain and read on the star and the ring. (e) The detuning weights of Theorem E (d): W_b(p),
/// the sum over the hops across bond b of 2|t|²ΔV²/Γ³ (Γ = 4, the rate of a hop's coherence; t read off H, ΔV the
/// change of the ZZ energy Δ·Σ_b z_a z_b across the hop, computed from the bond list at the exact value of the
/// double Δ, since the diagonal of H in floats is a float sum of terms ±Δ, which can round),
/// equal their closed forms, 4Δ²C(N − 4, p − 2) on a bond with two outer neighbours (chain bulk, every ring bond),
/// Δ²C(N − 2, p − 1)/2 on an end bond of the chain, Δ²(N − 2p)²C(N − 2, p − 1)/2 on every bond of the star.</para>
///
/// <para><b>What is read (floats, each beside its law).</b> δ_p, the detuning form's quotient on the lifted wave,
/// Σ_b W_b(p)(f_a − f_b)²/(‖f‖²·C(N − 2, p − 1)) with f the site Laplacian's lowest nonconstant eigenvector, and the
/// blocks where it is largest; on the Liouvillian, at x = 0.02 and 0.04, the slowest non-stationary rate of every
/// block, its leading coefficient against 2λ₁, the block(s) holding the slowest one, and the split, each block's
/// slowest level less block 1's over x⁴, against δ_p − δ₁, the deviation's ratio per doubling of x printed beside
/// its law 4, the next even order. Levels within 64·eps·‖B‖₁ of each other, the eigensolver's rounding, are read as
/// tied; where the read split is at that rounding no split is read, and where only its deviation from δ_p − δ₁ is,
/// as at small Δ, the split is read as predicted and its law is not read, no ratio being printed. That the hops' own correction is shared by every filling is Theorem E (d) on the chain, the winding
/// argument on the rings N ≥ 5 and, on the star, the gate's exact reading (row Z9); this object reads its
/// consequence, the predicted and the measured block agreeing or not.</para>
///
/// <para>Guards: chain and star N = 4..<see cref="MaxN"/>, ring N = 5..<see cref="MaxN"/> (N = 3 is left out, its
/// two blocks being X^⊗N partners that tie, with no filling to pick; on the ring N = 4 the hops' own correction is not shared, a path of four
/// hops winding around it, and block 2's slowest level is threefold).</para></summary>
public sealed class ZenoEndWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest N the live build admits (the (3,3) block at N = 7 is 1225-dimensional).</summary>
    public const int MaxN = 7;

    /// <summary>The two couplings x = J/γ of the Liouvillian reading.</summary>
    public static readonly double[] ReadingXs = { 0.02, 0.04 };

    public int N { get; }
    public string Topology { get; }
    public double Delta { get; }
    public IReadOnlyList<(int A, int B)> Bonds { get; }

    /// <summary>One block's exact checks and its readings.</summary>
    public sealed record BlockCheck(
        int P, int Dim, bool Ferromagnet, bool ZenoGeneratorIsEightLaplacian, bool LiftIntertwines,
        int CountBelowLow, int CountBelowHigh, bool DetuningClosedForm, double DetuningQuotient,
        IReadOnlyList<double> SlowestRates)
    {
        /// <summary>The slowest rate over x² at the first reading x; its limit is 2λ₁ in every block.</summary>
        public double LeadingCoefficient => SlowestRates[0] / (ReadingXs[0] * ReadingXs[0]);

        /// <summary>(rate/x² − 2λ₁) at the second x over the first: the next even order makes it 4.</summary>
        public double LeadingDeviationRatio(double lambda1) =>
            (SlowestRates[1] / (ReadingXs[1] * ReadingXs[1]) - 2 * lambda1) / (SlowestRates[0] / (ReadingXs[0] * ReadingXs[0]) - 2 * lambda1);
    }

    public IReadOnlyList<BlockCheck> Blocks { get; }
    public double Lambda1 { get; }
    public int Lambda1Multiplicity { get; }
    public BigRational QLow { get; }
    public BigRational QHigh { get; }

    public bool FerromagnetHolds => Blocks.All(b => b.Ferromagnet);
    public bool ZenoGeneratorHolds => Blocks.All(b => b.ZenoGeneratorIsEightLaplacian);
    public bool LiftHolds => Blocks.All(b => b.LiftIntertwines);
    public bool TokenGapHolds => Blocks.All(b => b.CountBelowLow == 1 && b.CountBelowHigh == 1 + Lambda1Multiplicity);
    public bool DetuningClosedFormHolds => Blocks.All(b => b.DetuningClosedForm);
    public bool IsTree => Topology != "ring";

    public IReadOnlyList<int> PredictedSlowestBlocks { get; }
    public IReadOnlyList<int> MeasuredSlowestBlocks { get; }
    /// <summary>The blocks read as slowest lie in the predicted set; a predicted set of more than one block is a tie at
    /// this order, which the next order may break.</summary>
    public bool SelectionAgrees => MeasuredSlowestBlocks.All(PredictedSlowestBlocks.Contains);

    /// <summary>The deviation of the split (each block's slowest level less block 1's, over x⁴) from δ_p − δ₁, max over
    /// p, at the second x over the first: the law is 4.</summary>
    public double DeviationRatio { get; }
    public double DeviationAtFirstX { get; }

    /// <summary>True where no ZZ term detunes a hop (every W_b(p) = 0): the predicted split vanishes identically.</summary>
    public bool PredictedSplitVanishes { get; }

    /// <summary>True where the read split itself, the largest |block p's slowest level less block 1's| at the first x, is
    /// within the eigensolver's rounding, 64·eps·‖B‖₁: no split read.</summary>
    public bool ReadSplitAtRounding { get; }

    /// <summary>True where the read split's deviation from δ_p − δ₁ is within that rounding at the first x (64·eps·‖B‖₁/x⁴
    /// in the split's units): the split is read as predicted, its law not read, and <see cref="DeviationRatio"/> is NaN.</summary>
    public bool DeviationAtRounding { get; }

    public ZenoEndWitness(int n = 5, string? topology = "chain", double delta = 1.0)
    {
        Topology = (topology ?? "chain").ToLowerInvariant();
        if (Topology is not ("chain" or "ring" or "star"))
            throw new ArgumentException($"topology must be chain, ring or star; got '{topology}'", nameof(topology));
        int minN = Topology == "ring" ? 5 : 4;
        if (n < minN || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"{Topology} needs {minN} ≤ N ≤ {MaxN}");
        if (double.IsNaN(delta) || double.IsInfinity(delta))
            throw new ArgumentException("Δ must be finite", nameof(delta));
        N = n;
        Delta = delta;
        Bonds = Topology switch
        {
            "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToList(),
            "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToList(),
            _ => Enumerable.Range(1, n - 1).Select(i => (0, i)).ToList(),
        };

        var H = Hamiltonian(n, Bonds, delta).ToMatrix();
        var Hheis = Hamiltonian(n, Bonds, 1.0).ToMatrix();
        BigRational deltaExact = ExactRational(delta);

        // the site Laplacian, its lowest nonconstant level and the rationals around it, verified exactly
        var site = SiteLaplacian(n, Bonds);
        var siteEvd = Matrix<double>.Build.DenseOfArray(ToDouble(site)).Evd();
        var levels = siteEvd.EigenValues.Select(z => z.Real).ToArray();
        var order = Enumerable.Range(0, n).OrderBy(i => levels[i]).ToArray();
        Lambda1 = levels[order[1]];
        Lambda1Multiplicity = order.Skip(1).Count(i => Math.Abs(levels[i] - Lambda1) <= 1e-9 * Math.Max(1.0, Lambda1));
        QLow = DyadicNear(Lambda1 * (1 - 1e-6), roundUp: false);
        QHigh = DyadicNear(Lambda1 * (1 + 1e-6), roundUp: true);
        if (CountBelow(site, QLow) != 1 || CountBelow(site, QHigh) != 1 + Lambda1Multiplicity)
            throw new InvalidOperationException($"the rationals around λ₁ = {Lambda1.ToString(Inv)} do not bracket it on the site Laplacian");
        var f = siteEvd.EigenVectors.Column(order[1]).ToArray();

        // the Liouvillian reading, per block and x
        var gamma = Enumerable.Repeat(1.0, n).ToArray();
        var slowest = new double[ReadingXs.Length, n - 1];
        var norm = new double[ReadingXs.Length];
        for (int k = 0; k < ReadingXs.Length; k++)
            for (int p = 1; p < n; p++)
            {
                var (rate, blockNorm) = SlowestRate(H, gamma, n, p, ReadingXs[k]);
                slowest[k, p - 1] = rate;
                norm[k] = Math.Max(norm[k], blockNorm);
            }

        var blocks = new List<BlockCheck>();
        bool allWeightsZero = true;
        for (int p = 1; p < n; p++)
        {
            var cf = Configurations(n, p);
            var L = ExclusionLaplacian(n, p, Bonds);
            bool ferro = FerromagnetIdentity(L, Hheis, cf, Bonds.Count);
            bool zeno = ZenoGeneratorIsEightLaplacian(L, H, cf);
            bool lift = LiftIntertwines(L, LiftMatrix(n, p), site);
            int lo = CountBelow(L, QLow) ?? throw new InvalidOperationException($"a zero leading minor at the lower rational, block {p}");
            int hi = CountBelow(L, QHigh) ?? throw new InvalidOperationException($"a zero leading minor at the upper rational, block {p}");
            var W = DetuningWeights(n, p, Bonds, H, deltaExact);
            bool closed = WeightsMatchClosedForm(Topology, n, p, Bonds, W, deltaExact);
            allWeightsZero &= W.Values.All(w => w.IsZero);
            double quotient = DetuningQuotient(n, p, Bonds, W, f);
            var rates = Enumerable.Range(0, ReadingXs.Length).Select(k => slowest[k, p - 1]).ToArray();
            blocks.Add(new BlockCheck(p, cf.Length, ferro, zeno, lift, lo, hi, closed, quotient, rates));
        }
        Blocks = blocks;

        double eps = Math.Pow(2, -52);
        double top = Blocks.Max(b => b.DetuningQuotient);
        // the quotients are double sums of n terms of exact weights: ties within 64 eps of their size
        PredictedSlowestBlocks = Blocks.Where(b => b.DetuningQuotient >= top - 64 * eps * Math.Max(1.0, Math.Abs(top))).Select(b => b.P).ToList();
        // the slowest rate is the smallest -Re; ties within the eigensolver's rounding, 64 eps times the block norm
        double least = Blocks.Min(b => b.SlowestRates[0]);
        MeasuredSlowestBlocks = Blocks.Where(b => b.SlowestRates[0] <= least + 64 * eps * norm[0]).Select(b => b.P).ToList();

        PredictedSplitVanishes = allWeightsZero;
        var quotients = Blocks.Select(b => b.DetuningQuotient).ToArray();
        var dev = SplitDeviations(slowest, quotients, ReadingXs);
        DeviationAtFirstX = dev[0];
        double readSplit = Blocks.Max(b => Math.Abs(b.SlowestRates[0] - Blocks[0].SlowestRates[0]));
        ReadSplitAtRounding = readSplit <= 64 * eps * norm[0];
        DeviationAtRounding = dev[0] <= 64 * eps * norm[0] / Math.Pow(ReadingXs[0], 4);
        DeviationRatio = DeviationAtRounding ? double.NaN : dev[1] / dev[0];
    }

    // ---------------------------------------------------------------- the pieces, public for the tests

    /// <summary>H = Σ_b (XX + YY + Δ·ZZ) from the repo's Pauli terms (Pauli book, J = 1).</summary>
    public static PauliHamiltonian Hamiltonian(int n, IReadOnlyList<(int A, int B)> bonds, double delta)
    {
        var terms = new List<PauliTerm>();
        foreach (var (a, b) in bonds)
        {
            terms.Add(PauliTerm.TwoSite(n, a, PauliLetter.X, b, PauliLetter.X, Complex.One));
            terms.Add(PauliTerm.TwoSite(n, a, PauliLetter.Y, b, PauliLetter.Y, Complex.One));
            if (delta != 0.0) terms.Add(PauliTerm.TwoSite(n, a, PauliLetter.Z, b, PauliLetter.Z, new Complex(delta, 0)));
        }
        return new PauliHamiltonian(n, terms);
    }

    /// <summary>The popcount-p configurations in ascending order; site l is bit N − 1 − l.</summary>
    public static int[] Configurations(int n, int p) =>
        Enumerable.Range(0, 1 << n).Where(x => BitOperations.PopCount((uint)x) == p).ToArray();

    private static int Occupied(int n, int x, int site) => (x >> (n - 1 - site)) & 1;

    /// <summary>The exclusion graph's Laplacian of block p: an edge of weight 1 per bond whose two sites disagree.</summary>
    public static BigInteger[,] ExclusionLaplacian(int n, int p, IReadOnlyList<(int A, int B)> bonds)
    {
        var cf = Configurations(n, p);
        var pos = new Dictionary<int, int>();
        for (int i = 0; i < cf.Length; i++) pos[cf[i]] = i;
        var L = new BigInteger[cf.Length, cf.Length];
        foreach (var (a, b) in bonds)
            foreach (int x in cf)
                if (Occupied(n, x, a) != Occupied(n, x, b))
                {
                    int y = x ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b));
                    L[pos[x], pos[x]] += 1;
                    L[pos[x], pos[y]] -= 1;
                }
        return L;
    }

    public static BigInteger[,] SiteLaplacian(int n, IReadOnlyList<(int A, int B)> bonds)
    {
        var L = new BigInteger[n, n];
        foreach (var (a, b) in bonds) { L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1; }
        return L;
    }

    /// <summary>F_p[S][i] = 1 if site i is occupied in configuration S.</summary>
    public static BigInteger[,] LiftMatrix(int n, int p)
    {
        var cf = Configurations(n, p);
        var F = new BigInteger[cf.Length, n];
        for (int r = 0; r < cf.Length; r++)
            for (int i = 0; i < n; i++) F[r, i] = Occupied(n, cf[r], i);
        return F;
    }

    /// <summary>2·L_p = n_B − (Σ_b XX + YY + ZZ)_p exactly, the right side from the given (Heisenberg) matrix.</summary>
    public static bool FerromagnetIdentity(BigInteger[,] L, ComplexMatrix heisenberg, int[] cf, int bondCount)
    {
        for (int i = 0; i < cf.Length; i++)
            for (int j = 0; j < cf.Length; j++)
            {
                BigRational rhs = (i == j ? new BigRational(bondCount) : new BigRational(0)) - ExactEntry(heisenberg[cf[i], cf[j]]);
                if (rhs != new BigRational(2 * L[i, j])) return false;
            }
        return true;
    }

    /// <summary>(A²)_PP = 8·L_p, with (A²)_PP[z][x] = 2·δ_zx·Σ_{y≠x}|H_yx|² − 2|H_zx|² (z ≠ x) from H's entries.</summary>
    public static bool ZenoGeneratorIsEightLaplacian(BigInteger[,] L, ComplexMatrix H, int[] cf)
    {
        int d = H.RowCount;
        for (int c = 0; c < cf.Length; c++)
        {
            int x = cf[c];
            BigRational diag = 0;
            for (int y = 0; y < d; y++)
                if (y != x) { var h = ExactEntry(H[y, x]); diag += 2 * h * h; }
            for (int r = 0; r < cf.Length; r++)
            {
                BigRational entry;
                if (r == c) entry = diag;
                else { var h = ExactEntry(H[cf[r], x]); entry = -2 * h * h; }
                if (entry != new BigRational(8 * L[r, c])) return false;
            }
        }
        return true;
    }

    /// <summary>L_p·F_p = F_p·L_site over the integers.</summary>
    public static bool LiftIntertwines(BigInteger[,] L, BigInteger[,] F, BigInteger[,] site)
    {
        int m = L.GetLength(0), n = site.GetLength(0);
        for (int r = 0; r < m; r++)
            for (int i = 0; i < n; i++)
            {
                BigInteger left = 0, right = 0;
                for (int k = 0; k < m; k++) left += L[r, k] * F[k, i];
                for (int k = 0; k < n; k++) right += F[r, k] * site[k, i];
                if (left != right) return false;
            }
        return true;
    }

    /// <summary>The number of eigenvalues of the symmetric integer matrix M below the rational q: Jacobi's rule on the
    /// leading principal minors of b·M − a·I (q = a/b), computed by Bareiss over the integers. Null when a leading
    /// minor vanishes, where the rule says nothing.</summary>
    public static int? CountBelow(BigInteger[,] M, BigRational q)
    {
        int n = M.GetLength(0);
        var a = new BigInteger[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i, j] = q.Denominator * M[i, j] - (i == j ? q.Numerator : BigInteger.Zero);
        BigInteger previous = BigInteger.One;
        int changes = 0, previousSign = 1;
        for (int k = 0; k < n; k++)
        {
            BigInteger pivot = a[k, k];                 // the (k+1)-th leading principal minor
            if (pivot.IsZero) return null;
            if (pivot.Sign != previousSign) changes++;
            previousSign = pivot.Sign;
            for (int i = k + 1; i < n; i++)
                for (int j = k + 1; j < n; j++)
                    a[i, j] = (pivot * a[i, j] - a[i, k] * a[k, j]) / previous;
            previous = pivot;
        }
        return changes;
    }

    /// <summary>W_b(p): over the hops e across bond b (each counted once), 2|t_e|²ΔV_e²/Γ³, Γ = 4 at γ = 1; t read off
    /// H (the hop amplitude, exact), ΔV the change of the ZZ energy Δ·Σ_b z_a z_b (z = 1 − 2·occupation) across the hop,
    /// exact in the given Δ.</summary>
    public static Dictionary<(int A, int B), BigRational> DetuningWeights(int n, int p, IReadOnlyList<(int A, int B)> bonds, ComplexMatrix H, BigRational delta)
    {
        var W = bonds.ToDictionary(b => b, _ => new BigRational(0));
        BigRational gamma3 = 64;
        foreach (int x in Configurations(n, p))
            foreach (var (a, b) in bonds)
                if (Occupied(n, x, a) == 1 && Occupied(n, x, b) == 0)
                {
                    int y = x ^ (1 << (n - 1 - a)) ^ (1 << (n - 1 - b));
                    var t = ExactEntry(H[y, x]);
                    var dV = ZzEnergy(n, y, bonds, delta) - ZzEnergy(n, x, bonds, delta);
                    W[(a, b)] += 2 * t * t * dV * dV / gamma3;
                }
        return W;
    }

    /// <summary>Δ·Σ_b z_a z_b of a configuration, z = 1 − 2·occupation, exact.</summary>
    public static BigRational ZzEnergy(int n, int x, IReadOnlyList<(int A, int B)> bonds, BigRational delta)
    {
        int sum = 0;
        foreach (var (a, b) in bonds) sum += (1 - 2 * Occupied(n, x, a)) * (1 - 2 * Occupied(n, x, b));
        return delta * sum;
    }

    /// <summary>Every W_b(p) equals its closed form: the comparison the constructor makes, open to a wrong input.</summary>
    public static bool WeightsMatchClosedForm(string topology, int n, int p, IReadOnlyList<(int A, int B)> bonds,
        Dictionary<(int A, int B), BigRational> W, BigRational delta) =>
        bonds.All(b => W[b] == ClosedFormWeight(topology, n, p, delta, b));

    /// <summary>Per reading x, the largest |(block p's slowest level less block 1's)/x⁴ − (δ_p − δ₁)| over p; slowest[k, p − 1]
    /// holds −Re of block p's slowest level at xs[k].</summary>
    public static double[] SplitDeviations(double[,] slowest, double[] quotients, double[] xs)
    {
        int blocks = quotients.Length;
        var dev = new double[xs.Length];
        for (int k = 0; k < xs.Length; k++)
        {
            double x4 = Math.Pow(xs[k], 4);
            for (int p = 0; p < blocks; p++)
            {
                double split = -(slowest[k, p] - slowest[k, 0]) / x4;
                dev[k] = Math.Max(dev[k], Math.Abs(split - (quotients[p] - quotients[0])));
            }
        }
        return dev;
    }

    /// <summary>The closed forms of W_b(p) at J = γ = 1 (Theorem E (d), "Which filling lasts longest").</summary>
    public static BigRational ClosedFormWeight(string topology, int n, int p, BigRational delta, (int A, int B) bond)
    {
        BigRational d2 = delta * delta;
        if (topology == "star")
            return d2 * (n - 2 * p) * (n - 2 * p) * Binomial(n - 2, p - 1) / 2;
        bool end = topology == "chain" && (bond.A == 0 || bond.B == n - 1);
        if (end) return d2 * Binomial(n - 2, p - 1) / 2;
        return 4 * d2 * (p >= 2 ? Binomial(n - 4, p - 2) : BigInteger.Zero);
    }

    /// <summary>δ_p = Σ_b W_b(p)(f_a − f_b)²/(‖f‖²·C(N − 2, p − 1)), f orthogonal to the constants.</summary>
    public static double DetuningQuotient(int n, int p, IReadOnlyList<(int A, int B)> bonds, Dictionary<(int A, int B), BigRational> W, double[] f)
    {
        double num = 0, norm = f.Sum(v => v * v);
        foreach (var b in bonds)
        {
            double diff = f[b.A] - f[b.B];
            num += ToDouble(W[b]) * diff * diff;
        }
        return num / (norm * (double)Binomial(n - 2, p - 1));
    }

    /// <summary>−Re of the slowest non-stationary eigenvalue of the repo's block (p, p) at x = J/γ, γ = 1, and the block's
    /// 1-norm, the scale of the eigensolver's rounding.</summary>
    public static (double Rate, double Norm) SlowestRate(ComplexMatrix H, double[] gamma, int n, int p, double x)
    {
        var B = PerBlockLiouvillianBuilder.BuildBlockZ(x * H, gamma, SectorBlock.SectorFlatIndices(n, p, p));
        double norm = 0;
        for (int c = 0; c < B.ColumnCount; c++)
        {
            double col = 0;
            for (int r = 0; r < B.RowCount; r++) col += B[r, c].Magnitude;
            norm = Math.Max(norm, col);
        }
        var ev = B.Evd().EigenValues;
        double best = double.PositiveInfinity;
        foreach (var z in ev)
            if (z.Magnitude > 1e-10 && -z.Real < best) best = -z.Real;
        return (best, norm);
    }

    // ---------------------------------------------------------------- exact numbers

    /// <summary>A double is a dyadic rational; this returns it exactly.</summary>
    public static BigRational ExactRational(double v)
    {
        if (double.IsNaN(v) || double.IsInfinity(v)) throw new ArgumentException("not finite", nameof(v));
        if (v == 0.0) return new BigRational(0);
        long bits = BitConverter.DoubleToInt64Bits(v);
        bool negative = bits < 0;
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xFFFFFFFFFFFFFL;
        if (exponent == 0) exponent = 1; else mantissa |= 1L << 52;
        exponent -= 1075;
        BigInteger num = mantissa, den = BigInteger.One;
        if (exponent > 0) num <<= exponent; else den <<= -exponent;
        return new BigRational(negative ? -num : num, den);
    }

    private static BigRational ExactEntry(Complex z)
    {
        if (z.Imaginary != 0.0) throw new InvalidOperationException($"a non-real entry {z} where the computational basis gives real ones");
        return ExactRational(z.Real);
    }

    private static BigRational DyadicNear(double v, bool roundUp)
    {
        double scaled = v * (1 << 24);
        double r = roundUp ? Math.Ceiling(scaled) : Math.Floor(scaled);
        return new BigRational(new BigInteger(r), new BigInteger(1 << 24));
    }

    private static BigInteger Binomial(int n, int k)
    {
        if (k < 0 || n < 0 || k > n) return BigInteger.Zero;
        BigInteger r = BigInteger.One;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }

    private static double ToDouble(BigRational q) => (double)q.Numerator / (double)q.Denominator;

    private static double[,] ToDouble(BigInteger[,] m)
    {
        var r = new double[m.GetLength(0), m.GetLength(1)];
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) r[i, j] = (double)m[i, j];
        return r;
    }

    // ---------------------------------------------------------------- IInspectable

    public string DisplayName => $"F50's Zeno end (Theorem E), N={N} {Topology}, Δ = {Delta.ToString(Inv)}";

    public string Summary =>
        $"ferromagnet {Mark(FerromagnetHolds)}, (A²)_PP = 8 L_p {Mark(ZenoGeneratorHolds)}, lift {Mark(LiftHolds)}, token gap = λ₁ = " +
        $"{Lambda1.ToString("0.000000", Inv)} in every block {Mark(TokenGapHolds)} ({GapRoute}), " +
        $"detuning weights closed form {Mark(DetuningClosedFormHolds)}; slowest at the Zeno end: blocks {{{string.Join(",", PredictedSlowestBlocks)}}} " +
        $"predicted, {{{string.Join(",", MeasuredSlowestBlocks)}}} read at x = {ReadingXs[0].ToString(Inv)}" +
        (ReadSplitAtRounding ? (PredictedSplitVanishes ? ", no split predicted and none read" : ", a split predicted and NONE READ")
         : DeviationAtRounding ? ", the split read as predicted to the eigensolver's rounding"
         : $", split deviation ×{DeviationRatio.ToString("0.00", Inv)} per doubling (law 4)");

    private string GapRoute => Topology switch
    {
        "chain" => "the gap by Theorem E (b), simple by (c)",
        "star" => "the gap by Theorem E (b) on this tree, its multiplicity read",
        _ => "read on the ring, where the proof cites the outside route",
    };

    private static string Mark(bool ok) => ok ? "exact" : "FAILS";

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the ferromagnet, exactly",
                summary: $"2·L_p = n_B − (Σ_b XX + YY + ZZ)_p on every block, the right side from the repo's Pauli matrices: " +
                         string.Join(", ", Blocks.Select(b => $"p={b.P} {Mark(b.Ferromagnet)}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the Zeno generator, exactly",
                summary: $"(A²)_PP = 8·L_p at Δ = {Delta.ToString(Inv)} on every block (the diagonal of H cancels in it by construction, " +
                         $"so this checks the construction; the leading coefficients below test the blindness): " +
                         string.Join(", ", Blocks.Select(b => $"p={b.P} {Mark(b.ZenoGeneratorIsEightLaplacian)}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the lift, exactly",
                summary: "L_p·F_p = F_p·L_site over the integers, F_p[S][i] = 1 if i ∈ S: " +
                         string.Join(", ", Blocks.Select(b => $"p={b.P} {Mark(b.LiftIntertwines)}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the token gap = λ₁, by exact counts",
                summary: $"λ₁ = {Lambda1.ToString("0.0000000000", Inv)}, multiplicity {Lambda1Multiplicity}, bracketed by {QLow.Numerator}/{QLow.Denominator} and " +
                         $"{QHigh.Numerator}/{QHigh.Denominator} (verified on the site Laplacian); eigenvalues of L_p below them, Jacobi's rule on Bareiss " +
                         "minors: " + string.Join(", ", Blocks.Select(b => $"p={b.P}: {b.CountBelowLow}, {b.CountBelowHigh}")) +
                         $" (expected 1 and {1 + Lambda1Multiplicity}): {(TokenGapHolds ? "the gap of every block is λ₁" : "A BLOCK HAS ANOTHER COUNT")}, " +
                         GapRoute,
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the detuning weights of Theorem E (d), exactly",
                summary: "W_b(p) = Σ over the hops across b of 2|t|²ΔV²/Γ³ against the closed forms: " +
                         string.Join(", ", Blocks.Select(b => $"p={b.P} {Mark(b.DetuningClosedForm)}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("which filling lasts longest (reading)",
                summary: "δ_p on the lifted wave: " + string.Join(", ", Blocks.Select(b => $"p={b.P} {b.DetuningQuotient.ToString("0.000000", Inv)}")) +
                         $"; predicted slowest {{{string.Join(",", PredictedSlowestBlocks)}}}; on the Liouvillian at x = {ReadingXs[0].ToString(Inv)} the slowest " +
                         $"rate over x² per block " + string.Join(", ", Blocks.Select(b => $"p={b.P} {b.LeadingCoefficient.ToString("0.000000", Inv)}")) +
                         $" against 2λ₁ = {(2 * Lambda1).ToString("0.000000", Inv)} (deviation ×" +
                         string.Join(", ×", Blocks.Select(b => b.LeadingDeviationRatio(Lambda1).ToString("0.00", Inv))) +
                         $" per doubling of x, law 4), slowest {{{string.Join(",", MeasuredSlowestBlocks)}}} " +
                         $"({(SelectionAgrees ? "inside the predicted set" : "OUTSIDE THE PREDICTED SET")}); the split (block p's slowest level less block 1's)/x⁴ off δ_p − δ₁ by {DeviationAtFirstX.ToString("0.0e0", Inv)} " +
                         (ReadSplitAtRounding
                             ? $"at x = {ReadingXs[0].ToString(Inv)}, the read split at the eigensolver's rounding: " +
                               (PredictedSplitVanishes ? "no split predicted (no ZZ term) and none read" : "a split predicted and NONE READ")
                             : DeviationAtRounding
                                 ? $"at x = {ReadingXs[0].ToString(Inv)}, within the eigensolver's rounding of the prediction, too close for its law to be read"
                                 : $"at x = {ReadingXs[0].ToString(Inv)}, ×{DeviationRatio.ToString("0.00", Inv)} at x = {ReadingXs[1].ToString(Inv)} (the law x² gives 4)"),
                provenance: NodeProvenance.Live);
        }
    }
}
