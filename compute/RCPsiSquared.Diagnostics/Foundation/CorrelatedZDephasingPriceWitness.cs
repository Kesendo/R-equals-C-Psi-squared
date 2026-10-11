using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for Z dephasing with correlated site fields
/// (<c>experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md</c>): the palindrome's price becomes the
/// operator {K, ·}, K = Σ_μ c_μ†c_μ = Σ_kl a_kl Z_l Z_k, so that for a real H commuting with X^⊗N the
/// residual M = Π·L·Π⁻¹ + L + 2σ·I (σ = tr a) is exactly −{K − σ, ·}, with
/// ‖M‖²_F = 8·4^N·Σ_{k&lt;l} a_kl².
///
/// <para>The page's verifier builds L from the Kossakowski matrix entry by entry on matrix units
/// and conjugates by Π_Z as a permutation. This witness takes the other road: it draws an integer
/// matrix B, realises a = BᵀB + I by the JUMPS c_μ = Σ_k B_μk Z_k plus one Z_k per site, builds L
/// with <see cref="LindbladianBuilder.Build"/>, and reads the residual in the Pauli-string basis with
/// the repository's canonical Π (<see cref="PalindromeResidual.Build"/>). Every entry along the way
/// is an integer or a dyadic rational, so every comparison is exact.</para>
///
/// <para>Children: the price identity (the assert, max|M + {K − σ, ·}| = 0, beside max|M| &gt; 0 for
/// the number 2σ alone); the residual's closed norm; the field control (a Z field on one site breaks
/// the identity through the same door); the shared field, independent jumps t_k·Z_k plus the
/// collective jump 2·Σ_k Z_k (κ = 4), whose generator differs from the independent one, the local
/// generator at the shifted rates γ_k − κ = t_k², by exactly −2κ(p − q)² on every joint-popcount
/// block, and whose K, averaged over a block's two popcounts, is the claim's
/// <see cref="CorrelatedZDephasingPriceClaim.BlockSigma"/> on every block; the edge blocks (0,1),
/// (1,0), (N−1,N), (N,N−1), where each coherence disagrees at one site l and decays at the rate 2·a_ll
/// whatever the covariances, while the image (1,N) of (0,1) does not; the four shapes, the popcounts on which K
/// takes one value read cell by cell (from the jumps for this witness's own a, from K = Σ_kl a_kl Z_l Z_k for a ring,
/// a star and one common covariance) against the claim's
/// <see cref="CorrelatedZDephasingPriceClaim.SingleValuePopcounts"/>, read from the covariances alone, with K there
/// the shared field's value at the mean covariance, compared in the integers; the third moment, Tr((L + σ)³) read from
/// L itself (2(L + σ) has integer entries, so the trace of its cube is an integer sum over its nonzero entries) against
/// the claim's <see cref="CorrelatedZDephasingPriceClaim.ThirdMoment"/>, for this witness's a and for a covariance −1
/// between the two end sites alone, where the summed covariance is negative, the GHZ argument says nothing, and the
/// third moment is nonzero; the far-end modes, a shared field of negative covariance realised by the jumps Z_k − Z_l
/// for every pair plus one Z_k per site (a_kk = N, a_kl = −1), on which P_p·X^⊗N is read live as an eigenvector of L
/// at −2K_p for every popcount p, the popcounts with K_p > σ (a rate above 2σ) against the claim's
/// <see cref="CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts"/>.</para>
///
/// <para>Guard: N ∈ {3, 4, 5} (4^N ≤ <see cref="MaxDim"/>; at N = 2 the image (1, N) is itself an
/// edge block). Typed home <see cref="CorrelatedZDephasingPriceClaim"/>.</para></summary>
public sealed class CorrelatedZDephasingPriceWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest residual dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public int Seed { get; }

    /// <summary>The Kossakowski matrix a = BᵀB + I over the Z_k, every off-diagonal entry nonzero.</summary>
    public int[,] Kossakowski { get; }

    /// <summary>σ = tr(a), F1's half-price.</summary>
    public double Sigma { get; }

    /// <summary>max over entries of |M + {K − σ, ·}| in the Pauli basis. The identity says 0.</summary>
    public double PriceIdentityResidual { get; }

    /// <summary>max |M|: the number 2σ alone as the shift. Under correlation this is not 0.</summary>
    public double ScalarShiftResidual { get; }

    /// <summary>‖M‖²_F, read off the Pauli-basis residual.</summary>
    public double ResidualNormSquared { get; }

    /// <summary>8·4^N·Σ_{k&lt;l} a_kl², the closed form.</summary>
    public double PredictedNormSquared { get; }

    /// <summary>max|M_f + {K − σ, ·}| with a Z field on site 0 added to H (H no longer commutes
    /// with X^⊗N). Must be nonzero: the control through the same door.</summary>
    public double FieldControlResidual { get; }

    /// <summary>The shared covariance of the shared-field construction (the collective jump
    /// 2·Σ_k Z_k carries κ = 4).</summary>
    public int SharedCovariance => 4;

    /// <summary>max over entries of |L_shared − L_independent − diag(−2κ(p − q)²)| on matrix units.</summary>
    public double SharedFieldBlockLawResidual { get; }

    /// <summary>max |L_shared − L_independent|: nonzero, so the block law is not read on an empty difference.</summary>
    public double SharedFieldShiftNorm { get; }

    /// <summary>max over the four edge blocks of |D(u,u) + 2·a_ll|, l the one site where ket and bra
    /// disagree, D the dissipator part of L (a decay at the rate 2·a_ll). The edge-block reading says 0.</summary>
    public double EdgeBlockResidual { get; }

    /// <summary>True iff some coherence of the image block (1, N) decays at a rate other than the
    /// local one, 2·Σ_{k disagreeing} a_kk.</summary>
    public bool ImageBlockCarriesCovariance { get; }

    /// <summary>max over every block (p, q) and every coherence in it of |(K(i) + K(j))/2 − σ(p, q)|
    /// for the shared field, σ(p, q) read from <see cref="CorrelatedZDephasingPriceClaim.BlockSigma"/>.
    /// The block law says 0.</summary>
    public double SharedFieldCentreResidual { get; }

    /// <summary>For each of the four patterns read (this witness's own a; and a ring, a star and one common covariance,
    /// those three with a covariance of 1 on each of their pairs), the popcounts on which K takes one value, read cell by
    /// cell, beside the claim's prediction.</summary>
    public IReadOnlyList<(string Pattern, int[] Read, int[] Predicted)> FourShapes { get; }

    /// <summary>True iff on every pattern the popcounts read cell by cell are the claim's
    /// <see cref="CorrelatedZDephasingPriceClaim.SingleValuePopcounts"/>.</summary>
    public bool FourShapesAgree => FourShapes.All(s => s.Read.SequenceEqual(s.Predicted));

    /// <summary>The number of distinct sets of popcounts among the patterns read, at the default seed 2 at N = 3, 4 at
    /// N = 4 and 3 at N = 5 (the N = 4 count needs the seed's own a to be generic), so the agreement is not read on one
    /// shape.</summary>
    public int DistinctShapes => FourShapes.Select(s => string.Join(",", s.Read)).Distinct().Count();

    /// <summary>True iff on every pattern and every popcount p that sees one value of K, C(N, 2)·K_p =
    /// C(N, 2)·σ − N·T + T·(N − 2p)², the shared field's value at the mean covariance T/C(N, 2), in the integers.</summary>
    public bool MeanCovarianceHolds { get; }

    /// <summary>Tr((2L + 2σ)³) read from this witness's L, real and imaginary part, exact.</summary>
    public (BigInteger Real, BigInteger Imaginary) CubeTrace { get; }

    /// <summary>Tr((L + σ)³) by the claim's closed form, from this witness's a and the chain's Tr(H²Z_kZ_l), Tr(H Z_k).</summary>
    public BigInteger ThirdMomentPredicted { get; }

    /// <summary>The noise term of the bracket for this witness's a, positive for a ⪰ 0 with a covariance.</summary>
    public BigInteger ThirdMomentNoiseTerm { get; }

    /// <summary>The control the GHZ argument cannot reach: jumps Z_0 − Z_{N−1} plus one Z_k per site, so the only
    /// covariance is −1, between the two end sites, and the summed covariance is negative.</summary>
    public long[,] EndCovarianceKossakowski { get; }

    /// <summary>Σ_{k≠l} a_kl·Tr(H²Z_kZ_l) for the end-site control: at most zero, so the break follows.</summary>
    public BigInteger EndCovarianceWeightedSum { get; }

    /// <summary>Tr((2L + 2σ)³) read from the end-site control's L.</summary>
    public (BigInteger Real, BigInteger Imaginary) EndCovarianceCubeTrace { get; }

    /// <summary>Tr((L + σ)³) by the closed form for the end-site control.</summary>
    public BigInteger EndCovarianceThirdMomentPredicted { get; }

    /// <summary>The far-end modes of a negatively shared field (jumps Z_k − Z_l for every pair plus one Z_k per site, so
    /// a_kk = N and a_kl = −1): for each popcount p, whether P_p·X^⊗N is an eigenvector of L at −2K_p exactly, and K_p.</summary>
    public IReadOnlyList<(int Popcount, bool Eigenvector, long K)> FarEndModes { get; }

    /// <summary>σ of the negatively shared field, N².</summary>
    public long FarEndSigma { get; }

    /// <summary>The popcounts whose far-end mode decays faster than 2σ, read from L, against the claim's law.</summary>
    public IReadOnlyList<int> FarEndRead { get; }
    public IReadOnlyList<int> FarEndPredicted { get; }

    /// <summary>Every P_p·X^⊗N is an eigenvector at −2K_p, the middle decays faster than 2σ, and the read popcounts
    /// are the claim's law.</summary>
    public bool FarEndHolds => FarEndModes.All(m => m.Eigenvector) && FarEndRead.SequenceEqual(FarEndPredicted)
                               && FarEndRead.Contains(N / 2);

    /// <summary>The far-end popcounts of the end-site control by the claim's law: empty, so no far-end mode (and, the
    /// summed covariance being negative, no GHZ coherence) decides there.</summary>
    public IReadOnlyList<int> EndCovarianceFarEndPopcounts { get; }

    /// <summary>Both readings equal 8 times the closed form with a zero imaginary part, and the end-site control's third
    /// moment is nonzero where the summed covariance is negative, no far-end mode decides, and the weighted sum is at
    /// most zero.</summary>
    public bool ThirdMomentHolds =>
        CubeTrace == (8 * ThirdMomentPredicted, BigInteger.Zero)
        && EndCovarianceCubeTrace == (8 * EndCovarianceThirdMomentPredicted, BigInteger.Zero)
        && ThirdMomentNoiseTerm > 0 && EndCovarianceFarEndPopcounts.Count == 0
        && EndCovarianceWeightedSum <= 0 && EndCovarianceThirdMomentPredicted != 0;

    public CorrelatedZDephasingPriceWitness(int n = 3, int seed = 7)
    {
        if (n < 3) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 3 (at N = 2 the image of an edge block is an edge block); got {n}");
        if (n > 5) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≤ 5 (4^N ≤ {MaxDim}); got {n}");
        N = n;
        Seed = seed;
        int d = 1 << n;

        var zs = Enumerable.Range(0, n)
            .Select(k => new PauliHamiltonian(n, new[] { PauliTerm.SingleSite(n, k, PauliLetter.Z, Complex.One) }).ToMatrix())
            .ToArray();
        var zval = new int[n, d];
        var pop = new int[d];
        for (int k = 0; k < n; k++)
            for (int i = 0; i < d; i++)
            {
                zval[k, i] = (int)Math.Round(zs[k][i, i].Real);
                if (zval[k, i] < 0) pop[i]++;
            }

        // The integer XXZ chain, real and commuting with X^⊗N.
        var terms = new List<PauliTerm>();
        for (int i = 0; i < n - 1; i++)
        {
            terms.Add(PauliTerm.TwoSite(n, i, PauliLetter.X, i + 1, PauliLetter.X, new Complex(5, 0)));
            terms.Add(PauliTerm.TwoSite(n, i, PauliLetter.Y, i + 1, PauliLetter.Y, new Complex(5, 0)));
            terms.Add(PauliTerm.TwoSite(n, i, PauliLetter.Z, i + 1, PauliLetter.Z, new Complex(2, 0)));
        }
        ComplexMatrix H = new PauliHamiltonian(n, terms).ToMatrix();

        // a = BᵀB + 1 with every covariance nonzero, realised by integer jumps.
        var rng = new Random(seed);
        int[,] B;
        int[,] a;
        for (int attempt = 0; ; attempt++)
        {
            if (attempt == 10000)
                throw new InvalidOperationException($"seed {seed} drew no Kossakowski matrix with every covariance nonzero in 10000 tries");
            B = new int[n, n];
            for (int m = 0; m < n; m++)
                for (int k = 0; k < n; k++)
                    B[m, k] = rng.Next(-2, 3);
            a = new int[n, n];
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++)
                {
                    int s = k == l ? 1 : 0;
                    for (int m = 0; m < n; m++) s += B[m, k] * B[m, l];
                    a[k, l] = s;
                }
            bool allNonzero = true;
            for (int k = 0; k < n && allNonzero; k++)
                for (int l = 0; l < n; l++)
                    if (k != l && a[k, l] == 0) { allNonzero = false; break; }
            if (allNonzero) break;
        }
        Kossakowski = a;

        var jumps = new List<ComplexMatrix>();
        for (int m = 0; m < n; m++)
        {
            var c = Matrix<Complex>.Build.Dense(d, d);
            for (int k = 0; k < n; k++)
                if (B[m, k] != 0) c += (Complex)B[m, k] * zs[k];
            jumps.Add(c);
        }
        for (int k = 0; k < n; k++) jumps.Add(zs[k]);

        double sigma = 0;
        for (int k = 0; k < n; k++) sigma += a[k, k];
        Sigma = sigma;

        ComplexMatrix L = LindbladianBuilder.Build(H, jumps);
        ComplexMatrix M = PalindromeResidual.Build(L, n, sigma);

        // −{K − σ, ·} in the Pauli basis, by the same transform PalindromeResidual uses.
        var K = Matrix<Complex>.Build.Dense(d, d);
        foreach (var c in jumps) K += c.ConjugateTranspose() * c;
        var Id = Matrix<Complex>.Build.DenseIdentity(d);
        var Ks = K - (Complex)sigma * Id;
        var antiVec = Ks.KroneckerProduct(Id) + Id.KroneckerProduct(Ks.Transpose());
        var T = PauliBasis.VecToPauliBasisTransform(n);
        var antiPauli = (T.ConjugateTranspose() * antiVec * T) / Math.Pow(2, n);

        PriceIdentityResidual = MaxAbs(M + antiPauli);
        ScalarShiftResidual = MaxAbs(M);
        ResidualNormSquared = SumSquares(M);
        double cov2 = 0;
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++)
                cov2 += (double)a[k, l] * a[k, l];
        PredictedNormSquared = CorrelatedZDephasingPriceClaim.ResidualNormSquared(n, cov2);

        // The control through the same door: a Z field on site 0.
        ComplexMatrix Lf = LindbladianBuilder.Build(H + (Complex)3.0 * zs[0], jumps);
        ComplexMatrix Mf = PalindromeResidual.Build(Lf, n, sigma);
        FieldControlResidual = MaxAbs(Mf + antiPauli);

        // The shared field: independent jumps t_k·Z_k (t_k = k + 2) plus the collective jump 2·Σ_k Z_k
        // (κ = 4); its Kossakowski matrix has a_kk = t_k² + κ and a_kl = κ.
        var indep = Enumerable.Range(0, n).Select(k => (Complex)(k + 2) * zs[k]).ToList();
        var collective = Matrix<Complex>.Build.Dense(d, d);
        for (int k = 0; k < n; k++) collective += (Complex)2.0 * zs[k];
        var shared = new List<ComplexMatrix>(indep) { collective };
        var diff = LindbladianBuilder.Build(H, shared) - LindbladianBuilder.Build(H, indep);
        SharedFieldShiftNorm = MaxAbs(diff);
        double blockLaw = 0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                int u = i * d + j;
                int dd = pop[i] - pop[j];
                double want = -2.0 * SharedCovariance * dd * dd;
                for (int v = 0; v < d * d; v++)
                {
                    double target = v == u ? want : 0.0;
                    blockLaw = Math.Max(blockLaw, Complex.Abs(diff[u, v] - target));
                }
            }
        SharedFieldBlockLawResidual = blockLaw;

        // The block centre, read live against the claim's closed form: (K(i) + K(j))/2 = σ(p, q).
        var Ksh = Matrix<Complex>.Build.Dense(d, d);
        foreach (var c in shared) Ksh += c.ConjugateTranspose() * c;
        double sigmaShared = Enumerable.Range(0, n).Sum(k => (double)((k + 2) * (k + 2) + SharedCovariance));
        double centre = 0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                double mean = (Ksh[i, i].Real + Ksh[j, j].Real) / 2.0;
                double want = CorrelatedZDephasingPriceClaim.BlockSigma(n, sigmaShared, SharedCovariance, pop[i], pop[j]);
                centre = Math.Max(centre, Math.Abs(mean - want) + Math.Abs(Ksh[i, i].Imaginary));
            }
        SharedFieldCentreResidual = centre;

        // The edge blocks, read on the dissipator part of L.
        var D = L - LindbladianBuilder.Build(H, Array.Empty<ComplexMatrix>());
        double edge = 0;
        bool imageCarries = false;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                int u = i * d + j;
                var disagree = Enumerable.Range(0, n).Where(k => zval[k, i] != zval[k, j]).ToList();
                bool isEdge = (pop[i], pop[j]) == (0, 1) || (pop[i], pop[j]) == (1, 0)
                              || (pop[i], pop[j]) == (n - 1, n) || (pop[i], pop[j]) == (n, n - 1);
                if (isEdge)
                    edge = Math.Max(edge, Complex.Abs(D[u, u] + 2.0 * a[disagree[0], disagree[0]]));
                if ((pop[i], pop[j]) == (1, n))
                {
                    double local = -2.0 * disagree.Sum(k => (double)a[k, k]);
                    if (D[u, u] != (Complex)local) imageCarries = true;
                }
            }
        EdgeBlockResidual = edge;
        ImageBlockCarriesCovariance = imageCarries;

        // The four shapes: the popcounts on which K takes one value, cell by cell, against the claim's closed form.
        long[,] Pattern(Func<int, int, bool> on)
        {
            var c = new long[n, n];
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++)
                    if (k != l && on(k, l)) c[k, l] = 1;
            return c;
        }
        var own = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) own[k, l] = a[k, l];
        var patterns = new List<(string Name, long[,] Cov, long[] Diagonal, Func<int, long> K)>
        {
            ("this witness's a", own, Enumerable.Range(0, n).Select(k => (long)a[k, k]).ToArray(), i =>
            {
                var entry = K[i, i];
                if (entry.Imaginary != 0.0 || entry.Real != Math.Round(entry.Real))
                    throw new InvalidOperationException($"K from integer jumps must be an integer on every configuration; got {entry}");
                return (long)entry.Real;
            }),
        };
        foreach (var (name, cov) in new[]
                 {
                     ("a ring", Pattern((k, l) => (k - l + n) % n == 1 || (l - k + n) % n == 1)),
                     ("a star", Pattern((k, l) => k == 0 || l == 0)),
                     ("one common covariance", Pattern((k, l) => true)),
                 })
        {
            var diag = Enumerable.Repeat(5L, n).ToArray();
            var cv = cov;
            patterns.Add((name, cv, diag, i =>
            {
                long e = 0;
                for (int k = 0; k < n; k++)
                    for (int l = 0; l < n; l++)
                        e += (k == l ? diag[k] : cv[k, l]) * zval[k, i] * zval[l, i];
                return e;
            }));
        }
        var shapes = new List<(string, int[], int[])>();
        bool meanHolds = true;
        long pairs = n * (n - 1) / 2;
        foreach (var (name, cov, diag, kOf) in patterns)
        {
            var read = new List<int>();
            long sigmaP = diag.Sum(), total = 0;
            for (int k = 0; k < n; k++)
                for (int l = k + 1; l < n; l++) total += cov[k, l];
            for (int p = 0; p <= n; p++)
            {
                var values = Enumerable.Range(0, d).Where(i => pop[i] == p).Select(kOf).Distinct().ToList();
                if (values.Count != 1) continue;
                read.Add(p);
                long np2 = n - 2 * p;
                if (pairs * values[0] != pairs * sigmaP - n * total + total * np2 * np2) meanHolds = false;
            }
            shapes.Add((name, read.ToArray(), CorrelatedZDephasingPriceClaim.SingleValuePopcounts(cov).ToArray()));
        }
        FourShapes = shapes;
        MeanCovarianceHolds = meanHolds;

        // The third moment: Tr(H²Z_kZ_l) and Tr(H Z_k) of the (traceless) chain, read off its diagonal.
        var H2 = H * H;
        var twoSite = new long[n, n];
        var single = new long[n];
        long ReadInteger(double x, string what)
        {
            if (x != Math.Round(x)) throw new InvalidOperationException($"{what} must be an integer; got {x}");
            return (long)x;
        }
        long traceH = 0;
        for (int i = 0; i < d; i++) traceH += ReadInteger(H[i, i].Real, "a diagonal entry of H");
        if (traceH != 0) throw new InvalidOperationException($"the closed form is written for a traceless H; Tr H = {traceH}");
        for (int k = 0; k < n; k++)
        {
            double t1 = 0;
            for (int i = 0; i < d; i++) t1 += H[i, i].Real * zval[k, i];
            single[k] = ReadInteger(t1, "Tr(H Z_k)");
            for (int l = 0; l < n; l++)
            {
                double t2 = 0;
                for (int i = 0; i < d; i++) t2 += H2[i, i].Real * zval[k, i] * zval[l, i];
                twoSite[k, l] = ReadInteger(t2, "Tr(H² Z_k Z_l)");
            }
        }
        var aLong = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++) aLong[k, l] = a[k, l];
        CubeTrace = CubeTraceOfTwoL(L, sigma);
        ThirdMomentPredicted = CorrelatedZDephasingPriceClaim.ThirdMoment(aLong, twoSite, single);
        ThirdMomentNoiseTerm = CorrelatedZDephasingPriceClaim.ThirdMomentNoiseTerm(aLong);

        var endJumps = new List<ComplexMatrix> { zs[0] - zs[n - 1] };
        for (int k = 0; k < n; k++) endJumps.Add(zs[k]);
        var aEnd = new long[n, n];
        for (int k = 0; k < n; k++) aEnd[k, k] = k == 0 || k == n - 1 ? 2 : 1;
        aEnd[0, n - 1] = aEnd[n - 1, 0] = -1;
        EndCovarianceKossakowski = aEnd;
        var endCov = new long[n, n];
        endCov[0, n - 1] = endCov[n - 1, 0] = -1;
        EndCovarianceFarEndPopcounts = CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(endCov);
        BigInteger weighted = 0;
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) weighted += (BigInteger)aEnd[k, l] * twoSite[k, l];
        EndCovarianceWeightedSum = weighted;
        double sigmaEnd = 0;
        for (int k = 0; k < n; k++) sigmaEnd += aEnd[k, k];
        EndCovarianceCubeTrace = CubeTraceOfTwoL(LindbladianBuilder.Build(H, endJumps), sigmaEnd);
        EndCovarianceThirdMomentPredicted = CorrelatedZDephasingPriceClaim.ThirdMoment(aEnd, twoSite, single);

        // The far-end modes of a negatively shared field: jumps Z_k − Z_l for every pair plus one Z_k per site.
        var antiJumps = new List<ComplexMatrix>();
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++) antiJumps.Add(zs[k] - zs[l]);
        for (int k = 0; k < n; k++) antiJumps.Add(zs[k]);
        ComplexMatrix La = LindbladianBuilder.Build(H, antiJumps);
        FarEndSigma = (long)n * n;
        var modes = new List<(int, bool, long)>();
        var farRead = new List<int>();
        for (int p = 0; p <= n; p++)
        {
            var v = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
            long kp = 0;
            for (int i = 0; i < d; i++)
            {
                if (pop[i] != p) continue;
                v[i * d + (d - 1 - i)] = Complex.One;
                long s = (long)n * n;                                    // K(i) = σ − Σ_{k≠l} z_k z_l
                for (int k = 0; k < n; k++)
                    for (int l = 0; l < n; l++)
                        if (k != l) s -= zval[k, i] * zval[l, i];
                kp = s;
            }
            var Lv = La * v;
            bool eig = true;
            for (int r = 0; r < d * d && eig; r++)
                if (Lv[r] != (Complex)(-2.0 * kp) * v[r]) eig = false;
            modes.Add((p, eig, kp));
            if (kp > FarEndSigma) farRead.Add(p);
        }
        FarEndModes = modes;
        FarEndRead = farRead;
        var antiCov = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) antiCov[k, l] = -1;
        FarEndPredicted = CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(antiCov);
    }

    /// <summary>Tr((2L + 2σ)³) exactly: the jumps are integer and the dissipator halves them, so 2(L + σ) has integer
    /// entries, and the trace of its cube is a sum of integer products over its nonzero entries, accumulated in
    /// BigInteger.</summary>
    private static (BigInteger Real, BigInteger Imaginary) CubeTraceOfTwoL(ComplexMatrix L, double sigma)
    {
        int dim = L.RowCount;
        var re = new long[dim, dim];
        var im = new long[dim, dim];
        var rows = new List<int>[dim];
        for (int r = 0; r < dim; r++)
        {
            rows[r] = new List<int>();
            for (int c = 0; c < dim; c++)
            {
                var z = L[r, c] + (r == c ? (Complex)sigma : Complex.Zero);
                double x = 2 * z.Real, y = 2 * z.Imaginary;
                if (x != Math.Round(x) || y != Math.Round(y) || Math.Abs(x) > 1e15 || Math.Abs(y) > 1e15)
                    throw new InvalidOperationException($"2(L + σ) must have integer entries; got {z} at ({r}, {c})");
                re[r, c] = (long)x;
                im[r, c] = (long)y;
                if (re[r, c] != 0 || im[r, c] != 0) rows[r].Add(c);
            }
        }
        BigInteger tRe = 0, tIm = 0;
        for (int i = 0; i < dim; i++)
            foreach (int j in rows[i])
                foreach (int k in rows[j])
                {
                    if (re[k, i] == 0 && im[k, i] == 0) continue;
                    BigInteger pr = (BigInteger)re[i, j] * re[j, k] - (BigInteger)im[i, j] * im[j, k];
                    BigInteger pi = (BigInteger)re[i, j] * im[j, k] + (BigInteger)im[i, j] * re[j, k];
                    tRe += pr * re[k, i] - pi * im[k, i];
                    tIm += pr * im[k, i] + pi * re[k, i];
                }
        return (tRe, tIm);
    }

    private static double MaxAbs(ComplexMatrix A)
    {
        double m = 0;
        for (int r = 0; r < A.RowCount; r++)
            for (int c = 0; c < A.ColumnCount; c++)
                m = Math.Max(m, Complex.Abs(A[r, c]));
        return m;
    }

    private static double SumSquares(ComplexMatrix A)
    {
        double s = 0;
        for (int r = 0; r < A.RowCount; r++)
            for (int c = 0; c < A.ColumnCount; c++)
            {
                var z = A[r, c];
                s += z.Real * z.Real + z.Imaginary * z.Imaginary;
            }
        return s;
    }

    public string DisplayName => $"CorrelatedZDephasingPriceWitness (N={N}, seed={Seed}, σ={Sigma.ToString("0", Inv)})";

    private bool PriceHolds => PriceIdentityResidual == 0.0;
    private bool NormHolds => ResidualNormSquared == PredictedNormSquared;
    private bool BlockLawHolds => SharedFieldBlockLawResidual == 0.0;
    private bool CentreHolds => SharedFieldCentreResidual == 0.0;
    private bool EdgeHolds => EdgeBlockResidual == 0.0;

    public string Summary =>
        $"N={N}: under Z dephasing whose site fields are correlated (every covariance nonzero), the palindrome " +
        $"residual of the XXZ chain {(PriceHolds ? "is exactly" : "is NOT")} −{{K − σ, ·}} (max deviation " +
        $"{PriceIdentityResidual.ToString("e2", Inv)}), the shift {(ScalarShiftResidual > 0.0 ? "not" : "possibly")} the number 2σ alone (max|M| = {ScalarShiftResidual.ToString("0.##", Inv)}), " +
        $"with ‖M‖² = {ResidualNormSquared.ToString("0", Inv)} {(NormHolds ? "=" : "≠")} 8·4^N·Σ_{{k<l}} a_kl²; against the " +
        $"local generator at the shifted rates, a shared field {(BlockLawHolds ? "shifts" : "does NOT shift")} every block " +
        $"by −2κ(p − q)² with K averaging to σ(p, q) on each block{(CentreHolds ? "" : " (that average FAILS)")}, and the four edge blocks {(EdgeHolds ? "decay at the local rates" : "do NOT decay at the local rates")} " +
        $"whatever the covariances{(ImageBlockCarriesCovariance ? "" : " (their image control FAILS)")}; the popcounts on which K takes one " +
        $"value {(FourShapesAgree ? "are" : "are NOT")} the claim's four-shapes law on {FourShapes.Count} patterns ({DistinctShapes} shapes){(MeanCovarianceHolds ? "" : ", the mean-covariance values FAIL")}; " +
        $"the third moment read from L {(ThirdMomentHolds ? "is" : "is NOT")} the claim's closed form, a negative covariance between the end sites alone " +
        $"{(EndCovarianceCubeTrace.Real != 0 ? "breaks" : "does NOT break")} the complex pairing, and under a negatively shared field the far-end modes " +
        $"{(FarEndHolds ? "are eigenvectors and the middle decays faster than 2σ, as the claim's law says" : "FAIL the claim's law")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the price identity, live (the assert)",
                summary: $"max|M + {{K − σ, ·}}| = {PriceIdentityResidual.ToString("e3", Inv)} in the Pauli basis, K = Σ c†c " +
                         $"for the integer jumps c_μ = Σ_k B_μk Z_k plus one Z_k per site ({(PriceHolds ? "the identity holds" : "the identity FAILS")}); " +
                         $"beside it max|M| = {ScalarShiftResidual.ToString("0.##", Inv)}, so the number 2σ alone " +
                         $"{(ScalarShiftResidual > 0.0 ? "does not serve" : "would serve")} as the shift.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the residual's closed norm",
                summary: $"‖M‖²_F = {ResidualNormSquared.ToString("0", Inv)} against 8·4^N·Σ_{{k<l}} a_kl² = " +
                         $"{PredictedNormSquared.ToString("0", Inv)}: {(ResidualNormSquared == PredictedNormSquared ? "equal" : "DIFFERENT")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the field control",
                summary: $"a Z field on site 0, odd under X^⊗N, through the same door: max|M + {{K − σ, ·}}| = " +
                         $"{FieldControlResidual.ToString("0.##", Inv)} ({(FieldControlResidual > 0.0 ? "the identity breaks, so the gate can fail" : "NO break, control compromised")}).",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the shared field, block by block",
                summary: $"independent jumps t_k·Z_k plus the collective jump 2·Σ_k Z_k (κ = {SharedCovariance}): against the " +
                         $"independent generator, the local one at the shifted rates, the difference {(BlockLawHolds ? "is" : "is NOT")} −2κ(p − q)² on every " +
                         $"joint-popcount block (max deviation {SharedFieldBlockLawResidual.ToString("e3", Inv)}, the shift itself up to " +
                         $"{SharedFieldShiftNorm.ToString("0", Inv)}), and (K(i) + K(j))/2 {(CentreHolds ? "is" : "is NOT")} σ(p, q) = σ − Nκ + κ(d² + m²) on every block " +
                         $"(max deviation {SharedFieldCentreResidual.ToString("e3", Inv)}, read against the claim's BlockSigma); " +
                         $"{(BlockLawHolds && CentreHolds && SharedFieldShiftNorm > 0.0 ? "both hold, on a nonzero shift" : "a law FAILS or the shift is empty")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the edge blocks",
                summary: $"on (0,1), (1,0), (N−1,N), (N,N−1) every coherence disagrees at one site l and {(EdgeHolds ? "decays" : "does NOT decay")} at the rate 2·a_ll " +
                         $"whatever the covariances (max deviation {EdgeBlockResidual.ToString("e3", Inv)}); the image (1,N) " +
                         $"of (0,1) {(ImageBlockCarriesCovariance ? "carries the covariances" : "does NOT carry them, gate compromised")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the four shapes",
                summary: string.Join("; ", FourShapes.Select(s =>
                             $"{s.Pattern}: one value of K on popcounts {{{string.Join(",", s.Read)}}} read cell by cell, " +
                             $"{{{string.Join(",", s.Predicted)}}} from the covariances alone")) +
                         $"; {(FourShapesAgree ? "the two routes meet" : "the routes DIFFER")} on {DistinctShapes} distinct shapes, and on every such " +
                         $"popcount K {(MeanCovarianceHolds ? "is" : "is NOT")} the shared field's value at the mean covariance.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the third moment",
                summary: $"Tr((2L + 2σ)³) read from L = {CubeTrace.Real} + {CubeTrace.Imaginary}i against 8 times the claim's ThirdMoment " +
                         $"{8 * ThirdMomentPredicted} for this witness's a (noise term {ThirdMomentNoiseTerm}, positive for a ⪰ 0 with a covariance); " +
                         $"for a covariance −1 between the two end sites alone (summed covariance −1, so the GHZ argument says nothing, and " +
                         $"no far-end mode decides: FarEndBreakPopcounts {{{string.Join(",", EndCovarianceFarEndPopcounts)}}}; " +
                         $"Σ_{{k≠l}} a_kl·Tr(H²Z_kZ_l) = {EndCovarianceWeightedSum}) read {EndCovarianceCubeTrace.Real} + {EndCovarianceCubeTrace.Imaginary}i " +
                         $"against {8 * EndCovarianceThirdMomentPredicted}: {(ThirdMomentHolds ? "both meet, and the end-site third moment is nonzero, so no single centre serves that spectrum's complex eigenvalues" : "a reading DIFFERS or the control is void")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the far-end modes",
                summary: $"a negatively shared field (jumps Z_k − Z_l for every pair plus one Z_k per site: a_kk = N, a_kl = −1, σ = {FarEndSigma}): " +
                         string.Join("; ", FarEndModes.Select(m => $"P_{m.Popcount}·X^⊗N {(m.Eigenvector ? "is" : "is NOT")} an eigenvector at −2K = {-2 * m.K}")) +
                         $"; K_p > σ, a rate above 2σ whose partner about −σ would grow, on popcounts {{{string.Join(",", FarEndRead)}}} read from L against " +
                         $"{{{string.Join(",", FarEndPredicted)}}} from the claim's FarEndBreakPopcounts: " +
                         $"{(FarEndHolds ? "they meet, so the summed covariance being negative, the middle breaks the single mirror, the rates included" : "they DIFFER")}.",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
