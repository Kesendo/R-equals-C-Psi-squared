using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The sideways-ladder chain walk of the <c>sideways_spin_ladder</c> open arc: the ONE
/// construction behind both the C# gate (<c>SidewaysSpinLadderGateTests</c>, Categories SIDEWAYS +
/// SLOW_SIDEWAYS) and the live witness (<see cref="SidewaysSpinLadderWitness"/>), so the gate's numbers
/// and the inspect-time numbers cannot drift apart. Port of <c>simulations/eta_ladder_chain.py</c>.
///
/// <para>What one <see cref="Run"/> computes, in the order that makes the later parts readable:
/// CONTROL — Φ (<see cref="SpectatorIntertwiner.BuildW"/>) and S⁺
/// (<see cref="SpectatorIntertwiner.BuildSPlus"/>) both intertwine L on every rung of the two chains
/// p+q̃ = N∓1, worst residual compared to 0.0 exactly (the Python gate measures exactly 0.0 here, and so
/// does this one; see <see cref="ExactResidual"/> for why the zero is order-independent). SHAPE — the interior blocks whose spectrum carries
/// the cross-fold partner −λ_A − 2N, against the two chain interiors. COUNT — fold + band sectors vs the
/// F125 orbit size 4N−8. LADDER — per chain, the transport norms ‖S⁺v‖ of the fold eigenvector against
/// the Clebsch-Gordan coefficients √(ℓ(ℓ+1) − m(m+1)), ℓ = (N−3)/2, plus the terminal step into the
/// boundary block as a RATIO to the interior norms (the vector comes from an eigensolver, so the exact
/// zero of the highest-weight death S⁺|ℓ,ℓ⟩ = 0 is not measurable, its floor ratio is) and the survivors
/// negative control (the terminal map is not the zero map).</para>
///
/// <para>The EVEN-N walk lives beside it: <see cref="RunN4"/> ports
/// <c>simulations/eta_ladder_chain_n4.py</c> (2026-08-10), walking the four real-q defective loci of
/// the N=4 (1,2) block (<see cref="SeedsN4"/>), where band and fold freight are a conjugate pair
/// sharing the confined 4-orbit.</para>
///
/// <para><b>Conventions</b> (identical to <c>SpectatorIntertwinerGateTests</c>): H = XYChain(N, 2q*) so
/// the hop amplitude is J = 2q*, γ = 1 per site, basis = <see cref="BlockBasis.PopcountStates"/>
/// ascending, flat = pIdx·Mq + qIdx, blocks via <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>,
/// eigensolver = LAPACK zgeev through <see cref="MklDirect"/>. The recorded λ_A is a 4-decimal value, so
/// the acceptance tolerance is DERIVED from the run (ten times the smallest offset actually found)
/// rather than set.</para></summary>
public static class SidewaysSpinLadderChain
{
    /// <summary>One interior rung of a chain: the fold eigenvector of (P,Q) transported to
    /// (P+1,Q−1) at norm <paramref name="Norm"/>, staying an eigenvector at the same λ up to
    /// <paramref name="EigResidual"/> (relative). <paramref name="TargetDim"/> is the target block's
    /// dimension, the scale of the eigensolver error model the gate normalizes the residual by.</summary>
    public readonly record struct Rung(int P, int Q, double Norm, double EigResidual, int TargetDim);

    /// <summary>One chain p+q̃ = Total: the interior rungs, the terminal norm into the boundary
    /// block, its ratio to the largest interior norm, and the survivors negative control
    /// (how many of the terminal source block's eigenvectors transport at norm &gt; 0.1;
    /// must be ≥ TerminalTargetDim for the death of ONE vector to mean anything).</summary>
    public readonly record struct Chain(
        int Total, IReadOnlyList<Rung> Rungs, double Terminal, double TerminalRatio,
        int Survivors, int TerminalSourceDim, int TerminalTargetDim);

    /// <summary>Everything one run at (N, q*, λ_A) produces; the gate asserts on these fields, the
    /// witness renders them.</summary>
    public sealed record LadderRun(
        int N, double QStar, double LambdaA, double Fold,
        double WorstControlResidual, double Tolerance,
        IReadOnlyList<(int P, int Q)> FoldSet, IReadOnlyList<(int P, int Q)> BandSet,
        IReadOnlyList<(int P, int Q)> PredictedSet,
        IReadOnlyList<double> CgPredicted, IReadOnlyList<Chain> Chains);

    /// <summary>The recorded seed used per N: the same census entries the Python gate hard-codes
    /// (N=5 locus 1 and the N=7 q*≈1.515 locus, both R-parity +1). Verified against
    /// <see cref="RealDefectiveSeeds.All"/> by the gate's seed test. N=4 deliberately throws too:
    /// its loci are not census SEEDS (real λ) but real-q defective LOCI with complex λ, walked by
    /// <see cref="RunN4"/> from <see cref="SeedsN4"/>.</summary>
    public static (double QStar, double LambdaA) Seed(int n) => n switch
    {
        5 => (0.620878, -4.6189),
        7 => (1.514833, -4.8846),
        _ => throw new ArgumentOutOfRangeException(nameof(n),
            $"no seed wired for N={n} (RealDefectiveSeeds lists odd N; the dense chain walk is seeded " +
            "at N=5 and N=7, the N=9 seed lives in SidewaysSpinLadderSparse.Seed, the N=4 complex-λ " +
            "loci in SeedsN4/RunN4)"),
    };

    /// <summary>The four real-q defective loci of the N=4 (1,2) block, the walk inputs of the even-N
    /// chain reading (octic-book q; J = 2q*). Same values the Python gate
    /// <c>simulations/eta_ladder_chain_n4.py</c> hard-codes; the gate's seed test verifies them against
    /// <see cref="ReferenceDefectiveLoci.For"/>(4) (the scanned, conjugation-closed defective set).
    /// Their λ are NOT recorded here: unlike the odd-N census seeds the N=4 λ* are complex, and
    /// <see cref="RunN4"/> measures the pair from the (1,2) spectrum itself.</summary>
    public static IReadOnlyList<double> SeedsN4 { get; } = new[] { 0.460212, 0.854438, 0.857458, 1.738181 };

    /// <summary>The Clebsch-Gordan transport norms √(ℓ(ℓ+1) − m(m+1)) for ℓ = (N−3)/2,
    /// m = −ℓ … ℓ−1 dropped one short (the interior rungs; length N−3).</summary>
    public static IReadOnlyList<double> CgCoefficients(int n)
    {
        double spin = (n - 3) / 2.0;
        var cg = new List<double>();
        for (double m = -spin; m < spin - 0.5; m += 1.0)
            cg.Add(Math.Sqrt(spin * (spin + 1) - m * (m + 1)));
        return cg;
    }

    /// <summary>The (p,q̃) joint-popcount block of L at real q: H = XYChain(n, 2q), γ = 1 per site,
    /// flat = pIdx·Mq + qIdx, exactly as <c>SpectatorIntertwinerGateTests.BuildBlockN</c>.</summary>
    public static ComplexMatrix BuildBlock(int n, int p, int qt, double q)
        => BuildBlockWithH(PauliHamiltonian.XYChain(n, 2.0 * q).ToMatrix(), n, p, qt);

    /// <summary>Same block, from a prebuilt H (one Run touches ~200 blocks; H is built once).</summary>
    private static ComplexMatrix BuildBlockWithH(ComplexMatrix h, int n, int p, int qt)
    {
        int d = 1 << n;
        var statesP = BlockBasis.PopcountStates(n, p);
        var statesQ = BlockBasis.PopcountStates(n, qt);
        var gamma = Enumerable.Repeat(1.0, n).ToList();
        var flat = new int[statesP.Count * statesQ.Count];
        for (int i = 0; i < statesP.Count; i++)
            for (int j = 0; j < statesQ.Count; j++)
                flat[i * statesQ.Count + j] = (int)(statesP[i] * d + statesQ[j]);
        return PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);
    }

    /// <summary>Max-entry magnitude of L_target·M − M·L_source (max-entry rather than the odd-N
    /// script's Frobenius; the N=4 sibling <c>eta_ladder_chain_n4.py</c> is max-entry too:
    /// identical as an exact-zero test, different as a diagnostic magnitude if it ever went nonzero). The
    /// comparison to 0.0 is the house case-1 rule: the operator identity cancels pairwise per entry, both
    /// sides of each entry draw the same multiset of addends (products of L entries with the ±1 ladder
    /// entries; γ = 1 real parts are integer multiples of 2, imaginary parts integer multiples of the one
    /// hop amplitude J), and the measured residual IS exactly 0.0, in numpy's zgemm and in this ordered
    /// walk alike, on every rung at N = 5 and 7. A nonzero here is a finding about the construction, not a
    /// tolerance case. Sparse column walk: every column of M has ≤ n nonzeros and every column of L has
    /// ≤ 2(n−1)+1, so the residual costs O(dim²·n), not O(dim³).</summary>
    public static double ExactResidual(Complex[,] lt, Complex[,] m, Complex[,] ls)
    {
        int tdim = lt.GetLength(0), sdim = ls.GetLength(0);
        double worst = 0.0;
        var colA = new Complex[tdim];
        var colB = new Complex[tdim];
        for (int j = 0; j < sdim; j++)
        {
            Array.Clear(colA);
            Array.Clear(colB);
            for (int k = 0; k < tdim; k++)          // (Lt·M)[:,j] = Σ_k M[k,j]·Lt[:,k], ascending k
            {
                var mkj = m[k, j];
                if (mkj == Complex.Zero) continue;
                for (int i = 0; i < tdim; i++)
                    if (lt[i, k] != Complex.Zero)
                        colA[i] += lt[i, k] * mkj;
            }
            for (int k = 0; k < sdim; k++)          // (M·Ls)[:,j] = Σ_k Ls[k,j]·M[:,k], ascending k
            {
                var lkj = ls[k, j];
                if (lkj == Complex.Zero) continue;
                for (int i = 0; i < tdim; i++)
                    if (m[i, k] != Complex.Zero)
                        colB[i] += m[i, k] * lkj;
            }
            for (int i = 0; i < tdim; i++)
                worst = Math.Max(worst, (colA[i] - colB[i]).Magnitude);
        }
        return worst;
    }

    /// <summary>CONTROL: both ladders (Φ and S⁺) intertwine L exactly on every rung of the two chains
    /// p+q̃ = N∓1; the worst residual is compared to 0.0 EXACTLY by the callers (see
    /// <see cref="ExactResidual"/> for why the zero is order-independent).</summary>
    private static double WorstControl(int n, ComplexMatrix h)
    {
        int[] chains = { n - 1, n + 1 };
        double worstControl = 0.0;
        foreach (bool phi in new[] { true, false })
        {
            int dq = phi ? 1 : -1;
            foreach (int total in chains)
                for (int p = 0; p <= n; p++)
                {
                    int q = total - p;
                    if (q < 0 || q > n || p + 1 > n || q + dq < 0 || q + dq > n) continue;
                    int sdim = Binom(n, p) * Binom(n, q);
                    int tdim = Binom(n, p + 1) * Binom(n, q + dq);
                    if (sdim == 0 || tdim == 0) continue;
                    var ls = BuildBlockWithH(h, n, p, q).ToArray();
                    var lt = BuildBlockWithH(h, n, p + 1, q + dq).ToArray();
                    var m = (phi ? SpectatorIntertwiner.BuildW(n, p, q)
                                 : SpectatorIntertwiner.BuildSPlus(n, p, q)).ToArray();
                    worstControl = Math.Max(worstControl, ExactResidual(lt, m, ls));
                }
        }
        return worstControl;
    }

    private static int Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }

    private static bool Interior(int n, int p, int q) => p >= 1 && p <= n - 1 && q >= 1 && q <= n - 1;

    private static double Abs2(Complex x) => x.Real * x.Real + x.Imaginary * x.Imaginary;

    private static Complex[] ToColMajor(Complex[,] m)
    {
        int rows = m.GetLength(0), cols = m.GetLength(1);
        var a = new Complex[(long)rows * cols];
        for (int j = 0; j < cols; j++)
            for (int i = 0; i < rows; i++)
                a[i + (long)j * rows] = m[i, j];
        return a;
    }

    /// <summary>The full chain walk at the recorded seed for <paramref name="n"/> (see
    /// <see cref="Seed"/>), or at an explicit (q*, λ_A).</summary>
    public static LadderRun Run(int n)
    {
        var (qStar, lamA) = Seed(n);
        return Run(n, qStar, lamA);
    }

    /// <summary>The full chain walk at an explicit (q*, λ_A): CONTROL residuals, the per-block spectra
    /// with the derived tolerance, the fold/band sector sets, and both chains' transport norms with
    /// terminal ratio and survivors control. The gate asserts on the returned fields; the witness
    /// renders them.</summary>
    public static LadderRun Run(int n, double qStar, double lamA)
    {
        double fold = -lamA - 2.0 * n;
        int[] chains = { n - 1, n + 1 };
        var h = PauliHamiltonian.XYChain(n, 2.0 * qStar).ToMatrix();   // built once, ~200 blocks reuse it

        double worstControl = WorstControl(n, h);

        // Spectra once per block, and a tolerance derived from the run (λ_A is 4-decimal).
        var spec = new Dictionary<(int P, int Q), Complex[]>();
        double bestOffset = double.PositiveInfinity;
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                int dim = Binom(n, p) * Binom(n, q);
                if (dim == 0) continue;
                var eigs = MklDirect.EigenvaluesOnlyRaw(ToColMajor(BuildBlockWithH(h, n, p, q).ToArray()), dim);
                spec[(p, q)] = eigs;
                bestOffset = Math.Min(bestOffset, eigs.Min(e => (e - fold).Magnitude));
                bestOffset = Math.Min(bestOffset, eigs.Min(e => (e - lamA).Magnitude));
            }
        double tol = Math.Max(bestOffset * 10, 1e-9);

        var foldSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.Min(e => (e - fold).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();
        var bandSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.Min(e => (e - lamA).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();
        var predSet = chains.SelectMany(total => Enumerable.Range(0, n + 1)
                                .Where(p => total - p >= 0 && total - p <= n && Interior(n, p, total - p))
                                .Select(p => (P: p, Q: total - p)))
                            .OrderBy(k => k).ToList();

        // LADDER: transport norms per chain against the Clebsch-Gordan coefficients.
        var chainResults = new List<Chain>();
        foreach (int total in chains)
        {
            var walk = Enumerable.Range(0, n + 1).Where(p => total - p >= 0 && total - p <= n)
                                 .Select(p => (P: p, Q: total - p)).ToList();
            var rungs = new List<Rung>();
            double terminal = double.NaN;
            Complex[,]? termLs = null, termS = null;
            int termTgtDim = 0, termSrcDim = 0;
            foreach (var (p, q) in walk.Take(walk.Count - 1))
            {
                int sdim = Binom(n, p) * Binom(n, q);
                int tdim = Binom(n, p + 1) * Binom(n, q - 1);
                if (sdim == 0 || tdim == 0) continue;
                var lsArr = BuildBlockWithH(h, n, p, q).ToArray();
                var (w, vr) = MklDirect.EigenvaluesAndVectorsDirectRaw(ToColMajor(lsArr), sdim);
                int jBest = 0;
                for (int j = 1; j < sdim; j++)
                    if ((w[j] - fold).Magnitude < (w[jBest] - fold).Magnitude) jBest = j;
                if ((w[jBest] - fold).Magnitude > tol) continue;   // source does not carry the fold
                var v = new Complex[sdim];
                double vn = 0.0;
                for (int i = 0; i < sdim; i++) { v[i] = vr[jBest * sdim + i]; vn += Abs2(v[i]); }
                vn = Math.Sqrt(vn);
                for (int i = 0; i < sdim; i++) v[i] /= vn;
                var s = SpectatorIntertwiner.BuildSPlus(n, p, q).ToArray();
                var img = new Complex[tdim];
                for (int k = 0; k < sdim; k++)
                {
                    if (v[k] == Complex.Zero) continue;
                    for (int i = 0; i < tdim; i++)
                        if (s[i, k] != Complex.Zero)
                            img[i] += s[i, k] * v[k];
                }
                double nrm = Math.Sqrt(img.Sum(Abs2));
                if (Interior(n, p + 1, q - 1))
                {
                    var ltArr = BuildBlockWithH(h, n, p + 1, q - 1).ToArray();
                    double res = 0.0;
                    for (int i = 0; i < tdim; i++)
                    {
                        Complex acc = -w[jBest] * img[i];
                        for (int k = 0; k < tdim; k++) acc += ltArr[i, k] * img[k];
                        res += Abs2(acc);
                    }
                    res = Math.Sqrt(res) / nrm;
                    rungs.Add(new Rung(p, q, nrm, res, tdim));
                }
                else
                {
                    terminal = nrm;
                    termLs = lsArr; termS = s; termTgtDim = tdim; termSrcDim = sdim;
                }
            }

            // The terminal step is NOT an exact zero and cannot be one: the intertwining residual IS
            // exactly 0.0 (operator identity, cancels pairwise), but v comes from an eigensolver, so the
            // meaningful quantity is the RATIO to the interior steps. The terminal check is weak on its
            // own (most of the source block is in the kernel at the last rung); the survivors negative
            // control is what makes the death readable.
            int survivors = 0;
            if (termLs is not null && termS is not null)
            {
                var (_, vrAll) = MklDirect.EigenvaluesAndVectorsDirectRaw(ToColMajor(termLs), termSrcDim);
                for (int j = 0; j < termSrcDim; j++)
                {
                    double cn = 0.0;
                    for (int i = 0; i < termSrcDim; i++) cn += Abs2(vrAll[j * termSrcDim + i]);
                    cn = Math.Sqrt(cn);
                    double outn = 0.0;
                    for (int i = 0; i < termTgtDim; i++)
                    {
                        Complex acc = Complex.Zero;
                        for (int k = 0; k < termSrcDim; k++)
                            if (termS[i, k] != Complex.Zero)
                                acc += termS[i, k] * vrAll[j * termSrcDim + k];
                        outn += Abs2(acc);
                    }
                    if (Math.Sqrt(outn) / cn > 0.1) survivors++;
                }
            }
            double ratio = rungs.Count > 0 && !double.IsNaN(terminal)
                ? terminal / rungs.Max(r => r.Norm)
                : double.NaN;
            chainResults.Add(new Chain(total, rungs, terminal, ratio, survivors, termSrcDim, termTgtDim));
        }

        return new LadderRun(n, qStar, lamA, fold, worstControl, tol,
            foldSet, bandSet, predSet, CgCoefficients(n), chainResults);
    }

    // ------------------------------------------------------------------------------- the N=4 walk

    /// <summary>One ℓ = 1/2 doublet of the N=4 walk: the chain B0 → B1 → B2 walked at ONE of the two
    /// conjugate values (λ* or its holomorphic fold −λ*−2N). The split pair's two members near that
    /// value (<paramref name="NormM1"/>/<paramref name="NormM2"/>) and their mix all transport at the
    /// CG norm 1; near the EP the members are almost parallel
    /// (<paramref name="OverlapDefect"/> = 1 − |⟨m1,m2⟩|, measured 2.6e-7..3.0e-6 across the loci), so
    /// the mix agreement is a smoke check here, not the N=9 mix-robustness. <paramref name="EigResidual"/>
    /// is the relative image-eigenvector residual of member 1 in B1 at its own eigenvalue,
    /// <paramref name="RungTargetDim"/> = dim(B1) the scale of the eigensolver error model the gate
    /// normalizes it by; <paramref name="TerminalRatio"/> the death of the transported vector at the
    /// step into the boundary block B2, as a ratio to the largest transport over ALL of B1's
    /// eigenvectors; <paramref name="Survivors"/> how many of those transport at norm &gt; 0.1 (must
    /// equal <paramref name="TerminalTargetDim"/>: the terminal map is not the zero map). The terminal
    /// denominator and the survivors count are computed once per CHAIN (B1's eigenbasis is shared), so
    /// across the four doublets they are two independent controls presented four times.</summary>
    public sealed record DoubletN4(
        (int P, int Q) B0, (int P, int Q) B1, (int P, int Q) B2, bool AtFoldValue,
        double OverlapDefect, double NormM1, double NormM2, double NormMix,
        double EigResidual, int RungTargetDim,
        double TerminalRatio, int Survivors, int TerminalTargetDim);

    /// <summary>Everything one N=4 run at q* produces; the gate asserts on these fields.
    /// <paramref name="LambdaStar"/> is MEASURED from the (1,2) spectrum (mean of the closest pair,
    /// split <paramref name="PairSplit"/>), not read from a census table; the tolerance is derived
    /// from the split, and the split itself is TRUNCATION-derived, not an EP floor: q* is recorded to
    /// 6 decimals, so the pair opens as ~√|δq| (measured 6.1e-4..3.9e-3), decades above an exact-EP
    /// √eps. <paramref name="FoldValue"/> is the HOLOMORPHIC fold −λ*−2N (§7's block-level identity),
    /// complex-valued here; §7's antiunitary partner map −conj(λ)−2N (the F125 fold in the F-registry)
    /// is the correct MAP but is the identity on its own fixed line, so as the twin-target formula it
    /// could never fail, which is why the holomorphic form is the gate-bearing choice.
    /// <paramref name="FoldVsConjDistance"/> = 2·<paramref name="SelfFoldPinResidual"/> exactly
    /// (algebra), both printed. READ the pin honestly: the (1,2) block is self-folded, its spectrum
    /// closed under σ: λ ↦ −conj(λ)−2N, and an isometry maps the closest pair to a closest pair, so
    /// wherever the closest pair is σ-STABLE its mean sits on Re λ = −4, machine-exact. That is not
    /// generic at N=4: where the closest pair is instead a non-stable pair, its σ-image ties its
    /// distance EXACTLY over the whole interval on which that pair is closest, the argmin pick is
    /// arbitrary, and the pin fails by O(1) — measured over ~27% of q ∈ [0.05, 3], contiguous
    /// intervals below q ≈ 0.86. So the pin is two-sidedly informative at 0.460212 (±1.5e-3),
    /// bounded only from BELOW at the 0.854/0.857 twins (the failure-region boundary q ≈ 0.852 sits
    /// just beneath them; upward it is free), and free at 1.738181 (downward to q ≈ 0.852, upward
    /// unbounded within the sweep box [0.05, 3]); it certifies
    /// σ-stability of the closest pair, never defectiveness. The loci LABELS are inherited from
    /// <see cref="ReferenceDefectiveLoci"/> (EpCharacter-classified there), and this walk gates
    /// transport structure at those q, exactly as the Python sibling states.</summary>
    public sealed record LadderRunN4(
        double QStar, Complex LambdaStar, double PairSplit, double Tolerance,
        double SelfFoldPinResidual, Complex FoldValue, double FoldVsConjDistance,
        IReadOnlyList<((int P, int Q) Block, double DistLambda, double DistFold)> OrbitCarries,
        IReadOnlyList<(int P, int Q)> FoldSet, IReadOnlyList<(int P, int Q)> BandSet,
        IReadOnlyList<(int P, int Q)> Orbit,
        double WorstControlResidual, IReadOnlyList<DoubletN4> Doublets);

    /// <summary>The N=4 chain walk at one real-q defective locus (<see cref="SeedsN4"/>): band and
    /// fold freight share the four confined-orbit sectors. Port of
    /// <c>simulations/eta_ladder_chain_n4.py</c>, with the arc's recorded port deltas: the seed is
    /// the complex pair (q*, λ*) with λ* measured from the (1,2) block; the fold formula stays the
    /// holomorphic −λ−2N, now complex-valued; the sector COUNT is the union |fold ∪ band| = 4N−12 = 4
    /// (the odd-N sum |fold| + |band| = 8 = 4N−8 still holds numerically at N=4 but double-counts the
    /// shared sectors, so the union is the count that carries the even-N accounting); and each block's
    /// split pair is read by the TWO closest eigenvalues per target value, so the conjugate twin is
    /// never silently dropped by a single argmin.
    ///
    /// <para>What is measurement and what is structure, stated once: S⁺(4,1,2) has singular spectrum
    /// {1×20, 2×4} and S⁺(4,2,1) is rank 4 with all singular values √3, so CG norm 1, terminal death
    /// and survivors = 4 hold for most eigenvectors — and, measured, the full LADDER passes at any
    /// generic q ≳ 0.58, so it is largely structure there (the sharpened form of the σ_min caveat).
    /// Below q ≈ 0.57 the closest-pair eigenvector transports at norm 2 instead and the LADDER
    /// fails; the locus 0.460212 sits in a ~1e-3 norm-1 island inside that failing region, so at
    /// THAT locus the CG reading is itself q-selective. The run-specific content is the VALUE-level
    /// λ-sharing: λ* in every orbit block at the derived tolerance, with fold set = band set = the
    /// confined 4-orbit (union gate). The conjugate twin's PRESENCE per block is a corollary, not a
    /// second measurement: block spectra are conjugation-closed (block-local T-conjugacy) and the pin
    /// makes fold = conj(λ*), so d(fold) tracks d(λ*) to ~1e-11. And the union gate's selectivity is
    /// calibrated, not absolute: away from a locus it fails by OVER-inclusion (tol balloons with the
    /// split), with measured pass-windows in q of ±1.3e-4 around the first two loci, ±3e-5 around the
    /// twin 0.857458, but −0.076/+0.101 around q* = 1.738181, whose nearest non-orbit interior block
    /// sits ~450× beyond tol. Gate evidence is walk structure AT the recorded loci, never a locus
    /// detector.</para></summary>
    public static LadderRunN4 RunN4(double qStar)
    {
        const int n = 4;
        var orbit = new[] { (P: 1, Q: 2), (P: 2, Q: 1), (P: 2, Q: 3), (P: 3, Q: 2) };
        var h = PauliHamiltonian.XYChain(n, 2.0 * qStar).ToMatrix();

        double worstControl = WorstControl(n, h);

        // λ*: the closest pair of the (1,2) spectrum, its mean; tolerance derived from the split.
        var w12 = MklDirect.EigenvaluesOnlyRaw(ToColMajor(BuildBlockWithH(h, n, 1, 2).ToArray()),
                                               Binom(n, 1) * Binom(n, 2));
        double split = double.PositiveInfinity;
        Complex lam = Complex.Zero;
        for (int i = 0; i < w12.Length; i++)
            for (int j = i + 1; j < w12.Length; j++)
                if ((w12[i] - w12[j]).Magnitude < split)
                {
                    split = (w12[i] - w12[j]).Magnitude;
                    lam = (w12[i] + w12[j]) / 2.0;
                }
        double tol = Math.Max(5.0 * split, 1e-9);   // Python: max(10·half, 1e-9), half = split/2

        Complex fold = -lam - 2.0 * n;              // holomorphic; = conj(λ*) exactly iff Re λ* = −4
        double pin = Math.Abs(lam.Real + n);
        double foldVsConj = (fold - Complex.Conjugate(lam)).Magnitude;

        // Every orbit block carries BOTH values, and fold/band sets over ALL interior blocks.
        var spec = new Dictionary<(int P, int Q), (Complex[] W, Complex[] Vr, int Dim)>();
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                int dim = Binom(n, p) * Binom(n, q);
                if (dim == 0) continue;
                var (w, vr) = MklDirect.EigenvaluesAndVectorsDirectRaw(
                    ToColMajor(BuildBlockWithH(h, n, p, q).ToArray()), dim);
                spec[(p, q)] = (w, vr, dim);
            }

        var orbitCarries = orbit.Select(b =>
        {
            var w = spec[b].W;
            return (Block: b, DistLambda: w.Min(e => (e - lam).Magnitude),
                              DistFold: w.Min(e => (e - fold).Magnitude));
        }).ToList();

        var foldSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.W.Min(e => (e - fold).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();
        var bandSet = spec.Where(kv => Interior(n, kv.Key.P, kv.Key.Q)
                                       && kv.Value.W.Min(e => (e - lam).Magnitude) < tol)
                          .Select(kv => kv.Key).OrderBy(k => k).ToList();

        // The four doublets: two chains × two conjugate values.
        var doublets = new List<DoubletN4>();
        foreach (var (b0, b1, b2) in new[]
                 {
                     ((P: 1, Q: 2), (P: 2, Q: 1), (P: 3, Q: 0)),
                     ((P: 2, Q: 3), (P: 3, Q: 2), (P: 4, Q: 1)),
                 })
        {
            var (ws, vs, sdim) = spec[b0];
            var (_, vt, tdimB1) = spec[b1];
            var s01 = SpectatorIntertwiner.BuildSPlus(n, b0.P, b0.Q).ToArray();
            var s12 = SpectatorIntertwiner.BuildSPlus(n, b1.P, b1.Q).ToArray();
            var ltArr = BuildBlockWithH(h, n, b1.P, b1.Q).ToArray();
            int termTgtDim = Binom(n, b2.P) * Binom(n, b2.Q);

            // transport norms of ALL of b1's eigenvectors under S12: the terminal denominator and
            // the survivors control, shared by both values of this chain.
            var normsAll = new double[tdimB1];
            for (int j = 0; j < tdimB1; j++)
            {
                double cn = 0.0;
                for (int i = 0; i < tdimB1; i++) cn += Abs2(vt[j * tdimB1 + i]);
                cn = Math.Sqrt(cn);
                double outn = 0.0;
                for (int i = 0; i < termTgtDim; i++)
                {
                    Complex acc = Complex.Zero;
                    for (int k = 0; k < tdimB1; k++)
                        if (s12[i, k] != Complex.Zero)
                            acc += s12[i, k] * vt[j * tdimB1 + k];
                    outn += Abs2(acc);
                }
                normsAll[j] = Math.Sqrt(outn) / cn;
            }
            double normsAllMax = normsAll.Max();
            int survivors = normsAll.Count(x => x > 0.1);

            foreach (var (atFold, val) in new[] { (false, lam), (true, fold) })
            {
                // the split pair near val: the TWO closest (conjugate-twin coverage; one argmin
                // per block would silently drop one of the two values' twins).
                var order = Enumerable.Range(0, sdim).OrderBy(j => (ws[j] - val).Magnitude).ToList();
                int k1 = order[0], k2 = order[1];

                Complex[] Col(int k)
                {
                    var v = new Complex[sdim];
                    double vn = 0.0;
                    for (int i = 0; i < sdim; i++) { v[i] = vs[k * sdim + i]; vn += Abs2(v[i]); }
                    vn = Math.Sqrt(vn);
                    for (int i = 0; i < sdim; i++) v[i] /= vn;
                    return v;
                }
                var v1 = Col(k1);
                var v2 = Col(k2);
                Complex ov = Complex.Zero;
                for (int i = 0; i < sdim; i++) ov += Complex.Conjugate(v1[i]) * v2[i];
                double overlapDefect = 1.0 - ov.Magnitude;

                var mix = new Complex[sdim];
                double mn = 0.0;
                for (int i = 0; i < sdim; i++) { mix[i] = v1[i] + v2[i]; mn += Abs2(mix[i]); }
                mn = Math.Sqrt(mn);
                for (int i = 0; i < sdim; i++) mix[i] /= mn;

                double Transport(Complex[] v, out Complex[] img)
                {
                    img = new Complex[tdimB1];
                    for (int k = 0; k < sdim; k++)
                    {
                        if (v[k] == Complex.Zero) continue;
                        for (int i = 0; i < tdimB1; i++)
                            if (s01[i, k] != Complex.Zero)
                                img[i] += s01[i, k] * v[k];
                    }
                    return Math.Sqrt(img.Sum(Abs2));
                }
                double n1 = Transport(v1, out var img1);
                double n2 = Transport(v2, out _);
                double nMix = Transport(mix, out _);

                // image of member 1 is an eigenvector of B1 at its own eigenvalue (relative residual)
                double res = 0.0;
                for (int i = 0; i < tdimB1; i++)
                {
                    Complex acc = -ws[k1] * img1[i];
                    for (int k = 0; k < tdimB1; k++) acc += ltArr[i, k] * img1[k];
                    res += Abs2(acc);
                }
                res = Math.Sqrt(res) / n1;

                // terminal: the transported vector dies at the step into the boundary block
                double outn2 = 0.0;
                for (int i = 0; i < termTgtDim; i++)
                {
                    Complex acc = Complex.Zero;
                    for (int k = 0; k < tdimB1; k++)
                        if (s12[i, k] != Complex.Zero)
                            acc += s12[i, k] * img1[k] / n1;
                    outn2 += Abs2(acc);
                }
                double terminalRatio = Math.Sqrt(outn2) / normsAllMax;

                doublets.Add(new DoubletN4(b0, b1, b2, atFold, overlapDefect, n1, n2, nMix,
                    res, tdimB1, terminalRatio, survivors, termTgtDim));
            }
        }

        return new LadderRunN4(qStar, lam, split, tol, pin, fold, foldVsConj,
            orbitCarries, foldSet, bandSet, orbit.OrderBy(k => k).ToList(), worstControl, doublets);
    }
}
