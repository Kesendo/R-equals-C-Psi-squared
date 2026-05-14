using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Probes;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F89 per-Bloch-mode sigma anatomy at the F86 c=2 stratum: extracts the
/// quadratic self-pairing amplitude σ_i = |c0_i|² · (R†·S·R)[i,i].Real for every
/// eigenmode of the uniform-J block-L. Filters F_a modes by Re(λ) ≈ −2γ₀ and
/// assigns each one a Bloch index n in the S_2-anti orbit
/// {2, 4, ..., 2·floor(N_block/2)} via Im(λ_i) / (J/2) ≈ y_n = 4·cos(πn/(N_block+1))
/// (F89-J = F86-J/2 per F90 bridge).
///
/// <para>Counterpart to <see cref="C2FullBlockEigenAnatomy"/>: same biorthogonal
/// eigendecomposition + Dicke probe + spectral coordinates, but projected through
/// the F86 spatial-sum kernel S = Σ_l 2·|w_l⟩⟨w_l| instead of the per-bond
/// Hellmann-Feynman matrix M_h_per_bond[b]. Verifies bit-exactly against the
/// analytical <see cref="F89UnifiedFaClosedFormClaim.Sigma"/> table for
/// path-3..9 (path-7 closed form derived 2026-05-13 from this anatomy; path-8 and
/// path-9 extracted via Vandermonde fit on the same anatomy at C2Block(9) and
/// C2Block(10) respectively); path-10+ open. Empirical rule: odd-part(D_path) = (odd-part(k))².</para>
///
/// <para>Tier outcome: Tier2Verified. Pure numerical anatomy at a fixed Q value;
/// the witnesses ARE the data, no derivation. Q defaults to 1.5 (away from Q_EP
/// to avoid Jordan-block degeneracy, matches the Python convention in
/// <c>simulations/_f89_pathk_survey.py</c>).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2) +
/// <see cref="F89UnifiedFaClosedFormClaim"/> +
/// <see cref="F89PathKAtLockMechanismClaim"/>.</para>
/// </summary>
public sealed class C2FullBlockSigmaAnatomy : Claim
{
    public const double DefaultQ = 1.5;

    public CoherenceBlock Block { get; }
    public double Q { get; }
    public IReadOnlyList<SigmaModeWitness> SigmaSpectrum { get; }

    /// <summary>True when this anatomy was produced by <see cref="BuildFaOnly"/> via
    /// targeted inverse iteration: <see cref="SigmaSpectrum"/> then holds ONLY the
    /// floor(N/2) F_a witnesses, not one witness per eigenmode. False for <see cref="Build"/>
    /// (full-spectrum zgeev).</summary>
    public bool FaModesOnly { get; }

    public static C2FullBlockSigmaAnatomy Build(CoherenceBlock block, double q = DefaultQ)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2FullBlockSigmaAnatomy applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        MathNetSetup.EnsureInitialized();
        var witnesses = ComputeAnatomy(block, q);
        return new C2FullBlockSigmaAnatomy(block, q, witnesses);
    }

    /// <summary>Build the F_a-mode sigma anatomy via targeted inverse iteration on the
    /// analytically-known F_a eigenvalues (Re(λ) = −2γ₀, Im(λ) = (J/2)·y_n), using dense
    /// LU instead of a full dense zgeev. Produces ONLY the floor(N/2) F_a witnesses — the
    /// rest of the spectrum is never computed.
    /// <para>Use at large N where a full zgeev is infeasible: zgeev's QR/Schur iteration
    /// (zhseqr) is sequential, so a single dim-17424 decomposition is hours. The inverse
    /// iteration here is dense LU (zgetrf), which is BLAS-3 and threads. For small N, or
    /// when the full spectrum is needed, use <see cref="Build"/>.</para>
    /// <para>Verified bit-exactly (1e-9) against the <see cref="Build"/> zgeev path for the
    /// F_a σ values at path-3..9; see <c>C2FullBlockSigmaAnatomyTests</c>.</para></summary>
    public static C2FullBlockSigmaAnatomy BuildFaOnly(CoherenceBlock block, double q = DefaultQ)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2FullBlockSigmaAnatomy applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        MathNetSetup.EnsureInitialized();
        var witnesses = ComputeFaAnatomyInverseIteration(block, q);
        return new C2FullBlockSigmaAnatomy(block, q, witnesses, faModesOnly: true);
    }

    private C2FullBlockSigmaAnatomy(
        CoherenceBlock block, double q,
        IReadOnlyList<SigmaModeWitness> witnesses,
        bool faModesOnly = false)
        : base($"c=2 full block-L sigma anatomy at Q={q:G4}: F89 per-Bloch-mode σ_n extraction via R†·S·R diagonal",
               Tier.Tier2Verified,
               Item1Anchors.Root + "; F89UnifiedFaClosedFormClaim")
    {
        Block = block;
        Q = q;
        SigmaSpectrum = witnesses;
        FaModesOnly = faModesOnly;
    }

    private static IReadOnlyList<SigmaModeWitness> ComputeAnatomy(
        CoherenceBlock block, double q)
        => ComputeAnatomyMklDirect(block, q);

    private static IReadOnlyList<SigmaModeWitness> ComputeAnatomyMklDirect(
        CoherenceBlock block, double q)
    {
        double j = q * block.GammaZero;
        var basis = block.Basis;
        int mTotal = basis.MTotal;

        // Build L as a column-major raw Complex[] (consumed + DESTROYED by MklDirect).
        Complex[] lRaw = BlockLDecomposition.BuildUniformLColumnMajorRaw(block, j);

        // Direct LAPACK zgeev_ via fixed pointers (no MathNet marshaller cap).
        var (values, leftVecs, rightVecs) =
            MklDirect.EigenvaluesLeftRightDirectRaw(lRaw, mTotal);
        // lRaw is now destroyed; do not reuse.

        // Probe: uniform Dicke-coherence vector, every entry identical.
        ComplexVector probeVec = DickeBlockProbe.Build(block);
        Complex probeConst = probeVec[0];

        // Per-site overlap index lists for the spatial-sum kernel.
        var siteOverlap = BuildSiteOverlapIndices(basis);

        int nBlock = block.N;
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        double faRateTarget = -2.0 * block.GammaZero;
        double faRateTolerance = 1e-3 * block.GammaZero;
        double faFreqTolerance = 1e-4;

        var witnesses = new List<SigmaModeWitness>(mTotal);

        for (int i = 0; i < mTotal; i++)
        {
            long colBase = (long)i * mTotal;

            // Biorthogonal c0[i] = (VL_i^H · probe) / (VL_i^H · VR_i).
            // probe is uniform: VL_i^H · probe = conj(Σ_r VL_i[r]) · probeConst.
            Complex vlSum = Complex.Zero;
            Complex vlHvr = Complex.Zero;
            for (int r = 0; r < mTotal; r++)
            {
                Complex vl = leftVecs[colBase + r];
                vlSum += vl;
                vlHvr += Complex.Conjugate(vl) * rightVecs[colBase + r];
            }
            Complex vlHprobe = Complex.Conjugate(vlSum) * probeConst;
            Complex c0 = vlHvr.Magnitude < 1e-300 ? Complex.Zero : vlHprobe / vlHvr;
            double overlapSq = c0.Real * c0.Real + c0.Imaginary * c0.Imaginary;

            // S-kernel diagonal: (R†·S·R)[i,i] = Σ_l 2·|Σ_{r in siteOverlap[l]} VR_i[r]|².
            // SpatialSumKernel.Build bakes in the factor 2; the loop below mirrors that.
            double sDiag = 0.0;
            foreach (int[] siteIndices in siteOverlap)
            {
                Complex proj = Complex.Zero;
                foreach (int r in siteIndices)
                    proj += rightVecs[colBase + r];
                sDiag += 2.0 * (proj.Real * proj.Real + proj.Imaginary * proj.Imaginary);
            }

            // SpatialSumKernel.Build returns Σ_l 2·|w_l⟩⟨w_l| (factor 2 baked in;
            // see compute/RCPsiSquared.Core/Probes/SpatialSumKernel.cs:51). The F89
            // σ convention omits this factor: per simulations/_f89_pathk_survey.py:87
            // (sig = sum |a|²) the per-site reduction stores the bare amplitude.
            // Half of (R†·S·R)[i,i] gives the correct match against F89UnifiedFaClosedForm.
            double sigma = 0.5 * overlapSq * sDiag;

            int? blochIndexN = null;
            if (Math.Abs(values[i].Real - faRateTarget) <= faRateTolerance)
            {
                // F90F86C2BridgeIdentity: full block-L uses J/2 as the effective
                // single-particle coupling, so F89's y_n = 4cos(πn/(N+1)) in units
                // of J_F89 maps to Im(λ) = (J_F86/2)·y_n in block-L eigenvalues.
                // Divide by 0.5·J to recover the dimensionless y_n correctly.
                double imOverHalfJ = values[i].Imaginary / (0.5 * j);
                foreach (int n in orbit)
                {
                    double yN = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
                    if (Math.Abs(imOverHalfJ - yN) <= faFreqTolerance)
                    {
                        blochIndexN = n;
                        break;
                    }
                }
            }

            witnesses.Add(new SigmaModeWitness(
                EigenIndexAtQ: i,
                EigenvalueReal: values[i].Real,
                EigenvalueImag: values[i].Imaginary,
                ProbeOverlapSquared: overlapSq,
                SKernelDiagonal: sDiag,
                Sigma: sigma,
                BlochIndexN: blochIndexN));
        }

        witnesses.Sort((a, b) => b.Sigma.CompareTo(a.Sigma));
        return witnesses;
    }

    private static IReadOnlyList<SigmaModeWitness> ComputeAnatomyMathNet(
        CoherenceBlock block, double q)
    {
        double j = q * block.GammaZero;
        // Use BuildUniformLAt: skips per-bond Mh storage (saves NumBonds × Mtot² memory),
        // critical at large N (e.g. nBlock=29 saves ~62 GB vs the full Decomposition path).
        // Equivalent to block.Decomposition.AssembleUniform(j) but does not allocate MhPerBond[].
        ComplexMatrix L = BlockLDecomposition.BuildUniformLAt(block, j);

        var evd = L.Evd();
        ComplexMatrix R = evd.EigenVectors;
        ComplexMatrix rInv = R.Inverse();
        var lambdas = evd.EigenValues;

        ComplexVector probe = DickeBlockProbe.Build(block);
        ComplexVector c0 = rInv * probe;

        ComplexMatrix sKernel = SpatialSumKernel.Build(block);
        ComplexMatrix sEigenBasis = R.ConjugateTranspose() * sKernel * R;

        int dim = block.Basis.MTotal;
        int nBlock = block.N;
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        double faRateTarget = -2.0 * block.GammaZero;
        double faRateTolerance = 1e-3 * block.GammaZero;
        double faFreqTolerance = 1e-4;

        var witnesses = new List<SigmaModeWitness>(dim);
        for (int i = 0; i < dim; i++)
        {
            double overlapSq = c0[i].Real * c0[i].Real + c0[i].Imaginary * c0[i].Imaginary;
            double sDiag = sEigenBasis[i, i].Real;
            double sigma = overlapSq * sDiag;

            int? blochIndexN = null;
            if (Math.Abs(lambdas[i].Real - faRateTarget) <= faRateTolerance)
            {
                // F90F86C2BridgeIdentity: full block-L uses J/2 as the effective
                // single-particle coupling, so F89's y_n = 4cos(πn/(N+1)) in units
                // of J_F89 maps to Im(λ) = (J_F86/2)·y_n in block-L eigenvalues.
                // Divide by 0.5·J to recover the dimensionless y_n correctly.
                double imOverJ = lambdas[i].Imaginary / (0.5 * j);
                foreach (int n in orbit)
                {
                    double yN = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
                    if (Math.Abs(imOverJ - yN) <= faFreqTolerance)
                    {
                        blochIndexN = n;
                        break;
                    }
                }
            }

            // SpatialSumKernel.Build returns Σ_l 2·|w_l⟩⟨w_l| (factor 2 baked in;
            // see compute/RCPsiSquared.Core/Probes/SpatialSumKernel.cs:51). The F89
            // σ convention omits this factor: per simulations/_f89_pathk_survey.py:87
            // (sig = sum |a|²) the per-site reduction stores the bare amplitude.
            // Half of (R†·S·R)[i,i] gives the correct match against F89UnifiedFaClosedForm.
            witnesses.Add(new SigmaModeWitness(
                EigenIndexAtQ: i,
                EigenvalueReal: lambdas[i].Real,
                EigenvalueImag: lambdas[i].Imaginary,
                ProbeOverlapSquared: overlapSq,
                SKernelDiagonal: sDiag,
                Sigma: 0.5 * sigma,
                BlochIndexN: blochIndexN));
        }

        witnesses.Sort((a, b) => b.Sigma.CompareTo(a.Sigma));
        return witnesses;
    }

    /// <summary>Per-site overlap index lists for the spatial-sum kernel: for each site l,
    /// the flat indices of basis pairs (p, q) with p_l = 0 and q = p | bit_l. Mirrors the
    /// construction inside <see cref="Probes.SpatialSumKernel.Build"/> but returns only the
    /// index lists, so the (R†·S·R) diagonal can be computed per-mode without materialising
    /// the full S matrix (which would itself hit the 2 GB cap at nBlock >= 29).</summary>
    private static List<int[]> BuildSiteOverlapIndices(BlockBasis basis)
    {
        int N = basis.N;
        int n = basis.LowerPopcount;
        var lists = new List<int[]>(N);
        for (int site = 0; site < N; site++)
        {
            var indices = new List<int>(64);
            long maskI = 1L << (N - 1 - site);
            foreach (long p in basis.StatesP)
            {
                if (((p >> (N - 1 - site)) & 1) != 0) continue;   // need p_site = 0
                long qState = p | maskI;
                if (System.Numerics.BitOperations.PopCount((ulong)qState) != n + 1) continue;
                indices.Add(basis.FlatIndex(p, qState));
            }
            lists.Add(indices.ToArray());
        }
        return lists;
    }

    /// <summary>F_a-mode σ extraction via targeted inverse iteration. For each Bloch index
    /// n in the S_2-anti orbit, shifts the block-L by the analytic eigenvalue
    /// λ_n = −2γ₀ + i·(J/2)·y_n, LU-factorizes (L − λ_n·I) once, and inverse-iterates for
    /// the right (TRANS='N') and left (TRANS='C') eigenvectors from the same factorization.
    /// σ_n is then the same biorthogonal quadratic form as in
    /// <see cref="ComputeAnatomyMklDirect"/>; the witness eigenvalue is the post-iteration
    /// Rayleigh quotient.</summary>
    private static IReadOnlyList<SigmaModeWitness> ComputeFaAnatomyInverseIteration(
        CoherenceBlock block, double q)
    {
        double j = q * block.GammaZero;
        var basis = block.Basis;
        int mTotal = basis.MTotal;
        int nBlock = block.N;

        // L as a column-major raw Complex[]; kept intact for per-mode cloning + Rayleigh matvec.
        Complex[] lRaw = BlockLDecomposition.BuildUniformLColumnMajorRaw(block, j);

        ComplexVector probeVec = DickeBlockProbe.Build(block);
        Complex probeConst = probeVec[0];   // uniform Dicke probe: every entry identical

        var siteOverlap = BuildSiteOverlapIndices(basis);
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);

        // 15 steps (was 3): at large N the c=2 block spectrum is dense, so fixed-shift
        // inverse iteration converges slower (the gap to the nearest non-target eigenvalue
        // shrinks). 15 steps drive the eigenvector residual well below the Vandermonde-
        // conditioning floor, so the k>=31 D_k integrality deviation is cleanly ISOLATED
        // as cond(V) amplification, not confounded with sigma-extraction error. It does
        // NOT make the k>=31 verification test pass: that ~1.5-2e-4 deviation is
        // Vandermonde extraction conditioning itself, a deliberate red signal (see
        // PredictDenominatorDeviationDiagnosticTests and docs/proofs/PROOF_F86B_OBSTRUCTION.md).
        // Each extra step is one cheap zgetrs on the already-computed LU factors.
        const int inverseIterationSteps = 15;
        const double rayleighResidualSanityTolerance = 1e-6;

        var witnesses = new List<SigmaModeWitness>();

        foreach (int n in orbit)
        {
            double yN = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
            // Analytic F_a eigenvalue: Re = −2γ₀, Im = (J/2)·y_n (F90 bridge: F89-J = F86-J/2).
            Complex shift = new Complex(-2.0 * block.GammaZero, 0.5 * j * yN);

            // (L − shift·I): clone L, subtract shift from the diagonal, LU-factorize in place.
            Complex[] shifted = (Complex[])lRaw.Clone();
            for (int r = 0; r < mTotal; r++)
                shifted[(long)r * mTotal + r] -= shift;
            var ipiv = new int[mTotal];
            MklDirect.LuFactorizeRaw(shifted, mTotal, ipiv);

            // Right eigenvector: inverse-iterate (L − shift·I) vr = vr.
            var vr = new Complex[mTotal];
            Array.Fill(vr, Complex.One);
            for (int step = 0; step < inverseIterationSteps; step++)
            {
                MklDirect.LuSolveRaw(shifted, mTotal, ipiv, vr, conjugateTranspose: false);
                NormalizeInPlace(vr);
            }

            // Left eigenvector: inverse-iterate (L − shift·I)^H vl = vl, same LU factors.
            var vl = new Complex[mTotal];
            Array.Fill(vl, Complex.One);
            for (int step = 0; step < inverseIterationSteps; step++)
            {
                MklDirect.LuSolveRaw(shifted, mTotal, ipiv, vl, conjugateTranspose: true);
                NormalizeInPlace(vl);
            }

            // Biorthogonal c0 = (vl^H · probe) / (vl^H · vr); probe is uniform.
            Complex vlSum = Complex.Zero;
            Complex vlHvr = Complex.Zero;
            for (int r = 0; r < mTotal; r++)
            {
                vlSum += vl[r];
                vlHvr += Complex.Conjugate(vl[r]) * vr[r];
            }
            Complex vlHprobe = Complex.Conjugate(vlSum) * probeConst;
            Complex c0 = vlHvr.Magnitude < 1e-300 ? Complex.Zero : vlHprobe / vlHvr;
            double overlapSq = c0.Real * c0.Real + c0.Imaginary * c0.Imaginary;

            // S-kernel diagonal: Σ_l 2·|Σ_{r ∈ siteOverlap[l]} vr[r]|² (the factor 2 mirrors
            // SpatialSumKernel.Build; see ComputeAnatomyMklDirect for the convention note).
            double sDiag = 0.0;
            foreach (int[] siteIndices in siteOverlap)
            {
                Complex proj = Complex.Zero;
                foreach (int r in siteIndices)
                    proj += vr[r];
                sDiag += 2.0 * (proj.Real * proj.Real + proj.Imaginary * proj.Imaginary);
            }
            double sigma = 0.5 * overlapSq * sDiag;

            // Rayleigh quotient λ = (vl^H · L · vr) / (vl^H · vr) for the witness eigenvalue,
            // plus a coarse residual sanity gate ‖L·vr − λ·vr‖ / ‖vr‖.
            Complex[] lvr = DenseMatVecColumnMajor(lRaw, vr, mTotal);
            Complex vlHlvr = Complex.Zero;
            for (int r = 0; r < mTotal; r++)
                vlHlvr += Complex.Conjugate(vl[r]) * lvr[r];
            Complex rayleigh = vlHvr.Magnitude < 1e-300 ? shift : vlHlvr / vlHvr;

            double residualSq = 0.0, vrNormSq = 0.0;
            for (int r = 0; r < mTotal; r++)
            {
                Complex diff = lvr[r] - rayleigh * vr[r];
                residualSq += diff.Real * diff.Real + diff.Imaginary * diff.Imaginary;
                vrNormSq += vr[r].Real * vr[r].Real + vr[r].Imaginary * vr[r].Imaginary;
            }
            double residual = Math.Sqrt(residualSq / vrNormSq);
            if (residual > rayleighResidualSanityTolerance)
                throw new InvalidOperationException(
                    $"Inverse iteration did not converge to a clean eigenvector for Bloch " +
                    $"index n={n} (N={nBlock}, Q={q:G4}): Rayleigh residual {residual:G6} " +
                    $"exceeds {rayleighResidualSanityTolerance:G1}. The analytic shift may " +
                    $"have hit a degenerate or non-isolated mode.");

            witnesses.Add(new SigmaModeWitness(
                EigenIndexAtQ: n,            // no global Evd index on this path; the Bloch index identifies the mode
                EigenvalueReal: rayleigh.Real,
                EigenvalueImag: rayleigh.Imaginary,
                ProbeOverlapSquared: overlapSq,
                SKernelDiagonal: sDiag,
                Sigma: sigma,
                BlochIndexN: n));
        }

        witnesses.Sort((a, b) => b.Sigma.CompareTo(a.Sigma));
        return witnesses;
    }

    /// <summary>2-norm normalization in place. No-op if the vector is numerically zero.</summary>
    private static void NormalizeInPlace(Complex[] v)
    {
        double sumSq = 0.0;
        for (int i = 0; i < v.Length; i++)
            sumSq += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
        if (sumSq < 1e-300) return;
        double inv = 1.0 / Math.Sqrt(sumSq);
        for (int i = 0; i < v.Length; i++)
            v[i] *= inv;
    }

    /// <summary>Dense matrix-vector product L·v for a column-major L of length n²:
    /// result[i] = Σ_j L[i,j]·v[j] = Σ_j lRaw[(long)j·n + i]·v[j].</summary>
    private static Complex[] DenseMatVecColumnMajor(Complex[] lRaw, Complex[] v, int n)
    {
        var result = new Complex[n];
        for (int col = 0; col < n; col++)
        {
            Complex vCol = v[col];
            long colBase = (long)col * n;
            for (int row = 0; row < n; row++)
                result[row] += lRaw[colBase + row] * vCol;
        }
        return result;
    }

    public override string DisplayName =>
        $"c=2 block-L sigma anatomy (N={Block.N}, Q={Q:G4}, dim={Block.Basis.MTotal}" +
        (FaModesOnly ? ", F_a-only via inverse iteration)" : ", full spectrum)");

    public override string Summary =>
        $"σ-spectrum over {SigmaSpectrum.Count} modes at Q={Q:G4} ({Tier.Label()})";

    /// <summary>Return the σ value of the F_a mode at Bloch index n, or <c>null</c>
    /// if no F_a mode was assigned that index. Convenience accessor for
    /// verification against <see cref="F89UnifiedFaClosedFormClaim.Sigma"/>.</summary>
    public double? SigmaForBlochIndex(int n)
    {
        foreach (var w in SigmaSpectrum)
            if (w.BlochIndexN == n) return w.Sigma;
        return null;
    }

    /// <summary>Extract the F_a polynomial coefficients (low-to-high degree) via
    /// Vandermonde fit through the F_a-mode σ values. Returns the polynomial
    /// P(y) such that σ_n = P(y_n) / [D_k · N² · (N-1)] up to a denominator D_k
    /// (D_k itself is NOT determined by this method; use
    /// <see cref="F89UnifiedFaClosedFormClaim.PredictDenominator"/> for that).
    /// The returned coefs are rationals as floating-point — multiply by the
    /// expected D_k and N²(N-1) and check integrality to verify the closed form.
    ///
    /// <para>F_a count must equal floor(N_block/2). Throws if fewer F_a modes are
    /// detected (could indicate an off-AT-lock case, e.g. degenerate Q).</para></summary>
    public double[] ExtractRawPolynomialCoefficients()
    {
        var fa = SigmaSpectrum
            .Where(w => w.BlochIndexN.HasValue)
            .OrderBy(w => w.BlochIndexN!.Value)
            .ToList();

        int expectedCount = Block.N / 2;
        if (fa.Count != expectedCount)
            throw new InvalidOperationException(
                $"Expected {expectedCount} F_a modes (floor(N/2) for N={Block.N}); got {fa.Count}.");

        int n = fa.Count;
        int chainN = Block.N;
        double[] yVals = new double[n];
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++)
        {
            int blochN = fa[i].BlochIndexN!.Value;
            yVals[i] = F89PathKAtLockMechanismClaim.BlochEigenvalueY(chainN, blochN);
            // P(y_i) = σ_n · D_k · N² · (N-1); we return the σ · N² · (N-1) part
            // (without D_k) so the caller can rationalize by D_k.
            rhs[i] = fa[i].Sigma * chainN * chainN * (chainN - 1);
        }

        // Vandermonde: V[i, j] = y_i^j (low-to-high powers)
        var V = Matrix<double>.Build.Dense(n, n,
            (i, j) => Math.Pow(yVals[i], j));
        var rhsVec = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(rhs);
        return V.Solve(rhsVec).ToArray();
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Q", Q, "G6");
            yield return new InspectableNode("dim (block-L)", summary: Block.Basis.MTotal.ToString());
        }
    }
}

/// <summary>Per-eigenmode sigma witness for <see cref="C2FullBlockSigmaAnatomy"/>.
/// Captures the spectral coordinate, biorthogonal probe overlap, S-kernel diagonal,
/// the composite F89 σ value, and the optional F_a Bloch-index assignment.</summary>
/// <param name="EigenIndexAtQ">Index in the L(Q) Evd output.</param>
/// <param name="EigenvalueReal">Re(λ_i). F_a modes sit at ≈ −2γ₀.</param>
/// <param name="EigenvalueImag">Im(λ_i). F_a modes match J·y_n on the S_2-anti orbit.</param>
/// <param name="ProbeOverlapSquared"><c>|c0_i|² = |R⁻¹·probe|_i²</c>.</param>
/// <param name="SKernelDiagonal"><c>(R†·S·R)[i,i].Real</c>. Stores the genuine S-kernel
/// diagonal (with the factor-2 that <see cref="Probes.SpatialSumKernel.Build"/> bakes in).
/// The F89 σ convention omits that factor, so the relationship is
/// <c>Sigma = 0.5 · ProbeOverlapSquared · SKernelDiagonal</c>.</param>
/// <param name="Sigma">The composite F89 σ value: <c>0.5 · ProbeOverlapSquared · SKernelDiagonal</c>.
/// The 0.5 undoes the factor-2 that <see cref="Probes.SpatialSumKernel.Build"/> bakes into
/// the spatial-sum kernel; the F89 σ convention (per <c>simulations/_f89_pathk_survey.py</c>)
/// omits that factor.</param>
/// <param name="BlochIndexN">Assigned S_2-anti Bloch index n if this is an F_a mode
/// (Re(λ) ≈ −2γ₀ and Im(λ)/J matches y_n); <c>null</c> otherwise.</param>
public sealed record SigmaModeWitness(
    int EigenIndexAtQ,
    double EigenvalueReal,
    double EigenvalueImag,
    double ProbeOverlapSquared,
    double SKernelDiagonal,
    double Sigma,
    int? BlochIndexN
) : IInspectable
{
    public string DisplayName =>
        $"λ_{EigenIndexAtQ} = {EigenvalueReal:+0.0000;-0.0000} {EigenvalueImag:+0.0000i;-0.0000i}" +
        (BlochIndexN.HasValue ? $" [F_a, n={BlochIndexN.Value}]" : "");

    public string Summary =>
        $"σ = {Sigma:G6}, |c|² = {ProbeOverlapSquared:G4}, (R†SR)_ii = {SKernelDiagonal:G4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("EigenvalueReal", EigenvalueReal, "F4");
            yield return InspectableNode.RealScalar("EigenvalueImag", EigenvalueImag, "F4");
            yield return InspectableNode.RealScalar("ProbeOverlapSquared", ProbeOverlapSquared, "G6");
            yield return InspectableNode.RealScalar("SKernelDiagonal", SKernelDiagonal, "G6");
            yield return InspectableNode.RealScalar("Sigma", Sigma, "G6");
            if (BlochIndexN.HasValue)
                yield return new InspectableNode("BlochIndexN", summary: BlochIndexN.Value.ToString());
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"Sigma (eigenmode {EigenIndexAtQ})", Sigma, "G6");
}
