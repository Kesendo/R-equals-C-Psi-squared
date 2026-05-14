using System.Numerics;
using System.Runtime.InteropServices;

namespace RCPsiSquared.Core.Numerics;

/// <summary>
/// Direct LAPACK P/Invoke with pinned arrays to bypass the .NET 2GB marshaling limit.
/// LP64 (int32) for N≤7, ILP64 (int64) for N=8+ where n² overflows int.MaxValue.
/// Fallback chain for zgeev_: libMathNetNumericsMKL → mkl_rt.2 → libopenblas.
/// </summary>
public static class MklDirect
{
    // Threshold: n*n > int.MaxValue requires ILP64. sqrt(2^31-1) ≈ 46340.
    private const int Ilp64Threshold = 46340;

    private enum ZgeevBackend { Unknown, MathNet, MklRt, OpenBlas }
    private static ZgeevBackend _resolvedBackend = ZgeevBackend.Unknown;
    private static bool _threadsConfigured;

    // --- MathNet wrapper (computes eigenvectors, needs 3x memory) ---
    [DllImport("libMathNetNumericsMKL", CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
    private static extern unsafe int z_eigen(
        [MarshalAs(UnmanagedType.Bool)] bool isSymmetric,
        int n,
        Complex* a,
        Complex* vectors,
        Complex* values,
        Complex* d);

    // ========== LP64 backends (int32 parameters, for n ≤ 46340) ==========

    [DllImport("libMathNetNumericsMKL", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "zgeev_", ExactSpelling = true)]
    private static extern unsafe void zgeev_mathnet(
        byte* jobvl, byte* jobvr, int* n,
        Complex* a, int* lda, Complex* w,
        Complex* vl, int* ldvl, Complex* vr, int* ldvr,
        Complex* work, int* lwork, double* rwork, int* info);

    [DllImport("mkl_rt.2", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "zgeev_", ExactSpelling = true)]
    private static extern unsafe void zgeev_mkl(
        byte* jobvl, byte* jobvr, int* n,
        Complex* a, int* lda, Complex* w,
        Complex* vl, int* ldvl, Complex* vr, int* ldvr,
        Complex* work, int* lwork, double* rwork, int* info);

    [DllImport("libopenblas", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "zgeev_", ExactSpelling = true)]
    private static extern unsafe void zgeev_openblas(
        byte* jobvl, byte* jobvr, int* n,
        Complex* a, int* lda, Complex* w,
        Complex* vl, int* ldvl, Complex* vr, int* ldvr,
        Complex* work, int* lwork, double* rwork, int* info);

    // zgesv: solve A x = b via LU factorization (LP64)
    [DllImport("libopenblas", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "zgesv_", ExactSpelling = true)]
    private static extern unsafe void zgesv_openblas(
        int* n, int* nrhs, Complex* a, int* lda,
        int* ipiv, Complex* b, int* ldb, int* info);

    // ========== OpenBLAS threading control ==========
    // Both LP64 and ILP64 builds need their own thread config.

    [DllImport("libopenblas", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "openblas_set_num_threads")]
    private static extern void openblas_set_num_threads_lp64(int num_threads);

    [DllImport("libopenblas", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "openblas_get_num_threads")]
    private static extern int openblas_get_num_threads_lp64();

    [DllImport("libopenblas64", CallingConvention = CallingConvention.Cdecl)]
    private static extern void openblas_set_num_threads(int num_threads);

    [DllImport("libopenblas64", CallingConvention = CallingConvention.Cdecl)]
    private static extern int openblas_get_num_threads();

    /// <summary>
    /// Configure threading for BOTH OpenBLAS builds (LP64 + ILP64).
    /// </summary>
    public static void ConfigureThreads(int numThreads, Action<string>? log = null)
    {
        // LP64 (libopenblas) - used by N<=7 eigendecomposition
        try
        {
            openblas_set_num_threads_lp64(numThreads);
            int actual = openblas_get_num_threads_lp64();
            log?.Invoke($"OpenBLAS LP64 threads: requested={numThreads}, actual={actual}");
        }
        catch (Exception ex) when (ex is DllNotFoundException or EntryPointNotFoundException)
        {
            log?.Invoke($"OpenBLAS LP64 thread config skipped ({ex.GetType().Name})");
        }

        // ILP64 (libopenblas64) - used by N=8
        try
        {
            openblas_set_num_threads(numThreads);
            int actual = openblas_get_num_threads();
            log?.Invoke($"OpenBLAS ILP64 threads: requested={numThreads}, actual={actual}");
        }
        catch (Exception ex) when (ex is DllNotFoundException or EntryPointNotFoundException)
        {
            log?.Invoke($"OpenBLAS ILP64 thread config skipped ({ex.GetType().Name})");
        }
    }

    /// <summary>Lazily configure OpenBLAS to use all logical processors, once, on the first
    /// zgeev dispatch. When <see cref="MklDirect"/> was lifted out of RCPsiSquared.Compute
    /// (whose startup called <see cref="ConfigureThreads"/> explicitly) into Core, callers
    /// like C2FullBlockSigmaAnatomy had no equivalent startup hook, leaving the direct-LAPACK
    /// path at OpenBLAS's default thread count. This makes the path self-sufficient. No-op for
    /// the MKL backends (their threading is governed separately) and harmless if libopenblas
    /// is absent. Note: this only accelerates the BLAS-3 phases of zgeev (Hessenberg reduction,
    /// eigenvector back-transform); the QR/Schur iteration is sequential by nature.</summary>
    private static void EnsureThreadsConfigured()
    {
        if (_threadsConfigured) return;
        ConfigureThreads(Environment.ProcessorCount);
        _threadsConfigured = true;
    }

    // ========== ILP64 backend (int64 parameters, for n > 46340) ==========

    [DllImport("libopenblas64", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "zgeev_", ExactSpelling = true)]
    private static extern unsafe void zgeev_openblas64(
        byte* jobvl, byte* jobvr, long* n,
        Complex* a, long* lda, Complex* w,
        Complex* vl, long* ldvl, Complex* vr, long* ldvr,
        Complex* work, long* lwork, double* rwork, long* info);

    // ========== LP64 dispatch (unchanged from before) ==========

    private static unsafe void CallZgeev(
        byte* jobvl, byte* jobvr, int* n,
        Complex* a, int* lda, Complex* w,
        Complex* vl, int* ldvl, Complex* vr, int* ldvr,
        Complex* work, int* lwork, double* rwork, int* info,
        Action<string>? log)
    {
        EnsureThreadsConfigured();
        if (_resolvedBackend != ZgeevBackend.Unknown)
        {
            DispatchZgeev(_resolvedBackend, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            return;
        }

        try
        {
            zgeev_mathnet(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            _resolvedBackend = ZgeevBackend.MathNet;
            log?.Invoke("zgeev_ backend: libMathNetNumericsMKL");
            return;
        }
        catch (EntryPointNotFoundException) { log?.Invoke("zgeev_ not in libMathNetNumericsMKL"); }
        catch (DllNotFoundException) { log?.Invoke("libMathNetNumericsMKL not found"); }

        try
        {
            zgeev_mkl(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            _resolvedBackend = ZgeevBackend.MklRt;
            log?.Invoke("zgeev_ backend: mkl_rt.2");
            return;
        }
        catch (EntryPointNotFoundException) { log?.Invoke("zgeev_ not in mkl_rt.2"); }
        catch (DllNotFoundException) { log?.Invoke("mkl_rt.2 not found"); }

        try
        {
            zgeev_openblas(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
            _resolvedBackend = ZgeevBackend.OpenBlas;
            log?.Invoke("zgeev_ backend: libopenblas (OpenBLAS LP64)");
            return;
        }
        catch (EntryPointNotFoundException) { log?.Invoke("zgeev_ not in libopenblas"); }
        catch (DllNotFoundException) { log?.Invoke("libopenblas not found"); }

        throw new InvalidOperationException(
            "Cannot find LAPACK zgeev_ in any available library. " +
            "Tried: libMathNetNumericsMKL, mkl_rt.2, libopenblas.");
    }

    private static unsafe void DispatchZgeev(ZgeevBackend backend,
        byte* jobvl, byte* jobvr, int* n,
        Complex* a, int* lda, Complex* w,
        Complex* vl, int* ldvl, Complex* vr, int* ldvr,
        Complex* work, int* lwork, double* rwork, int* info)
    {
        switch (backend)
        {
            case ZgeevBackend.MathNet:
                zgeev_mathnet(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
                break;
            case ZgeevBackend.MklRt:
                zgeev_mkl(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
                break;
            case ZgeevBackend.OpenBlas:
                zgeev_openblas(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
                break;
            default:
                throw new InvalidOperationException($"Unknown zgeev backend: {backend}");
        }
    }

    // ========== Public API ==========

    /// <summary>
    /// Eigenvalues WITH eigenvectors from column-major Complex[] array.
    /// Uses MathNet z_eigen wrapper. Needs 3 × n×n arrays.
    /// Good for N ≤ 7.
    /// </summary>
    public static unsafe Complex[] EigenvaluesRaw(Complex[] a, int n)
    {
        var vectors = new Complex[(long)n * n];
        var values = new Complex[n];
        var d = new Complex[(long)n * n];

        fixed (Complex* pA = a)
        fixed (Complex* pVectors = vectors)
        fixed (Complex* pValues = values)
        fixed (Complex* pD = d)
        {
            int result = z_eigen(false, n, pA, pVectors, pValues, pD);
            if (result != 0)
                throw new InvalidOperationException($"MKL z_eigen failed with code {result}");
        }

        return values;
    }

    /// <summary>
    /// Eigenvalues AND right eigenvectors from column-major Complex[] array.
    /// vectors[j * n + i] = i-th component of the j-th eigenvector.
    /// Uses MathNet z_eigen wrapper. Needs 3 × n×n arrays.
    /// Good for N ≤ 7.
    /// </summary>
    public static unsafe (Complex[] values, Complex[] vectors) EigenvaluesAndVectorsRaw(Complex[] a, int n)
    {
        var vectors = new Complex[(long)n * n];
        var values = new Complex[n];
        var d = new Complex[(long)n * n];

        fixed (Complex* pA = a)
        fixed (Complex* pVectors = vectors)
        fixed (Complex* pValues = values)
        fixed (Complex* pD = d)
        {
            int result = z_eigen(false, n, pA, pVectors, pValues, pD);
            if (result != 0)
                throw new InvalidOperationException($"MKL z_eigen failed with code {result}");
        }

        return (values, vectors);
    }

    /// <summary>
    /// Eigenvalues ONLY from column-major Complex[] managed array (LP64).
    /// For N ≤ 7 where n² fits in int.MaxValue.
    /// WARNING: input array 'a' is destroyed by LAPACK.
    /// </summary>
    public static unsafe Complex[] EigenvaluesOnlyRaw(Complex[] a, int n, Action<string>? log = null)
    {
        var w = new Complex[n];
        var rwork = new double[2 * n];
        var workQuery = new Complex[1];
        int lwork = -1;
        int info = 0;
        int ldvl = 1, ldvr = 1;
        byte jobN = (byte)'N';
        var vl = new Complex[1];
        var vr = new Complex[1];

        log?.Invoke("LAPACK zgeev workspace query...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = workQuery)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobN, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, log);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev workspace query failed: info={info}");

        int optimalWork = (int)workQuery[0].Real;
        log?.Invoke($"Optimal workspace: {optimalWork} complex values ({optimalWork * 16.0 / 1e6:F1} MB)");
        var work = new Complex[optimalWork];
        lwork = optimalWork;

        log?.Invoke($"Running zgeev eigenvalue-only on {n}x{n} matrix...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = work)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobN, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, null);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev failed: info={info}");

        log?.Invoke($"zgeev complete. {n} eigenvalues computed.");
        return w;
    }

    /// <summary>
    /// Eigenvalues AND right eigenvectors via direct LAPACK zgeev (LP64).
    /// Bypasses MathNet z_eigen wrapper which has internal allocation issues at large n.
    /// WARNING: input array 'a' is destroyed by LAPACK.
    /// </summary>
    public static unsafe (Complex[] values, Complex[] vectors) EigenvaluesAndVectorsDirectRaw(
        Complex[] a, int n, Action<string>? log = null)
    {
        var w = new Complex[n];
        var rwork = new double[2 * n];
        var workQuery = new Complex[1];
        int lwork = -1;
        int info = 0;
        int ldvl = 1;
        int ldvr = n;
        byte jobN = (byte)'N';
        byte jobV = (byte)'V';
        var vl = new Complex[1];
        var vr = new Complex[(long)n * n];

        log?.Invoke("LAPACK zgeev (with eigenvectors) workspace query...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = workQuery)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobV, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, log);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev workspace query failed: info={info}");

        int optimalWork = (int)workQuery[0].Real;
        log?.Invoke($"Optimal workspace: {optimalWork} complex values ({optimalWork * 16.0 / 1e6:F1} MB)");
        var work = new Complex[optimalWork];
        lwork = optimalWork;

        log?.Invoke($"Running zgeev with eigenvectors on {n}x{n} matrix...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = work)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobV, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, null);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev failed: info={info}");

        log?.Invoke($"zgeev complete. {n} eigenvalues + eigenvectors computed.");
        return (w, vr);
    }

    /// <summary>
    /// Eigenvalues AND both left+right eigenvectors via direct LAPACK zgeev (LP64).
    /// JOBVL='V', JOBVR='V'. Returns (values, leftVectors, rightVectors).
    /// Left eigenvectors u satisfy u^H * A = λ * u^H (stored column-major as rows of VL^H).
    /// Right eigenvectors v satisfy A * v = λ * v (stored column-major).
    /// WARNING: input array 'a' is destroyed by LAPACK. Needs ~4 × n² memory for n=16384.
    /// </summary>
    public static unsafe (Complex[] values, Complex[] leftVectors, Complex[] rightVectors)
        EigenvaluesLeftRightDirectRaw(Complex[] a, int n, Action<string>? log = null)
    {
        var w = new Complex[n];
        var rwork = new double[2 * n];
        var workQuery = new Complex[1];
        int lwork = -1;
        int info = 0;
        int ldvl = n;
        int ldvr = n;
        byte jobV = (byte)'V';
        var vl = new Complex[(long)n * n];
        var vr = new Complex[(long)n * n];

        log?.Invoke("LAPACK zgeev (left+right eigenvectors) workspace query...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = workQuery)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobV, &jobV, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, log);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev workspace query failed: info={info}");

        int optimalWork = (int)workQuery[0].Real;
        log?.Invoke($"Optimal workspace: {optimalWork} complex values ({optimalWork * 16.0 / 1e6:F1} MB)");
        var work = new Complex[optimalWork];
        lwork = optimalWork;

        log?.Invoke($"Running zgeev with left+right eigenvectors on {n}x{n} matrix...");

        fixed (Complex* pA = a)
        fixed (Complex* pW = w)
        fixed (Complex* pVL = vl)
        fixed (Complex* pVR = vr)
        fixed (Complex* pWork = work)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobV, &jobV, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, null);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev failed: info={info}");

        log?.Invoke($"zgeev complete. {n} eigenvalues + left + right eigenvectors computed.");
        return (w, vl, vr);
    }

    /// <summary>
    /// Solve R^T y = e_k directly via LAPACK zgesv (LP64).
    /// R is column-major in vectors[]. We build R^T column-major by transposing indices.
    /// WARNING: the transposed copy is destroyed by LAPACK.
    /// </summary>
    public static unsafe Complex[] SolveRTranspose(
        Complex[] vectors, int d2, int targetIndex, Action<string>? log = null)
    {
        int n = d2; // local copy for fixed block (can't take address of parameter used in lambda)

        // Build R^T in column-major: RT[col*n + row] = R^T[row,col] = R[col,row] = vectors[row*n + col]
        log?.Invoke($"Building R^T ({n}x{n}, {(long)n * n * 16 / 1e9:F1} GB)...");
        var rt = new Complex[(long)n * n];
        int nc = n; // capture for lambda
        Parallel.For(0, n, col =>
        {
            long colOff = (long)col * nc;
            for (int row = 0; row < nc; row++)
                rt[colOff + row] = vectors[(long)row * nc + col];
        });

        // RHS: e_k
        var b = new Complex[n];
        b[targetIndex] = Complex.One;

        var ipiv = new int[n];
        int nrhs = 1;
        int info = 0;

        log?.Invoke($"LAPACK zgesv ({n}x{n})...");
        fixed (Complex* pA = rt)
        fixed (Complex* pB = b)
        fixed (int* pIpiv = ipiv)
        {
            zgesv_openblas(&n, &nrhs, pA, &n, pIpiv, pB, &n, &info);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgesv failed: info={info}");

        log?.Invoke("zgesv complete.");
        return b;
    }

    /// <summary>
    /// Eigenvalues ONLY from native pointer (unmanaged memory).
    /// Automatically selects LP64 or ILP64 based on matrix dimension.
    /// ILP64 required when n > 46340 (n² overflows int32).
    /// </summary>
    public static unsafe Complex[] EigenvaluesOnlyNative(IntPtr matrixPtr, int n, Action<string>? log = null)
    {
        if (n > Ilp64Threshold)
            return EigenvaluesOnlyNativeIlp64(matrixPtr, n, log);

        // LP64 path (same as before)
        var w = new Complex[n];
        var rwork = new double[2 * n];
        var workQuery = new Complex[1];
        int lwork = -1;
        int info = 0;
        int ldvl = 1, ldvr = 1;
        byte jobN = (byte)'N';
        var vl = new Complex[1];
        var vr = new Complex[1];
        var pA = (Complex*)matrixPtr;

        log?.Invoke($"LAPACK zgeev LP64 workspace query (n={n})...");

        fixed (Complex* pW = w, pVL = vl, pVR = vr, pWork = workQuery)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobN, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, log);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev workspace query failed: info={info}");

        int optimalWork = (int)workQuery[0].Real;
        log?.Invoke($"Optimal workspace: {optimalWork} ({optimalWork * 16.0 / 1e6:F1} MB)");
        var work = new Complex[optimalWork];
        lwork = optimalWork;

        log?.Invoke($"Running zgeev eigenvalue-only on {n}x{n} matrix...");

        fixed (Complex* pW = w, pVL = vl, pVR = vr, pWork = work)
        fixed (double* pRwork = rwork)
        {
            CallZgeev(&jobN, &jobN, &n, pA, &n, pW, pVL, &ldvl, pVR, &ldvr,
                       pWork, &lwork, pRwork, &info, null);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev failed: info={info}");

        log?.Invoke($"zgeev complete. {n} eigenvalues computed.");
        return w;
    }

    /// <summary>
    /// ILP64 path: all LAPACK integer parameters are 64-bit.
    /// Required when n > 46340 (n² > int.MaxValue).
    /// Uses libopenblas64.dll (OpenBLAS built with USE64BITINT).
    /// </summary>
    public static unsafe Complex[] EigenvaluesOnlyNativeIlp64(IntPtr matrixPtr, int n, Action<string>? log = null)
    {
        EnsureThreadsConfigured();
        log?.Invoke($"Using ILP64 (64-bit int) path: n={n}, n²={(long)n * n:N0} > int.MaxValue");

        var w = new Complex[n];
        var rwork = new double[2 * n];
        var workQuery = new Complex[1];
        long lwork = -1;
        long info = 0;
        long n64 = n;
        long ldvl = 1, ldvr = 1;
        byte jobN = (byte)'N';
        var vl = new Complex[1];
        var vr = new Complex[1];
        var pA = (Complex*)matrixPtr;

        log?.Invoke("LAPACK zgeev ILP64 workspace query...");

        fixed (Complex* pW = w, pVL = vl, pVR = vr, pWork = workQuery)
        fixed (double* pRwork = rwork)
        {
            zgeev_openblas64(&jobN, &jobN, &n64, pA, &n64, pW,
                             pVL, &ldvl, pVR, &ldvr,
                             pWork, &lwork, pRwork, &info);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev ILP64 workspace query failed: info={info}");

        long optimalWork = (long)workQuery[0].Real;
        log?.Invoke($"Optimal workspace: {optimalWork} complex values ({optimalWork * 16.0 / 1e6:F1} MB)");
        var work = new Complex[optimalWork];
        lwork = optimalWork;

        log?.Invoke($"Running zgeev ILP64 eigenvalue-only on {n}x{n} matrix...");

        fixed (Complex* pW = w, pVL = vl, pVR = vr, pWork = work)
        fixed (double* pRwork = rwork)
        {
            zgeev_openblas64(&jobN, &jobN, &n64, pA, &n64, pW,
                             pVL, &ldvl, pVR, &ldvr,
                             pWork, &lwork, pRwork, &info);
        }

        if (info != 0)
            throw new InvalidOperationException($"zgeev ILP64 failed: info={info}");

        log?.Invoke($"zgeev ILP64 complete. {n} eigenvalues computed.");
        return w;
    }
}
