using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The spectator-intertwiner gate (codim-1-by-additivity plan, Task B1): verify from trusted
/// builders that the site-summed spectator W(ρ) = Σ_l c_l†ρc_l (JW strings included) EXACTLY intertwines
/// both parts of the per-block Liouvillian pencil L(q) = A + q·C, block-shifting (p,q̃) → (p+1,q̃+1);
/// record the documented refutation of the single-mode reading V_k(ρ) = η_k†ρη_k (C-part exact, A-part
/// O(1) for every k); and measure the kernel-death boundary numbers at the N=5 defective locus
/// q* = 0.620878 (does W kill the transported Jordan pair at the diamond-boundary step (3,3)→(4,4)?).
///
/// <para><b>Conventions</b> (all read from source): blocks built exactly as
/// <c>SectorBraidModeGeometry.BuildBlock</c> — H = XYChain(N, 2q) so J = 2q, γ = 1 per site,
/// basis = <see cref="BlockBasis.PopcountStates"/> ascending, flat index = pIdx·Mq + qIdx, Liouville flat
/// = ket·2^N + bra via <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>. Site indexing is big-endian
/// (site 0 = MSB, site l = bit N−1−l), matching <see cref="JwBlockBasis"/>/<see
/// cref="XyJordanWignerModes"/>; the JW string sign is s_l(a) = (−1)^(# occupied sites m &lt; l in a),
/// i.e. popcount of the bits strictly above bit N−1−l. Chain bonds connect site-adjacent (= bit-adjacent)
/// sites, so H is quadratic in these site fermions and Lemma 1 applies.</para>
///
/// <para><b>Normalization</b> (scale-honest, per review): residual(X) = ‖X₂W − WX₁‖_F / (‖X₁‖₂·‖W‖_F)
/// with ‖X₁‖₂ the SPECTRAL norm (largest singular value, via <see cref="EpCharacter.SpectralNorm"/>).</para>
///
/// <para>Items 4a/4b are printed measurements without asserts (first run; the genuinely open number is
/// whether the (3,3) Jordan pair survives W into the census-braid-free (4,4) block).</para></summary>
public class SpectatorIntertwinerGateTests
{
    private const int N = 5;
    private const double QDefective = 0.620878;   // real defective locus 1 of the (1,2) block, N=5
    private const double QGeneric = 0.37;
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly ITestOutputHelper _out;
    public SpectatorIntertwinerGateTests(ITestOutputHelper o) => _out = o;

    // ---------------------------------------------------------------- shared construction helpers
    // The W construction and its bit/JW-sign conventions live in the shared Diagnostics helper
    // SpectatorIntertwiner (one construction for the gate AND the live SectorBraidWitness node 2);
    // the thin wrappers below just fix N.

    /// <summary>Site l occupied in computational-basis state (big-endian: site 0 = MSB = bit N−1).</summary>
    private static bool Occupied(long state, int site) => SpectatorIntertwiner.Occupied(N, state, site);

    private static long SiteMask(int site) => SpectatorIntertwiner.SiteMask(N, site);

    /// <summary>JW string sign s_l(a) = (−1)^(# occupied sites strictly before site l); see
    /// <see cref="SpectatorIntertwiner.JwSign"/> for the bit-order convention.</summary>
    private static int JwSign(long state, int site) => SpectatorIntertwiner.JwSign(N, state, site);

    /// <summary>Build the (p,q̃) joint-popcount block of L at real q, exactly as
    /// <c>SectorBraidModeGeometry.BuildBlock</c> does for (1,2): H = XYChain(N, 2q) (or an explicit
    /// override for the Δ-gate), γ=1 per site, flat = pIdx·Mq + qIdx, Liouville flat = ket·d + bra.</summary>
    private static ComplexMatrix BuildBlock(int p, int qt, double q, ComplexMatrix? hOverride = null)
    {
        int d = 1 << N;
        var statesP = BlockBasis.PopcountStates(N, p);
        var statesQ = BlockBasis.PopcountStates(N, qt);
        var H = hOverride ?? PauliHamiltonian.XYChain(N, 2.0 * q).ToMatrix();
        var gamma = Enumerable.Repeat(1.0, N).ToList();
        var flat = new int[statesP.Count * statesQ.Count];
        for (int i = 0; i < statesP.Count; i++)
            for (int j = 0; j < statesQ.Count; j++)
                flat[i * statesQ.Count + j] = (int)(statesP[i] * d + statesQ[j]);
        return PerBlockLiouvillianBuilder.BuildBlockZ(H, gamma, flat);
    }

    /// <summary>The site-summed spectator W: block(p,q̃) → block(p+1,q̃+1). Delegates to the shared
    /// <see cref="SpectatorIntertwiner.BuildW"/> (W(ρ) = Σ_l c_l†ρc_l restricted to the block; the
    /// bra-side JW sign is real).</summary>
    private static ComplexMatrix BuildW(int p, int qt) => SpectatorIntertwiner.BuildW(N, p, qt);

    /// <summary>The single-mode spectator V_k(ρ) = η_k†ρη_k with η_k† = Σ_l ψ_k(l)c_l†:
    /// V_k(|a⟩⟨b|) = Σ_{l∉a, m∉b} ψ_k(l)ψ_k(m)·s_l(a)s_m(b)·|a+e_l⟩⟨b+e_m| (independent ket/bra sites;
    /// ψ real so no conjugation issue). The provably-refuted single-mode reading.</summary>
    private static ComplexMatrix BuildVk(int p, int qt, int k, XyJordanWignerModes modes)
    {
        var inP = BlockBasis.PopcountStates(N, p);
        var inQ = BlockBasis.PopcountStates(N, qt);
        var outP = BlockBasis.PopcountStates(N, p + 1);
        var outQ = BlockBasis.PopcountStates(N, qt + 1);
        var outPIdx = new Dictionary<long, int>();
        for (int i = 0; i < outP.Count; i++) outPIdx[outP[i]] = i;
        var outQIdx = new Dictionary<long, int>();
        for (int i = 0; i < outQ.Count; i++) outQIdx[outQ[i]] = i;

        var v = Matrix<Complex>.Build.Dense(outP.Count * outQ.Count, inP.Count * inQ.Count);
        for (int i = 0; i < inP.Count; i++)
            for (int j = 0; j < inQ.Count; j++)
            {
                long a = inP[i], b = inQ[j];
                int col = i * inQ.Count + j;
                for (int l = 0; l < N; l++)
                {
                    if (Occupied(a, l)) continue;
                    for (int m = 0; m < N; m++)
                    {
                        if (Occupied(b, m)) continue;
                        int row = outPIdx[a | SiteMask(l)] * outQ.Count + outQIdx[b | SiteMask(m)];
                        v[row, col] += modes.SineMode(k, l) * modes.SineMode(k, m) * JwSign(a, l) * JwSign(b, m);
                    }
                }
            }
        return v;
    }

    /// <summary>Scale-honest intertwining residual ‖X₂W − WX₁‖_F / (‖X₁‖₂·‖W‖_F), ‖·‖₂ = spectral norm.</summary>
    private static double Residual(ComplexMatrix x2, ComplexMatrix w, ComplexMatrix x1)
        => (x2 * w - w * x1).FrobeniusNorm() / (EpCharacter.SpectralNorm(x1) * w.FrobeniusNorm());

    /// <summary>The two eigenvalues nearest a target center. NOT the global minimal-gap pair: at N≥5 the
    /// blocks carry exactly degenerate SEMISIMPLE pairs (e.g. λ = −6 in (1,2) at q*, gap ~1e-15) that mask
    /// the defective EP under global closest-pair clustering — the recorded trap
    /// (reference_at_masking_defective_anywhere_probe); the defective center must be seeded, exactly as
    /// <c>SectorBraidBranchSweepTests</c> seeds λ_start = −4.6189 at q* = 0.620878.</summary>
    private static (Complex E1, Complex E2, double Gap) NearestPair(ComplexVector eigs, Complex target)
    {
        var idx = Enumerable.Range(0, eigs.Count).OrderBy(i => (eigs[i] - target).Magnitude).Take(2).ToArray();
        return (eigs[idx[0]], eigs[idx[1]], (eigs[idx[0]] - eigs[idx[1]]).Magnitude);
    }

    /// <summary>The unit-norm smallest right singular vector of (L − λI) — the numerical eigenvector at
    /// a defective (near-)coalescence — plus the two smallest singular values (printed for honesty).</summary>
    private static (ComplexVector X1, double SMin, double S2) DefectiveEigenvector(ComplexMatrix l, Complex lambda)
    {
        var shifted = l - lambda * Matrix<Complex>.Build.DenseIdentity(l.RowCount);
        var svd = shifted.Svd(true);
        int n = svd.S.Count;
        var x1 = svd.VT.Row(n - 1).Conjugate();
        x1 = x1.Multiply(new Complex(1.0 / x1.L2Norm(), 0));
        return (x1, svd.S[n - 1].Real, svd.S[n - 2].Real);
    }

    /// <summary>Least-squares generalized vector: minimal-norm x₂ with (L − λ)x₂ ≈ x₁, computed by SVD
    /// pseudo-inverse dropping the smallest singular direction (the near-kernel = x₁'s own direction)
    /// plus anything below 1e-12·σ_max; then projected orthogonal to x₁ and unit-normalized.</summary>
    private static ComplexVector GeneralizedVector(ComplexMatrix l, Complex lambda, ComplexVector x1)
    {
        var shifted = l - lambda * Matrix<Complex>.Build.DenseIdentity(l.RowCount);
        var svd = shifted.Svd(true);
        int n = svd.S.Count;
        double sMax = svd.S[0].Real;
        var y = svd.U.ConjugateTranspose() * x1;
        var z = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(n);
        for (int i = 0; i < n - 1; i++)                       // drop index n−1 = the near-kernel
            if (svd.S[i].Real > 1e-12 * sMax)
                z[i] = y[i] / svd.S[i];
        var x2 = svd.VT.ConjugateTranspose() * z;
        Complex ip = Complex.Zero;                            // ⟨x₁|x₂⟩ with x₁ unit-norm
        for (int i = 0; i < x2.Count; i++) ip += Complex.Conjugate(x1[i]) * x2[i];
        x2 -= x1.Multiply(ip);
        return x2.Multiply(new Complex(1.0 / x2.L2Norm(), 0));
    }

    private static string F(double v) => v.ToString("E3", Inv);

    // ------------------------------------------------------------------------------------- item 1 + 3

    [Fact(DisplayName = "B1 items 1+3: W exactly intertwines A and C on (1,2)→(2,3); σ_min(W) = √2")]
    [Trait("Category", "SLOW_MSM")]
    public void W_Intertwines_BothPencilParts_And_IsInjective()
    {
        var w = BuildW(1, 2);
        var a1 = BuildBlock(1, 2, 0.0);                       // A = L(0) = −2·diag(n_diff)
        var a2 = BuildBlock(2, 3, 0.0);

        double resA = Residual(a2, w, a1);
        _out.WriteLine($"[item 1] W:(1,2)→(2,3), N={N}. Normalization: ‖X₂W−WX₁‖_F/(‖X₁‖₂·‖W‖_F), ‖·‖₂ = spectral norm.");
        _out.WriteLine($"  A-part residual (q-independent)      : {F(resA)}");

        ComplexMatrix? cPrev = null;
        foreach (double q in new[] { QDefective, QGeneric })
        {
            var c1 = (BuildBlock(1, 2, q) - a1).Divide(new Complex(q, 0));   // C = (L(q) − L(0))/q
            var c2 = (BuildBlock(2, 3, q) - a2).Divide(new Complex(q, 0));
            double resC = Residual(c2, w, c1);
            _out.WriteLine($"  C-part residual at q={q.ToString("0.######", Inv)}     : {F(resC)}");
            Assert.True(resC <= 1e-10, $"C-part residual {resC:E3} at q={q} exceeds 1e-10");
            if (cPrev is not null)
                _out.WriteLine($"  linearity witness ‖C(q₁)−C(q₂)‖_F    : {F((c1 - cPrev).FrobeniusNorm())}");
            cPrev = c1;
        }
        Assert.True(resA <= 1e-10, $"A-part residual {resA:E3} exceeds 1e-10");

        // item 3: injectivity on the climbing rung.
        var s = BuildW(1, 2).Svd().S;
        double sMin = s[s.Count - 1].Real;
        _out.WriteLine($"[item 3] σ_min(W) on (1,2)→(2,3)       : {sMin.ToString("0.###############", Inv)}  (predicted √2 = {Math.Sqrt(2.0).ToString("0.###############", Inv)})");
        Assert.True(Math.Abs(sMin - Math.Sqrt(2.0)) <= 1e-10, $"σ_min(W) = {sMin:R} deviates from √2");
    }

    // ----------------------------------------------------------------------------------------- item 2

    [Fact(DisplayName = "B1 item 2: single-mode V_k is NOT an intertwiner — C-part exact, A-part O(1), every k")]
    [Trait("Category", "SLOW_MSM")]
    public void SingleMode_Vk_Refuted()
    {
        var modes = XyJordanWignerModes.Build(N, 1.0);        // ψ_k only; ε_k never enters V_k
        var a1 = BuildBlock(1, 2, 0.0);
        var a2 = BuildBlock(2, 3, 0.0);
        var c1Def = (BuildBlock(1, 2, QDefective) - a1).Divide(new Complex(QDefective, 0));
        var c2Def = (BuildBlock(2, 3, QDefective) - a2).Divide(new Complex(QDefective, 0));
        var c1Gen = (BuildBlock(1, 2, QGeneric) - a1).Divide(new Complex(QGeneric, 0));
        var c2Gen = (BuildBlock(2, 3, QGeneric) - a2).Divide(new Complex(QGeneric, 0));

        _out.WriteLine($"[item 2] V_k on (1,2)→(2,3), N={N}: per-k residuals (same normalization as item 1).");
        _out.WriteLine("  k | A-part      | C-part (q=0.620878) | C-part (q=0.37)");
        for (int k = 1; k <= N; k++)
        {
            var vk = BuildVk(1, 2, k, modes);
            double rA = Residual(a2, vk, a1);
            double rCd = Residual(c2Def, vk, c1Def);
            double rCg = Residual(c2Gen, vk, c1Gen);
            _out.WriteLine($"  {k} | {F(rA)}  | {F(rCd)}          | {F(rCg)}");
            Assert.True(rCd <= 1e-10 && rCg <= 1e-10, $"V_{k} C-part residual not machine zero: {rCd:E3}, {rCg:E3}");
            Assert.True(rA >= 0.01, $"V_{k} A-part residual unexpectedly small: {rA:E3} (refutation requires O(1))");
        }
    }

    // ----------------------------------------------------------------------------------------- item 4

    [Fact(DisplayName = "B1 item 4: kernel-death boundary check at q* = 0.620878 (printed measurement) + edge echo")]
    [Trait("Category", "SLOW_MSM")]
    public void KernelDeath_At_DiamondBoundary()
    {
        // (4b setup first: λ_A from the (1,2) block's defective pair at q*, seeded on the known
        // defective center −4.6189 — the SectorBraidBranchSweepTests seed. Global minimal-gap
        // clustering is NOT used: it locks onto the exactly degenerate semisimple pair at λ = −6.)
        var seed = new Complex(-4.6189, 0.0);
        var l12 = BuildBlock(1, 2, QDefective);
        var (e1A, e2A, gapA) = NearestPair(l12.Evd().EigenValues, seed);
        Complex lambdaA = (e1A + e2A) / 2.0;
        _out.WriteLine($"[item 4] q* = {QDefective.ToString("0.######", Inv)}, N={N}, defective seed λ = {seed.Real.ToString("0.####", Inv)}.");
        _out.WriteLine($"  (1,2) defective pair: λ = {e1A.ToString("G6", Inv)}, {e2A.ToString("G6", Inv)}; gap = {F(gapA)}; center λ_A = {lambdaA.ToString("G6", Inv)}");

        // 4a: the boundary step (3,3)→(4,4). λ_B = −conj(λ_A) − 2N (F89d cross-fold shift).
        Complex lambdaBPred = -Complex.Conjugate(lambdaA) - 2.0 * N;
        var l33 = BuildBlock(3, 3, QDefective);
        var (e1B, e2B, gapB) = NearestPair(l33.Evd().EigenValues, lambdaBPred);
        Complex lambdaB = (e1B + e2B) / 2.0;
        _out.WriteLine($"  λ_B predicted = −conj(λ_A) − 2N = {lambdaBPred.ToString("G6", Inv)}");
        _out.WriteLine($"  (3,3) pair near λ_B: {e1B.ToString("G6", Inv)}, {e2B.ToString("G6", Inv)}; gap = {F(gapB)}; |mean − pred| = {F((lambdaB - lambdaBPred).Magnitude)}");

        var (x1B, sMinB, s2B) = DefectiveEigenvector(l33, lambdaB);
        var x2B = GeneralizedVector(l33, lambdaB, x1B);
        var w34 = BuildW(3, 3);                               // (3,3) → (4,4), 25×100
        double wx1B = (w34 * x1B).L2Norm();
        double wx2B = (w34 * x2B).L2Norm();
        _out.WriteLine($"  (3,3)−λ_B singular floor: σ_min = {F(sMinB)}, σ_2 = {F(s2B)}");
        _out.WriteLine($"  [4a boundary (3,3)→(4,4)] ‖W·x₁‖/‖x₁‖ = {F(wx1B)}   ‖W·x₂‖/‖x₂‖ = {F(wx2B)}   (printed, no assert; census says (4,4) braid-free ⟹ NOT BOTH may survive)");

        // 4b: the interior step (1,2)→(2,3) with the (1,2) defective pair.
        var (x1A, sMinA, s2A) = DefectiveEigenvector(l12, lambdaA);
        var x2A = GeneralizedVector(l12, lambdaA, x1A);
        var w12 = BuildW(1, 2);
        double wx1A = (w12 * x1A).L2Norm();
        double wx2A = (w12 * x2A).L2Norm();
        _out.WriteLine($"  (1,2)−λ_A singular floor: σ_min = {F(sMinA)}, σ_2 = {F(s2A)}");
        _out.WriteLine($"  [4b interior (1,2)→(2,3)] ‖W·x₁‖/‖x₁‖ = {F(wx1A)}   ‖W·x₂‖/‖x₂‖ = {F(wx2A)}   (predicted both O(1): the chain transports)");

        // 4c: edge echo — W:(3,4)→(4,5) rank, and the edge block (4,5) has n_diff ≡ 1 ⟹ A = −2I.
        var w45 = BuildW(3, 4);                               // 5×50
        var s45 = w45.Svd().S;
        double sMax45 = s45[0].Real;
        int rank = s45.Count(sv => sv.Real > 1e-10 * sMax45);
        var a45 = BuildBlock(4, 5, 0.0);
        double maxDev = 0.0;
        for (int i = 0; i < a45.RowCount; i++)
            for (int j = 0; j < a45.ColumnCount; j++)
                maxDev = Math.Max(maxDev, (a45[i, j] - (i == j ? new Complex(-2, 0) : Complex.Zero)).Magnitude);
        _out.WriteLine($"  [4c edge echo] rank(W:(3,4)→(4,5)) = {rank} of min(5,50)=5; A₍₄,₅₎ vs −2I max deviation = {F(maxDev)}");
        Assert.True(maxDev <= 1e-12, $"A(4,5) deviates from −2I by {maxDev:E3}");
    }

    // ------------------------------------------------------------------------------- item 5 (Δ-gate)

    [Fact(DisplayName = "B1 item 5 (Δ-gate): under XXZ Δ=0.3 the H-part of the W identity dies, the D-part survives")]
    [Trait("Category", "SLOW_MSM")]
    public void DeltaGate_XXZ_HPartDies_DPartSurvives()
    {
        const double delta = 0.3;
        double q = QGeneric;
        var terms = PauliHamiltonian.XYChain(N, 2.0 * q).Terms.ToList();
        for (int b = 0; b < N - 1; b++)
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Z, new Complex(delta, 0)));
        var hXxz = new PauliHamiltonian(N, terms).ToMatrix(); // H = q·Σ(XX+YY) + Δ·Σ ZZ

        var a1 = BuildBlock(1, 2, 0.0);
        var a2 = BuildBlock(2, 3, 0.0);
        var h1 = BuildBlock(1, 2, 0.0, hXxz) - a1;            // the full −i[H_XXZ,·] part of the block
        var h2 = BuildBlock(2, 3, 0.0, hXxz) - a2;

        var w = BuildW(1, 2);
        double resH = Residual(h2, w, h1);
        double resD = Residual(a2, w, a1);
        _out.WriteLine($"[item 5] XXZ Δ={delta.ToString("0.##", Inv)}, q={q.ToString("0.##", Inv)}, W:(1,2)→(2,3):");
        _out.WriteLine($"  H-part residual = {F(resH)}   (predicted O(1): additivity dies with Δ)");
        _out.WriteLine($"  D-part residual = {F(resD)}   (predicted machine zero: Lemma 2 never sees H)");
        // Printed measurement per plan; no hard assert on the O(1) side (first run).
    }

    // --------------------------------------------------------- item 6 (§6 sl(2) structure of the kernel)

    /// <summary>C(n,k), the block dimension factor.</summary>
    private static int Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return (int)r;
    }

    [Fact(DisplayName = "B1 item 6 (§6): [W,W†] = N̂_bra+N̂_ket−N closes an sl(2); Lefschetz kernel dims (N=3,4,5); [L,W†]=0")]
    [Trait("Category", "SLOW_MSM")]
    public void W_ClosesSl2_Cartan_And_Lefschetz()
    {
        // The kernel death of PROOF_CODIM1_BY_ADDITIVITY §6 is highest-weight annihilation: W is the raising
        // operator of an sl(2) with Cartan H₀ = N̂_bra+N̂_ket−N (weight m = p+q̃−N per block). Verified here
        // block-by-block from the trusted BuildW: on block (p,q̃),
        //   [W,W†] = W[p−1,q̃−1]·W[p−1,q̃−1]† − W[p,q̃]†·W[p,q̃]   should equal (p+q̃−N)·I
        // (missing neighbour at a boundary contributes the zero map). And the Lefschetz count:
        //   dim ker W[p,q̃] = 0 for m<0 (injective below the anti-diagonal), else dim(p,q̃)−dim(p+1,q̃+1).
        foreach (int n in new[] { 3, 4, 5 })
        {
            for (int p = 0; p <= n; p++)
                for (int qt = 0; qt <= n; qt++)
                {
                    int dpq = Binom(n, p) * Binom(n, qt);
                    if (dpq == 0) continue;

                    var comm = Matrix<Complex>.Build.Dense(dpq, dpq);
                    if (p >= 1 && qt >= 1)
                    {
                        var wBelow = SpectatorIntertwiner.BuildW(n, p - 1, qt - 1);   // (p−1,q̃−1)→(p,q̃)
                        comm += wBelow * wBelow.ConjugateTranspose();
                    }
                    ComplexMatrix? wHere = null;
                    if (p <= n - 1 && qt <= n - 1)
                    {
                        wHere = SpectatorIntertwiner.BuildW(n, p, qt);                // (p,q̃)→(p+1,q̃+1)
                        comm -= wHere.ConjugateTranspose() * wHere;
                    }
                    double cartan = p + qt - n;
                    double dev = 0.0;
                    for (int i = 0; i < dpq; i++)
                        for (int j = 0; j < dpq; j++)
                            dev = Math.Max(dev, (comm[i, j] - (i == j ? new Complex(cartan, 0) : Complex.Zero)).Magnitude);
                    Assert.True(dev <= 1e-10, $"N={n} block ({p},{qt}): [W,W†] − (p+q̃−N)I max dev {dev:E3}");

                    if (wHere is not null)
                    {
                        int m = p + qt - n;
                        int kerPred = m < 0 ? 0 : dpq - Binom(n, p + 1) * Binom(n, qt + 1);
                        var s = wHere.Svd().S;
                        double sMax = s[0].Real;
                        int rank = s.Count(sv => sv.Real > 1e-9 * sMax);
                        int kerAct = dpq - rank;
                        Assert.True(kerAct == kerPred,
                            $"N={n} block ({p},{qt}) m={m}: ker dim {kerAct} != Lefschetz pred {kerPred}");
                    }
                }
        }
        _out.WriteLine("[item 6] sl(2) closure [W,W†] = N̂_bra+N̂_ket−N and Lefschetz kernel dims: exact/machine-zero at every block, N=3,4,5.");
        _out.WriteLine($"         (N=5 core (3,3): m=+1, ker W = {100 - 25} = dim(3,3)−dim(4,4); climbing rung (1,2): m=−2, ker = 0, injective.)");

        // L commutes with W† (reverse spectator) too, completing L ↔ the whole sl(2). Rung (2,3)→(1,2), N=5.
        var a1 = BuildBlock(1, 2, 0.0);
        var a2 = BuildBlock(2, 3, 0.0);
        var wd = BuildW(1, 2).ConjugateTranspose();                     // (2,3)→(1,2) = W†
        double resAd = (a1 * wd - wd * a2).FrobeniusNorm() / (EpCharacter.SpectralNorm(a2) * wd.FrobeniusNorm());
        var c1 = (BuildBlock(1, 2, QDefective) - a1).Divide(new Complex(QDefective, 0));
        var c2 = (BuildBlock(2, 3, QDefective) - a2).Divide(new Complex(QDefective, 0));
        double resCd = (c1 * wd - wd * c2).FrobeniusNorm() / (EpCharacter.SpectralNorm(c2) * wd.FrobeniusNorm());
        _out.WriteLine($"[L,W†] residuals (2,3)→(1,2): A-part {F(resAd)}, C-part {F(resCd)}");
        Assert.True(resAd <= 1e-10 && resCd <= 1e-10, $"[L,W†] not machine zero: A {resAd:E3}, C {resCd:E3}");
    }
}
