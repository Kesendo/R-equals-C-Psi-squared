using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F121, the qudit partial palindrome (2026-06-11). Where the typed
/// claims <see cref="QuditPartialPalindromeCeiling"/> and <see cref="QuditProductMirrorCap"/>
/// carry the closed forms as integer arithmetic, this witness BUILDS the actual full-Cartan
/// dephasing dissipator at inspect time and reads its spectrum: it materialises the
/// d^{2N}×d^{2N} Liouvillian L_D in the computational coherence basis (diagonal, rate of
/// |i⟩⟨j| = −2γ·Hamming(i, j) per the Absorption Theorem one dimension up), counts the modes
/// that actually pair under the palindrome reflection λ ↦ −2(Nγ) − λ about the physical
/// center −Nγ, and asserts that the live count equals the closed-form ceiling
/// Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)}. It then reads the product-mirror cap (2d)^N and the
/// non-product remainder off the same spectrum.
///
/// <para>The dissipator is diagonal in this basis, so the spectrum is the matrix diagonal
/// (eigenvalues of a diagonal matrix), no O(n³) eigendecomposition needed; the matrix is built
/// in full so the witness is a genuine compute, not a formula lookup. Guard:
/// d^{2N} ≤ <see cref="MaxDim"/> (= 1024), admitting (d, N) ∈ {(2,2):16, (2,3):64, (3,2):81,
/// (4,2):256, (3,3):729} and excluding (4,3):4096 and beyond.</para>
///
/// <para>Children: one node per Hamming/disagreement rung k with (multiplicity c_k, tilt
/// (d−1)^{min(k,N−k)}, mirror rung N−k, paired count), plus a cap node splitting the live
/// paired count into the product-attained (2d)^N and the non-product remainder. Summary, e.g.
/// "d=3 N=2: paired 54/81 (ceiling met), cap 36, non-product 18; full iff d=2".</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> +
/// <c>simulations/qutrit_partial_palindrome.py</c> +
/// <c>simulations/qudit_product_mirror_cap.py</c>.</para></summary>
public sealed class QuditPartialPalindromeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest Liouville dimension d^{2N} the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int D { get; }
    public int N { get; }
    public double Gamma { get; }

    /// <summary>Liouville dimension d^{2N}.</summary>
    public int Dim { get; }

    /// <summary>The physical palindrome center −Nγ. Pairing reflects λ ↦ −2(Nγ) − λ.</summary>
    public double Center => -N * Gamma;

    private readonly double[] _spectrum;     // the live diagonal of L_D, length Dim
    private readonly int[] _hamming;         // Hamming(i, j) per coherence index, length Dim

    public QuditPartialPalindromeWitness(int d = 3, int n = 2, double gamma = 0.05)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");

        long dim = QuditPartialPalindromeCeiling.IntPow(d, 2 * n);
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(d),
                $"d^(2N) = {dim} exceeds the live-build guard {MaxDim} for (d={d}, N={n}); " +
                "pick a smaller (d, N), e.g. (3,2)=81, (4,2)=256, (3,3)=729.");

        D = d;
        N = n;
        Gamma = gamma;
        Dim = (int)dim;

        var L = BuildDissipator(d, n, gamma, out _hamming);
        _spectrum = new double[Dim];
        for (int x = 0; x < Dim; x++) _spectrum[x] = L[x, x].Real;   // diagonal = spectrum
    }

    /// <summary>Builds the full-Cartan dephasing dissipator L_D as a d^{2N}×d^{2N} matrix in
    /// the computational coherence basis. The basis index of |i⟩⟨j| is i·d^N + j (site 0 the
    /// most significant base-d digit, matching the typed claims). Under full-Cartan dephasing
    /// the d levels are equidistant, so L_D is diagonal with entry −2γ·Hamming(i, j); we build
    /// it as a dense matrix (not just the diagonal) so the witness is a real compute. Also
    /// returns the per-index Hamming distance.</summary>
    public static ComplexMatrix BuildDissipator(int d, int n, double gamma, out int[] hamming)
    {
        int dN = (int)QuditPartialPalindromeCeiling.IntPow(d, n);
        int dim = dN * dN;
        hamming = new int[dim];
        var L = Matrix<Complex>.Build.Dense(dim, dim);
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
            {
                int idx = i * dN + j;
                int h = HammingDistance(i, j, d, n);
                hamming[idx] = h;
                L[idx, idx] = new Complex(-2.0 * gamma * h, 0.0);
            }
        return L;
    }

    /// <summary>Per-digit disagreement count between the base-d representations of i and j.</summary>
    private static int HammingDistance(int i, int j, int d, int n)
    {
        int a = i, b = j, h = 0;
        for (int l = 0; l < n; l++)
        {
            if (a % d != b % d) h++;
            a /= d;
            b /= d;
        }
        return h;
    }

    /// <summary>The live paired count: greedily matches each eigenvalue with a partner under
    /// λ ↦ −2(Nγ) − λ (the reflection about the physical center −Nγ), counting how many modes
    /// find a mirror within <paramref name="tol"/>. A self-mirroring mode (λ on the center)
    /// counts 1; a matched pair counts 2. This is the spectral analogue of the Python
    /// <c>palindrome_pairs</c>.</summary>
    public int PairedCount(double tol = 1e-9)
    {
        double twoCenter = 2.0 * Center;   // reflection target = twoCenter − λ
        var used = new bool[Dim];
        int paired = 0;
        for (int k = 0; k < Dim; k++)
        {
            if (used[k]) continue;
            double target = twoCenter - _spectrum[k];
            int best = -1;
            double bestDiff = double.PositiveInfinity;
            for (int u = 0; u < Dim; u++)
            {
                if (used[u]) continue;
                double diff = Math.Abs(_spectrum[u] - target);
                if (diff < bestDiff) { bestDiff = diff; best = u; }
            }
            if (best >= 0 && bestDiff < tol)
            {
                paired += (k == best) ? 1 : 2;
                used[k] = true;
                used[best] = true;
            }
            else
            {
                used[k] = true;   // k has no partner; leave it unpaired
            }
        }
        return paired;
    }

    /// <summary>Live multiplicity per Hamming rung, read off the built spectrum (not the
    /// formula): how many coherence indices sit at each Hamming distance k.</summary>
    public long[] LiveRungCounts()
    {
        var c = new long[N + 1];
        foreach (int h in _hamming) c[h]++;
        return c;
    }

    public string DisplayName =>
        $"QuditPartialPalindromeWitness (F121 live lab, d={D}, N={N}, γ={Gamma.ToString("0.###", Inv)}, dim={Dim})";

    public string Summary
    {
        get
        {
            int paired = PairedCount();
            long ceil = QuditPartialPalindromeCeiling.Ceiling(D, N);
            long cap = QuditProductMirrorCap.ProductCap(D, N);
            long nonProduct = QuditProductMirrorCap.NonProductPart(D, N);
            bool met = paired == ceil;
            bool full = paired == Dim;
            return $"d={D} N={N}: paired {paired}/{Dim} " +
                   $"({(met ? "ceiling met" : $"ceiling {ceil} NOT met")}), cap {cap}, " +
                   $"non-product {nonProduct}; full {(full ? "yes" : "no")} (full iff d=2). " +
                   $"Live spectrum built ({Dim}×{Dim} L_D), counted about center −Nγ={Center.ToString("0.###", Inv)}.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            int paired = PairedCount();
            long ceil = QuditPartialPalindromeCeiling.Ceiling(D, N);
            long cap = QuditProductMirrorCap.ProductCap(D, N);
            long nonProduct = QuditProductMirrorCap.NonProductPart(D, N);
            var live = LiveRungCounts();

            // 1. The live spectrum vs the closed-form ceiling (the assert made visible).
            yield return new InspectableNode(
                displayName: "the live count vs the closed form",
                summary: $"built the {Dim}×{Dim} dissipator L_D, counted modes pairing under λ ↦ −2(Nγ)−λ " +
                         $"about the center −Nγ={Center.ToString("0.###", Inv)}: paired = {paired}; " +
                         $"closed-form ceiling Σ_k d^N·C(N,k)·(d−1)^min(k,N−k) = {ceil}. " +
                         (paired == ceil ? "MATCH (live = closed form)." : "MISMATCH."));

            // 2. One node per Hamming/disagreement rung k.
            for (int k = 0; k <= N; k++)
            {
                int m = N - k;
                long ck = live[k];
                long formCk = QuditPartialPalindromeCeiling.Multiplicity(D, N, k);
                long tilt = QuditPartialPalindromeCeiling.IntPow(D - 1, Math.Min(k, m));
                long rungPaired = (k < m)
                    ? 2L * Math.Min(live[k], live[m])
                    : (k == m ? live[k] : 0);   // k > m already counted from the mirror rung
                string note = k == m ? "self-mirrors (middle rung)"
                                     : (k < m ? $"pairs into rung {m}" : $"counted from rung {m}");
                yield return new InspectableNode(
                    displayName: $"rung k={k} (rate −2γ·{k}={(-2.0 * Gamma * k).ToString("0.###", Inv)})",
                    summary: $"multiplicity c_{k} = {ck} (live) = {formCk} (formula d^N·C(N,k)·(d−1)^{k}); " +
                             $"mirror rung N−k = {m}; tilt (d−1)^min(k,N−k) = {tilt}; " +
                             $"paired here = {rungPaired}; {note}.");
            }

            // 3. The cap node: product-attained (2d)^N vs non-product remainder of the live count.
            yield return new InspectableNode(
                displayName: "the cap split (product vs non-product)",
                summary: $"the live paired {paired} splits as product-mirror cap (2d)^N = {cap} " +
                         $"(reached by Π_d(ρ)=ρᵀ·Shift^⊗N) + non-product remainder ceiling−(2d)^N = {nonProduct} " +
                         $"(needs a global non-product isometry / translation-invariant mirror). " +
                         (D == 2 ? "At d=2 the cap is full, non-product is 0: the qubit magic." : "Non-product > 0: d>2."));

            // 4. The d=2 anchor: full iff d=2 (the d²−2d=0 column).
            yield return new InspectableNode(
                displayName: "full iff d=2 (the d²−2d=0 column)",
                summary: $"the live paired count equals the total dim {Dim} ⟺ d=2 " +
                         $"(here paired {paired}, total {Dim}: {(paired == Dim ? "FULL" : "partial")}); " +
                         "d=2 is the unique fully-paired column, the (d−1)^k tilt base = 1 only there.");

            // 5. The spectrum as a payload curve (sorted real parts), for --draw.
            var sorted = (double[])_spectrum.Clone();
            Array.Sort(sorted);
            var idxAxis = new double[Dim];
            for (int x = 0; x < Dim; x++) idxAxis[x] = x;
            yield return new InspectableNode(
                displayName: "the live spectrum (sorted Re λ)",
                summary: $"the {Dim} eigenvalues of L_D, all real (rungs at 0, −2γ, …, −2Nγ); " +
                         "palindromic overlap about −Nγ is the paired count.",
                payload: new InspectablePayload.Curve("sorted Re λ", idxAxis, sorted, "index", "Re λ"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
