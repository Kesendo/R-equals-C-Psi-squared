using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Multi-Sector Monodromy Census REFERENCE loci: the (1,2) symmetric-residual octic's
/// DEFECTIVE exceptional-point q-values, conjugation-closed, with the silent DIABOLIC node EXCLUDED. These
/// are the q-values at which the census (a later task) probes every OTHER joint-popcount sector's raw block.
///
/// <para><b>What a locus is.</b> As the coupling q = J/γ moves in the complex plane the (SE,DE)=(1,2)
/// residual eigenvalues braid. A DEFECTIVE EP is a genuine √-branch point: two eigenvalues coalesce AND
/// their eigenvectors coalesce (a Jordan block), so a small q-loop about it carries a transposition (it is
/// LOUD, Petermann factor → ∞). These are the simple zeros of the P₁₀ factor of the octic λ-discriminant
/// (disc(F₈) = const·q²⁴·(3q⁴+q²−1)²·P₁₀(q²)), the S₈-generating branch points (F89Path3OcticEpClaim). A
/// DIABOLIC point is the squared-factor degeneracy (3q⁴+q²−1) at q≈0.659: the eigenvalues coalesce but the
/// eigenvectors stay independent (semisimple), so the loop is the identity, silent. Only the DEFECTIVE loci
/// are returned.</para>
///
/// <para><b>Footing (why the octic).</b> The reference is deliberately the octic's (symmetric-residual)
/// defective EPs, the S₈-carrying object, even though the census probes other sectors' RAW blocks at these
/// q-values: the reference q-values are the same either way, and <see cref="PathKMonodromyScout.FindDiabolicsExact"/>
/// (the trusted exact complement-compression scanner, no tracking-flood, EpCharacter classifying each
/// coalescence Diabolic vs Defective) is the correct instrument. <b>DEFECTIVE = !IsSemisimple.</b></para>
///
/// <para><b>N=4 result (validated by <c>ReferenceDefectiveLociTests</c>).</b> Four REAL defective EPs at
/// q≈0.460212, 0.854438, 0.857458, 1.738181, plus complex EPs (0.509819±0.130259i, 0.673890±0.053808i,
/// 2.312084±1.248236i): the complete right-half-plane defective set (10 of the 20 total, the q→−q mirror
/// giving the left half). The 0.854/0.857 pair is a ~0.003-split near-degenerate twin that the 0.05-cell
/// monodromy lasso (<c>gmscan</c>) encloses TOGETHER; the fine <see cref="Cell"/> below RESOLVES it. Sibling
/// live probe of one such point: <see cref="SectorEpProbe"/>.</para></summary>
public static class ReferenceDefectiveLoci
{
    // Scan box for the path-(N-1) (1,2) symmetric-residual defective-EP hunt. re[0.2,3.0]×im[-1.5,1.5] is the
    // documented widened F89 branch-locus window restricted to the RIGHT half-plane: reLo=0.2 masks off both
    // the negative-Re q→−q mirror and the pure-imaginary ±0.876i DIABOLIC nodes (which sit at Re=0); reHi/imHi
    // comfortably enclose the outermost EPs (the real 1.738 and the remote quartet 2.312±1.248i).
    private const double ReLo = 0.2, ReHi = 3.0, ImLo = -1.5, ImHi = 1.5;

    // Grid cell for FindDiabolicsExact. Fine enough that its gap-field seeder AND its 2·cell dedup keep the two
    // near-degenerate real twins q≈0.854438 / 0.857458 apart (separation ≈0.00302 > 2·cell = 0.002): at
    // cell≥0.0015 they merge into a single detection (the lasso's blind spot). Validated at N=4; a higher N
    // (denser F_d residual) may need a finer/retuned scan, so this is a from-N=4 default, not a proof.
    private const double Cell = 0.001;

    // Merge tolerance in q for dedup, conjugate-closure, and diabolic exclusion. Well below the smallest true
    // separation (the 0.003 twin split) and well above the EVD/refine noise on a single locus (~1e-8).
    private const double Tol = 1e-4;

    private static readonly ConcurrentDictionary<int, IReadOnlyList<Complex>> Cache = new();

    /// <summary>The (1,2) octic's DEFECTIVE (P₁₀) EP q-values at chain length <paramref name="N"/>
    /// (path-k = N−1), conjugation-closed, with the silent diabolic excluded. Memoized per N (the scan is
    /// deterministic). N=4 returns the 10 right-half-plane defective EPs (4 real + 3 conjugate pairs).</summary>
    /// <param name="N">Chain length; the (1,2) octic is the path-(N−1) symmetric residual (N ≥ 2).</param>
    public static IReadOnlyList<Complex> For(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 2 (path-k = N-1 >= 1).");
        return Cache.GetOrAdd(N, Compute);
    }

    private static IReadOnlyList<Complex> Compute(int N)
    {
        // The trusted exact-residual coalescence scanner: it already classifies each coalescence
        // Diabolic (IsSemisimple) vs Defective (!IsSemisimple) via EpCharacter. No re-derivation, no hand
        // re-classification.
        var pts = PathKMonodromyScout.FindDiabolicsExact(
            k: N - 1, reLo: ReLo, reHi: ReHi, imLo: ImLo, imHi: ImHi, cell: Cell);

        // The DIABOLIC nodes found in the SAME scan (e.g. q≈0.659): the exclusion set. Closed under conjugation
        // so a defective locus is dropped whether it lands near a diabolic or near its mirror.
        var diabolics = pts.Where(d => d.IsSemisimple)
                           .SelectMany(d => new[] { d.QValue, Complex.Conjugate(d.QValue) })
                           .ToList();

        var defective = pts.Where(d => !d.IsSemisimple).Select(d => d.QValue).ToList();

        bool NearDiabolic(Complex q) => diabolics.Any(z => (z - q).Magnitude < Tol);

        // Conjugation-close (add each mirror), then dedup within Tol, dropping anything near a diabolic. A real
        // locus and its own (≈identical) conjugate collapse to one representative that still answers the
        // conjugation-closed test (it lies within Tol of its own conjugate); a complex locus keeps both signs.
        var closed = new List<Complex>();
        foreach (var q in defective.Concat(defective.Select(Complex.Conjugate)))
        {
            if (NearDiabolic(q)) continue;
            if (closed.Any(r => (r - q).Magnitude < Tol)) continue;
            closed.Add(q);
        }
        return closed;
    }
}
