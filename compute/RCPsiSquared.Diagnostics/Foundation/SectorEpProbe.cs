using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The core Multi-Sector Monodromy Census primitive: for a joint-popcount sector (p, q̃) at a
/// coupling q0, build the RAW block, find the coalescing eigenvalue pair, and classify it DEFECTIVE
/// (a Jordan block, eigenvectors coalesced, a braid-carrying √-EP) vs DIABOLIC (a semisimple repeated
/// eigenvalue, eigenvectors independent, silent) via the trusted artifact-free
/// <see cref="EpCharacter"/>.
///
/// <para><b>Why it exists.</b> It PINS THE CONVENTION end-to-end: the Task-1 <see cref="SectorBlock"/>
/// materialiser (L(q)/γ = −i·q·[Ĥ,·] + D, Ĥ = Σ_bonds (XX + YY), Δ=0, γ=1, q injected as q·Ĥ) meets the
/// Task-2 EP-character reader (Riesz spectral-projector ‖P‖ / departure-from-normality / geometric-vs-
/// algebraic multiplicity — never an <c>eig</c> eigenvector). The load-bearing validation is that this
/// composition reproduces the exact-math from-below result on the N=4 (1,2) block: q≈0.658983 (a
/// (3q⁴+q²−1) root) is SEMISIMPLE DIABOLIC (SVD nullity 2, the pair merges at λ≈−4+1.318i), q≈0.460212
/// (a P₁₀ real root) is a DEFECTIVE EP (SVD nullity 1).</para>
///
/// <para><b>The contour discipline.</b> <see cref="EpCharacter.Characterize"/> draws a Riesz circle that
/// must ENCLOSE the coalescing pair and EXCLUDE every other eigenvalue of the (possibly large) raw
/// block. The radius is chosen exactly as in <c>EpCharacterWitness.CharacterizeAtQ</c>: a fraction of the
/// distance from the pair midpoint to the THIRD-nearest eigenvalue (0.40·third), floored so it comfortably
/// encloses the split pair (5·gap) and capped so it never reaches the third eigenvalue (0.49·third). At a
/// genuine coalescence the pair's gap → 0, so it is the global-min-gap pair and the two nearest to the
/// midpoint; the third-nearest is the first intruder the radius must stay inside of.</para>
///
/// <para><b>RAW footing.</b> The block is the genuine principal submatrix of the full L on the sector's
/// flat indices, with NO spatial-symmetry reduction. If a third eigenvalue (e.g. from the block's
/// antisymmetric part) sits within any workable radius of the pair, EpCharacter would read Algebraic ≥ 3
/// and the defective/diabolic split would be corrupted; the observed N=4 (1,2) layout keeps the pair
/// cleanly isolated (reported in the census scout). Sibling live lab (the coherence-horizon SE block):
/// <c>EpCharacterWitness</c> (<c>inspect --root epcharacter</c>).</para></summary>
public static class SectorEpProbe
{
    /// <summary>One probe reading: the EP character <paramref name="Kind"/> of the coalescing pair,
    /// its <paramref name="MinGap"/> = |λ_i − λ_j| (→ 0 at a genuine coalescence), the coupling
    /// <paramref name="QValue"/> = q0 probed, the pair midpoint <paramref name="Center"/> the Riesz
    /// contour is drawn about, and the supplementary Riesz <paramref name="ProjectorNorm"/> ‖P‖₂.</summary>
    public readonly record struct ProbeReading(
        EpCharacter.EpKind Kind,
        double MinGap,
        Complex QValue,
        Complex Center,
        double ProjectorNorm);

    /// <summary>The AT-aware "defective-anywhere" reading (<see cref="ProbeDefectiveAnywhere"/>): does ANY
    /// near-coalescing eigenvalue cluster of the (p, q̃) block carry a genuine √-EP (a Jordan block,
    /// <see cref="EpCharacter.EpKind.Defective"/>), even when a PERMANENT semisimple (AT) degeneracy is the
    /// globally-closest pair and hides it from <see cref="Probe"/>. <paramref name="HasDefective"/> is true iff at
    /// least one multi-member cluster reads Defective; <paramref name="MinDefectiveGap"/> is that cluster's intra
    /// gap |λ_i − λ_j| (the smallest over all defective clusters; +∞ if none); <paramref name="DefectiveCenter"/>
    /// is its midpoint λ, the Riesz-contour center (default if none); <paramref name="DefectiveProjectorNorm"/>
    /// the supplementary ‖P‖₂ there; <paramref name="MultiMemberClusterCount"/> the number of near-coalescing
    /// (size ≥ 2) clusters examined.</summary>
    public readonly record struct DefectiveAnywhereReading(
        bool HasDefective,
        double MinDefectiveGap,
        Complex DefectiveCenter,
        double DefectiveProjectorNorm,
        int MultiMemberClusterCount);

    /// <summary>Fraction of the midpoint-to-third-eigenvalue distance used as the nominal contour radius
    /// (the <c>EpCharacterWitness</c> convention: enclose the pair, stay well inside the third).</summary>
    private const double NominalRadiusFraction = 0.40;

    /// <summary>Hard cap on the radius as a fraction of the third-eigenvalue distance: the contour must
    /// NEVER reach the third eigenvalue (which would pull it into the enclosed cluster and read alg ≥ 3).</summary>
    private const double MaxRadiusFraction = 0.49;

    /// <summary>Floor multiple of the pair gap: the contour must comfortably enclose the (split) pair.</summary>
    private const double GapEnclosureFactor = 5.0;

    /// <summary>Fallback third-eigenvalue distance for a block with fewer than three eigenvalues (a closed
    /// tiny block where the contour trivially encloses everything). Large so the radius is unconstrained.</summary>
    private const double LoneClusterFallback = 10.0;

    /// <summary>Single-linkage gap threshold for the <see cref="ProbeDefectiveAnywhere"/> cluster search: two
    /// eigenvalues join a cluster when |λ_i − λ_j| is below this. Chosen (from the N=4 census review) well ABOVE
    /// the permanent semisimple AT pair (gap ~1e-15) AND the residual √-EP split (gap ~5e-4), and well BELOW the
    /// generic eigenvalue spacing (~0.1–0.6), so the permanent AT pair and the residual EP land in SEPARATE
    /// clusters and each is characterized on its own contour.</summary>
    private const double ClusterLinkThreshold = 1e-3;

    /// <summary>Probe the (p, q̃) sector of the N-site XY (Δ=0) + Z-dephasing Liouvillian at coupling
    /// <paramref name="q0"/>: build the RAW block via <see cref="SectorBlock.Build"/>, find the closest
    /// eigenvalue pair (the coalescing pair at a genuine EP/diabolic q0), draw a Riesz contour that
    /// isolates ONLY that pair, and read its DEFECTIVE-vs-DIABOLIC character off
    /// <see cref="EpCharacter.Characterize"/>.</summary>
    /// <param name="N">Number of sites (N ≥ 1).</param>
    /// <param name="p">Column popcount (0 ≤ p ≤ N).</param>
    /// <param name="qTilde">Row popcount (0 ≤ q̃ ≤ N).</param>
    /// <param name="q0">The complex coupling q = J/γ the sector is probed at.</param>
    /// <returns>The <see cref="ProbeReading"/> for the coalescing pair at <paramref name="q0"/>.</returns>
    public static ProbeReading Probe(int N, int p, int qTilde, Complex q0)
    {
        var raw = SectorBlock.Build(N, p, qTilde, q0);
        var m = Matrix<Complex>.Build.DenseOfArray(raw);

        // Eigenvalues only, via the same MathNet dense EVD path used across the block-spectrum layer.
        var spectrum = m.Evd().EigenValues.ToArray();

        // A size-1 sector has no pair to coalesce: the corners (0,0), (0,N), (N,0), (N,N) are C(N,0)² = 1
        // at every N, so a full census sweep hits them. Report "no coalescence here" (MinGap = +∞, Normal,
        // center = the sole eigenvalue) so the sweep skips it, rather than indexing a second eigenvalue
        // that does not exist (or drawing a Riesz contour around a whole 1×1 block).
        if (spectrum.Length < 2)
        {
            Complex sole = spectrum.Length == 1 ? spectrum[0] : default;
            return new ProbeReading(EpCharacter.EpKind.Normal, double.PositiveInfinity, q0, sole, 0.0);
        }

        var (iBest, jBest, minGap) = ClosestPair(spectrum);
        Complex center = 0.5 * (spectrum[iBest] + spectrum[jBest]);
        double radius = IsolatingRadius(spectrum, center, minGap);

        var reading = EpCharacter.Characterize(m, center, radius);
        return new ProbeReading(reading.Kind, minGap, q0, center, reading.ProjectorNorm);
    }

    /// <summary>The AT-aware "defective-anywhere" probe: build the RAW (p, q̃) block via
    /// <see cref="SectorBlock.Build"/>, take its full spectrum, single-linkage-cluster ALL near-coalescences
    /// below <see cref="ClusterLinkThreshold"/>, and run <see cref="EpCharacter.Characterize"/> on EACH
    /// multi-member cluster (each on a Riesz contour that encloses only that cluster). Reports whether ANY
    /// cluster is <see cref="EpCharacter.EpKind.Defective"/> (a genuine √-EP / Jordan block).
    ///
    /// <para><b>Why it exists (the mask <see cref="Probe"/> cannot see through).</b> <see cref="Probe"/> keys on
    /// the GLOBAL closest pair only. At N ≥ 5 many RAW joint-popcount blocks carry a PERMANENT semisimple (AT)
    /// degeneracy (a gap ~1e-15 eigenvalue pair at EVERY q, on the dissipative AT rate lines); that permanent pair
    /// IS the global-closest pair, so <see cref="Probe"/> reads <see cref="EpCharacter.EpKind.Diabolic"/> and the
    /// residual braid-carrying √-EP (gap ~5e-4, at a DIFFERENT λ) is masked. This probe examines every
    /// near-coalescence, so it reaches the residual EP: the permanent semisimple clusters read Diabolic and are
    /// skipped for the verdict; a Defective cluster sets <paramref name="HasDefective"/>. This re-enables the
    /// census's N ≥ 5 cross-fold gate ((1, N−2) shares the (1,2) braid at the conjugate locus, F89d), which the
    /// global-closest probe under-reports as all-Diabolic.</para>
    ///
    /// <para><b>Strict improvement over <see cref="Probe"/>.</b> Where a block has NO permanent pair (the N=4
    /// (1,2) Klein orbit: the residual EP IS the global-closest pair), both agree — this returns
    /// <paramref name="HasDefective"/> = true exactly where <see cref="Probe"/> returns Defective. Where the only
    /// coalescence is a permanent semisimple pair (the N=4 node sectors), both agree it is not a braid
    /// (<paramref name="HasDefective"/> = false). It differs ONLY by ALSO seeing a defective cluster that is not
    /// the global-closest pair: the masked residual √-EP.</para>
    ///
    /// <para>The cluster contour discipline mirrors <see cref="Probe"/> / <c>EpCharacterWitness</c>, generalized
    /// from a pair to a cluster: interpolate the radius <see cref="NominalRadiusFraction"/> of the way from the
    /// cluster's outer edge to the nearest OUTSIDE eigenvalue, floored to comfortably enclose the split cluster
    /// (<see cref="GapEnclosureFactor"/>·edge) and capped at <see cref="MaxRadiusFraction"/> of the
    /// nearest-outsider distance so it never reaches it.</para></summary>
    /// <param name="N">Number of sites (N ≥ 1).</param>
    /// <param name="p">Column popcount (0 ≤ p ≤ N).</param>
    /// <param name="qTilde">Row popcount (0 ≤ q̃ ≤ N).</param>
    /// <param name="q0">The complex coupling q = J/γ the sector is probed at.</param>
    /// <returns>The <see cref="DefectiveAnywhereReading"/> for <paramref name="q0"/>.</returns>
    public static DefectiveAnywhereReading ProbeDefectiveAnywhere(int N, int p, int qTilde, Complex q0)
    {
        var raw = SectorBlock.Build(N, p, qTilde, q0);
        var m = Matrix<Complex>.Build.DenseOfArray(raw);
        var spectrum = m.Evd().EigenValues.ToArray();
        int d = spectrum.Length;

        // Fewer than two eigenvalues (the size-1 corners C(N,0)² = 1 at every N): no pair can coalesce.
        if (d < 2)
        {
            Complex sole = d == 1 ? spectrum[0] : default;
            return new DefectiveAnywhereReading(false, double.PositiveInfinity, sole, 0.0, 0);
        }

        bool anyDefective = false;
        double minDefectiveGap = double.PositiveInfinity;
        Complex defectiveCenter = default;
        double defectiveProjectorNorm = 0.0;
        int multiMemberClusters = 0;

        foreach (var cluster in ClusterByLinkage(spectrum, ClusterLinkThreshold))
        {
            if (cluster.Count < 2) continue;   // a lone eigenvalue cannot be a coalescence
            multiMemberClusters++;

            var members = new HashSet<int>(cluster);
            Complex center = new(cluster.Average(k => spectrum[k].Real), cluster.Average(k => spectrum[k].Imaginary));
            double maxIntra = cluster.Max(k => (spectrum[k] - center).Magnitude);
            double minInter = Enumerable.Range(0, d).Where(k => !members.Contains(k))
                                        .Select(k => (spectrum[k] - center).Magnitude)
                                        .DefaultIfEmpty(LoneClusterFallback).Min();

            // Enclose the whole cluster, stay well inside the nearest outsider (the pair discipline, generalized).
            double radius = maxIntra + NominalRadiusFraction * (minInter - maxIntra);
            radius = Math.Max(radius, GapEnclosureFactor * Math.Max(maxIntra, 1e-13));
            radius = Math.Min(radius, MaxRadiusFraction * minInter);

            var reading = EpCharacter.Characterize(m, center, radius);
            if (reading.Kind == EpCharacter.EpKind.Defective)
            {
                anyDefective = true;
                double intraGap = ClusterMinGap(spectrum, cluster);
                if (intraGap < minDefectiveGap)
                {
                    minDefectiveGap = intraGap;
                    defectiveCenter = center;
                    defectiveProjectorNorm = reading.ProjectorNorm;
                }
            }
        }

        return new DefectiveAnywhereReading(anyDefective, minDefectiveGap, defectiveCenter,
            defectiveProjectorNorm, multiMemberClusters);
    }

    /// <summary>The closest eigenvalue pair (indices and gap): the (i, j) minimising |λ_i − λ_j|. At a
    /// genuine coalescence the coalescing pair has gap → 0, so it IS the global-min-gap pair.</summary>
    private static (int I, int J, double Gap) ClosestPair(Complex[] spectrum)
    {
        int bestI = 0, bestJ = 1;
        double bestGap = double.PositiveInfinity;
        for (int i = 0; i < spectrum.Length; i++)
            for (int j = i + 1; j < spectrum.Length; j++)
            {
                double gap = (spectrum[i] - spectrum[j]).Magnitude;
                if (gap < bestGap) { bestGap = gap; bestI = i; bestJ = j; }
            }
        return (bestI, bestJ, bestGap);
    }

    /// <summary>A radius that encloses ONLY the coalescing pair: 0.40·(distance from <paramref name="center"/>
    /// to the third-nearest eigenvalue), floored at 5·<paramref name="gap"/> (comfortably enclose the split
    /// pair) and capped at 0.49·third (never reach the third eigenvalue). Mirrors the proven contour
    /// discipline in <c>EpCharacterWitness.CharacterizeAtQ</c>.</summary>
    private static double IsolatingRadius(Complex[] spectrum, Complex center, double gap)
    {
        // Distances from the pair midpoint to every eigenvalue, ascending: [0],[1] are the pair, [2] the third.
        var dist = spectrum.Select(e => (e - center).Magnitude).OrderBy(d => d).ToArray();
        double third = dist.Length > 2 ? dist[2] : LoneClusterFallback;

        double r = NominalRadiusFraction * third;
        r = Math.Max(r, GapEnclosureFactor * gap);   // enclose the (split) pair comfortably
        r = Math.Min(r, MaxRadiusFraction * third);  // ...but never reach the third eigenvalue
        return r;
    }

    /// <summary>Single-linkage clustering of the spectrum by the gap <paramref name="threshold"/>: two
    /// eigenvalues share a cluster iff a chain of pairwise gaps all below <paramref name="threshold"/> connects
    /// them (union-find). Returns every cluster (including singletons); the caller keeps the size ≥ 2 ones. This
    /// is what lets a permanent semisimple pair (gap ~1e-15) and a residual √-EP (gap ~5e-4) at a DIFFERENT λ sit
    /// in SEPARATE clusters, so each is characterized on its own contour.</summary>
    private static List<List<int>> ClusterByLinkage(Complex[] spectrum, double threshold)
    {
        int d = spectrum.Length;
        var parent = new int[d];
        for (int i = 0; i < d; i++) parent[i] = i;

        int Find(int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
        void Union(int a, int b) { parent[Find(a)] = Find(b); }

        for (int i = 0; i < d; i++)
            for (int j = i + 1; j < d; j++)
                if ((spectrum[i] - spectrum[j]).Magnitude < threshold) Union(i, j);

        var groups = new Dictionary<int, List<int>>();
        for (int i = 0; i < d; i++)
        {
            int root = Find(i);
            if (!groups.TryGetValue(root, out var list)) { list = new List<int>(); groups[root] = list; }
            list.Add(i);
        }
        return groups.Values.ToList();
    }

    /// <summary>The smallest pairwise gap within a cluster: min over (a, b) ∈ cluster of |λ_a − λ_b| (→ 0 at a
    /// genuine coalescence). Reported as the defective cluster's
    /// <see cref="DefectiveAnywhereReading.MinDefectiveGap"/>.</summary>
    private static double ClusterMinGap(Complex[] spectrum, List<int> cluster)
    {
        double gap = double.PositiveInfinity;
        for (int a = 0; a < cluster.Count; a++)
            for (int b = a + 1; b < cluster.Count; b++)
                gap = Math.Min(gap, (spectrum[cluster[a]] - spectrum[cluster[b]]).Magnitude);
        return gap;
    }
}
