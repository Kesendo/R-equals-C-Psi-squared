using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

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

        var (iBest, jBest, minGap) = ClosestPair(spectrum);
        Complex center = 0.5 * (spectrum[iBest] + spectrum[jBest]);
        double radius = IsolatingRadius(spectrum, center, minGap);

        var reading = EpCharacter.Characterize(m, center, radius);
        return new ProbeReading(reading.Kind, minGap, q0, center, reading.ProjectorNorm);
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
}
