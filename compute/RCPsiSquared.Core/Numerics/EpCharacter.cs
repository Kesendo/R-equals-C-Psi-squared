using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The artifact-free EP-character diagnostic: tells a genuine DEFECTIVE exceptional point
/// (a Jordan block, two eigenvectors coalesced) from a DIABOLIC degeneracy (a normal repeated
/// eigenvalue, eigenvectors independent) from a NEAR-EP (a simple isolated eigenvalue with a
/// non-orthogonal neighbour). It is the non-<c>eig</c>-eigenvector SIBLING of
/// <see cref="PhaseRigidity"/>.
///
/// <para><b>Why it exists.</b> <see cref="PhaseRigidity"/> reads the Petermann factor r = 1/√K off a
/// single <c>Evd(L)</c> — it reads a raw eigenvector pairing, and that family MISFIRED in the F86a
/// retraction (it reported r → 0 / K huge on a merely near-degenerate, NON-defective spectrum,
/// because <c>eig</c> returns an ill-conditioned, near-parallel eigenvector basis on a near-degenerate
/// cluster). None of the three measures here reads an <c>eig</c> eigenvector:</para>
///
/// <list type="number">
/// <item><b>Riesz spectral-projector</b> P = (1/2πi)∮(zI−L)⁻¹ dz on a small circle enclosing only the
/// target cluster, by M-point midpoint quadrature of the resolvent. ‖P‖₂ ≫ 1 ⟺ the enclosed eigenspace
/// is non-orthogonal to the rest (non-normality leaking in from OUTSIDE the contour). NOTE: on a CLOSED
/// small block (the contour encloses the WHOLE matrix, e.g. a bare 2×2) ‖P‖₂ ≈ 1 trivially and CANNOT
/// gate — so it is supplementary, never load-bearing. The discriminators are 2 and 3.</item>
/// <item><b>Departure-from-normality</b> of the compression A = VᴴLV (V = an orthonormal basis of
/// range(P), from the leading left singular vectors of P): dep = √(max(0, ‖A‖_F² − Σ|λ(A)|²)). For a
/// normal A, ‖A‖_F² = Σ|λ|² exactly, so dep = 0; a Jordan block leaves the coupling in the off-diagonal,
/// so dep > 0. Basis-independent (a unitary change of V leaves both terms invariant).</item>
/// <item><b>Geometric vs algebraic multiplicity.</b> alg = #eigenvalues enclosed = round(Re tr P);
/// geo = nullity of (A − λ̄·I) at the compression's mean eigenvalue λ̄, by SVD with an ABSOLUTE floor
/// (the F86a/reference lesson: a near-zero (A − λ̄I) must not be read full-rank because a relative
/// tolerance scales to 0). geo &lt; alg ⟹ DEFECTIVE (a Jordan block); geo = alg ⟹ DIABOLIC.</item>
/// </list>
///
/// <para><b>Verdict.</b> alg ≥ 2 with geo &lt; alg and dep above tolerance ⟹ <see cref="EpKind.Defective"/>;
/// alg ≥ 2 with geo = alg and dep ≈ 0 ⟹ <see cref="EpKind.Diabolic"/>; exactly one eigenvalue enclosed
/// with a large ‖P‖ ⟹ <see cref="EpKind.NearEp"/> (an off-axis EP felt as non-orthogonality);
/// otherwise <see cref="EpKind.Normal"/>.</para>
///
/// <para>Port of the artifact-free machinery in <c>simulations/review_f86a_diabolic_vs_defective.py</c>
/// and <c>simulations/review_coherence_horizon_ep.py</c>; the live lab is
/// <c>EpCharacterWitness</c> (<c>inspect --root epcharacter</c>).</para></summary>
public static class EpCharacter
{
    /// <summary>The four character verdicts.</summary>
    public enum EpKind
    {
        /// <summary>A genuine 2nd-order (or higher) exceptional point: a Jordan block, geo &lt; alg,
        /// departure-from-normality bounded away from 0. The eigenvectors have coalesced.</summary>
        Defective,

        /// <summary>A diabolic point: a normal repeated eigenvalue, geo = alg, departure ≈ 0. The
        /// eigenvectors stay independent (the eigenvalues merely touch).</summary>
        Diabolic,

        /// <summary>A simple isolated eigenvalue (alg = 1) whose enclosing projector norm is large:
        /// an off-axis EP felt as non-orthogonality with neighbours, but not coalesced on this
        /// contour.</summary>
        NearEp,

        /// <summary>No coalescence and no non-normality: a normal, well-separated eigenvalue.</summary>
        Normal,
    }

    /// <summary>The three artifact-free measures plus the verdict. <paramref name="ProjectorNorm"/> is the
    /// supplementary ‖P‖₂; <paramref name="Algebraic"/> = round(Re tr P); <paramref name="Geometric"/> =
    /// SVD nullity of (A − λ̄I); <paramref name="Departure"/> = departure-from-normality of A;
    /// <paramref name="CompressionEigenvalues"/> are the eigenvalues of the compression A (its dimension
    /// is <paramref name="Algebraic"/>); <paramref name="EigenvectorMergeCos"/> is |cos∠| of the two
    /// compression eigenvectors when alg = 2 (→ 1 as a Jordan pair coalesces), else NaN.</summary>
    public readonly record struct Reading(
        EpKind Kind,
        double ProjectorNorm,
        double AlgebraicRaw,
        int Algebraic,
        int Geometric,
        double Departure,
        double CompressionNorm,
        Complex[] CompressionEigenvalues,
        double EigenvectorMergeCos);

    /// <summary>Quadrature points on the Riesz contour. 600 midpoints match the Python reference; the
    /// resolvent is smooth on a circle that cleanly separates the cluster, so this is plenty.</summary>
    public const int DefaultQuadraturePoints = 600;

    private const double DepartureRelTolerance = 1e-2;   // dep/‖A‖ below this ⟹ "normal" compression
    private const double NearEpProjectorNorm = 10.0;     // ‖P‖₂ above this on a simple eigenvalue ⟹ near-EP

    /// <summary>The Riesz spectral projector P = (1/2πi)∮(zI − L)⁻¹ dz on the circle of radius
    /// <paramref name="radius"/> about <paramref name="center"/>, by <paramref name="quadPoints"/>-point
    /// midpoint quadrature. A pure resolvent contour integral — never an <c>eig</c> eigenvector. P is the
    /// spectral projector onto the eigenvalues strictly inside the contour; tr P = Σ (algebraic
    /// multiplicities enclosed), and ‖P‖₂ ≥ 1 with equality iff the enclosed eigenspace is orthogonal to
    /// its complement.</summary>
    public static ComplexMatrix RieszProjector(ComplexMatrix l, Complex center, double radius,
        int quadPoints = DefaultQuadraturePoints)
    {
        int d = l.RowCount;
        var id = Matrix<Complex>.Build.DenseIdentity(d);
        var p = Matrix<Complex>.Build.Dense(d, d);
        double twoPi = 2.0 * Math.PI;
        for (int k = 0; k < quadPoints; k++)
        {
            double th = twoPi * (k + 0.5) / quadPoints;
            Complex z = center + radius * Complex.FromPolarCoordinates(1.0, th);
            // dz = i·r·e^{iθ}·(2π/M); the integrand is (zI − L)⁻¹.
            Complex dz = Complex.ImaginaryOne * radius * Complex.FromPolarCoordinates(1.0, th) * (twoPi / quadPoints);
            var resolvent = (z * id - l).Inverse();
            p += resolvent.Multiply(dz);
        }
        return p.Divide(new Complex(0, twoPi));   // /(2πi)
    }

    /// <summary>An orthonormal basis of range(P): the first <paramref name="rank"/> left singular vectors
    /// of P (its dominant <paramref name="rank"/>-dimensional column space). Returned as a d×rank matrix
    /// V with VᴴV = I.</summary>
    private static ComplexMatrix RangeBasis(ComplexMatrix p, int rank)
    {
        var svd = p.Svd();
        return svd.U.SubMatrix(0, p.RowCount, 0, rank);
    }

    /// <summary>Departure-from-normality of A: √(max(0, ‖A‖_F² − Σ|λ(A)|²)). Zero for a normal matrix,
    /// positive for a defective one (the Frobenius norm counts the Jordan off-diagonal that the
    /// eigenvalues miss). Schur-free: uses only the Frobenius norm and the eigenvalue magnitudes, both
    /// unitarily invariant.</summary>
    public static double DepartureFromNormality(ComplexMatrix a)
    {
        double froSq = a.FrobeniusNorm();
        froSq *= froSq;
        double eigSumSq = a.Evd().EigenValues.Sum(e => e.Magnitude * e.Magnitude);
        return Math.Sqrt(Math.Max(0.0, froSq - eigSumSq));
    }

    /// <summary>The SVD nullity of M = A − shift·I with an absolute-and-relative threshold: #singular
    /// values below max(relTol·σ_max, absTol·max(‖A‖₂, 1)). The absolute floor is the F86a lesson — a
    /// matrix that is exactly or nearly the zero matrix (a diabolic A − λ̄I where every eigenvalue equals
    /// λ̄) must read as fully null, but a purely relative tolerance scales the threshold to 0 and reads it
    /// full-rank. Returned nullity is clamped to ≥ 1 (a shift at an eigenvalue has at least one null
    /// direction).</summary>
    public static int Nullity(ComplexMatrix a, Complex shift, double relTol = 1e-6, double absTol = 1e-7)
    {
        var m = a - shift * Matrix<Complex>.Build.DenseIdentity(a.RowCount);
        // For a complex matrix MathNet returns the singular values as a Vector<Complex> (non-negative
        // reals stored in the real part); take their real magnitudes.
        var s = m.Svd().S.Select(sv => sv.Real).ToArray();
        double aNorm = Math.Max(SpectralNorm(a), 1.0);
        double sMax = s.Length > 0 ? s.Max() : 0.0;
        double thr = Math.Max(relTol * sMax, absTol * aNorm);
        int nullity = s.Count(sv => sv < thr);
        return Math.Max(nullity, 1);
    }

    /// <summary>The spectral norm ‖A‖₂ = the largest singular value of A.</summary>
    public static double SpectralNorm(ComplexMatrix a) => a.Svd().S.Select(sv => sv.Real).Max();

    /// <summary>|cos∠| between the two eigenvectors of a 2×2 (normalised). → 1 as a Jordan pair
    /// coalesces (the corroborating eigenvector-merge measure); bounded away from 1 for a diabolic pair.
    /// Returns NaN if A is not 2×2.</summary>
    public static double EigenvectorMergeCosine(ComplexMatrix a)
    {
        if (a.RowCount != 2) return double.NaN;
        var evd = a.Evd();
        var v0 = evd.EigenVectors.Column(0);
        var v1 = evd.EigenVectors.Column(1);
        v0 = v0.Divide(v0.L2Norm());
        v1 = v1.Divide(v1.L2Norm());
        // |⟨v0, v1⟩| with the Hermitian inner product.
        Complex ip = v0.ConjugateDotProduct(v1);
        return ip.Magnitude;
    }

    /// <summary>Characterise the eigenvalue cluster of L enclosed by a circle of radius
    /// <paramref name="radius"/> about <paramref name="center"/>, via the three artifact-free measures.
    /// <paramref name="center"/> should be the cluster midpoint and <paramref name="radius"/> chosen to
    /// enclose ONLY that cluster (outside any split, inside the next eigenvalue).
    ///
    /// <para>geo is read from the compression A at its OWN mean eigenvalue λ̄ = mean λ(A), not from
    /// (L − center·I): when a Jordan pair is split to ~√(machine-ε), the contour midpoint is not an exact
    /// eigenvalue of L, so the full-matrix nullity at <paramref name="center"/> would read 0. The
    /// compression localises the question to the enclosed 2×2/3×3, where λ̄ is the genuine repeated
    /// value.</para></summary>
    public static Reading Characterize(ComplexMatrix l, Complex center, double radius,
        int quadPoints = DefaultQuadraturePoints)
    {
        var p = RieszProjector(l, center, radius, quadPoints);
        double projNorm = SpectralNorm(p);                // ‖P‖₂ = largest singular value
        double algRaw = p.Trace().Real;
        int alg = Math.Max((int)Math.Round(algRaw), 1);

        var v = RangeBasis(p, alg);
        var a = v.ConjugateTranspose() * l * v;           // the compression onto range(P)
        var eigA = a.Evd().EigenValues.ToArray();
        double dep = DepartureFromNormality(a);
        double aNorm = a.FrobeniusNorm();

        Complex lamBar = eigA.Length > 0
            ? new Complex(eigA.Average(e => e.Real), eigA.Average(e => e.Imaginary))
            : center;
        int geo = Nullity(a, lamBar);
        double mergeCos = EigenvectorMergeCosine(a);

        double relDep = dep / Math.Max(1.0, aNorm);
        EpKind kind = Classify(alg, geo, relDep, projNorm);

        return new Reading(kind, projNorm, algRaw, alg, geo, dep, aNorm, eigA, mergeCos);
    }

    /// <summary>The verdict from the three measures. Mirrors the Python <c>classify</c>: the primary
    /// discriminators are geo-vs-alg and the RELATIVE departure (both basis-independent, neither an
    /// <c>eig</c> pairing); ‖P‖ is supplementary (informative only when eigenvalues OUTSIDE the contour
    /// have non-orthogonal eigenspaces — trivially 1 on a closed block — so it gates only the alg = 1
    /// near-EP case, never the defective/diabolic split).</summary>
    private static EpKind Classify(int alg, int geo, double relDep, double projNorm)
    {
        if (alg >= 2 && geo < alg && relDep > DepartureRelTolerance) return EpKind.Defective;
        if (alg >= 2 && geo == alg && relDep < DepartureRelTolerance) return EpKind.Diabolic;
        if (alg < 2 && projNorm > NearEpProjectorNorm) return EpKind.NearEp;
        // alg >= 2 but the two signals disagree (e.g. dep large yet geo = alg): non-normal but not a
        // clean Jordan on this contour — report NearEp (a felt non-orthogonality) rather than Normal.
        if (alg >= 2 && relDep > DepartureRelTolerance) return EpKind.NearEp;
        return EpKind.Normal;
    }
}
