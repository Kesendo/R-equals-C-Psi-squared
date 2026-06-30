using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One N's reading of the reflection-sector conjugacy structure of the full (SE,DE) block:
/// the dimensions of the R-even and R-odd sectors, the count of reflection-fixed singletons, and how
/// close each sector's spectrum is to being self-conjugate vs the two sectors being conjugate partners.</summary>
public sealed record ReflectionSectorReading(
    int N, int FullDim, int EvenDim, int OddDim, int FixedSingletons,
    double SelfConjEven, double SelfConjOdd, double CrossConj);

/// <summary>One N's reading of the REAL R-even residual population over the physical real-q window: the residual
/// degree F_d, how many of the AT-stripped residual strands are real at a generic q and at most over the window,
/// and the closest two real strands ever get (the real-real approach). The within-odd grounding: even N carries
/// NO real residual strands at all (MaxRealCount = 0, the realness is entirely cross-sector), so it cannot host
/// a real-q diabolic in R-even; odd N does (population growing with N), and the real-q diabolics live among these
/// real coalescences. A small GlobalApproach (≈ 0) flags a real-real CROSSING (the geometry where two real
/// strands meet, e.g. N=7); odd N can also host a real-q diabolic as a conjugate pair tangent to the axis (e.g.
/// N=9), which leaves GlobalApproach bounded — so GlobalApproach detects one geometry, not the diabolic count.</summary>
public sealed record RealResidualReading(
    int N, int ResidualDegree, int RealCountGenericQ, int MaxRealCount, double GlobalApproach);

/// <summary>The from-below grounding of the odd-N real-q diabolic onset: the dimension-mismatch / sector-swap.
///
/// <para>The site reflection R: i ↦ nBlock−1−i commutes with the full (SE,DE) block
/// (<see cref="F89PathKSeDeBlock.BuildFullBlock"/> + <see cref="F89PathKSeDeBlock.ReflectionPermutation"/>)
/// and splits it into an R-even and an R-odd sector. This witness builds the block at inspect time, splits
/// it, eigendecomposes each sector, and measures whether the realness antiunitarity (the chiral Σ L Σ = L†,
/// which forces the spectrum real-or-conjugate-paired) acts WITHIN each sector (each self-conjugate) or
/// ACROSS them (the two sectors conjugate partners):</para>
///
/// <para>At EVEN nBlock the two sectors have equal dimension and σ(R-even) = conj(σ(R-odd)) exactly
/// (cross-conjugacy defect ~10⁻¹³, each sector's self-conjugacy defect O(20)): the realness lives across
/// the sectors, so a real-axis collision inside R-even (where the diabolic scout works) is the generic
/// pseudo-Hermitian DEFECTIVE EP. At ODD nBlock the reflection-fixed central site makes
/// dim(R-even) − dim(R-odd) = (nBlock−1)/2 ≠ 0; the dimension mismatch forbids the cross-pairing, forcing
/// the antiunitarity to act within each sector (each self-conjugate, defect ~10⁻¹⁴), so R-even carries
/// genuine real eigenvalues, two of which can cross at real q as a SEMISIMPLE (diabolic) coalescence.
/// The witness also reproduces the C# scout's N=7 (λ=−4.942) and N=9 (λ=−5.424) diabolics in the R-even
/// sector, from the full block.</para>
///
/// <para>Anchor: <c>experiments/F89_PATH_K_DIABOLIC.md</c> (the "mechanism, grounded from below" section);
/// the precedent is <c>experiments/SLOW_MODE_R_PARITY.md</c> (the reflection-fixed odd-N JW band-centre
/// zero mode). The persistent evidence for the <c>diabolic_over_higher_n</c> arc's mechanism.</para></summary>
public sealed class DiabolicReflectionParityWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static readonly int[] SweepN = { 5, 6, 7, 8, 9 };
    private const double GenericQ = 1.0;          // a generic real q, away from any diabolic

    /// <summary>Build the full (SE,DE) block at real coupling q, split by R, eigendecompose each sector, and
    /// measure the self- vs cross-sector conjugacy structure.</summary>
    public ReflectionSectorReading Read(int nBlock, double q)
    {
        var L = F89PathKSeDeBlock.BuildFullBlock(nBlock, new Complex(q, 0));
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        var (even, odd, fixedCount) = ParityColumns(perm);
        var se = SectorEigenvalues(L, even);
        var so = SectorEigenvalues(L, odd);
        double selfE = MatchDist(se, se.Select(Complex.Conjugate).ToArray());
        double selfO = MatchDist(so, so.Select(Complex.Conjugate).ToArray());
        double cross = even.Count == odd.Count
            ? MatchDist(se.Select(Complex.Conjugate).ToArray(), so)     // is conj(σ_even) == σ_odd?
            : double.NaN;                                              // dims differ ⟹ no cross-pairing
        return new ReflectionSectorReading(
            nBlock, L.GetLength(0), even.Count, odd.Count, fixedCount, selfE, selfO, cross);
    }

    /// <summary>The two R-even eigenvalues nearest a target real λ at real coupling qRe, and their gap.
    /// At the scout's real-q diabolic the gap collapses (the coalescence), reproduced from the full block.</summary>
    public (double Lambda1, double Lambda2, double Gap) ReproduceDiabolic(int nBlock, double qRe, double target)
    {
        var L = F89PathKSeDeBlock.BuildFullBlock(nBlock, new Complex(qRe, 0));
        var perm = F89PathKSeDeBlock.ReflectionPermutation(nBlock);
        var (even, _, _) = ParityColumns(perm);
        var se = SectorEigenvalues(L, even);
        var nearest = se.OrderBy(z => Math.Abs(z.Real - target)).Take(2).ToArray();
        return (nearest[0].Real, nearest[1].Real, (nearest[0] - nearest[1]).Magnitude);
    }

    private const double QLo = 0.2, QHi = 3.0, QStep = 0.01;   // the physical real-q window (matched to the diabolic scan)
    private const double RealTol = 1e-6;                       // |Im λ| < RealTol ⟹ the residual root is real

    /// <summary>The sorted real parts of the R-even residual roots that are REAL (|Im λ| &lt; <see cref="RealTol"/>)
    /// at real coupling q, read off the exact AT-complement compression
    /// (<see cref="PathKMonodromyScout.ResidualRootsExact"/>, the same AT-free F_d strands the diabolic scout
    /// uses; <c>ResidualRootsExact</c> compresses onto the complement of the q-independent AT invariant subspace,
    /// and that subspace lives in the S₂-symmetric block = the R-even sector). At odd N the R-even sector is
    /// self-conjugate, so it carries genuine real eigenvalues; at even N it is conjugate-paired with R-odd, so
    /// this is empty (a real residual root would need a real partner in R-odd, which the cross-pairing forbids).</summary>
    private static double[] RealResidualParts(int k, double q)
    {
        var roots = PathKMonodromyScout.ResidualRootsExact(k, new Complex(q, 0));
        var reals = new List<double>();
        foreach (var z in roots) if (Math.Abs(z.Imaginary) < RealTol) reals.Add(z.Real);
        reals.Sort();
        return reals.ToArray();
    }

    /// <summary>The within-odd grounding: scan the physical real-q window and measure the REAL R-even residual
    /// population. The parity mechanism makes R-even self-conjugate (so it carries real eigenvalues) only at ODD
    /// N; at even N R-even is conjugate-paired with R-odd, so the AT-stripped residual carries NO real eigenvalues
    /// anywhere (MaxRealCount = 0) and cannot host a real-q diabolic. Odd N carries a real residual population
    /// that grows with N, and the scout's real-q diabolics are real-λ coalescences among these. Also reports the
    /// closest two real strands ever get (GlobalApproach): ≈ 0 flags a real-real CROSSING (e.g. N=7 at q≈2.628);
    /// bounded away does NOT mean no diabolic, since an odd-N real-q diabolic can also be a conjugate pair tangent
    /// to the axis (e.g. N=9 at q≈0.4755, where the pair's Im → 0 only at q*) — that geometry is the scout's to
    /// count, not this population reading's.</summary>
    public RealResidualReading ReadRealResidual(int nBlock)
    {
        int k = nBlock - 1;
        int degree = PathKMonodromyScout.ResidualRootsExact(k, new Complex(QLo, 0)).Length;
        int realGeneric = RealResidualParts(k, GenericQ).Length;

        int steps = (int)Math.Round((QHi - QLo) / QStep);
        int maxReal = 0;
        double globalApproach = double.PositiveInfinity;
        for (int t = 0; t <= steps; t++)
        {
            var r = RealResidualParts(k, QLo + t * QStep);
            maxReal = Math.Max(maxReal, r.Length);
            for (int i = 0; i + 1 < r.Length; i++)
                globalApproach = Math.Min(globalApproach, r[i + 1] - r[i]);
        }
        return new RealResidualReading(nBlock, degree, realGeneric, maxReal, globalApproach);
    }

    /// <summary>Orthonormal R-even / R-odd basis columns of the reflection permutation, as sparse
    /// (index, weight) lists: a singleton (perm[t]==t, reflection-fixed) is the R-even unit e_t; a 2-cycle
    /// {t, t2} gives the R-even (e_t + e_t2)/√2 and the R-odd (e_t − e_t2)/√2.</summary>
    private static (List<(int Idx, double W)[]> Even, List<(int Idx, double W)[]> Odd, int Fixed) ParityColumns(int[] perm)
    {
        var even = new List<(int, double)[]>();
        var odd = new List<(int, double)[]>();
        var seen = new bool[perm.Length];
        int fixedCount = 0;
        double inv2 = 1.0 / Math.Sqrt(2.0);
        for (int t = 0; t < perm.Length; t++)
        {
            if (seen[t]) continue;
            int t2 = perm[t];
            if (t2 == t)
            {
                even.Add(new[] { (t, 1.0) });
                fixedCount++;
                seen[t] = true;
            }
            else
            {
                even.Add(new[] { (t, inv2), (t2, inv2) });
                odd.Add(new[] { (t, inv2), (t2, -inv2) });
                seen[t] = seen[t2] = true;
            }
        }
        return (even, odd, fixedCount);
    }

    /// <summary>Eigenvalues of the block restricted to a sector with orthonormal (sparse) columns:
    /// M[r,s] = colᵣᵀ · L · colₛ (the columns are real, so no conjugation), then a dense complex EVD.</summary>
    private static Complex[] SectorEigenvalues(Complex[,] L, List<(int Idx, double W)[]> cols)
    {
        int d = cols.Count;
        var M = Matrix<Complex>.Build.Dense(d, d);
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
            {
                Complex sum = Complex.Zero;
                foreach (var (a, wa) in cols[r])
                    foreach (var (b, wb) in cols[s])
                        sum += wa * L[a, b] * wb;
                M[r, s] = sum;
            }
        return M.Evd().EigenValues.ToArray();
    }

    /// <summary>Symmetric nearest-match distance between two complex multisets (sort by (Re, Im), compare
    /// element-wise). For a self-conjugate set, conj(σ) equals σ as a multiset, so the sorted comparison is
    /// ~0; otherwise it is O(spectral width). NaN if the multisets differ in size.</summary>
    private static double MatchDist(Complex[] a, Complex[] b)
    {
        if (a.Length != b.Length) return double.NaN;
        var sa = a.OrderBy(z => Math.Round(z.Real, 9)).ThenBy(z => Math.Round(z.Imaginary, 9)).ToArray();
        var sb = b.OrderBy(z => Math.Round(z.Real, 9)).ThenBy(z => Math.Round(z.Imaginary, 9)).ToArray();
        double m = 0;
        for (int i = 0; i < sa.Length; i++) m = Math.Max(m, (sa[i] - sb[i]).Magnitude);
        return m;
    }

    public string DisplayName =>
        "DiabolicReflectionParityWitness (the odd-N real-q diabolic onset, grounded: dimension-mismatch / sector-swap)";

    public string Summary
    {
        get
        {
            var odd = Read(7, GenericQ);
            var even = Read(8, GenericQ);
            return $"odd N (e.g. N={odd.N}): dim(even)−dim(odd) = {odd.EvenDim - odd.OddDim} = (N−1)/2 fixed singletons, " +
                   $"both sectors SELF-conjugate (defect {odd.SelfConjEven.ToString("E1", Inv)}/{odd.SelfConjOdd.ToString("E1", Inv)}) " +
                   $"⟹ R-even carries real eigenvalues that cross = the real-q diabolic. " +
                   $"even N (e.g. N={even.N}): balanced, σ_even = conj σ_odd (cross-conj {even.CrossConj.ToString("E1", Inv)}), " +
                   $"neither sector self-conjugate (defect {even.SelfConjEven.ToString("E1", Inv)}) ⟹ real collisions DEFECTIVE. " +
                   "The dimension-mismatch / sector-swap grounding of the odd-N ≥ 7 onset.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (int n in SweepN)
            {
                var r = Read(n, GenericQ);
                bool oddN = n % 2 == 1;
                string verdict = oddN
                    ? $"dim(even)−dim(odd) = {r.EvenDim - r.OddDim} = (N−1)/2 reflection-fixed singletons; " +
                      $"both sectors SELF-conjugate (defect {r.SelfConjEven.ToString("E2", Inv)} / {r.SelfConjOdd.ToString("E2", Inv)}); " +
                      "cross-pairing impossible (dims differ). R-even carries real eigenvalues ⟹ a real-q semisimple (diabolic) crossing is possible."
                    : $"dim(even) = dim(odd) = {r.EvenDim}; σ_even = conj σ_odd (cross-conj {r.CrossConj.ToString("E2", Inv)}); " +
                      $"neither sector self-conjugate (defect {r.SelfConjEven.ToString("E2", Inv)}). A real-axis collision in R-even is the generic pseudo-Hermitian DEFECTIVE EP.";
                yield return new InspectableNode(
                    displayName: $"N={n} ({(oddN ? "odd" : "even")}), q={GenericQ.ToString("0.##", Inv)}, full dim {r.FullDim}",
                    summary: verdict,
                    provenance: NodeProvenance.Live);
            }

            foreach (var (n, q, lam) in new[] { (7, 1.1264, -4.942), (9, 0.4755, -5.424) })
            {
                var (l1, l2, gap) = ReproduceDiabolic(n, q, lam);
                yield return new InspectableNode(
                    displayName: $"gate: N={n} real-q diabolic reproduced in R-even (from the full block)",
                    summary: $"at q={q.ToString("0.####", Inv)} the two R-even eigenvalues nearest the scout's λ={lam.ToString("0.###", Inv)} " +
                             $"are {l1.ToString("0.####", Inv)} and {l2.ToString("0.####", Inv)}, gap {gap.ToString("E2", Inv)} " +
                             "(a real-q semisimple coalescence, the C# scout's value recovered from below).",
                    provenance: NodeProvenance.Live);
            }

            // The within-odd grounding: the REAL R-even residual population vs N (even N: none; odd N: present, growing).
            foreach (int n in SweepN)
            {
                var rr = ReadRealResidual(n);
                bool oddN = n % 2 == 1;
                string verdict = rr.MaxRealCount == 0
                    ? $"R-even residual (degree F_{rr.ResidualDegree}) carries NO real eigenvalues anywhere in q∈[{QLo.ToString("0.#", Inv)},{QHi.ToString("0.#", Inv)}] " +
                      "(even N: the realness lives across the sectors) ⟹ no real strands ⟹ cannot host a real-q diabolic in R-even."
                    : $"R-even residual (degree F_{rr.ResidualDegree}) carries real eigenvalues (up to {rr.MaxRealCount} over the window, " +
                      $"{rr.RealCountGenericQ} at q={GenericQ.ToString("0.#", Inv)}); the scout's real-q diabolics are real-λ coalescences among these. " +
                      (rr.GlobalApproach < 1e-3
                          ? $"Two real strands CROSS in the window (closest approach {rr.GlobalApproach.ToString("E1", Inv)}, a geometry-A real-real crossing)."
                          : $"The closest two real strands approach to {rr.GlobalApproach.ToString("0.###", Inv)} (no real-real crossing here; an odd-N real-q diabolic can still be a conjugate-pair tangency).");
                yield return new InspectableNode(
                    displayName: $"within-odd grounding: N={n} ({(oddN ? "odd" : "even")}) real R-even residual population",
                    summary: verdict,
                    provenance: NodeProvenance.Live);
            }
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
