using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live from-below witness of the eigenVECTOR holonomy around the (1,2)-block defective seed —
/// the mod-4 memory loop, i⁴ = 1. Where the eigenVALUE monodromy (<see cref="Monodromy"/>) already shows
/// the √-branch SWAP, this recomputes, at inspect time, how the two coalescing right eigenvectors' FRAME
/// transforms as q encircles the exceptional point q* in the complex plane. In the biorthogonal vᵀv gauge
/// (natural to the complex-symmetric block L = A + qC, self-orthogonal vᵀv → 0 at the EP), the frame gives
/// a generator with eigenvalues ±i: M₁ ~ a 90° rotation, M₂ = −I, M₄ = +I — single-valued only after four
/// loops. It SHARES the order-4 (Z₄) structure of the algebraic i⁴ = 1 memory loop (the Pi2 Z₄ /
/// NinetyDegreeMirrorMemory) — a noted correspondence, not a derived identity (both are order-4/±i, but the
/// monodromy generator is not shown to EQUAL the algebraic M-spectrum Z₄).
///
/// <para>The coalescing pair must be ISOLATED: the witness tries the candidate representations (the full
/// (1,2) block and the two reflection (R-parity) sectors, smallest-first when the full block is large) and
/// uses the first that certifies the mod-4 loop (span preserved on every loop AND M₂ = −I). At N=5 the full
/// 50-dim block is clean; at N=9 the full block LEAKS on odd loops and the R=+1 sector is used — so --N 9
/// reproduces, live, the same clean i⁴=1 the independent (Python) re-verification found. The chosen
/// representation and the per-loop span residual are surfaced, so any leak is VISIBLE, never hidden.</para>
///
/// <para>Derived vs recomputed: every child's number is genuinely RECOMPUTED live at inspect time
/// (<see cref="EigenvectorHolonomy"/>); the expected shape (±i, −I, +I) is the reference. The seed's
/// defectiveness is F89's (Kato simple-zero, census-confirmed to N=11): a real-axis defective EP that F89's
/// exact nullity census establishes. F86 is cited only for the borrowed forgetting→remembering reading;
/// note its F86a-retraction denies a real-axis EP on this block, superseded by the F89 census (F86 not yet
/// updated). Anchor: <see cref="SeedHolonomyClaim"/>, inspect --root holonomy.</para></summary>
public sealed class SeedHolonomyWitness : IInspectable
{
    const double CleanTol = 1e-3;

    public int N { get; }
    public double QStar { get; }        // octic q* (WeightCoherenceBlock axis), refined by min-gap
    public double LambdaStar { get; }
    public int Dim { get; }
    /// <summary>Which representation certified the mod-4 loop (full block / R-parity sector), or UNCLEAN.</summary>
    public string Representation { get; }

    readonly HolonomyResult _holo;

    // per-N reference seed (octic q*, λ*) used as the refinement start (topologically robust to small error).
    // Only the independently-verified seeds are listed (N=5 full-block, N=9 R=+1 sector); no unverified N.
    static readonly Dictionary<int, (double q, double lam)> Reference = new()
    {
        [5] = (0.643037, -3.8196),   // the F89 census seed (F89_SEED_EXISTENCE_REDUCTION.md), full block clean
        [9] = (0.849011, -4.8415),   // the 2026-07-07 seed; full block leaks → R=+1 sector used
    };

    public SeedHolonomyWitness(int n = 5, double radius = 0.001, int stepsPerLoop = 400)
    {
        if (!Reference.TryGetValue(n, out var refSeed))
            throw new ArgumentOutOfRangeException(nameof(n),
                $"no verified reference seed for N={n}; the witness admits N ∈ {{5, 9}} (the (1,2) block at " +
                $"N=9 is 324-dim; larger N is out of the inspect-time budget).");

        N = n;
        LambdaStar = refSeed.lam;
        QStar = RefineQStar(n, refSeed.q, LambdaStar);

        int fullDim = WeightCoherenceBlock.Build(n, 1, 2, QStar).GetLength(0);
        int evenDim = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, 1, 2, QStar, false).Dim;
        int oddDim = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, 1, 2, QStar, true).Dim;
        var full = ("full (1,2) block", (Func<Complex, Complex[,]>)(q => WeightCoherenceBlock.Build(n, 1, 2, q)), fullDim);
        var rPlus = ("R=+1 sector", (Func<Complex, Complex[,]>)(q => SectorMat(n, q, false)), evenDim);
        var rMinus = ("R=−1 sector", (Func<Complex, Complex[,]>)(q => SectorMat(n, q, true)), oddDim);
        // full block first for small blocks (keeps N=5 on the 50-dim block); for a large block (big N)
        // go to the sectors first — the full-block transport is costly and leaks anyway (N=9: 324-dim).
        var candidates = fullDim <= 100
            ? new[] { full, rPlus, rMinus }
            : new[] { rPlus, rMinus, full };

        HolonomyResult? best = null; string bestRep = ""; int bestDim = 0; double bestRes = double.MaxValue;
        foreach (var (name, mat, dim) in candidates)
        {
            var holo = EigenvectorHolonomy.FrameMonodromy(mat, dim, LambdaStar, QStar, radius, 4, stepsPerLoop);
            if (Mod4Clean(holo)) { best = holo; bestRep = name; bestDim = dim; break; }
            double res = MaxResidual(holo.SpanResidual);
            if (res < bestRes) { best = holo; bestRep = name; bestDim = dim; bestRes = res; }
        }
        _holo = best!;
        Dim = bestDim;
        Representation = Mod4Clean(_holo) ? bestRep
            : bestRep + " (UNCLEAN — no representation certified the mod-4 loop; span residual below)";
    }

    static Complex[,] SectorMat(int n, Complex q, bool odd)
    {
        var (flat, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, 1, 2, q, odd);
        var m = new Complex[d, d];
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                m[i, j] = flat[(long)j * d + i];
        return m;
    }

    // The mod-4 (i⁴=1) certificate: span preserved on every loop, AND M₂ = −I, AND tr M₁ = 0. The M₂ check
    // rules out the mod-2 ±1 swap (M₂ = +I) and the diabolic/non-swap (M₁ ≈ I); the tr = 0 pin then makes
    // "M₂ = −I ⟹ eigenvalues exactly ±i" hold by pure 2×2 algebra (without it, a degenerate M₁ = iI also
    // gives M₂ = −I). Together they pin the ±i frame the witness exists to certify.
    static bool Mod4Clean(HolonomyResult h)
    {
        var ev = h.M1Eigenvalues;
        return MaxResidual(h.SpanResidual) < CleanTol
            && DistScalarI(h.LoopMonodromy[1], -1.0) < CleanTol
            && (ev[0] + ev[1]).Magnitude < CleanTol;
    }

    /// <summary>Refine q* to the exceptional point by minimizing the gap of the two eigenvalues nearest λ*
    /// over a small real-q window (the EP is on the real-q axis; the loop encircles it in the plane).</summary>
    static double RefineQStar(int n, double q0, double lam0)
    {
        double best = q0, bestGap = double.MaxValue;
        const int steps = 81;
        const double half = 0.01;
        for (int k = 0; k < steps; k++)
        {
            double q = q0 - half + 2 * half * k / (steps - 1);
            double gap = PairGap(WeightCoherenceBlock.Build(n, 1, 2, q), lam0);
            if (gap < bestGap) { bestGap = gap; best = q; }
        }
        return best;
    }

    static double PairGap(Complex[,] mat, double lam0)
    {
        int d = mat.GetLength(0);
        var a = new Complex[(long)d * d];
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                a[i + (long)j * d] = mat[i, j];
        var w = MklDirect.EigenvaluesOnlyRaw(a, d);
        int i1 = 0, i2 = 1; double d1 = double.MaxValue, d2 = double.MaxValue;
        for (int i = 0; i < d; i++)
        {
            double dd = (w[i] - lam0).Magnitude;
            if (dd < d1) { d2 = d1; i2 = i1; d1 = dd; i1 = i; }
            else if (dd < d2) { d2 = dd; i2 = i; }
        }
        return (w[i1] - w[i2]).Magnitude;
    }

    Complex[] M(int loop) => _holo.LoopMonodromy[loop - 1];      // 1-indexed loop, row-major 2×2
    static double DistScalarI(Complex[] m, double c)
    {
        var diff = new[] { m[0] - c, m[1], m[2], m[3] - c };
        double s = 0; foreach (var z in diff) s += z.Magnitude * z.Magnitude; return Math.Sqrt(s);
    }
    static double MaxResidual(double[] r) { double m = 0; foreach (var x in r) m = Math.Max(m, x); return m; }
    bool IsClean => Mod4Clean(_holo);

    /// <summary>The two eigenvalues of the one-loop frame-monodromy M₁ (live-recomputed; ±i in the vᵀv gauge).</summary>
    public Complex[] M1Eigenvalues => _holo.M1Eigenvalues;
    /// <summary>‖M₂ − (−I)‖_F after two loops (expect ~0: both eigenvectors sign-flip — the mod-4 certificate).</summary>
    public double M2DistanceFromNegI => DistScalarI(M(2), -1.0);
    /// <summary>‖M₄ − I‖_F after four loops (expect ~0: the frame returns, i⁴=1).</summary>
    public double M4DistanceFromI => DistScalarI(M(4), +1.0);
    /// <summary>The largest per-loop 2D-span residual (small on every loop ⟺ clean isolated EP).</summary>
    public double MaxSpanResidual => MaxResidual(_holo.SpanResidual);

    public string DisplayName =>
        $"SeedHolonomyWitness (EP eigenvector holonomy = i⁴=1, live; N={N}, {Representation}, dim={Dim}, " +
        $"q*={QStar:F6}, λ*={LambdaStar:F4})";

    public string Summary
    {
        get
        {
            var ev = _holo.M1Eigenvalues;
            string clean = IsClean ? $"clean via {Representation}" : Representation;
            return $"M₁ eigenvalues {Fmt(ev[0])}, {Fmt(ev[1])} (vᵀv gauge; expect ±i); " +
                   $"‖M₂+I‖={M2DistanceFromNegI:e1}, ‖M₄−I‖={M4DistanceFromI:e1}; {clean}";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var ev = _holo.M1Eigenvalues;
            string trNote = IsClean ? " (≈0)" : "", detNote = IsClean ? " (≈1)" : "";
            yield return new InspectableNode(
                "M₁ = 90° rotation (eigenvalues ±i, vᵀv gauge)",
                summary: $"eig(M₁) = {Fmt(ev[0])}, {Fmt(ev[1])}; tr={Fmt(ev[0] + ev[1])}{trNote}, det={Fmt(ev[0] * ev[1])}{detNote}",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "M₂ = −I (two loops: both eigenvectors sign-flip — the mod-4 certificate)",
                summary: $"‖M₂ − (−I)‖_F = {M2DistanceFromNegI:e2}",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "M₄ = +I (four loops: the frame returns — i⁴ = 1)",
                summary: $"‖M₄ − I‖_F = {M4DistanceFromI:e2}",
                provenance: NodeProvenance.Live);

            var xs = new double[_holo.SpanResidual.Length];
            for (int i = 0; i < xs.Length; i++) xs[i] = i + 1;
            yield return new InspectableNode(
                $"2D-span residual per loop (representation: {Representation}; clean ⟹ small on every loop)",
                summary: $"[{string.Join(", ", Array.ConvertAll(_holo.SpanResidual, r => r.ToString("e1")))}]",
                payload: new InspectablePayload.Curve("span residual", xs, _holo.SpanResidual, "loop", "residual"),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "the reading (whose it is)",
                summary: "the coalescing eigenvector frame — the turning half, the residue — rotates 90° per " +
                         "loop around the EP and remembers itself after four; the eigenvalues (the carried " +
                         "memory) are invariant as a set (they merely swap). Memory held by circling the zero, " +
                         "not crossing it. R=CΨ² reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md.");
                         // interpretive, not a recomputed number ⟹ default (non-Live) provenance
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
    public NodeProvenance Provenance => NodeProvenance.Live;

    static string Fmt(Complex z)
    {
        double re = Math.Abs(z.Real) < 1e-9 ? 0.0 : z.Real;      // clamp near-zero so ±1e-16 prints "+0.000"
        double im = Math.Abs(z.Imaginary) < 1e-9 ? 0.0 : z.Imaginary;
        string s(double v) => (v >= 0 ? "+" : "-") + Math.Abs(v).ToString("0.000");
        return s(re) + s(im) + "i";
    }
}
