using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The (SE,DE) = (weight-1 ket, weight-2 bra) coherence block of the XXZ-chain Liouvillian,
/// parameterized by q = J/γ AND the ZZ-anisotropy Δ, for the F89 path-4 (N=5) integrability test
/// (docs/superpowers/plans/2026-06-27-f89-path4-delta-test.md). Convention (bare Pauli, matching
/// F89Path3OcticBlock / f89_zz_break_gate.py): H(Δ) = J·Σ_b(X_bX_{b+1}+Y_bY_{b+1}) + J·Δ·Σ_b Z_bZ_{b+1},
/// γ=1, q=J. L = −i[H,·] + Z-dephasing:
/// <list type="bullet">
/// <item>dephasing diagonal −2·n_diff (Absorption Theorem 2γ·HammingDistance(ket,bra), ∈{−2,−6} here);</item>
/// <item>ket excitation hops −2qi, bra excitation hops +2qi (the −i(H⊗I − I⊗Hᵀ) split), NN, exclusion;</item>
/// <item>the Δ·ZZ term is DIAGONAL in the computational basis and Hermitian, so the AT rate is UNTOUCHED
/// (Re λ = −2γ⟨n_XY⟩ preserved); it adds only the frequency −i·qΔ·(zz(ket) − zz(bra)), with
/// zz(c) = Σ_bonds(+1 if the two sites are equal, −1 if they differ) = Σ_b ⟨c|Z_bZ_{b+1}|c⟩.</item>
/// </list>
/// The S₂ site-reflection (s → N−1−s) R=+1 symmetric sector carries the diabolics (same sector as
/// F89Path3OcticBlock.BuildSeDeSymBlock at N=4 / the path-k residual at N≥5). At Δ=0 this reproduces the
/// XY (SE,DE) block exactly (the trusted anchor). The Δ-test object: does a path-4 diabolic flip defective
/// or lift as Δ turns on, generalizing DIABOLIC_BY_INTEGRABILITY's N=4 gate off N=4.</summary>
public static class XxzCoherenceBlock
{
    // all bitmasks on n sites with exactly w set bits (excitations), ascending.
    private static List<int> Weight(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (System.Numerics.BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }

    // zz(c) = Σ_{bond (b,b+1)} ⟨c|Z_bZ_{b+1}|c⟩ = Σ_b (+1 if bits b,b+1 equal, −1 if differ).
    private static int Zz(int n, int c)
    {
        int s = 0;
        for (int b = 0; b < n - 1; b++)
            s += (((c >> b) & 1) == ((c >> (b + 1)) & 1)) ? 1 : -1;
        return s;
    }

    // reverse the n-bit string (the site reflection s → n−1−s).
    private static int Reflect(int n, int m)
    {
        int r = 0;
        for (int s = 0; s < n; s++)
            if ((m & (1 << s)) != 0) r |= 1 << (n - 1 - s);
        return r;
    }

    private static (List<(int ket, int bra)> basis, Dictionary<(int, int), int> index) Basis(int n)
    {
        var kets = Weight(n, 1);
        var bras = Weight(n, 2);
        var basis = new List<(int, int)>();
        var index = new Dictionary<(int, int), int>();
        foreach (var k in kets)
            foreach (var b in bras) { index[(k, b)] = basis.Count; basis.Add((k, b)); }
        return (basis, index);
    }

    /// <summary>The full (SE,DE) coherence block at complex q and real Δ.</summary>
    public static Complex[,] BuildFull(int n, Complex q, double delta)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        var l = new Complex[d, d];
        for (int col = 0; col < d; col++)
        {
            var (kc, bc) = basis[col];
            int nDiff = System.Numerics.BitOperations.PopCount((uint)(kc ^ bc));
            // AT dephasing rate (−2·n_diff) + the Δ·ZZ frequency (−i·qΔ·(zz_ket − zz_bra)), q complex.
            l[col, col] += new Complex(-2.0 * nDiff, 0)
                         + (-Complex.ImaginaryOne) * q * (delta * (Zz(n, kc) - Zz(n, bc)));
            for (int s = 0; s < n; s++)                                  // ket excitation hops −2qi
                if ((kc & (1 << s)) != 0)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] += new Complex(0, -2) * q;
            for (int s = 0; s < n; s++)                                  // bra excitation hops +2qi
                if ((bc & (1 << s)) != 0)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] += new Complex(0, 2) * q;
        }
        return l;
    }

    // fieldEnergy(config) = Σ_k w[k]·z_k, with z_k = −1 if site k is excited (bit set), +1 otherwise.
    private static double FieldEnergy(int n, double[] w, int config)
    {
        double e = 0;
        for (int k = 0; k < n; k++) e += w[k] * (((config >> k) & 1) == 1 ? -1.0 : 1.0);
        return e;
    }

    /// <summary>The full (SE,DE) block at (q, Δ) plus a per-site Z-field Σ_k w_k Z_k. The field is diagonal
    /// and Hermitian, so it leaves the absorption-theorem real rate untouched and only shifts the imaginary
    /// frequency by −i·q·(fieldEnergy(ket) − fieldEnergy(bra)), with fieldEnergy(c) = Σ_k w_k·z_k (z_k = −1
    /// if site k excited, +1 else). Field strength is dimensionless (scaled by q like the hopping). A random
    /// w breaks integrability, the S₂ reflection, AND conjugation symmetry (Stage 2 uses OffReal + BuildFull).
    /// w=null reproduces BuildFull(n,q,Δ).</summary>
    public static Complex[,] BuildFullWithField(int n, Complex q, double delta, double[] w)
    {
        var l = BuildFull(n, q, delta);
        if (w == null) return l;
        var (basis, _) = Basis(n);
        for (int col = 0; col < basis.Count; col++)
        {
            var (kc, bc) = basis[col];
            double fe = FieldEnergy(n, w, kc) - FieldEnergy(n, w, bc);
            l[col, col] += (-Complex.ImaginaryOne) * q * fe;
        }
        return l;
    }

    /// <summary>The R=+1 (site reflection s→N−1−s) symmetric sector of the full block. Same construction as
    /// F89Path3OcticBlock.BuildSeDeSymBlock, generalized to any N (the reflection commutes with XX+YY, ZZ and
    /// uniform dephasing, so the sector is invariant including the Δ·ZZ term).</summary>
    public static Matrix<Complex> BuildSym(int n, Complex q, double delta)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        var full = Matrix<Complex>.Build.DenseOfArray(BuildFull(n, q, delta));

        var cols = new List<Complex[]>();
        var handled = new HashSet<int>();
        for (int col = 0; col < d; col++)
        {
            if (handled.Contains(col)) continue;
            var (kc, bc) = basis[col];
            int mcol = index[(Reflect(n, kc), Reflect(n, bc))];
            var v = new Complex[d];
            if (mcol == col) v[col] = Complex.One;                      // reflection-fixed coherence
            else { double s = 1.0 / Math.Sqrt(2); v[col] = s; v[mcol] = s; handled.Add(mcol); }
            handled.Add(col);
            cols.Add(v);
        }
        var p = Matrix<Complex>.Build.Dense(d, cols.Count, (r, c) => cols[c][r]);
        return p.ConjugateTranspose() * full * p;
    }

    /// <summary>The spectrum of the R=+1 symmetric (SE,DE) sector at (q, Δ). At Δ=0, N=4 this reproduces
    /// F89Path3OcticBlock.BuildSeDeSymBlock(q, 1)'s spectrum (the trusted XY anchor).</summary>
    public static Complex[] SeDeSymSpectrum(int n, Complex q, double delta)
        => BuildSym(n, q, delta).Evd().EigenValues.ToArray();

    /// <summary>The residual (H_B-mixed, non-AT) strands of the XXZ (q,Δ) coherence block at (q, Δ), tracked by
    /// nearest-neighbour continuity from the base (q0=2, Δ=0). There the residual SET is identified exactly via
    /// <see cref="PathKMonodromyScout.ResidualIndices"/> (the XXZ block equals the F89 block at Δ=0). The split
    /// is Δ-stable: the ZZ term is Hermitian, so the AT rate Re λ = −2γ⟨n_XY⟩ is Δ-independent and the AT-locked
    /// half never changes membership. AT-free by construction — the AT-locked exact degeneracies that flood the
    /// full-block Δ-test box scan at N≥6 are simply not in this set.</summary>
    public static Complex[] ResidualRootsTrackedXxz(int k, Complex q, double delta, int trackSteps = 160)
    {
        int n = k + 1;
        var q0 = new Complex(2, 0);
        var r0 = SeDeSymSpectrum(n, q0, 0.0);
        var (residual, _) = PathKMonodromyScout.ResidualIndices(k, r0);
        // track the full spectrum along (q0,0) -> (q,Δ); cur[i] = position of the strand that started at index i.
        var cur = (Complex[])r0.Clone();
        for (int s = 1; s <= trackSteps; s++)
        {
            double t = (double)s / trackSteps;
            var next = SeDeSymSpectrum(n, q0 + (q - q0) * t, delta * t);
            var used = new bool[next.Length];
            var moved = new Complex[cur.Length];
            for (int i = 0; i < cur.Length; i++)                 // each strand -> its nearest unused next eigenvalue
            {
                int best = -1; double bd = double.PositiveInfinity;
                for (int j = 0; j < next.Length; j++)
                    if (!used[j]) { double dd = (next[j] - cur[i]).Magnitude; if (dd < bd) { bd = dd; best = j; } }
                used[best] = true; moved[i] = next[best];
            }
            cur = moved;
        }
        return residual.Select(i => cur[i]).ToArray();
    }

    /// <summary>The character (diabolic / defective / …) of a coalescence in the symmetric (SE,DE) sector at
    /// (q, Δ, λ): the load-bearing semisimplicity discriminant (EpCharacter Riesz projector). Geometric=
    /// Algebraic ∧ Departure≈0 ⟹ Diabolic (semisimple); Geometric&lt;Algebraic ⟹ Defective (Jordan).</summary>
    public static EpCharacter.Reading CharacterAt(int n, Complex q, double delta, Complex lambda, double radius)
        => EpCharacter.Characterize(BuildSym(n, q, delta), lambda, radius);

    /// <summary>One Δ-track datum: the character of the coalescence nearest qSeed at this Δ, with the
    /// refined q*, the merge λ*, and the min gap. <paramref name="gap"/> &gt; liftTol ⟹ the degeneracy has
    /// LIFTED (no EP nearby); otherwise Kind/Geometric/Departure read defective-vs-diabolic.</summary>
    public readonly record struct DeltaFlipReading(
        EpCharacter.EpKind Kind, int Algebraic, int Geometric, double Departure, double ProjectorNorm,
        Complex QStar, Complex LambdaStar, double Gap);

    /// <summary>Re-locate the coalescence nearest qSeed at this Δ (GapRefine on the symmetric-sector min-gap,
    /// since Δ shifts the diabolic) and characterize it with an adaptive radius (enclosing only the pair).
    /// The Δ&gt;0 datum of the flip test: at the N=4 q_EP this reads DIABOLIC at Δ=0 and DEFECTIVE at Δ&gt;0
    /// (the committed f89_zz_break_gate.py table). <paramref name="lambdaSeed"/> is the expected merge λ
    /// (advisory; the refiner finds the coalescence by gap).</summary>
    public static DeltaFlipReading CharacterAtDiabolicNear(
        int n, double delta, Complex qSeed, Complex lambdaSeed, double cell = 0.01)
    {
        Func<Complex, Complex[]> roots = qq => SeDeSymSpectrum(n, qq, delta);
        var qd = PathKMonodromyScout.GapRefine(roots, qSeed, cell);
        var spec = roots(qd);
        int ai = 0, bi = 1; double best = double.PositiveInfinity;
        for (int i = 0; i < spec.Length; i++)
            for (int j = i + 1; j < spec.Length; j++)
            {
                double g = (spec[i] - spec[j]).Magnitude;
                if (g < best) { best = g; ai = i; bi = j; }
            }
        var mid = (spec[ai] + spec[bi]) / 2;
        double distOther = double.PositiveInfinity;
        for (int i = 0; i < spec.Length; i++)
            if (i != ai && i != bi) distOther = Math.Min(distOther, (spec[i] - mid).Magnitude);
        double radius = Math.Clamp(0.4 * distOther, 0.05, 0.5);
        var r = EpCharacter.Characterize(BuildSym(n, qd, delta), mid, radius);
        return new DeltaFlipReading(r.Kind, r.Algebraic, r.Geometric, r.Departure, r.ProjectorNorm, qd, mid, best);
    }

    // The residualOnly Δ-track (N>=6) with LOCAL continuity, so the box scan + descent never re-track from the
    // base per probe (the nested O(boxScan x trackSteps) cost). Anchor the residual roots at qSeed once (one
    // global track), then identify the residual subset at each probe by sub-stepped nearest continuity from the
    // running point. AT strands are in the full spectrum but never adopted (a residual strand's continuation is
    // its own nearest root), so no AT capture. geo/alg via EpCharacter on the full block at the residual pair's
    // midpoint with an AT-aware radius (nearest non-pair eigenvalue over the FULL block).
    private static DeltaTrackResult TrackDiabolicUnderDeltaResidual(
        int n, Complex qSeed, double delta, double boxHalf, double boxCell, double coalesceTol, double depTol, int trackSteps)
    {
        int k = n - 1;
        var anchor = ResidualRootsTrackedXxz(k, qSeed, delta, trackSteps);
        Complex[] Local(Complex[] from, Complex qFrom, Complex qTo)
        {
            int sub = Math.Max(1, (int)Math.Ceiling((qTo - qFrom).Magnitude / 0.01));
            var pos = from;
            for (int s = 1; s <= sub; s++)
            {
                var full = SeDeSymSpectrum(n, qFrom + (qTo - qFrom) * ((double)s / sub), delta);
                var used = new bool[full.Length]; var next = new Complex[pos.Length];
                for (int i = 0; i < pos.Length; i++)
                {
                    int best = -1; double bd = double.PositiveInfinity;
                    for (int j = 0; j < full.Length; j++)
                        if (!used[j]) { double dd = (full[j] - pos[i]).Magnitude; if (dd < bd) { bd = dd; best = j; } }
                    used[best] = true; next[i] = full[best];
                }
                pos = next;
            }
            return pos;
        }
        // box scan (local from the anchor) re-finds the coalescence region under the Δ-shift.
        double boxMin = double.PositiveInfinity; Complex boxArg = qSeed; var boxRoots = anchor;
        int steps = Math.Max(1, (int)Math.Round(2 * boxHalf / boxCell));
        for (int ir = 0; ir <= steps; ir++)
            for (int ii = 0; ii <= steps; ii++)
            {
                var q = new Complex(qSeed.Real - boxHalf + ir * boxCell, qSeed.Imaginary - boxHalf + ii * boxCell);
                var rr = Local(anchor, qSeed, q);
                double g = PathKMonodromyScout.MinGap(rr);
                if (g < boxMin) { boxMin = g; boxArg = q; boxRoots = rr; }
            }
        // local descent on the residual min-gap from the box-min.
        Complex center = boxArg; var cur = boxRoots; double half = boxCell / 2, curGap = PathKMonodromyScout.MinGap(cur);
        for (int it = 0; it < 40 && half > 1e-8; it++)
        {
            Complex best = center; double bestGap = curGap; var bestRoots = cur; bool moved = false;
            foreach (var off in new[] { new Complex(half, 0), new Complex(-half, 0), new Complex(0, half), new Complex(0, -half),
                                        new Complex(half, half), new Complex(half, -half), new Complex(-half, half), new Complex(-half, -half) })
            {
                var probe = Local(cur, center, center + off);
                double g = PathKMonodromyScout.MinGap(probe);
                if (g < bestGap) { bestGap = g; best = center + off; bestRoots = probe; moved = true; }
            }
            if (moved) { center = best; curGap = bestGap; cur = bestRoots; } else half /= 2;
        }
        double refined = PathKMonodromyScout.MinGap(cur);
        if (refined > coalesceTol)
            return new DeltaTrackResult(DeltaFlipVerdict.Lifted, 0, 0, double.NaN, center, Complex.Zero, refined);
        // character: closest residual pair -> midpoint; AT-aware radius over the FULL block; EpCharacter (geo/alg).
        int ai = 0, bi = 1; double bb = double.PositiveInfinity;
        for (int i = 0; i < cur.Length; i++)
            for (int j = i + 1; j < cur.Length; j++)
            { double g = (cur[i] - cur[j]).Magnitude; if (g < bb) { bb = g; ai = i; bi = j; } }
        var mid = (cur[ai] + cur[bi]) / 2;
        var ds = SeDeSymSpectrum(n, center, delta).Select(z => (z - mid).Magnitude).OrderBy(x => x).ToArray();
        double radius = ds.Length > 2 ? Math.Clamp(0.4 * ds[2], 0.05, 0.5) : 0.1;
        var rr2 = EpCharacter.Characterize(BuildSym(n, center, delta), mid, radius);
        var verdict = (rr2.Geometric == rr2.Algebraic && rr2.Departure < depTol)
            ? DeltaFlipVerdict.Diabolic : DeltaFlipVerdict.Defective;
        return new DeltaTrackResult(verdict, rr2.Algebraic, rr2.Geometric, rr2.Departure, center, mid, bb);
    }

    /// <summary>The verdict of a Δ-track step. DIABOLIC = the coalescence survives semisimply (geo=alg,
    /// dep≈0); DEFECTIVE = it persists as a Jordan EP (geo&lt;alg); LIFTED = the degeneracy is gone (no
    /// coalescence in the local q-box). For the integrability test, DEFECTIVE or LIFTED at Δ&gt;0 confirms
    /// the diabolic was integrability-protected; DIABOLIC surviving at Δ&gt;0 would refute it.</summary>
    public enum DeltaFlipVerdict { Diabolic, Defective, Lifted }

    public sealed record DeltaTrackResult(
        DeltaFlipVerdict Verdict, int Algebraic, int Geometric, double Departure,
        Complex QStar, Complex LambdaStar, double Gap)
    {
        public bool Survived => Verdict == DeltaFlipVerdict.Diabolic;
    }

    /// <summary>Track the coalescence near qSeed at this Δ and classify it. R-4: a 2D q-BOX scan re-finds the
    /// region (not GapRefine-from-previous, which assumes continuity), THEN GapRefine from the box-minimum to
    /// resolve the cusp — this is essential because a DEFECTIVE EP is a √-branch (a SHARP cusp a coarse grid
    /// overshoots, so the box-grid-min alone overestimates the gap and reads a false LIFT; the descent drives
    /// a real coalescence's gap → 0, while a genuine LIFT stays at a shallow min above coalesceTol). Then read
    /// geo vs alg (the load-bearing discriminant, R-3): refined gap &gt; coalesceTol ⟹ LIFTED (degeneracy
    /// gone); else geo=alg ∧ dep≈0 ⟹ DIABOLIC (survives), geo&lt;alg ⟹ DEFECTIVE. Box ±boxHalf (the path-4
    /// diabolics' nearest neighbour is ~0.1 away, so ±0.04 is safe).</summary>
    public static DeltaTrackResult TrackDiabolicUnderDelta(int n, Complex qSeed, Complex lambdaSeed, double delta,
        double boxHalf = 0.04, double boxCell = 0.008, double coalesceTol = 1e-3, double depTol = 1e-6,
        bool residualOnly = false, int trackSteps = 160)
    {
        // residualOnly (N>=6): scan the residual strands only, so the box scan + refine cannot be captured by the
        // AT-locked exact degeneracies that crowd the full sym block at N>=6 (the local-tracking variant below).
        if (residualOnly)
            return TrackDiabolicUnderDeltaResidual(n, qSeed, delta, boxHalf, boxCell, coalesceTol, depTol, trackSteps);

        Func<Complex, Complex[]> roots = qq => SeDeSymSpectrum(n, qq, delta);
        double boxMin = double.PositiveInfinity; Complex boxArg = qSeed;
        int steps = Math.Max(1, (int)Math.Round(2 * boxHalf / boxCell));
        for (int ir = 0; ir <= steps; ir++)
            for (int ii = 0; ii <= steps; ii++)
            {
                var q = new Complex(qSeed.Real - boxHalf + ir * boxCell, qSeed.Imaginary - boxHalf + ii * boxCell);
                double g = PathKMonodromyScout.MinGap(roots(q));
                if (g < boxMin) { boxMin = g; boxArg = q; }
            }
        // GapRefine from the box-minimum resolves the √-cusp; the refined gap (not the box-grid-min) decides LIFT.
        var qd = PathKMonodromyScout.GapRefine(roots, boxArg, boxCell);
        double refined = PathKMonodromyScout.MinGap(roots(qd));
        if (refined > coalesceTol)
            return new DeltaTrackResult(DeltaFlipVerdict.Lifted, 0, 0, double.NaN, qd, Complex.Zero, refined);

        var r = CharacterAtDiabolicNear(n, delta, qd, lambdaSeed, boxCell);
        var verdict = (r.Geometric == r.Algebraic && r.Departure < depTol)
            ? DeltaFlipVerdict.Diabolic
            : DeltaFlipVerdict.Defective;
        return new DeltaTrackResult(verdict, r.Algebraic, r.Geometric, r.Departure, r.QStar, r.LambdaStar, r.Gap);
    }
}
