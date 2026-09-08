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
/// <item>the Δ·ZZ term is DIAGONAL in the computational basis and leaves the dissipator unchanged.
/// At real q the AT identity Re λ = −2γ⟨n_XY⟩ still holds, but the eigenmodes and their rates can change.
/// The added diagonal frequency is −i·qΔ·(zz(ket) − zz(bra)), with
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

    /// <summary>The full (SE,DE) block plus the diagonal term -i*q*(fieldEnergy(ket)-fieldEnergy(bra)),
    /// with fieldEnergy(c) = Σ_k w_k·z_k (z_k = −1 if site k excited, +1 else).
    /// Field strength is dimensionless and scaled by q; for complex q this need not be purely imaginary.
    /// The unchanged dissipator entries do not imply unchanged eigenmode real parts.
    /// A generic w can break the S₂ reflection and sector conjugation symmetry; Stage 2 uses OffReal + BuildFull.
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

    /// <summary>The residual labels are chosen exactly at the base (q0=2, Delta=0) via
    /// <see cref="PathKMonodromyScout.ResidualIndices"/> (the XXZ block equals the F89 block at Delta=0).
    /// The corresponding full-spectrum strands are tracked by nearest-neighbour continuation to (q, Delta).
    /// This is finite-path labeling, not a proof of a globally Delta-invariant AT/residual split.
    /// Full-block character remains required.</summary>
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

    /// <summary>Propose a nearby coalescence, then use the same strict full-block certificate as
    /// <see cref="TrackDiabolicUnderDelta"/>. The caller's lambdaSeed must isolate a pair at Delta=0;
    /// the candidate must remain in that isolation disk, coincide within the independent 1e-6 bound,
    /// and have isolated algebraic-2 Diabolic/Defective character. Otherwise the result is Uncertified,
    /// with unavailable character and unknown survival. A search null never certifies lifting.
    /// This full-spectrum locator is unsupported at N&gt;=7 and returns Uncertified there.</summary>
    public static DeltaTrackResult CharacterAtDiabolicNear(
        int n, double delta, Complex qSeed, Complex lambdaSeed, double cell = 0.01)
        => TrackDiabolicUnderDelta(n, qSeed, lambdaSeed, delta, boxHalf: 0, boxCell: cell);

    // The residualOnly Δ-track (N>=6) with LOCAL continuity, so the box scan + descent never re-track from the
    // base per probe (the nested O(boxScan x trackSteps) cost). Anchor the residual roots at qSeed once (one
    // global track), then identify the residual subset at each probe by sub-stepped nearest continuity from the
    // running point. Residual-labeled candidates are carried by greedy nearest-unused local matching;
    // proposal tracking can exchange identity near dense or degenerate encounters.
    // geo/alg via EpCharacter on the full block at the candidate pair's
    // midpoint with an AT-aware radius (nearest non-pair eigenvalue over the FULL block).
    private static DeltaTrackResult TrackDiabolicUnderDeltaResidual(
        int n, Complex qSeed, Complex lambdaSeed, double seedRadius, double delta, double boxHalf, double boxCell, int trackSteps)
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
        // The search proposes q only. Select the pair by the supplied seed, never an unrelated global minimum.
        var pair = cur.OrderBy(z => (z - lambdaSeed).Magnitude).Take(2).ToArray();
        var mid = (pair[0] + pair[1]) / 2;
        return CertifyFullBlockProposal(BuildSym(n, center, delta), center, mid, lambdaSeed, seedRadius);
    }

    // Compression supplies proposals only: its fixed Delta=0 AT complement need not remain invariant.
    // A small compressed gap never certifies a full-block coalescence or a defect/lift verdict.
    private static DeltaTrackResult TrackDiabolicUnderDeltaCompressed(
        int n, Complex qSeed, Complex lambdaSeed, double seedRadius, double delta, double boxHalf, double boxCell)
    {
        int k = n - 1;
        Func<Complex, Complex[]> resRoots = q => PathKMonodromyScout.ResidualRootsCompressedXxz(k, q, delta);

        // Box scan and refinement minimize the compressed proposal gap, not a full-block gap.
        double boxMin = double.PositiveInfinity; Complex boxArg = qSeed;
        int steps = Math.Max(1, (int)Math.Round(2 * boxHalf / boxCell));
        for (int ir = 0; ir <= steps; ir++)
            for (int ii = 0; ii <= steps; ii++)
            {
                var q = new Complex(qSeed.Real - boxHalf + ir * boxCell, qSeed.Imaginary - boxHalf + ii * boxCell);
                double g = PathKMonodromyScout.MinGap(resRoots(q));
                if (g < boxMin) { boxMin = g; boxArg = q; }
            }
        var qd = PathKMonodromyScout.GapRefine(resRoots, boxArg, boxCell);
        var cur = resRoots(qd);

        // Select the full-block pair nearest the proposal midpoint, not the full spectrum's global minimum.
        int ai = 0, bi = 1; double bb = double.PositiveInfinity;
        for (int i = 0; i < cur.Length; i++)
            for (int j = i + 1; j < cur.Length; j++)
            { double g = (cur[i] - cur[j]).Magnitude; if (g < bb) { bb = g; ai = i; bi = j; } }
        var mid = (cur[ai] + cur[bi]) / 2;
        return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);
    }

    // Numerical coincidence/correspondence floor in the gamma=1 generator's eigenvalue units.
    // This is independent of caller-controlled proposal-search tolerances, not an exact rank proof.
    internal const double FullBlockCoincidenceTolerance = 1e-6;

    internal static DeltaFlipVerdict CertifiedCharacterVerdict(EpCharacter.EpKind kind) => kind switch
    {
        EpCharacter.EpKind.Diabolic => DeltaFlipVerdict.Diabolic,
        EpCharacter.EpKind.Defective => DeltaFlipVerdict.Defective,
        _ => DeltaFlipVerdict.Uncertified
    };

    internal static DeltaTrackResult CertifyFullBlockProposal(Matrix<Complex> block, Complex qd, Complex mid,
        Complex? lambdaSeed = null, double seedRadius = double.PositiveInfinity)
    {
        var full = block.Evd().EigenValues.OrderBy(z => (z - mid).Magnitude).ToArray();
        double fullGap = (full[0] - full[1]).Magnitude;
        // Both roots must correspond to the proposed midpoint; an unrelated coincident AT pair is no certificate.
        if (lambdaSeed.HasValue && !((mid - lambdaSeed.Value).Magnitude < seedRadius)
            || !double.IsFinite(fullGap) || fullGap > FullBlockCoincidenceTolerance
            || !((full[0] - mid).Magnitude <= FullBlockCoincidenceTolerance)
            || !((full[1] - mid).Magnitude <= FullBlockCoincidenceTolerance))
            return new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qd, mid, fullGap);
        mid = (full[0] + full[1]) / 2;
        var ds = full.Select(z => (z - mid).Magnitude).OrderBy(x => x).ToArray();
        double radius = ds.Length > 2 ? Math.Min(0.4 * ds[2], 0.5) : 0.1;
        if (!(radius > 2 * ds[1]))
            return new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qd, mid, fullGap);
        var rr = EpCharacter.Characterize(block, mid, radius);
        var verdict = rr.Algebraic == 2 ? CertifiedCharacterVerdict(rr.Kind) : DeltaFlipVerdict.Uncertified;
        if (verdict == DeltaFlipVerdict.Uncertified)
            return new DeltaTrackResult(verdict, 0, 0, double.NaN, qd, mid, fullGap);
        return new DeltaTrackResult(verdict, rr.Algebraic, rr.Geometric, rr.Departure, qd, mid, fullGap);
    }

    /// <summary>The verdict of a Δ-track step. DIABOLIC = the coalescence survives semisimply (geo=alg,
    /// dep≈0); DEFECTIVE = it persists as a Jordan EP (geo&lt;alg); LIFTED = the degeneracy is gone (no
    /// coalescence in the local q-box, requiring a separate exclusion certificate; this tracker does not
    /// issue LIFTED from a search null). This is a finite-N Delta response, compared with a defective
    /// control and consistent with the conditional residual mechanism. DEFECTIVE or LIFTED at a sampled
    /// nonzero Δ does not alone prove causality or all-N protection. DIABOLIC survival would falsify
    /// the defect-or-lift prediction at that sampled locus and Δ.
    /// UNCERTIFIED means no corresponding isolated full-block coincident pair or no definite
    /// Diabolic/Defective character was certified numerically. It is neither survival nor death.
    /// Algebraic/Geometric=0 and Departure=NaN then denote unavailable character, not measured multiplicities.</summary>
    public enum DeltaFlipVerdict { Diabolic, Defective, Lifted, Uncertified }

    public sealed record DeltaTrackResult(
        DeltaFlipVerdict Verdict, int Algebraic, int Geometric, double Departure,
        Complex QCandidate, Complex LambdaCandidate, double Gap)
    {
        /// <summary>True only for a certified Diabolic character; false is not a death verdict.</summary>
        public bool IsCertifiedDiabolic => Verdict == DeltaFlipVerdict.Diabolic;

        /// <summary>True for Diabolic, false for Defective/Lifted, null for Uncertified (unknown).</summary>
        public bool? Survived => Verdict switch
        {
            DeltaFlipVerdict.Diabolic => true,
            DeltaFlipVerdict.Defective or DeltaFlipVerdict.Lifted => false,
            _ => null
        };
    }

    /// <summary>Search near qSeed at this Δ, then certify an isolated full-block pair, not a search null.
    /// For every proposal path, lambdaSeed selects the candidate within its initial
    /// Delta=0 isolation disk (half the distance from lambdaSeed to the third eigenvalue).
    /// Leaving that disk is Uncertified, not a claim of identity transport or non-existence.
    /// Default/residual modes at N&gt;=7 are Uncertified because the dense AT spectrum invalidates their locator.
    /// The legacy exact flag selects compressed proposals, not an invariant residual restriction.
    /// All paths ignore legacy coalesceTol/depTol for certification: both full-block roots must lie within
    /// 1e-6 of the proposal midpoint and each other (gamma=1 units), an isolated pair must be enclosed,
    /// and EpCharacter must return Diabolic or Defective; otherwise Uncertified, with unknown survival.</summary>
    public static DeltaTrackResult TrackDiabolicUnderDelta(int n, Complex qSeed, Complex lambdaSeed, double delta,
        double boxHalf = 0.04, double boxCell = 0.008, double coalesceTol = 1e-3, double depTol = 1e-6,
        bool residualOnly = false, int trackSteps = 160, bool exact = false)
    {
        var seedPair = SeDeSymSpectrum(n, qSeed, delta).OrderBy(z => (z - lambdaSeed).Magnitude).Take(2).ToArray();
        if (n >= 7 && !exact)
            return new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qSeed,
                lambdaSeed, (seedPair[0] - seedPair[1]).Magnitude);
        var seedDistances = SeDeSymSpectrum(n, qSeed, 0).Select(z => (z - lambdaSeed).Magnitude).OrderBy(x => x).ToArray();
        double seedRadius = seedDistances.Length > 2 ? .5 * seedDistances[2] : double.PositiveInfinity;
        if (!(seedDistances[1] < seedRadius))
            return new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qSeed,
                lambdaSeed, (seedPair[0] - seedPair[1]).Magnitude);

        // Compressed proposals use the same caller-seed isolation disk as every other path.
        // A coincident pair outside it is not evidence about the requested seed, even at Delta=0.
        if (exact)
            return TrackDiabolicUnderDeltaCompressed(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell);

        // Local residual labels propose q; they do not establish an invariant subset or exclude AT capture.
        if (residualOnly)
            return TrackDiabolicUnderDeltaResidual(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell, trackSteps);

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
        // Refinement only proposes q; a bounded search cannot certify lifting or seed identity.
        var qd = PathKMonodromyScout.GapRefine(roots, boxArg, boxCell);
        var proposedPair = roots(qd).OrderBy(z => (z - lambdaSeed).Magnitude).Take(2).ToArray();
        var mid = (proposedPair[0] + proposedPair[1]) / 2;
        return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);
    }
}
