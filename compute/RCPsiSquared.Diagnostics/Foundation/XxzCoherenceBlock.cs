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
/// XY (SE,DE) block exactly (the trusted anchor).
///
/// <para>The Δ-test object and what it reads. The certifier (<see cref="CertifyCoalescenceNear"/>,
/// <see cref="CertifySplitUnderDelta"/>, and the refinement inside <see cref="TrackDiabolicUnderDelta"/>) follows
/// one eigenvalue pair and its discriminant f(q) = (λ_a − λ_b)². A semisimple crossing is a DOUBLE zero of f (the
/// gap closes linearly in q); a SIMPLE zero forces a Jordan EP2 (the gap closes as √|q − q*|). The kinds are told
/// apart by f together with geo against alg, not by f alone. At every sampled Δ=0 seed (N = 4, 5, 6, 7) the pair's
/// zero is double and the compression reads geo = alg = 2; at the sampled Δ &gt; 0 no crossing is certified, and
/// every zero the certifier claims inside the seed's isolation disk is simple and reads alg = 2, geo = 1 (N=4:
/// both zeros at real q on Re λ = −4, departures 0.0176 and 0.0222 at Δ = 0.02).
/// The N=5 defective control is a simple zero at Δ=0 and stays one. This carries DIABOLIC_BY_INTEGRABILITY's
/// N=4 reading to N = 5, 6, 7 as a finite-N Delta response; it does not by itself isolate integrability as the
/// cause.</para></summary>
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

    /// <summary>The full (SE,DE) coherence block at complex q and real Δ. <paramref name="bondWeights"/>
    /// scales the hopping bond by bond (length n−1, in chain order); null is the uniform chain and
    /// reproduces the plain block exactly. A non-uniform profile breaks the site reflection unless it is
    /// itself palindromic, so the R-sector builders below stay meaningful only on a palindromic profile.</summary>
    public static Complex[,] BuildFull(int n, Complex q, double delta, double[] bondWeights = null)
    {
        if (bondWeights is not null && bondWeights.Length != n - 1)
            throw new ArgumentException($"expected {n - 1} bond weights", nameof(bondWeights));
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
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] +=
                                new Complex(0, -2) * q * (bondWeights?[Math.Min(s, s2)] ?? 1.0);
            for (int s = 0; s < n; s++)                                  // bra excitation hops +2qi
                if ((bc & (1 << s)) != 0)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] +=
                                new Complex(0, 2) * q * (bondWeights?[Math.Min(s, s2)] ?? 1.0);
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

    /// <summary>Columns spanning the R = −1 (site reflection) sector of the full block. Reflection-fixed
    /// coherences carry no odd component and are dropped; each free pair contributes e_col − e_mirror.
    /// <paramref name="integerBasis"/> leaves the entries at ±1, so U<sup>T</sup>·M·U is an exact
    /// rearrangement with no 1/√2 rounding in it and residuals there can be compared to 0.0; the default
    /// normalises to an orthonormal basis. The two differ by the exact factor 2 in U<sup>T</sup>U.</summary>
    public static Matrix<Complex> BuildOddColumns(int n, bool integerBasis = false)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        double scale = integerBasis ? 1.0 : 1.0 / Math.Sqrt(2);
        var cols = new List<Complex[]>();
        var handled = new HashSet<int>();
        for (int col = 0; col < d; col++)
        {
            if (handled.Contains(col)) continue;
            var (kc, bc) = basis[col];
            int mirror = index[(Reflect(n, kc), Reflect(n, bc))];
            handled.Add(col);
            if (mirror == col) continue;                    // reflection-fixed: purely even
            handled.Add(mirror);
            var v = new Complex[d];
            v[col] = scale; v[mirror] = -scale;
            cols.Add(v);
        }
        return Matrix<Complex>.Build.Dense(d, cols.Count, (r, c) => cols[c][r]);
    }

    /// <summary>The site reflection s → n−1−s as a permutation matrix on the full block. Exact: both sides
    /// of R·L = L·R are the same float values rearranged, so the commutator is compared to 0.0.</summary>
    public static Matrix<Complex> ReflectionPermutation(int n)
    {
        var (basis, index) = Basis(n);
        int d = basis.Count;
        var r = Matrix<Complex>.Build.Dense(d, d);
        for (int i = 0; i < d; i++)
        {
            var (kc, bc) = basis[i];
            r[index[(Reflect(n, kc), Reflect(n, bc))], i] = Complex.One;
        }
        return r;
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

    /// <summary>The Riesz-projector reading of a coalescence in the symmetric (SE,DE) sector at (q, Δ, λ)
    /// (<see cref="EpCharacter"/>). The load-bearing discriminant is Geometric vs Algebraic on the compression
    /// (R-3): Geometric = Algebraic ⟹ semisimple (diabolic), Geometric &lt; Algebraic ⟹ a Jordan block
    /// (defective). The reading's own Kind does not decide it; see <see cref="CertifiedCharacterVerdict"/>.</summary>
    public static EpCharacter.Reading CharacterAt(int n, Complex q, double delta, Complex lambda, double radius)
        => EpCharacter.Characterize(BuildSym(n, q, delta), lambda, radius);

    /// <summary>Propose a nearby coalescence, refine it by the discriminant Newton and certify it with the same
    /// full-block certificate as <see cref="TrackDiabolicUnderDelta"/>. The caller's lambdaSeed must isolate a
    /// pair at Delta=0; the located pair must remain in that isolation disk, coincide within the independent 1e-6
    /// bound, and read an isolated algebraic-2 character. Otherwise the result is Uncertified, with unavailable
    /// character and unknown survival. A search null never certifies lifting. This box-scan locator is
    /// unsupported at N&gt;=7 and returns Uncertified there; the seed-following <see cref="CertifyCoalescenceNear"/>
    /// has no such limit.</summary>
    public static DeltaTrackResult CharacterAtDiabolicNear(
        int n, double delta, Complex qSeed, Complex lambdaSeed, double cell = 0.01)
        => TrackDiabolicUnderDelta(n, qSeed, lambdaSeed, delta, boxHalf: 0, boxCell: cell);

    // The residualOnly Δ-track (N>=6) with LOCAL continuity, so the box scan + descent never re-track from the
    // base per probe (the nested O(boxScan x trackSteps) cost). Anchor the residual roots at qSeed once (one
    // global track), then identify the residual subset at each probe by sub-stepped nearest continuity from the
    // running point. Residual-labeled candidates are carried by greedy nearest-unused local matching;
    // proposal tracking can exchange identity near dense or degenerate encounters. The proposal is then refined
    // by the discriminant Newton on the full R-even block and certified there (CertifyCoalescence).
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
        return CertifyCoalescence(n, delta, center, mid, lambdaSeed, seedRadius);
    }

    // Compression supplies proposals only: its fixed Delta=0 AT complement need not remain invariant.
    // A small compressed gap never certifies a full-block coalescence; the proposal is refined by the
    // discriminant Newton on the full R-even block and certified there (CertifyCoalescence).
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
        return CertifyCoalescence(n, delta, qd, mid, lambdaSeed, seedRadius);
    }

    // Coincidence/correspondence bound in the gamma=1 generator's eigenvalue units, independent of the caller's
    // proposal-search tolerances. Its law: after the discriminant Newton the residual pair gap sits at the
    // eigensolver's rounding floor, about u·‖M‖ times the eigenvector condition at a semisimple crossing, and
    // 2·√(departure·u·‖M‖) at a Jordan EP2 (the square-root branch turns a backward error u·‖M‖ into that split).
    // The bound sits above both while departure·‖M‖ stays below about 2e3. A pattern-search proposal alone, whose
    // gap at an EP2 stays near 1e-4 at |δq| ~ 1e-8, cannot pass it; that is why the Newton refinement runs first.
    internal const double FullBlockCoincidenceTolerance = 1e-6;

    // The floor a Jordan departure must clear, as a fraction of max(‖A‖_F, 1). Its law: EpCharacter computes the
    // departure as √(‖A‖_F² − Σ|λ|²), a difference of two numbers near ‖A‖_F², so on a semisimple compression it
    // reads √(k·u)·‖A‖_F rather than 0; over 4000 random semisimple 2×2 compressions (|λ| up to 11) k ≤ 17, i.e.
    // ≤ 4.3e-8·‖A‖_F. The floor sits 23 times above that and 31 times below the smallest sampled Jordan
    // departure relative to ‖A‖_F (1.98e-4 at ‖A‖_F ≈ 6.4, the N=7 q = 0.6788 EP2s at Δ = 1e-4).
    internal const double DepartureRoundingFloor = 1e-6;

    /// <summary>Character from the compression's multiplicities, the R-3 discriminant: alg = 2 and geo = 2 ⟹
    /// Diabolic; alg = 2, geo = 1 and a departure above the rounding floor ⟹ Defective (a Jordan block); anything
    /// else Uncertified.
    ///
    /// <para>EpCharacter's Kind is not consulted. It labels a small-departure Jordan block as Normal: it compares
    /// the departure with 1e-2·max(1, ‖A‖_F), and ‖A‖_F ≈ √2·|λ| for a 2×2 compression at λ, so its
    /// Defective/Normal split depends on where the eigenvalue sits rather than on the Jordan coupling. Both N=4
    /// Δ = 0.02 EP2s (departures 0.0176 and 0.0222 at |λ| ≈ 4.2) and the N=7 Δ = 0.02 EP2 with departure 0.0015
    /// at |λ| ≈ 4.94 read Kind = Normal, while alg = 2 and geo = 1 at all three.</para>
    ///
    /// <para>The geo read has its own law. At a located EP2 with residual gap g the matrix A − λ̄I has singular
    /// values ≈ departure and ≈ g²/(4·departure); EpCharacter.Nullity counts a singular value as zero below
    /// t = max(1e-6·σ_max, 1e-7·max(‖A‖₂, 1)), so geo reads 1 when g²/(4·departure) &lt; t &lt; departure. At a
    /// semisimple crossing both singular values are ≈ g/2, so geo = 2 needs g below 2t, which the Newton-located
    /// crossings meet with margin. A split crossing whose g is too large reads geo = 1 with a departure at its
    /// rounding level, below the departure floor, and is Uncertified rather than Defective.</para>
    ///
    /// <para>What Diabolic resolves. A Jordan coupling below t reads geo = 2, and two zeros closer than the winding
    /// radius read as one double zero, so Diabolic means: no split resolvable at radius 1e-7 in q and no coupling
    /// above about 1e-7·‖A‖. At the sampled Δ = 0 crossings the located gaps are 1e-15 to 1e-11, far inside that
    /// resolution. At a small nonzero Δ it is a statement about resolution only: at N=7 q = 0.6788 the coupling (≈ 2Δ) and the zero separation (≈ 0.56Δ) both stay
    /// under it below Δ ≈ 1.5e-7, where the certificate reads Diabolic although the crossing has split.</para></summary>
    internal static DeltaFlipVerdict CertifiedCharacterVerdict(int algebraic, int geometric, double departure,
        double departureFloor)
    {
        if (algebraic != 2) return DeltaFlipVerdict.Uncertified;
        if (geometric == 2) return DeltaFlipVerdict.Diabolic;
        if (geometric == 1 && departure > departureFloor) return DeltaFlipVerdict.Defective;
        return DeltaFlipVerdict.Uncertified;
    }

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
        var verdict = CertifiedCharacterVerdict(rr.Algebraic, rr.Geometric, rr.Departure,
            DepartureRoundingFloor * Math.Max(rr.CompressionNorm, 1.0));
        if (verdict == DeltaFlipVerdict.Uncertified)
            return new DeltaTrackResult(verdict, 0, 0, double.NaN, qd, mid, fullGap);
        return new DeltaTrackResult(verdict, rr.Algebraic, rr.Geometric, rr.Departure, qd, mid, fullGap);
    }

    // ---- The discriminant-Newton certifier ----

    // Central-difference step for f' and f''. f is holomorphic in q, so the difference along the real direction
    // is f'(q) itself.
    internal const double NewtonDifferenceStep = 1e-6;
    // One step never moves q by more than this; a zero farther away is reached in several steps, not one jump.
    internal const double NewtonStepCap = 0.05;
    // Stop once a step is below this fraction of max(1, |q|). Convergence is quadratic at both zero orders, so the
    // step after such a step is already at the double-rounding floor of q.
    internal const double NewtonStopRelative = 1e-12;
    internal const int NewtonMaxIterations = 60;

    // The zero-order circle about a located zero: radius and sampling of the discriminant's winding number.
    internal const double WindingRadius = 1e-7;
    internal const int WindingPoints = 64;
    // The circle closes on its opening value, so the phase increments telescope to 2π times an integer up to the
    // rounding of each principal argument, a few u per step: about WindingPoints·4u/2π ≈ 5e-15 turns (measured
    // ≤ 2e-15 from 20 Newton starts per sampled zero). The slack sits 200 times above the estimate. An integer residue is a check on the
    // arithmetic only; a closed sum is an integer whatever happens on the circle.
    internal const double WindingIntegerSlack = 1e-12;
    // The continuity guard, separate from the residue: on the circle a simple zero turns f by 2π/64 ≈ 0.098 per
    // step and a double zero by 0.196. An increment of π/4 or more is not a continuous step of the same pair (a
    // partner switch, or another zero within reach of the circle), and the winding is then not read.
    internal const double WindingMaxIncrement = Math.PI / 4;

    /// <summary>The pair discriminant at (q, Δ): f = (λ_a − λ_b)² for the two R-even eigenvalues nearest
    /// <paramref name="target"/>, divided by (q − deflate) when a known zero is deflated. Returns f, the pair
    /// midpoint and the undeflated pair gap |λ_a − λ_b|.</summary>
    private static (Complex F, Complex Mid, double Gap) PairDiscriminant(int n, double delta, Complex q,
        Complex target, Complex? deflate)
    {
        var pair = SeDeSymSpectrum(n, q, delta).OrderBy(z => (z - target).Magnitude).Take(2).ToArray();
        var split = pair[0] - pair[1];
        var f = split * split;
        if (deflate.HasValue) f /= q - deflate.Value;
        return (f, (pair[0] + pair[1]) / 2, split.Magnitude);
    }

    /// <summary>The discriminant-Newton locator. f(q) = (λ_a − λ_b)² is a symmetric function of the pair, so it is
    /// holomorphic in q wherever the pair stays isolated, the coalescence included: a semisimple crossing is a
    /// DOUBLE zero of f, a Jordan EP2 a SIMPLE zero. Newton is run on u = f/f', which has only simple zeros (step
    /// −u/u', u' = 1 − f·f''/f'²), so it converges quadratically at both without a multiplicity guess; f' and f''
    /// by central differences. The pair followed is the one nearest the running midpoint, starting from
    /// <paramref name="lambdaStart"/>. With <paramref name="deflate"/> the known zero is divided out of f, so the
    /// iteration cannot return to it. Returns the last iterate, its pair midpoint and gap, and the step count.
    /// The reached gap is the rounding floor of the eigensolver: about 1e-8 at the sampled EP2s (the square-root
    /// law), 1e-15 to 1e-11 at the sampled crossings.</summary>
    public static (Complex Q, Complex Lambda, double Gap, int Iterations) LocateCoalescence(
        int n, double delta, Complex qStart, Complex lambdaStart, Complex? deflate = null)
    {
        var q = qStart;
        var (f, mid, gap) = PairDiscriminant(n, delta, q, lambdaStart, deflate);
        for (int it = 0; it < NewtonMaxIterations; it++)
        {
            if (f == Complex.Zero) return (q, mid, gap, it);
            double h = NewtonDifferenceStep;
            var fPlus = PairDiscriminant(n, delta, q + h, mid, deflate).F;
            var fMinus = PairDiscriminant(n, delta, q - h, mid, deflate).F;
            var d1 = (fPlus - fMinus) / (2 * h);
            var d2 = (fPlus - 2 * f + fMinus) / (h * h);
            if (d1 == Complex.Zero) return (q, mid, gap, it);
            var u = f / d1;
            var du = Complex.One - f * d2 / (d1 * d1);
            var step = du == Complex.Zero ? -u : -u / du;
            if (step.Magnitude > NewtonStepCap) step *= NewtonStepCap / step.Magnitude;
            q += step;
            (f, mid, gap) = PairDiscriminant(n, delta, q, mid, deflate);
            if (step.Magnitude <= NewtonStopRelative * Math.Max(1.0, q.Magnitude)) return (q, mid, gap, it + 1);
        }
        return (q, mid, gap, NewtonMaxIterations);
    }

    /// <summary>The winding number of the pair discriminant f about <paramref name="qCenter"/> on a circle of
    /// radius <paramref name="radius"/>: the number of zeros of f inside, with multiplicity, while the pair stays
    /// isolated on the circle. A simple zero reads 1, a double zero 2 (a semisimple crossing is one; the
    /// character is told by geo against alg, not by f alone). The circle closes on its opening value, so the sum
    /// is exact up to rounding. Returned unrounded; NaN when a step fails the continuity guard
    /// (<see cref="WindingMaxIncrement"/>).</summary>
    public static double DiscriminantWinding(int n, double delta, Complex qCenter, Complex lambda,
        double radius = WindingRadius, int points = WindingPoints)
    {
        double turns = 0;
        var first = PairDiscriminant(n, delta, qCenter + radius, lambda, null).F;
        var previous = first;
        for (int k = 1; k <= points; k++)
        {
            var f = k == points
                ? first
                : PairDiscriminant(n, delta, qCenter + Complex.FromPolarCoordinates(radius, 2 * Math.PI * k / points),
                    lambda, null).F;
            double step = (f / previous).Phase;
            if (!(Math.Abs(step) < WindingMaxIncrement)) return double.NaN;
            turns += step;
            previous = f;
        }
        return turns / (2 * Math.PI);
    }

    // The seed's Delta=0 isolation disk: half the distance from lambdaSeed to the third-nearest eigenvalue of the
    // Delta=0 R-even block at qSeed (infinite with fewer than three); the two nearest must lie inside it.
    private static (double Radius, bool Isolated) SeedDisk(int n, Complex qSeed, Complex lambdaSeed)
    {
        var d = SeDeSymSpectrum(n, qSeed, 0).Select(z => (z - lambdaSeed).Magnitude).OrderBy(x => x).ToArray();
        double radius = d.Length > 2 ? .5 * d[2] : double.PositiveInfinity;
        return (radius, d[1] < radius);
    }

    /// <summary>Refine (q, mid) by <see cref="LocateCoalescence"/>, certify the reached pair with
    /// <see cref="CertifyFullBlockProposal"/>, and read the order of the discriminant's zero there
    /// (<see cref="DiscriminantWinding"/>). A certified character also needs an exact zero inside the winding
    /// circle (the existence certificate the gap alone does not give), with the order the character implies: a
    /// Diabolic reading needs order 2, since a semisimple double eigenvalue makes the discriminant vanish to
    /// second order; a Defective reading needs order 1, the certified object being an EP2. Anything else (an
    /// order-2 Jordan collision included) is Uncertified.</summary>
    private static DeltaTrackResult CertifyCoalescence(int n, double delta, Complex q, Complex mid,
        Complex lambdaSeed, double seedRadius)
    {
        var (qd, lambda, _, _) = LocateCoalescence(n, delta, q, mid);
        var result = CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, lambda, lambdaSeed, seedRadius);
        if (result.Verdict == DeltaFlipVerdict.Uncertified) return result;
        double winding = DiscriminantWinding(n, delta, qd, result.LambdaCandidate);
        int order = double.IsFinite(winding) ? (int)Math.Round(winding) : 0;
        bool consistent = double.IsFinite(winding) && Math.Abs(winding - order) < WindingIntegerSlack
                          && order == (result.Verdict == DeltaFlipVerdict.Diabolic ? 2 : 1);
        return consistent
            ? result with { DiscriminantZeroOrder = order }
            : new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qd, result.LambdaCandidate, result.Gap);
    }

    /// <summary>Follow the pair nearest lambdaSeed from qStart by the discriminant Newton at this Δ and certify the
    /// coalescence it reaches (<see cref="CertifyCoalescence"/>). The pair must be isolated in the Δ=0 disk of
    /// (qStart, lambdaSeed) and the reached midpoint must stay in it. No box scan is involved: the pair is followed,
    /// not picked from the smallest gap of a dense spectrum, so there is no N limit.</summary>
    public static DeltaTrackResult CertifyCoalescenceNear(int n, double delta, Complex qStart, Complex lambdaSeed)
    {
        var pair = SeDeSymSpectrum(n, qStart, delta).OrderBy(z => (z - lambdaSeed).Magnitude).Take(2).ToArray();
        var (seedRadius, isolated) = SeedDisk(n, qStart, lambdaSeed);
        if (!isolated)
            return new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qStart, lambdaSeed,
                (pair[0] - pair[1]).Magnitude);
        return CertifyCoalescence(n, delta, qStart, (pair[0] + pair[1]) / 2, lambdaSeed, seedRadius);
    }

    /// <summary>The two simple zeros a Δ=0 double zero of the pair discriminant splits into at this Δ. The Δ=0
    /// crossing q₀ is located from the seed first; the first zero at Δ is followed from q₀, the second from the
    /// reflection 2q₀ − z₁ with z₁ divided out of f, so it cannot return to z₁. Each is certified as in
    /// <see cref="CertifyCoalescenceNear"/>, both against the seed's Δ=0 disk. A zero that is not reached (from the
    /// N=7 q = 2.628 seed at Δ = 0.10), or whose midpoint has left the disk (the N=7 q = 1.1264 zeros at Δ = 0.10),
    /// comes back Uncertified; that says nothing about its existence elsewhere. Δ must be nonzero (at Δ = 0 there
    /// is no split, and both results are Uncertified); below Δ ≈ 1.5e-7 the split is under the certificate's
    /// resolution and First can read Diabolic (<see cref="CertifiedCharacterVerdict"/>).</summary>
    public static (DeltaTrackResult First, DeltaTrackResult Second) CertifySplitUnderDelta(
        int n, Complex qSeed, Complex lambdaSeed, double delta)
    {
        var (seedRadius, isolated) = SeedDisk(n, qSeed, lambdaSeed);
        if (!isolated || delta == 0.0)
        {
            var none = new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, qSeed, lambdaSeed, double.NaN);
            return (none, none);
        }
        var (q0, lambda0, _, _) = LocateCoalescence(n, 0.0, qSeed, lambdaSeed);
        var (q1, lambda1, _, _) = LocateCoalescence(n, delta, q0, lambda0);
        var first = CertifyCoalescence(n, delta, q1, lambda1, lambdaSeed, seedRadius);
        var (q2, lambda2, _, _) = LocateCoalescence(n, delta, 2 * q0 - q1, lambda0, deflate: q1);
        var second = CertifyCoalescence(n, delta, q2, lambda2, lambdaSeed, seedRadius);
        // Two zeros are distinct only if their winding circles do not overlap.
        if ((second.QCandidate - first.QCandidate).Magnitude <= 2 * WindingRadius)
            second = new DeltaTrackResult(DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, second.QCandidate,
                second.LambdaCandidate, second.Gap);
        return (first, second);
    }

    /// <summary>The verdict of a Δ-track step. DIABOLIC = a semisimple coalescence (geo = alg = 2, the
    /// discriminant vanishing to second order); DEFECTIVE = a Jordan EP (geo &lt; alg), at the sampled loci a
    /// simple zero of the discriminant, one of the two a Δ=0 double zero splits into; LIFTED = the degeneracy is gone
    /// (requiring a separate exclusion certificate; this tracker does not issue LIFTED from a search null). This is a
    /// finite-N Delta response, compared with a defective control and consistent with the conditional residual
    /// mechanism. DEFECTIVE at a sampled nonzero Δ does not alone prove causality or all-N protection; a diabolic
    /// double zero unfolds into EP2s under a generic perturbation. DIABOLIC at a sampled nonzero Δ would mean no
    /// split resolvable at the certificate's resolution (radius 1e-7 in q, coupling about 1e-7·‖A‖; see
    /// <see cref="CertifiedCharacterVerdict"/>), which below Δ ≈ 1.5e-7 a split crossing also gives.
    /// UNCERTIFIED means no corresponding isolated full-block coincident pair or no definite
    /// Diabolic/Defective character was certified numerically. It is neither survival nor death.
    /// Algebraic/Geometric=0 and Departure=NaN then denote unavailable character, not measured multiplicities.</summary>
    public enum DeltaFlipVerdict { Diabolic, Defective, Lifted, Uncertified }

    /// <summary>One certified (or Uncertified) reading. <paramref name="DiscriminantZeroOrder"/> is the winding of
    /// the pair discriminant about QCandidate: 1 = a simple zero (the square-root branch of an EP2), 2 = a double
    /// zero (a linear crossing); 0 = not read (Uncertified, or the matrix-level certificate, which has no q
    /// family).</summary>
    public sealed record DeltaTrackResult(
        DeltaFlipVerdict Verdict, int Algebraic, int Geometric, double Departure,
        Complex QCandidate, Complex LambdaCandidate, double Gap, int DiscriminantZeroOrder = 0)
    {
        /// <summary>True only for a certified Diabolic character; false is not a death verdict.</summary>
        public bool IsCertifiedDiabolic => Verdict == DeltaFlipVerdict.Diabolic;

        /// <summary>True for Diabolic (no split at the certificate's resolution), false for Defective/Lifted, null
        /// for Uncertified (unknown).</summary>
        public bool? Survived => Verdict switch
        {
            DeltaFlipVerdict.Diabolic => true,
            DeltaFlipVerdict.Defective or DeltaFlipVerdict.Lifted => false,
            _ => null
        };
    }

    /// <summary>Search near qSeed at this Δ, refine the proposal by the discriminant Newton, then certify an
    /// isolated full-block pair, not a search null. For every proposal path, lambdaSeed selects the candidate within
    /// its initial Delta=0 isolation disk (half the distance from lambdaSeed to the third eigenvalue).
    /// Leaving that disk is Uncertified, not a claim of identity transport or non-existence.
    /// Default/residual modes at N&gt;=7 are Uncertified because the dense AT spectrum invalidates their box-scan
    /// locator. The legacy exact flag selects compressed proposals, not an invariant residual restriction.
    /// All paths ignore legacy coalesceTol/depTol for certification: both full-block roots must lie within
    /// 1e-6 of the refined midpoint and each other (gamma=1 units), an isolated pair must be enclosed, the
    /// compression must read alg = 2 with geo = 2 (Diabolic) or geo = 1 above the departure rounding floor (Defective), and
    /// the discriminant must have a zero of matching order inside the winding circle; otherwise Uncertified, with
    /// unknown survival.</summary>
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
        // The pattern search only proposes q; the discriminant Newton then locates the coalescence itself.
        var qd = PathKMonodromyScout.GapRefine(roots, boxArg, boxCell);
        var proposedPair = roots(qd).OrderBy(z => (z - lambdaSeed).Magnitude).Take(2).ToArray();
        var mid = (proposedPair[0] + proposedPair[1]) / 2;
        return CertifyCoalescence(n, delta, qd, mid, lambdaSeed, seedRadius);
    }
}
