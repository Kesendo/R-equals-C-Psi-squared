using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One probed cell of the shell census: the (p,w) block × one shift of the s-invariant pair
/// {λ_A, μ = −λ_A − 2N}. Method: "lu-rparity" (probed, both R-parity sectors), "window" (analytically
/// excluded by the Bendixson rate window, margin recorded, no LU), "control" (window-excluded but probed
/// anyway to watch the margin sharpness), "deferred" (a sector exceeds the managed LP64 dim wall — named,
/// never silently skipped). σ_min semantics: an UPPER-bound-converging estimate (census evidence);
/// σ_min > m is the exclusion direction, membership verdicts lean on the containment corollary.</summary>
public sealed record ShellCensusEntry(
    int P, int W, long Dim, int DimEven, int DimOdd,
    string Shift, Complex ShiftValue,
    bool Probed, string Method,
    double SigmaMinEven, double SigmaMinOdd, double SigmaMin,
    bool Converged,
    double WindowMargin, int IterationsEven, int IterationsOdd, double Seconds);

/// <summary>The per-seed census result: the refined locus (two-stage: q* by golden-section gap
/// minimization on the (1,2) block, then λ_A = the near-defective pair midpoint), the member-noise
/// floor (the pair gap), the adaptive member threshold, and one entry per (block × shift) on the
/// fundamental-domain strip.</summary>
public sealed record ShellCensusResult(
    RealSeed Seed, double QRefined, Complex RefinedLambdaA, double PairGap,
    double MemberTol, bool SeedUsable,
    IReadOnlyList<ShellCensusEntry> Entries, TimeSpan Elapsed, SectorShellCensus.Options Opts)
{
    public ShellCensusSummary Summarize() => SectorShellCensus.Summarize(this);

    /// <summary>Tab-separated, InvariantCulture (the repo CSV convention).</summary>
    public void WriteCsv(string path)
    {
        var inv = CultureInfo.InvariantCulture;
        using var w = new StreamWriter(path);
        w.WriteLine($"# shell census N={Seed.N} qStar={Seed.QStar.ToString(inv)} qRefined={QRefined.ToString("F9", inv)} " +
                    $"lambdaA={RefinedLambdaA.Real.ToString("F9", inv)} pairGap={PairGap.ToString("E3", inv)} " +
                    $"memberTol={MemberTol.ToString("E3", inv)} rParity={Seed.RParity} usable={SeedUsable}");
        w.WriteLine("p\tw\tdim\tdimEven\tdimOdd\tshift\tshiftRe\tprobed\tmethod\tsigmaEven\tsigmaOdd\tsigmaMin\tconverged\twindowMargin\titersEven\titersOdd\tseconds");
        foreach (var e in Entries)
            w.WriteLine(string.Join("\t",
                e.P.ToString(inv), e.W.ToString(inv), e.Dim.ToString(inv), e.DimEven.ToString(inv), e.DimOdd.ToString(inv),
                e.Shift, e.ShiftValue.Real.ToString("F9", inv), e.Probed ? "1" : "0", e.Method,
                e.SigmaMinEven.ToString("E6", inv), e.SigmaMinOdd.ToString("E6", inv), e.SigmaMin.ToString("E6", inv),
                e.Converged ? "1" : "0", e.WindowMargin.ToString("E6", inv),
                e.IterationsEven.ToString(inv), e.IterationsOdd.ToString(inv), e.Seconds.ToString("F2", inv)));
    }
}

/// <summary>The three-valued census verdict (a clean partial run is PARTIAL, never DISAGREE): PASS =
/// found members equal the probeable expected set, nothing deferred, nothing ambiguous, all probed
/// non-members clear the exclusion threshold; PARTIAL = the same but some expected members sit behind
/// the LP64 wall (listed); DISAGREE = anything else.</summary>
public sealed record ShellCensusSummary(
    string Verdict,
    IReadOnlySet<(int P, int W, string Shift)> Expected,
    IReadOnlySet<(int P, int W, string Shift)> ExpectedProbeable,
    IReadOnlySet<(int P, int W, string Shift)> FoundMembers,
    IReadOnlyList<(int P, int W, string Shift)> DeferredMembers,
    IReadOnlyList<(int P, int W, string Shift)> Ambiguous,
    IReadOnlyList<(int P, int W, string Shift)> Unresolved,
    double WorstNonMemberSigma, double MemberTol, double PairGap);

/// <summary>The step-3 shell-census engine (the sectorbraid large-N exclusion program,
/// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6-§7): at a real defective seed locus of the (1,2) block,
/// probe every fundamental-domain block's σ_min(L(p,w) − s) for both shifts s ∈ {λ_A, μ = −λ_A − 2N},
/// LU-free of full spectra (ShiftedSigmaMin), split by R-parity (the ~¼ LU-cost lever; members carry
/// their value in the SEED's parity sector since W, Klein, fold and transpose all commute with R).
/// Probe policy: LU only where the Bendixson window contains Re s (elsewhere the window-shell lemma
/// already excludes analytically — margin recorded); sectors built and probed strictly sequentially;
/// sector dims above the managed LP64 wall are DEFERRED by name. Expected members come from the
/// containment corollary (<see cref="ExpectedMembers"/>), computed, not hardcoded.</summary>
public static class SectorShellCensus
{
    public sealed class Options
    {
        /// <summary>The managed-array LP64 wall: flat Complex[] caps at 46340² elements.</summary>
        public int MaxSectorDim { get; set; } = 46000;

        /// <summary>The member cut in units of the refined pair gap: MemberTol = max(1e-9,
        /// MemberCutFactor · pairGap). FIRST-LIGHT CALIBRATION (N=9, 2026-07-04): a member's σ_min is
        /// bounded by the SHIFT UNCERTAINTY (≤ ~1.5·pairGap; at a defective locus it reads
        /// quadratically deep, ~(gap/2)², the Jordan pseudospectrum signature — measured 5.5e-13),
        /// while the nearest DENSITY NEIGHBOR (a distinct eigenvalue, itself W-nested along the
        /// diagonal chain) sat at 2.7e-4 at the densest low-q seed — so a single adaptive cut at
        /// 10·pairGap separates by orders of magnitude. An absolute floor (the earlier 1e-3) is the
        /// calibration bug this replaces: it swallowed density neighbors as members.</summary>
        public double MemberCutFactor { get; set; } = 10.0;

        /// <summary>The gray band around the member cut: σ_min ∈ [MemberTol/3, 3·MemberTol] is
        /// AMBIGUOUS (flagged, never silently classified).</summary>
        public double AmbiguousBandFactor { get; set; } = 3.0;

        /// <summary>Pure REPORTING threshold: probed non-members with σ_min below this are listed as
        /// spectrally NEAR (local density information). NOT verdict-blocking — exclusion in the census
        /// is the byte-identical-sharing sense, decided by the member cut.</summary>
        public double NearTol { get; set; } = 1e-2;

        /// <summary>A seed whose refined pair gap exceeds this is UNUSABLE (refinement failed).</summary>
        public double MaxUsablePairGap { get; set; } = 1e-3;

        /// <summary>Extra (window-excluded) cells to probe anyway, watching margin sharpness.</summary>
        public List<(int P, int W, string Shift)> Controls { get; } = new();

        /// <summary>Route sectors past the LP64 wall (<see cref="MaxSectorDim"/>) through the sparse
        /// instruments instead of emitting a bare "deferred" row: MEMBER cells to the from-above
        /// W-transport witness (<see cref="SectorWitnessTransport.MemberUpperBounds"/>, method
        /// "sparse-witness"), NON-MEMBER / control cells to the sparse inverse-power estimator
        /// (<see cref="SparseShiftedSigmaMin"/>, method "sparse-invit"). Off restores the pure-dense
        /// landed behavior (the bare "deferred" rows). Calibrated at N=9 (SparseCensusCalibrationTests).</summary>
        public bool SparseForDeferred { get; set; } = true;

        /// <summary>Options for the sparse inverse-power estimator on the deferred non-member / control
        /// cells (and a member cell's secondary other-parity confirmation).</summary>
        public SparseShiftedSigmaMin.Options Sparse { get; set; } = SparseShiftedSigmaMin.Options.Default;

        public Action<string>? Log { get; set; }
    }

    /// <summary>The containment corollary's member set on the fundamental-domain strip: λ_A on the band
    /// blocks (p,p+1) ∩ FD, μ on the transpose-normalized fold images (p, N−p−1). The loop bound
    /// p ∈ [1, N−2] already keeps every fold image's popcounts in [1, N−2] and inside the FD.</summary>
    public static HashSet<(int P, int W, string Shift)> ExpectedMembers(int n)
    {
        var set = new HashSet<(int, int, string)>();
        for (int p = 1; p + 1 <= n - 1; p++)
        {
            if (BlockLattice.InFundamentalDomain(n, p, p + 1)) set.Add((p, p + 1, "lambdaA"));
            (int fp, int fw) = (p, n - (p + 1));
            var rep = fp <= fw ? (fp, fw) : (fw, fp);
            set.Add((rep.Item1, rep.Item2, "mu"));
        }
        return set;
    }

    /// <summary>Two-stage seed refinement: (a) golden-section minimization of the (1,2) block's
    /// near-defective pair gap over q ∈ [q*−5e-4, q*+5e-4] (the EP's √|q−q*| law makes the recorded
    /// few-decimal q* a ~1e-3 member-noise floor; this pushes it to ~1e-5); (b) at the refined q*, the
    /// minimal-gap pair within |λ − recorded λ_A| &lt; 0.05, midpoint = the probe shift. The pair search
    /// runs in the SEED'S R-PARITY SECTOR of (1,2) (the recorded parity — the full-block search can lock
    /// onto a wrong-parity density pair at low q and slide the bracket; the N=9 R-odd seed 0.511958 was
    /// the case that forced this). Exact semisimple degeneracies (gap &lt; 1e-12) are ignored — the
    /// at_masking trap.</summary>
    public static (double QRefined, Complex Lambda, double PairGap) RefineSeed(RealSeed seed)
    {
        double Gap(double q, out Complex mid)
        {
            var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(
                seed.N, 1, 2, new Complex(q, 0), odd: seed.RParity < 0);
            var m = Matrix<Complex>.Build.Dense(d, d);
            for (int c = 0; c < d; c++)
                for (int r = 0; r < d; r++)
                    m[r, c] = a[(long)c * d + r];
            var ev = m.Evd().EigenValues.Enumerate()
                .Where(z => Math.Abs(z.Real - seed.LambdaA) < 0.05).ToList();
            double best = double.PositiveInfinity;
            mid = new Complex(seed.LambdaA, 0);
            for (int i = 0; i < ev.Count; i++)
                for (int j = i + 1; j < ev.Count; j++)
                {
                    double g = (ev[i] - ev[j]).Magnitude;
                    if (g > 1e-12 && g < best) { best = g; mid = (ev[i] + ev[j]) / 2.0; }
                }
            return best;
        }

        double lo = seed.QStar - 5e-4, hi = seed.QStar + 5e-4;
        const double phi = 0.6180339887498949;
        double x1 = hi - phi * (hi - lo), x2 = lo + phi * (hi - lo);
        double f1 = Gap(x1, out var m1), f2 = Gap(x2, out var m2);
        for (int it = 0; it < 40; it++)
        {
            if (f1 <= f2) { hi = x2; x2 = x1; f2 = f1; m2 = m1; x1 = hi - phi * (hi - lo); f1 = Gap(x1, out m1); }
            else { lo = x1; x1 = x2; f1 = f2; m1 = m2; x2 = lo + phi * (hi - lo); f2 = Gap(x2, out m2); }
        }
        return f1 <= f2 ? (x1, m1, f1) : (x2, m2, f2);
    }

    private static int SafeInt(long v) => v > int.MaxValue ? int.MaxValue : (int)v;

    public static ShellCensusResult Run(RealSeed seed, Options opts)
    {
        int n = seed.N;
        var sw = Stopwatch.StartNew();
        var (qR, lambda, pairGap) = RefineSeed(seed);
        double memberTol = Math.Max(1e-9, opts.MemberCutFactor * pairGap);
        bool usable = pairGap < opts.MaxUsablePairGap;
        var entries = new List<ShellCensusEntry>();
        opts.Log?.Invoke($"seed N={n} q*={qR:F9} lambda_A={lambda.Real:F9} pairGap={pairGap:E3} memberTol={memberTol:E3} usable={usable}");
        if (!usable)
            return new ShellCensusResult(seed, qR, lambda, pairGap, memberTol, false, entries, sw.Elapsed, opts);

        var q = new Complex(qR, 0);
        var shifts = new (string Kind, Complex Value)[] { ("lambdaA", lambda), ("mu", -lambda - 2.0 * n) };
        var controls = new HashSet<(int, int, string)>(opts.Controls);
        if (n - 3 >= 1) controls.Add((1, n - 3, "lambdaA"));   // the default margin-sharpness control

        // The from-above member witness for the sparse deferred path, computed ONCE and only on the
        // first cell that actually crosses the wall (null while nothing defers, so the default-wall
        // dense runs pay nothing).
        IReadOnlyDictionary<(int P, int W, string Shift), (double Bound, bool CarriedOdd)>? memberBounds = null;

        foreach (var (p, w) in BlockLattice.FundamentalDomain(n).OrderBy(b => BlockLattice.Dim(n, b.P, b.W)))
        {
            long dim = BlockLattice.Dim(n, p, w);
            foreach (var (kind, s) in shifts)
            {
                double margin = BlockLattice.WindowDistance(n, p, w, s.Real);
                bool isControl = margin > 0 && controls.Contains((p, w, kind));
                if (margin > 0 && !isControl)
                {
                    entries.Add(new ShellCensusEntry(p, w, dim, -1, -1, kind, s, false, "window",
                        double.NaN, double.NaN, double.NaN, true, margin, 0, 0, 0.0));
                    continue;
                }

                // parity dims without building anything big
                var perm = WeightCoherenceBlock.ReflectionPermutation(n, p, w);
                int fixedCount = 0;
                for (int i = 0; i < perm.Length; i++) if (perm[i] == i) fixedCount++;
                int pairs = (perm.Length - fixedCount) / 2;
                int dimEven = fixedCount + pairs, dimOdd = pairs;
                if (Math.Max(dimEven, dimOdd) > opts.MaxSectorDim)
                {
                    if (!opts.SparseForDeferred)
                    {
                        entries.Add(new ShellCensusEntry(p, w, dim, dimEven, dimOdd, kind, s, false, "deferred",
                            double.NaN, double.NaN, double.NaN, true, margin, 0, 0, 0.0));
                        opts.Log?.Invoke($"({p},{w}) x {kind}: DEFERRED, sector dim {Math.Max(dimEven, dimOdd)} > {opts.MaxSectorDim}");
                        continue;
                    }

                    // Past the LP64 wall, but the sparse instruments reach it (calibrated cell-for-cell
                    // against the dense census at N=9, SparseCensusCalibrationTests). The witness is built
                    // once, lazily, on the first deferral.
                    memberBounds ??= SectorWitnessTransport.MemberUpperBounds(seed, opts.Log);
                    var swSparse = Stopwatch.StartNew();

                    if (memberBounds.TryGetValue((p, w, kind), out var wb))
                    {
                        // MEMBER → the from-above W-transport witness fills the CARRIED-parity σ column and
                        // IS the classification value. NEVER Math.Min across the two columns [F4]:
                        // Math.Min(e-13, NaN) = NaN would flip a witnessed member to a false DISAGREE when
                        // the other-parity invit stalls. The other parity is secondary confirmation only.
                        double carried = wb.Bound, other = double.NaN;
                        long otherIters = 0;
                        var csrOther = WeightCoherenceSectorCsr.BuildReflectionSector(n, p, w, q, odd: !wb.CarriedOdd);
                        if (csrOther.Dim > 0)
                        {
                            var ro = SparseShiftedSigmaMin.Estimate(csrOther, s, opts.Sparse);
                            otherIters = ro.InnerIterations;
                            if (ro.Converged) other = ro.SigmaMin;
                        }
                        // the other-parity invit ran on the !CarriedOdd sector, so its iterations belong to
                        // that column (CarriedOdd → the EVEN column carries the invit; the carried column is 0)
                        entries.Add(new ShellCensusEntry(p, w, dim, dimEven, dimOdd, kind, s, true, "sparse-witness",
                            wb.CarriedOdd ? other : carried, wb.CarriedOdd ? carried : other, carried,
                            true, margin,
                            wb.CarriedOdd ? SafeInt(otherIters) : 0, wb.CarriedOdd ? 0 : SafeInt(otherIters),
                            swSparse.Elapsed.TotalSeconds));
                        opts.Log?.Invoke($"({p},{w}) x {kind}: sparse-witness bound={carried:E3} carriedOdd={wb.CarriedOdd} other={other:E3}");
                    }
                    else
                    {
                        // NON-MEMBER / control → sparse inverse-power σ_min per R-parity sector. Converged
                        // reflects the estimator honestly; a stalled sector keeps NaN σ (never a fabricated
                        // exclusion) and forces the verdict to PARTIAL via the [F1] unresolved rule.
                        double sigE = double.PositiveInfinity, sigO = double.PositiveInfinity;
                        int itE = 0, itO = 0;
                        bool sconv = true;
                        foreach (bool odd in new[] { false, true })
                        {
                            var csr = WeightCoherenceSectorCsr.BuildReflectionSector(n, p, w, q, odd);
                            if (csr.Dim == 0) continue;
                            var rr = SparseShiftedSigmaMin.Estimate(csr, s, opts.Sparse);
                            double v = rr.Converged ? rr.SigmaMin : double.NaN;
                            if (odd) { sigO = v; itO = SafeInt(rr.InnerIterations); } else { sigE = v; itE = SafeInt(rr.InnerIterations); }
                            sconv &= rr.Converged;
                        }
                        double ssig = Math.Min(sigE, sigO);   // NaN if either parity stalled → unresolved [F1]
                        entries.Add(new ShellCensusEntry(p, w, dim, dimEven, dimOdd, kind, s, true, "sparse-invit",
                            sigE, sigO, ssig, sconv, margin, itE, itO, swSparse.Elapsed.TotalSeconds));
                        opts.Log?.Invoke($"({p},{w}) x {kind}: sparse-invit sigma_min={ssig:E3} (even {sigE:E3} / odd {sigO:E3}) conv={sconv}");
                    }
                    continue;
                }

                var swBlock = Stopwatch.StartNew();
                double sigEven = double.PositiveInfinity, sigOdd = double.PositiveInfinity;
                int itEven = 0, itOdd = 0;
                bool conv = true;
                foreach (bool odd in new[] { false, true })    // strictly sequential: build, probe, release
                {
                    // THE convention pin: the FD block labeled (p,w) IS Build(n, p, w) = wKet:p, wBra:w
                    // (the L(1,2) naming of WeightCoherenceBlockTests; a consistent swap is unitary-silent).
                    var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, p, w, q, odd);
                    if (d == 0) continue;                      // no odd sector on all-palindromic blocks
                    for (int i = 0; i < d; i++) a[(long)i * d + i] -= s;
                    var r = ShiftedSigmaMin.EstimateColumnMajor(a, d);
                    if (odd) { sigOdd = r.SigmaMin; itOdd = r.Iterations; } else { sigEven = r.SigmaMin; itEven = r.Iterations; }
                    conv &= r.Converged;
                }
                double sig = Math.Min(sigEven, sigOdd);
                entries.Add(new ShellCensusEntry(p, w, dim, dimEven, dimOdd, kind, s, true,
                    isControl ? "control" : "lu-rparity", sigEven, sigOdd, sig, conv, margin,
                    itEven, itOdd, swBlock.Elapsed.TotalSeconds));
                opts.Log?.Invoke($"({p},{w}) x {kind}: sigma_min={sig:E3} (even {sigEven:E3} / odd {sigOdd:E3}) margin={margin:F4} {swBlock.Elapsed.TotalSeconds:F1}s");
            }
        }
        return new ShellCensusResult(seed, qR, lambda, pairGap, memberTol, true, entries, sw.Elapsed, opts);
    }

    public static ShellCensusSummary Summarize(ShellCensusResult r)
    {
        var expected = ExpectedMembers(r.Seed.N);
        var probed = r.Entries.Where(e => e.Probed).ToList();
        var found = probed.Where(e => e.SigmaMin < r.MemberTol)
            .Select(e => (e.P, e.W, e.Shift)).ToHashSet();
        var deferredMembers = r.Entries
            .Where(e => e.Method == "deferred" && expected.Contains((e.P, e.W, e.Shift)))
            .Select(e => (e.P, e.W, e.Shift)).ToList();
        var expectedProbeable = expected.Where(m => !deferredMembers.Contains(m)).ToHashSet();
        // the gray band around the member cut — anything here is flagged, never silently classified
        var ambiguous = probed
            .Where(e => e.SigmaMin >= r.MemberTol / r.Opts.AmbiguousBandFactor
                     && e.SigmaMin <= r.MemberTol * r.Opts.AmbiguousBandFactor)
            .Select(e => (e.P, e.W, e.Shift)).ToList();
        var nonMemberSigmas = probed
            .Where(e => !expected.Contains((e.P, e.W, e.Shift)))
            .Select(e => e.SigmaMin).ToList();
        double worstNonMember = nonMemberSigmas.Count > 0 ? nonMemberSigmas.Min() : double.PositiveInfinity;

        // [F1] the exclusion half is the claim's promotion blocker: ANY probed cell whose estimator did
        // not converge is UNRESOLVED (its σ is NaN, never a fabricated exclusion), forcing PARTIAL at best
        // and named — converged=0 can NEVER contribute to a PASS. A sparse-witness member carries a
        // certified from-above bound and stays Converged=true, so it never poisons the verdict here.
        var unresolved = probed.Where(e => !e.Converged)
            .Select(e => (e.P, e.W, e.Shift)).ToList();
        // [F6] a PASS in which any found member rests on a witness reading is the witness-assisted twin of
        // the blind dense PASS, distinct in the ledger.
        bool witnessAssisted = probed.Any(e => e.Method == "sparse-witness" && found.Contains((e.P, e.W, e.Shift)));

        bool classificationClean = found.SetEquals(expectedProbeable) && ambiguous.Count == 0 && r.SeedUsable;
        string verdict = !classificationClean ? "DISAGREE"
            : deferredMembers.Count > 0 || unresolved.Count > 0 ? "PARTIAL"
            : witnessAssisted ? "PASS (witness-assisted)" : "PASS";
        return new ShellCensusSummary(verdict, expected, expectedProbeable, found,
            deferredMembers, ambiguous, unresolved, worstNonMember, r.MemberTol, r.PairGap);
    }
}
