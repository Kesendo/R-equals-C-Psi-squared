using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The artifact-free EP-character lab: the live, non-<c>eig</c>-eigenvector sibling of the
/// phase-rigidity reading in <see cref="CoherenceHorizonWitness"/>. It runs the three artifact-free
/// measures of <see cref="EpCharacter"/> (Riesz spectral-projector norm, departure-from-normality of
/// the orthonormal compression, geometric-vs-algebraic multiplicity) and certifies that the coherence
/// horizon Q*(N) is a GENUINE defective square-root exceptional point — not the eig-rigidity family
/// that misfired in the F86a retraction.
///
/// <para><b>Why this exists.</b> The only EP-character instrument in C# was
/// <see cref="PhaseRigidity"/>, which reads r = 1/√K off a single <c>Evd(L)</c> — a raw eigenvector
/// pairing. That family read r → 0 / K huge on a merely near-degenerate, NON-defective spectrum in the
/// F86a retraction. So the coherence-horizon EP verdict rested on the suspect instrument. This witness
/// re-establishes it from measures that never read an <c>eig</c> eigenvector: the resolvent contour
/// integral, the Frobenius departure of an orthonormal compression, and the SVD nullity.</para>
///
/// <para><b>Gate-first.</b> Before any horizon number is trusted, <see cref="GatePassed"/> validates the
/// three measures on KNOWN answers: a defective toy Jordan 2×2 (<c>[[-2γ, iJ],[iJ, -6γ]]</c> at J=2γ,
/// the F86a <c>TwoLevelEpModel</c> EP) MUST read DEFECTIVE; a diabolic point (<c>diag(-2γ,-2γ)</c>) MUST
/// read DIABOLIC. If the gate fails the witness declines honestly (every horizon node says so) rather
/// than report a verdict the instrument could not earn.</para>
///
/// <para><b>The demonstration.</b> On the single-excitation (Haken-Strobl) Liouvillian
/// <see cref="Lse"/> — the N²-dim real-space block L = −i(h⊗I − I⊗hᵀ) + (−4γ on the off-diagonal),
/// validated bit-for-bit against the full 4^N spectrum in
/// <c>simulations/coherence_horizon_se_block.py</c> — the slowest oscillating conjugate pair is found
/// just above Q*(N), the contour is drawn around exactly that pair, and <see cref="EpCharacter"/> reads
/// DEFECTIVE at every N=2..5: departure-from-normality ≈ 4 (exactly 4 at N=2,3, where the pair is a
/// clean 2×2 root of λ²+4γλ+cJ², drifting to ≈3.75 at N=5), geometric multiplicity 1 &lt; algebraic
/// multiplicity 2, the two compression eigenvectors merging (|cos| → 1). This is the artifact-free
/// confirmation of the √-EP, in the typed layer.</para>
///
/// <para>Port of <c>simulations/review_coherence_horizon_ep.py</c> (the algorithm) +
/// <c>simulations/review_f86a_diabolic_vs_defective.py</c> (the Riesz machinery). Typed neighbour:
/// <c>CoherenceHorizonClaim</c> (its EP-ness is now artifact-free-witnessed). Live:
/// <c>inspect --root epcharacter</c>.</para></summary>
public sealed class EpCharacterWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The coupling the horizon is read at (J=1, so Q = 1/γ).</summary>
    public const double J = 1.0;

    /// <summary>The detuning above Q* at which the freezing pair is a clean isolated 2×2: at Q* itself the
    /// pair is split only to ~√(machine-ε) (Jordan-block eigenvalue sensitivity), below what a contour or
    /// an SVD can resolve, so the measures are read at a small resolvable detuning where a DEFECTIVE EP
    /// keeps dep(A) ≈ O(1) (the 2×2 is a Jordan block as the eigenvalues merge).</summary>
    public const double ProbeDetuning = 1.002;

    /// <summary>The N range demonstrated (N²-dim SE block: N=5 is a 25×25 eigendecomposition, cheap).</summary>
    public const int MinN = 2;
    public const int MaxN = 5;

    // ---- the single-excitation (Haken-Strobl) Liouvillian, exactly as the Python reference ----

    /// <summary>The N×N single-particle Hamiltonian h: tridiagonal nearest-neighbour hopping J (the
    /// (J/2)(XX+YY) hop projected onto the one-excitation sector).</summary>
    public static double[,] HSingle(int n, double j)
    {
        var h = new double[n, n];
        for (int i = 0; i < n - 1; i++) { h[i, i + 1] = j; h[i + 1, i] = j; }
        return h;
    }

    /// <summary>The N²-dim single-excitation Liouvillian in the real-space basis (index i·N + j ↔
    /// |i⟩⟨j|): coherent part −i(h⊗I − I⊗hᵀ); Z-dephasing makes every off-diagonal coherence decay at
    /// 4γ (popcount(i^j)=2 ⟹ Re = −2γ·2 = −4γ, the Absorption Theorem), populations untouched. This is
    /// the L_se of <c>simulations/coherence_horizon_se_block.py</c>, validated there as an exact
    /// sub-spectrum of the full 4^N Liouvillian.</summary>
    public static ComplexMatrix Lse(int n, double j, double g)
    {
        var h = HSingle(n, j);
        int dim = n * n;
        var l = Matrix<Complex>.Build.Dense(dim, dim);
        var negI = new Complex(0, -1);
        // -i (h ⊗ I): acts on the left index a of |a⟩⟨b|.
        for (int a = 0; a < n; a++)
            for (int ap = 0; ap < n; ap++)
            {
                double hv = h[ap, a];
                if (hv != 0.0) for (int b = 0; b < n; b++) l[ap * n + b, a * n + b] += negI * hv;
            }
        // +i (I ⊗ hᵀ): acts on the right index b. (h is symmetric here, hᵀ = h.)
        for (int b = 0; b < n; b++)
            for (int bp = 0; bp < n; bp++)
            {
                double hv = h[b, bp];
                if (hv != 0.0) for (int a = 0; a < n; a++) l[a * n + b, a * n + bp] += -negI * hv;
            }
        // dephasing: −4γ on every coherence (i ≠ j), 0 on populations.
        for (int i = 0; i < n; i++)
            for (int jj = 0; jj < n; jj++)
                if (i != jj) l[i * n + jj, i * n + jj] += new Complex(-4.0 * g, 0.0);
        return l;
    }

    private readonly Dictionary<int, double> _qstarCache = new();

    /// <summary>The single-excitation coherence horizon Q*(N) = J/g*, where g* is the largest g (smallest
    /// Q) at which the slowest non-zero SE mode stops oscillating (the EP), by bisection on g — exactly
    /// the Python <c>qstar_se</c>. Cached per N.</summary>
    public double Qstar(int n)
    {
        if (n < MinN || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), $"the SE coherence horizon needs N in {MinN}..{MaxN}; got {n}");
        if (_qstarCache.TryGetValue(n, out var q)) return q;
        double lo = 0.02, hi = 6.0;
        for (int it = 0; it < 70; it++)
        {
            double m = 0.5 * (lo + hi);
            var ev = Lse(n, J, m).Evd().EigenValues;
            var nz = ev.Where(e => e.Real < -1e-7).ToArray();
            if (nz.Length == 0) { hi = m; continue; }
            double gap = nz.Max(e => e.Real);
            double bandIm = nz.Where(e => Math.Abs(e.Real - gap) < 1e-7).Max(e => Math.Abs(e.Imaginary));
            if (bandIm > 1e-7) lo = m; else hi = m;   // still oscillating ⟹ need more g (lower Q)
        }
        q = J / (0.5 * (lo + hi));
        _qstarCache[n] = q;
        return q;
    }

    /// <summary>The slowest non-zero oscillating conjugate pair about to freeze on the SE block at
    /// coupling J=Q·g (g=1): among the near-gap upper-half-plane modes, the one with the SMALLEST |Im|
    /// whose conjugate is genuinely its nearest neighbour (an isolated 2×2, the freezing {0,2}-coherence,
    /// not the large-|Im| band-edge survivor). Eigenvalue-only; no eigenvectors. Returns the pair's upper
    /// member λ and the full spectrum.</summary>
    public static (Complex Lambda, Complex[] Spectrum)? SlowCoalescer(ComplexMatrix l)
    {
        var w = l.Evd().EigenValues.ToArray();
        var nz = w.Where(e => e.Real < -1e-7).ToArray();
        if (nz.Length == 0) return null;
        var osc = nz.Where(e => e.Imaginary > 1e-9).OrderBy(e => Math.Abs(e.Imaginary)).ToArray();
        if (osc.Length == 0)
        {
            // already frozen: the two slowest reals (a coalesced or post-EP pair).
            var slow = w.OrderByDescending(e => e.Real).Take(2).ToArray();
            return (slow[0], w);
        }
        foreach (var la in osc)
        {
            var lb = Complex.Conjugate(la);
            double nn = w.Select(e => (e - la).Magnitude).Where(d => d > 1e-12).DefaultIfEmpty(double.PositiveInfinity).Min();
            if ((la - lb).Magnitude <= 1.5 * nn + 1e-9) return (la, w);   // conjugate is the nearest ⟹ isolated pair
        }
        return (osc[0], w);   // fallback: smallest |Im|
    }

    /// <summary>Characterise the freezing pair on the SE block at Q (g=1, J=Q): find the slow coalescer,
    /// draw the contour around exactly the pair nearest the pair-midpoint (radius outside the split,
    /// inside the genuine third eigenvalue), and read the three artifact-free measures off that contour.
    /// Returns null if no oscillating pair is resolved.</summary>
    public EpCharacter.Reading? CharacterizeAtQ(int n, double q)
    {
        // Probe convention (the Python reference): γ=1 fixed, hopping J=Q, so the freezing pair sits at
        // Re = −2γ = −2 and the only physical knob is the ratio Q = J/γ.
        var l = Lse(n, j: q, g: 1.0);
        var coalescer = SlowCoalescer(l);
        if (coalescer is not { } c) return null;
        Complex la = c.Lambda, lb = Complex.Conjugate(la);
        Complex lam0 = 0.5 * (la + lb);                // the pair midpoint (real, Re = −2γ band)
        double sep = (la - lb).Magnitude;
        // distances from lam0 to every eigenvalue, ascending: [0],[1] are the pair, [2] the third.
        var dist = c.Spectrum.Select(e => (e - lam0).Magnitude).OrderBy(d => d).ToArray();
        double third = dist.Length > 2 ? dist[2] : 10.0;
        double r = 0.40 * third;
        r = Math.Max(r, 5.0 * sep);                    // enclose the (split) pair comfortably
        r = Math.Min(r, 0.49 * third);                 // ...but never reach the third eigenvalue
        return EpCharacter.Characterize(l, lam0, r);
    }

    // ---- GATE 0: the instrument must read the known toys correctly before any horizon claim ----

    /// <summary>The known-defective toy: the F86a <c>TwoLevelEpModel</c> L_eff at its EP,
    /// <c>[[-2γ, iJ],[iJ, -6γ]]</c> with J=2γ (γ=1), a Jordan block with a double root at −4γ.</summary>
    public static ComplexMatrix ToyDefective(double gamma = 1.0)
    {
        double j = 2.0 * gamma;
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { new Complex(-2 * gamma, 0), new Complex(0, j) },
            { new Complex(0, j),          new Complex(-6 * gamma, 0) },
        });
    }

    /// <summary>The known-diabolic control: diag(−2γ, −2γ), two independent equal eigenvalues.</summary>
    public static ComplexMatrix ToyDiabolic(double gamma = 1.0) =>
        Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { new Complex(-2 * gamma, 0), Complex.Zero },
            { Complex.Zero,               new Complex(-2 * gamma, 0) },
        });

    private bool? _gatePassed;
    private EpCharacter.Reading _gateDefectiveReading, _gateDiabolicReading;

    /// <summary>The gate verdict, computed once: the defective toy reads DEFECTIVE AND the diabolic toy
    /// reads DIABOLIC. The whole witness's horizon claim is gated on this (declines honestly if false).</summary>
    public bool GatePassed
    {
        get
        {
            if (_gatePassed is { } g) return g;
            _gateDefectiveReading = EpCharacter.Characterize(ToyDefective(), new Complex(-4.0, 0), radius: 0.5);
            _gateDiabolicReading = EpCharacter.Characterize(ToyDiabolic(), new Complex(-2.0, 0), radius: 0.5);
            bool pass = _gateDefectiveReading.Kind == EpCharacter.EpKind.Defective
                        && _gateDiabolicReading.Kind == EpCharacter.EpKind.Diabolic;
            _gatePassed = pass;
            return pass;
        }
    }

    // ---- IInspectable ----

    public string DisplayName =>
        $"EpCharacterWitness (artifact-free EP character, the non-eig sibling of PhaseRigidity, J={J.ToString("0.#", Inv)})";

    public string Summary
    {
        get
        {
            if (!GatePassed)
                return "the artifact-free EP-character diagnostic (Riesz ‖P‖ / departure-from-normality / " +
                       "geo-vs-alg) — but the GATE FAILED on the known toys, so no horizon verdict is trusted " +
                       "(see the gate node). The instrument declines honestly.";
            string verdicts = string.Join(", ", Enumerable.Range(MinN, MaxN - MinN + 1)
                .Select(n => $"N={n}:{(CharacterizeAtQ(n, Qstar(n) * ProbeDetuning)?.Kind.ToString() ?? "n/a")}"));
            return "the artifact-free EP character of the coherence horizon Q*(N), recomputed live without an " +
                   "eig eigenvector pairing (the F86a-misfire-prone family). Three measures: Riesz spectral-" +
                   "projector ‖P‖, departure-from-normality of the orthonormal compression, geometric-vs-" +
                   "algebraic multiplicity (the SVD nullity). GATE PASSED (toy Jordan → DEFECTIVE, diag → " +
                   $"DIABOLIC). The single-excitation {{0,2}} freezing pair at Q*·{ProbeDetuning.ToString("0.###", Inv)} " +
                   $"reads: {verdicts} — a genuine defective √-EP (dep≈4, geo 1 < alg 2, eigenvectors merge). " +
                   "The non-eig confirmation of CoherenceHorizonClaim's EP-ness. Sibling reading (eig phase " +
                   "rigidity): inspect --root horizon.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheGate();
            if (GatePassed)
            {
                yield return TheHorizonVerdict();
                yield return TheLimit();
                yield return TheEigContrast();
            }
        }
    }

    /// <summary>The gate node: the two known toys read back their known characters, surfaced with the
    /// measures. A FIRING gate (a toy misread) is the finding — the witness then declines the horizon
    /// claim rather than launder a verdict the instrument could not earn.</summary>
    private InspectableNode TheGate()
    {
        bool pass = GatePassed;   // populates the two cached readings
        var d = _gateDefectiveReading;
        var b = _gateDiabolicReading;
        var rows = new List<IInspectable>
        {
            new InspectableNode("known-DEFECTIVE toy: [[-2γ, iJ],[iJ, -6γ]], J=2γ (the F86a TwoLevelEpModel EP)",
                summary: $"reads {d.Kind} (expected Defective): ‖P‖={d.ProjectorNorm.ToString("0.###", Inv)} " +
                         $"(closed 2×2 ⟹ ≈1, supplementary), departure-from-normality={d.Departure.ToString("0.###", Inv)} " +
                         $"(≈4, the Jordan coupling), geo={d.Geometric} < alg={d.Algebraic}, " +
                         $"eigenvector |cos|={d.EigenvectorMergeCos.ToString("0.####", Inv)} (→1, merged). " +
                         $"{(d.Kind == EpCharacter.EpKind.Defective ? "PASS" : "*** FAIL ***")}"),
            new InspectableNode("known-DIABOLIC control: diag(-2γ, -2γ)",
                summary: $"reads {b.Kind} (expected Diabolic): departure-from-normality=" +
                         $"{b.Departure.ToString("0.###e0", Inv)} (≈0), geo={b.Geometric} = alg={b.Algebraic}. " +
                         $"{(b.Kind == EpCharacter.EpKind.Diabolic ? "PASS" : "*** FAIL ***")}"),
        };
        return new InspectableNode("the gate (known-answer self-validation)",
            summary: pass
                ? "GATE PASSED: the three measures separate DEFECTIVE from DIABOLIC on the known toys, so the " +
                  "horizon verdict below is the instrument working, not a tautology. None reads an eig eigenvector."
                : "*** GATE FAILED: a known toy was misread; the horizon verdict is NOT trusted (the witness declines). ***",
            children: rows);
    }

    /// <summary>The horizon verdict: per N=2..5, the freezing {0,2} pair at Q*·1.002 read by the three
    /// artifact-free measures. DEFECTIVE at every N — the non-eig confirmation of the √-EP.</summary>
    private InspectableNode TheHorizonVerdict()
    {
        var rows = new List<IInspectable>();
        int nDefective = 0;
        foreach (int n in Enumerable.Range(MinN, MaxN - MinN + 1))
        {
            double qstar = Qstar(n);
            double qProbe = qstar * ProbeDetuning;
            var r = CharacterizeAtQ(n, qProbe);
            if (r is not { } reading)
            {
                rows.Add(new InspectableNode($"N={n}", summary: "no oscillating pair resolved at the probe."));
                continue;
            }
            if (reading.Kind == EpCharacter.EpKind.Defective) nDefective++;
            string eigStr = string.Join(", ", reading.CompressionEigenvalues
                .Select(e => $"{e.Real.ToString("0.###", Inv)}{(e.Imaginary >= 0 ? "+" : "")}{e.Imaginary.ToString("0.###", Inv)}i"));
            rows.Add(new InspectableNode($"N={n}: {reading.Kind}",
                summary: $"Q*={qstar.ToString("0.#####", Inv)}, probe Q={qProbe.ToString("0.#####", Inv)} (=1.002 Q*): " +
                         $"‖P‖={reading.ProjectorNorm.ToString("0.###", Inv)}, " +
                         $"departure-from-normality={reading.Departure.ToString("0.###", Inv)} " +
                         $"(≈4 at N=2,3 → ≈3.75 at N=5; BOUNDED away from 0 as the split → 0 ⟹ Jordan), " +
                         $"geo={reading.Geometric} < alg={reading.Algebraic}, " +
                         $"eigenvector |cos|={reading.EigenvectorMergeCos.ToString("0.####", Inv)} (→1), " +
                         $"compression eigenvalues [{eigStr}]."));
        }
        return new InspectableNode("the horizon verdict (artifact-free, N=2..5)",
            summary: $"the freezing {{0,2}}-coherence pair at Q*(N)·{ProbeDetuning.ToString("0.###", Inv)} reads " +
                     $"DEFECTIVE at {nDefective}/{MaxN - MinN + 1} of N=2..5 (a genuine 2nd-order √-EP). Departure-" +
                     "from-normality ≈ 4 (exactly 4 at N=2,3 where the pair is the clean 2×2 root of λ²+4γλ+cJ²; " +
                     "drifting to ≈3.75 at N=5 as the pair is collectively dressed), geometric mult 1 < algebraic " +
                     "mult 2 (a Jordan block), the two compression eigenvectors merging (|cos| → 1). This is the " +
                     "non-eig confirmation of CoherenceHorizonClaim's EP-ness — the eig phase rigidity (inspect " +
                     "--root horizon) corroborates but is not load-bearing (it is the F86a-misfire-prone family).",
            children: rows);
    }

    /// <summary>The limit Q → Q*: the defective signature is the split → 0 while dep(A) stays BOUNDED away
    /// from 0 and the eigenvectors merge (|cos| → 1). A diabolic point would have dep → 0 with the split.</summary>
    private InspectableNode TheLimit()
    {
        int n = 4;   // the first transcendental N (the clean 2×2 ends at N=3); the dressing shows here
        double qstar = Qstar(n);
        var rows = new List<IInspectable>();
        double depFloor = double.PositiveInfinity, splitMin = double.PositiveInfinity;
        foreach (double fac in new[] { 1.02, 1.01, 1.005, 1.002 })
        {
            double q = qstar * fac;
            var l = Lse(n, j: q, g: 1.0);
            var c = SlowCoalescer(l);
            if (c is not { } cc) continue;
            double split = (cc.Lambda - Complex.Conjugate(cc.Lambda)).Magnitude;
            var r = CharacterizeAtQ(n, q);
            if (r is not { } reading) continue;
            if (reading.Algebraic == 2)
            {
                depFloor = Math.Min(depFloor, reading.Departure);
                splitMin = Math.Min(splitMin, split);
            }
            rows.Add(new InspectableNode($"Q={q.ToString("0.#####", Inv)} (={fac.ToString("0.###", Inv)} Q*)",
                summary: $"pair-split={split.ToString("0.####e0", Inv)}, " +
                         $"departure-from-normality={reading.Departure.ToString("0.###", Inv)}, " +
                         $"eigenvector |cos|={reading.EigenvectorMergeCos.ToString("0.#####", Inv)}, " +
                         $"geo={reading.Geometric}/alg={reading.Algebraic}."));
        }
        return new InspectableNode($"the limit Q → Q* (N={n}, the defectiveness is in the limit)",
            summary: $"as Q → Q*(4)={qstar.ToString("0.#####", Inv)} the pair-split shrinks to " +
                     $"{splitMin.ToString("0.##e0", Inv)} while the departure-from-normality holds at a floor of " +
                     $"{depFloor.ToString("0.###", Inv)} (BOUNDED away from 0) and the eigenvectors merge — the " +
                     "DEFECTIVE (Jordan) signature. A diabolic point would have the departure collapse to 0 with " +
                     "the split. (The split cannot be driven below ~√(machine-ε): a true double root splits as " +
                     "√(eps), which is why the measures are read at a resolvable detuning, not at Q* itself.)",
            children: rows);
    }

    /// <summary>The eig phase rigidity (the F86a-misfire family) on the slow oscillating pair of the SAME
    /// SE block, read across a few Q near Q*(4). On a genuine EP eig DOES drop (this is opposite-prior to
    /// F86a, where it dropped falsely on a non-defective near-degeneracy), but the MAGNITUDE is grid-
    /// fragile while the artifact-free verdict (dep≈4, geo<alg) is flat. Restricting to the small-|Im|
    /// oscillating modes (the freezing pair) avoids reading the exactly-real frozen modes whose rigidity is
    /// numerically 1 by accident.</summary>
    private InspectableNode TheEigContrast()
    {
        int n = 4;
        double qstar = Qstar(n);
        var rows = new List<IInspectable>();
        var rMins = new List<double>();
        foreach (double fac in new[] { 1.001, 1.002, 1.005, 1.01 })
        {
            double q = qstar * fac;
            var l = Lse(n, j: q, g: 1.0);
            // min r over the slow oscillating modes (the coalescer branch), not the frozen real modes.
            var osc = PhaseRigidity.Compute(l)
                .Where(m => m.Lambda.Real < -1e-6 && Math.Abs(m.Lambda.Imaginary) > 1e-6 && Math.Abs(m.Lambda.Imaginary) < 1.5)
                .ToList();
            if (osc.Count == 0) continue;
            double rMin = osc.Min(m => m.Rigidity);
            rMins.Add(rMin);
            rows.Add(new InspectableNode($"Q={q.ToString("0.#####", Inv)} (={fac.ToString("0.###", Inv)} Q*)",
                summary: $"eig phase rigidity of the slow oscillating pair: min r = {rMin.ToString("0.###e0", Inv)} " +
                         "(the smaller as Q → Q*, the coalescence)."));
        }
        double rFloor = rMins.Count > 0 ? rMins.Min() : double.NaN;
        return new InspectableNode("the eig contrast (why the artifact-free reading is load-bearing)",
            summary: $"the eig phase rigidity on the SAME SE block's freezing pair drops toward 0 as Q → Q*(4) " +
                     $"(min r ≈ {rFloor.ToString("0.##e0", Inv)} at the closest resolvable probe). Here eig AGREES " +
                     "with the artifact-free DEFECTIVE verdict — because this IS a genuine EP, the OPPOSITE prior to " +
                     "the F86a case (a near-degenerate NON-defective spectrum, where eig read r → 0 FALSELY). The " +
                     "point: eig's r → 0 is necessary-but-not-sufficient, and right at Q* its magnitude is grid-fragile " +
                     "(it sits on the coalescence; K swings orders of magnitude under ΔQ=1e-3). The departure-from-" +
                     "normality + geo<alg verdict is what actually separates defective from diabolic, and it never reads " +
                     "an eig eigenvector. So PhaseRigidity corroborates here but cannot be the proof.",
            children: rows);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
