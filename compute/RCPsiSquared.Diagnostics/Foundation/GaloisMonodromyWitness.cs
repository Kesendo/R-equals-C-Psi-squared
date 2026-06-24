using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Where the Galois structure of the F89 path-3 octic lives spectrally: the q-parametric
/// MONODROMY. The fixed-q geometry was a null (--root galoischaos); but as q=J/γ loops the complex
/// plane, the 8 octic roots braid, and that braiding IS the Galois group (monodromy = Galois over the
/// function field). This witness tracks the 8 roots (the 12×12 S_2-symmetric block A + q·C, with the
/// four closed-form AT-locked roots −2+iJ·{α,β}, −6+iJ·{α,β} removed) around a loop, via the validated
/// Monodromy tracker.
///
/// <para>Gate G2 (live): a loop around the diabolic point q_EP = √((−1+√13)/6) ≈ 0.659 returns the
/// IDENTITY permutation. The two octic roots that coalesce there do NOT braid, because the discriminant
/// has a DOUBLE zero (the (3q⁴+q²−1)² factor): a transversal crossing of two analytic sheets, not a
/// √-branch point. This confirms f89octic's "diabolic, semisimple, NOT defective" by an INDEPENDENT
/// route (monodromy, not the Riesz projector / departure-from-normality). The gateway to the
/// S_8-generation gate: the genuine branch points (the simple zeros of P_10) carry transpositions that
/// generate Gal(F_8) = S_8 from below. The sequel to --root galoischaos and the sibling of
/// --root f89octic.</para></summary>
public sealed class GaloisMonodromyWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // the diabolic point q_EP from the discriminant factor 3q⁴+q²−1 = 0 (the squared factor's real root).
    public static readonly double QEp = Math.Sqrt((-1 + Math.Sqrt(13)) / 6);

    // AT-locked single-particle frequencies (the F_a/F_b roots: λ = −2γ + iJ·α etc.).
    private static readonly double Alpha = -1 + Math.Sqrt(5);
    private static readonly double Beta = -1 - Math.Sqrt(5);

    // The ×2-cleared 12×12 S_2-symmetric (SE,DE) block is linear in q: 2M(q) = A + q·C, eigenvalues 2λ_k.
    // A and C are extracted ONCE from the trusted integer builder at q0 = 0 (diagonal only) and q0 = 1.
    private static readonly (Complex[,] A, Complex[,] C) Lin = BuildLinear();

    private static Complex[,] ToComplex(GaussianInteger[,] g)
    {
        int d = g.GetLength(0);
        var c = new Complex[d, d];
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
                c[r, s] = new Complex((double)g[r, s].Re, (double)g[r, s].Im);
        return c;
    }

    private static (Complex[,], Complex[,]) BuildLinear()
    {
        var a = ToComplex(F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 0, nBlock: 4));
        var m1 = ToComplex(F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 1, nBlock: 4));
        int d = a.GetLength(0);
        var c = new Complex[d, d];
        for (int r = 0; r < d; r++)
            for (int s = 0; s < d; s++)
                c[r, s] = m1[r, s] - a[r, s];
        return (a, c);
    }

    /// <summary>The 8 octic roots λ_k(q) at an arbitrary complex q: eigenvalues of (A + q·C)/2 with the
    /// four closed-form AT-locked roots removed (nearest-match in the 2λ spectrum).</summary>
    public static Complex[] OcticRootsAt(Complex q)
    {
        int d = Lin.A.GetLength(0);                                   // 12
        var m = Matrix<Complex>.Build.Dense(d, d, (r, s) => Lin.A[r, s] + q * Lin.C[r, s]);
        var twoLambda = m.Evd().EigenValues.ToArray();               // 12 values = 2λ
        var vals = new List<Complex>(twoLambda);

        var i2q = 2 * Complex.ImaginaryOne * q;
        var atTargets = new[]                                         // AT-locked roots in the 2λ scale
        {
            new Complex(-4, 0) + i2q * Alpha,  new Complex(-4, 0) + i2q * Beta,
            new Complex(-12, 0) + i2q * Alpha, new Complex(-12, 0) + i2q * Beta,
        };
        foreach (var t in atTargets)
        {
            int idx = 0; double bd = double.PositiveInfinity;
            for (int i = 0; i < vals.Count; i++)
            {
                double dd = (vals[i] - t).Magnitude;
                if (dd < bd) { bd = dd; idx = i; }
            }
            vals.RemoveAt(idx);
        }
        return vals.Select(v => v / 2).ToArray();                    // 8 octic λ
    }

    /// <summary>All 12 roots λ_k(q) of the S_2-symmetric block at complex q (no AT removal). The 4
    /// AT-locked roots are closed-form single-valued (no branch ⟹ never braid); tracking all 12 avoids
    /// the per-q removal mis-assignment that corrupts the set when an octic root grazes the AT lines on
    /// the complex loop. The octic monodromy is the restriction to the 8 octic strands.</summary>
    public static Complex[] AllRootsAt(Complex q)
    {
        int d = Lin.A.GetLength(0);
        var m = Matrix<Complex>.Build.Dense(d, d, (r, s) => Lin.A[r, s] + q * Lin.C[r, s]);
        return m.Evd().EigenValues.Select(v => v / 2).ToArray();
    }

    /// <summary>G2: the monodromy permutation of all 12 roots around a loop centred on q_EP. The 4 AT
    /// roots stay fixed (single-valued); the octic strands carry the verdict.</summary>
    public static int[] MonodromyAroundDiabolic(double radius, int steps)
        => Monodromy.Permutation(AllRootsAt, new Complex(QEp, 0), radius, steps);

    /// <summary>G1 (grounding): the q on the real axis where two octic roots come closest, and their
    /// midpoint. The minimum is the diabolic point q_EP, λ_EP = −4γ + 2iJ.</summary>
    public static (double qMin, Complex lamMid, double gap) ClosestApproachOnRealAxis()
    {
        double bestQ = 0, bestGap = double.PositiveInfinity;
        Complex bestMid = Complex.Zero;
        const int n = 600;
        for (int s = 0; s <= n; s++)
        {
            double q = 0.3 + (1.3 - 0.3) * s / n;
            var r = OcticRootsAt(new Complex(q, 0));
            for (int i = 0; i < r.Length; i++)
                for (int j = i + 1; j < r.Length; j++)
                {
                    double g = (r[i] - r[j]).Magnitude;
                    if (g < bestGap) { bestGap = g; bestQ = q; bestMid = (r[i] + r[j]) / 2; }
                }
        }
        return (bestQ, bestMid, bestGap);
    }

    // ---- G3 exploration: map the octic branch points (where two octic roots collide in complex q) ----

    public readonly record struct BranchPoint(Complex Q, double Gap, int Moved);

    private static double OcticGap(Complex q)
    {
        var r = OcticRootsAt(q);
        double m = double.PositiveInfinity;
        for (int i = 0; i < r.Length; i++)
            for (int j = i + 1; j < r.Length; j++)
                m = Math.Min(m, (r[i] - r[j]).Magnitude);
        return m;
    }

    private static int Moved(int[] perm) => Enumerable.Range(0, perm.Length).Count(i => perm[i] != i);

    /// <summary>Find the genuine branch points (defective EPs) TOPOLOGICALLY: an EP is a √-branch (the
    /// gap ~ √|q−q*| is a cusp the gap-scan jumps over), but a small loop around any cell containing it
    /// returns a non-identity monodromy regardless of cusp sharpness. The diabolic point is INVISIBLE here
    /// (its loop is the identity), which is correct: it is not a branch point in the monodromy sense. Tiles
    /// a box, loops each cell, refines the firing cells by quartering. Returns the EPs with their local
    /// transposition (moved strands). Restricted to a box around q_EP for speed; EPs come in conjugate pairs
    /// (the box spans both halves) plus the spatial structure of P_10(q²).</summary>
    public static List<BranchPoint> FindBranchPoints(
        double reLo = 0.45, double reHi = 0.85, double imLo = -0.18, double imHi = 0.18, double cell = 0.02)
    {
        var found = new List<BranchPoint>();
        for (double re = reLo; re < reHi; re += cell)
            for (double im = imLo; im < imHi; im += cell)
            {
                var c = new Complex(re + cell / 2, im + cell / 2);
                int moved = Moved(Monodromy.Permutation(AllRootsAt, c, cell * 0.62, 80));
                if (moved == 0) continue;
                if (found.Any(b => (b.Q - c).Magnitude < 1.5 * cell)) continue;   // one cell per EP cluster
                var refined = QuarterRefine(c, cell);
                found.Add(new BranchPoint(refined, OcticGap(refined),
                    Moved(Monodromy.Permutation(AllRootsAt, refined, 0.6 * cell, 200))));
            }
        return found;
    }

    // box bisection: keep halving toward the quadrant whose enclosing loop still braids, until the box is
    // small. Converges geometrically to the EP so a tight lasso circle reliably encloses it.
    private static Complex QuarterRefine(Complex c, double cell)
    {
        double half = cell / 2;
        for (int it = 0; it < 9 && half > 3e-4; it++)
        {
            double q = half / 2;                                  // quadrant-centre offset
            Complex pick = c; bool found = false;
            foreach (var off in new[] { new Complex(q, q), new Complex(-q, q), new Complex(q, -q), new Complex(-q, -q) })
            {
                if (Moved(Monodromy.Permutation(AllRootsAt, c + off, half * 0.95, 100)) > 0) { pick = c + off; found = true; break; }
            }
            c = pick; half = q;
            if (!found) { /* EP sits near the centre, keep shrinking around c */ }
        }
        return c;
    }

    // ---- G3 gate: assemble the EPs' transpositions into the monodromy = Galois group ----

    // the 8 octic strand indices among the 12 roots at base point q0 (the other 4 are the AT-locked roots).
    private static int[] OcticIndices(Complex q0, Complex[] r0)
    {
        var at = new[]
        {
            new Complex(-2, 0) + Complex.ImaginaryOne * q0 * Alpha, new Complex(-2, 0) + Complex.ImaginaryOne * q0 * Beta,
            new Complex(-6, 0) + Complex.ImaginaryOne * q0 * Alpha, new Complex(-6, 0) + Complex.ImaginaryOne * q0 * Beta,
        };
        var atIdx = new HashSet<int>();
        foreach (var t in at)
        {
            int best = -1; double bd = double.PositiveInfinity;
            for (int i = 0; i < r0.Length; i++)
                if (!atIdx.Contains(i)) { double dd = (r0[i] - t).Magnitude; if (dd < bd) { bd = dd; best = i; } }
            atIdx.Add(best);
        }
        return Enumerable.Range(0, r0.Length).Where(i => !atIdx.Contains(i)).ToArray();
    }

    /// <summary>The G3 gate: lasso every EP from a common base point, read its transposition on the 8 octic
    /// strands in one consistent labelling, and assemble. Transpositions generate S_8 ⟺ their graph on the
    /// 8 strands is connected (one component). That is Gal(F_8) = S_8 reconstructed from below, purely by
    /// tracking eigenvalue braids: monodromy = Galois, made live.</summary>
    public static (bool connected, int components, int nEps, int nClean, List<(int a, int b)> edges) GeneratesS8()
    {
        var q0 = new Complex(2.0, 0);                             // REAL base: AT roots sit at Re −2/−6 exactly,
        var r0 = AllRootsAt(q0);                                  // octic roots strictly between ⟹ clean labelling
        var octic = OcticIndices(q0, r0);
        var pos = new Dictionary<int, int>();
        for (int p = 0; p < octic.Length; p++) pos[octic[p]] = p;  // strand index → 0..7

        var eps = FindBranchPoints(reLo: 0.10, reHi: 1.70, imLo: -0.55, imHi: 0.55, cell: 0.04);
        var edges = new List<(int, int)>();
        var uf = Enumerable.Range(0, 8).ToArray();
        int Find(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

        int clean = 0;
        foreach (var ep in eps)
        {
            var lasso = Monodromy.Lasso(q0, ep.Q, radius: 0.02);
            var perm = Monodromy.PermutationAlongPath(AllRootsAt, lasso);
            var moved = Enumerable.Range(0, perm.Length).Where(i => perm[i] != i).ToList();
            if (moved.Count != 2 || !moved.All(pos.ContainsKey)) continue;   // accept only clean octic transpositions
            clean++;
            int a = pos[moved[0]], b = pos[moved[1]];
            edges.Add((a, b));
            uf[Find(a)] = Find(b);
        }
        int comps = Enumerable.Range(0, 8).Select(Find).Distinct().Count();
        return (comps == 1, comps, eps.Count, clean, edges);
    }

    private static string Cycles(int[] perm)
    {
        var seen = new bool[perm.Length];
        var parts = new List<string>();
        for (int i = 0; i < perm.Length; i++)
        {
            if (seen[i]) continue;
            var cyc = new List<int>();
            int j = i;
            while (!seen[j]) { seen[j] = true; cyc.Add(j); j = perm[j]; }
            if (cyc.Count > 1) parts.Add("(" + string.Join(" ", cyc) + ")");
        }
        return parts.Count == 0 ? "identity" : string.Join("", parts);
    }

    // ---- IInspectable ----

    public string DisplayName => "Galois monodromy of the path-3 octic (live: the diabolic loop is identity)";

    public string Summary =>
        $"where the octic's Galois structure lives spectrally: the q-parametric monodromy. As q=J/γ loops " +
        $"the complex plane the 8 roots braid (monodromy = Galois). A loop around the diabolic point " +
        $"q_EP≈{QEp.ToString("0.000", Inv)} returns the identity (the coalescing pair does not braid: a " +
        "double discriminant zero, transversal crossing), confirming --root f89octic semisimple-not-defective " +
        "by an independent route. Built on the trusted 12×12 block + the validated Monodromy tracker.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (qMin, lamMid, gap) = ClosestApproachOnRealAxis();
            yield return new InspectableNode("G1 grounding: the two octic roots coalesce at q_EP",
                summary: $"on the real q axis the closest approach of any two octic roots is at " +
                         $"q={qMin.ToString("0.000", Inv)} (q_EP={QEp.ToString("0.000", Inv)}), gap={gap.ToString("0.0000", Inv)}, " +
                         $"midpoint λ={lamMid.Real.ToString("0.00", Inv)}{lamMid.Imaginary.ToString("+0.00;-0.00", Inv)}i " +
                         "(λ_EP = −4γ + 2iJ). The avoided-crossing structure of the octic in q.");

            var sweep = new[] { 0.1, 0.05, 0.02, 0.005 }
                .Select(r => $"r={r.ToString("0.###", Inv)}:{Cycles(MonodromyAroundDiabolic(r, 600))}");
            int[] diabolic = MonodromyAroundDiabolic(0.02, 600);
            bool id = diabolic.Select((p, i) => p == i).All(x => x);
            yield return new InspectableNode("G2: the diabolic loop monodromy is the identity",
                summary: $"radius sweep around q_EP (all 12 strands): {string.Join("  ", sweep)} " +
                         $"({(id ? "identity at the isolating radius ✓" : "NOT identity — investigate")}). " +
                         "The diabolic point is a double discriminant zero (transversal crossing of two analytic sheets), " +
                         "so the two coalescing roots do NOT swap: f89octic's semisimple-not-defective verdict reached by " +
                         "monodromy, independent of the Riesz-projector route. The wider r=0.1 loop additionally encloses " +
                         "genuine simple branch points (P_10 zeros = defective EPs) flanking q_EP, each a transposition: " +
                         "the silent diabolic vs the braiding EPs, side by side, and the seed of the S_8 gate.");

            var g3 = GeneratesS8();
            string edgeList = string.Join(" ", g3.edges.Select(e => $"({e.a} {e.b})"));
            yield return new InspectableNode(
                $"G3: monodromy = Galois, S_8 from below ({(g3.connected ? "CONNECTED ✓" : $"{g3.components} components")})",
                summary: $"every EP lassoed from a common base, its transposition read on the 8 octic strands in one " +
                         $"labelling: {g3.nEps} branch points found, {g3.nClean} clean octic transpositions {edgeList}. " +
                         $"The transposition graph on the 8 strands has {g3.components} component" +
                         $"{(g3.components == 1 ? "" : "s")} ⟹ {(g3.connected ? "CONNECTED, so the transpositions generate the " +
                         "FULL symmetric group: Gal(F_8) = S_8, reconstructed purely from eigenvalue braids (monodromy = " +
                         "Galois, from below, the independent route to the algebraic Frobenius certificate)" : "not yet connected " +
                         "(widen the EP search). Transpositions generate S_8 iff their graph is connected.")}. " +
                         "The diabolic q_EP contributes nothing (silent). Complete-graph contrast (solvable ⟹ small disconnected braids) is G4.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
