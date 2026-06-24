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

            yield return new InspectableNode("the sequel: the S_8-generation gate",
                summary: "the genuine branch points are the simple zeros of P_10(q²) in disc(F_8) = const·q²⁴·" +
                         "(3q⁴+q²−1)²·P_10(q²); each carries a transposition, and they generate Gal(F_8) = S_8 from " +
                         "below. That is the monodromy = Galois gate (G3), the q-direction realisation of the writability " +
                         "result, and the complete-graph contrast (solvable ⟹ disconnected small braids) is G4.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
