using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F117 Pascal-Gram positivity, recomputed live (the C# witness for the Python verifier
/// <c>simulations/f87_pascal_gram_positivity.py</c>). For a non-bipartite windowed diagonal-cell
/// pair, M = A + γQ has first nonvanishing odd power-sum at m* = 2ℓ + d, and the #Q = d class
/// coefficient is the equal-leg-total Pascal-Gram sum of squares
///
///   P_{m*,d} = (m*/d)·Σ_{l⃗∈[N]^d} Σ_{k⃗∈[0..ℓ]^d} |U^{(l⃗)}_{k⃗}|² ,
///   U^{(l⃗)}_{k⃗} = Σ_{|α⃗|=ℓ, d parts} ∏_i C(α_i, k_i)·T^{(l⃗)}_{α⃗} ,
///   T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁} Z_{l₂}H^{α₂} ··· Z_{l_d}H^{α_d}) ,
///
/// manifestly ≥ 0. At least one class is positive ⟹ p_{m*}(γ) > 0 for every γ > 0 ⟹ hard at one γ
/// is hard at ALL γ (the R-sign residual closed). The witness recomputes P_{m*,d} from H for the
/// five canonical branch representatives (d = 1, 3, 5) and surfaces the positivity verdict + the
/// mod-4 selection rule (d ≡ m*−2 (mod 4), d ≤ m*−2ℓ; singleness DERIVED for deg ≤ 3).
///
/// <para>The live complement to the closed-form <see cref="WindowedConverseAllGammaClaim"/>; reuses
/// the same k-body chain builder (<see cref="PauliKBodyChainExtensions.ChainKBody"/>) and per-site
/// Z (<see cref="PauliString.SiteOp"/>) the girth ladder uses. Inspect via <c>--root pascalgram</c>.</para></summary>
public sealed class PascalGramPositivityWitness : IInspectable
{
    /// <summary>A canonical F117 branch representative: a windowed pair (k-body templates summed into
    /// H on an N-site chain), its odd-girth ℓ, the firing #Q class d, m* = 2ℓ+d, and the exact
    /// integer coefficient the Python verifier pins (the cross-language ground truth).</summary>
    public sealed record Case(string Name, string[] Templates, int N, int Ell, int D, int MStar, double Expected);

    /// <summary>The five branch representatives of <c>f87_pascal_gram_positivity.py</c> BLOCK 3
    /// (d = 1: IXXZ+XIXZ; d = 3: K3, FLUX, MULTIZ; d = 5: IIXY+ZXZY).</summary>
    public static readonly IReadOnlyList<Case> Cases = new[]
    {
        new Case("IXXZ+XIXZ", new[] { "IXXZ", "XIXZ" }, N: 5, Ell: 3, D: 1, MStar: 7, Expected: 573440),
        new Case("K3",        new[] { "XXZ", "XZX" },   N: 4, Ell: 3, D: 3, MStar: 9, Expected: 2064384),
        new Case("FLUX",      new[] { "IXY", "XIY" },   N: 4, Ell: 3, D: 3, MStar: 9, Expected: 589824),
        new Case("MULTIZ",    new[] { "XXZ", "ZZZ" },   N: 4, Ell: 1, D: 3, MStar: 5, Expected: 61440),
        new Case("IIXY+ZXZY", new[] { "IIXY", "ZXZY" }, N: 5, Ell: 3, D: 5, MStar: 11, Expected: 86507520),
    };

    private readonly double[] _coefficients;   // cached live recompute, one per case

    public PascalGramPositivityWitness()
    {
        _coefficients = Cases.Select(ComputeCoefficient).ToArray();
    }

    public string DisplayName => "F117 Pascal-Gram positivity (live: 5 canonical cases)";

    public string Summary
    {
        get
        {
            int ok = 0;
            for (int i = 0; i < Cases.Count; i++)
                if (Math.Abs(_coefficients[i] - Cases[i].Expected) <= 1e-6 * Cases[i].Expected) ok++;
            return $"every #Q class at m* is a sum of squares ⟹ p_{{m*}}(γ) > 0 ∀γ>0 (hard at all γ); " +
                   $"{ok}/{Cases.Count} cases reproduce the exact coefficient from the SOS form live";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            for (int i = 0; i < Cases.Count; i++)
            {
                var c = Cases[i];
                double got = _coefficients[i];
                int[] allowed = AllowedClasses(c.MStar, c.Ell);
                string single = allowed.Length == 1 ? "monomiality DERIVED" : $"allowed {{{string.Join(",", allowed)}}}, measured single";
                yield return new InspectableNode(
                    $"{c.Name} (d={c.D}, N={c.N}, ℓ={c.Ell}, m*={c.MStar})",
                    summary: $"P_{{{c.MStar},{c.D}}} = (m*/d)·Σ_l⃗ Σ_k⃗ |U|² = {got:0} " +
                             $"(exact {c.Expected:0}, Δ={Math.Abs(got - c.Expected):0.0e+00}); > 0 ⟹ hard at every γ. " +
                             $"selection rule: d≡m*−2 (mod 4), d≤m*−2ℓ ⟹ {single}.",
                    payload: new InspectablePayload.Real($"P_{{{c.MStar},{c.D}}}", got, "0"));
            }
            yield return new InspectableNode("Live witness ⟷ WindowedConverseAllGammaClaim (F117)",
                summary: "the SOS coefficient is recomputed from H at inspect time; the closed-form theorem and the " +
                         "exact CRT cross-check live in WindowedConverseAllGammaClaim + f87_pascal_gram_positivity.py.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>Live readout: the per-site SOS coefficient P_{m*,d} = (m*/d)·Σ|U|² for one case,
    /// recomputed from H (no stored number). The cross-language ground truth is <see cref="Case.Expected"/>.</summary>
    public double Coefficient(int caseIndex) => _coefficients[caseIndex];

    /// <summary>The Pascal-Gram value P_{m*,d} for a case, from the d-leg moments T^{(l⃗)}_{α⃗}
    /// (|α⃗| = ℓ) through the binomial Vandermonde U and the sum of squares. Static so tests can
    /// call it without constructing the witness.</summary>
    public static double ComputeCoefficient(Case c)
    {
        int n = c.N, ell = c.Ell, d = c.D;
        int dim = 1 << n;
        var H = BuildHamiltonian(c.Templates, n);

        // P[l][a] = Z_l · H^a for l ∈ [0,n), a ∈ [0,ℓ].
        var hPow = new ComplexMatrix[ell + 1];
        hPow[0] = ComplexMatrix.Build.DenseIdentity(dim);
        for (int a = 1; a <= ell; a++) hPow[a] = hPow[a - 1] * H;
        var p = new ComplexMatrix[n][];
        for (int l = 0; l < n; l++)
        {
            var zl = PauliString.SiteOp(n, l, PauliLetter.Z);
            p[l] = new ComplexMatrix[ell + 1];
            for (int a = 0; a <= ell; a++) p[l][a] = zl * hPow[a];
        }

        // compositions of ℓ into d parts (the equal-leg total), indexed for the Vandermonde columns.
        var comps = Compositions(ell, d);
        var cidx = new Dictionary<string, int>();
        for (int i = 0; i < comps.Count; i++) cidx[Key(comps[i])] = i;

        // all d-leg moments T^{(l⃗)}_{α⃗}: dict (l-tuple) -> vector over compositions.
        var moments = AllLegMoments(p, n, d, ell, comps, cidx);

        // Vandermonde rows: all k⃗ ∈ [0..ℓ]^d; W[k⃗, α⃗] = ∏_i C(α_i, k_i). total = Σ_l⃗ Σ_k⃗ |U_k⃗|².
        var ks = Tuples(ell + 1, d);
        double total = 0.0;
        foreach (var vec in moments.Values)
        {
            if (vec.All(z => z.Magnitude < 1e-12)) continue;
            foreach (var k in ks)
            {
                Complex u = Complex.Zero;
                for (int cc = 0; cc < comps.Count; cc++)
                {
                    double w = 1.0;
                    for (int i = 0; i < d; i++) w *= Binomial(comps[cc][i], k[i]);
                    if (w != 0.0) u += w * vec[cc];
                }
                total += u.Real * u.Real + u.Imaginary * u.Imaginary;
            }
        }
        return (double)c.MStar / d * total;
    }

    /// <summary>The allowed #Q classes at (m*, ℓ): odd d ≤ m*−2ℓ with d ≡ m*−2 (mod 4). For
    /// deg = m*−2ℓ ∈ {1,3} this is a single class (monomiality derived).</summary>
    public static int[] AllowedClasses(int mStar, int ell)
    {
        var result = new List<int>();
        for (int d = 1; d <= mStar - 2 * ell; d += 2)
            if ((mStar - d) % 4 == 2) result.Add(d);
        return result.ToArray();
    }

    // ---- the d-leg moment recursion (port of all_leg_moments) ----
    private static Dictionary<string, Complex[]> AllLegMoments(
        ComplexMatrix[][] p, int n, int d, int ell, List<int[]> comps, Dictionary<string, int> cidx)
    {
        var moments = new Dictionary<string, Complex[]>();
        if (d == 1)
        {
            for (int l = 0; l < n; l++)
                moments[Key(new[] { l })] = new[] { p[l][ell].Trace() };
            return moments;
        }

        void Rec(int depth, List<int> ls, List<int> used, ComplexMatrix m)
        {
            int usedSum = used.Sum();
            if (depth == d - 1)
            {
                int aLast = ell - usedSum;
                for (int l = 0; l < n; l++)
                {
                    Complex t = (m * p[l][aLast]).Trace();
                    var key = Key(ls.Append(l).ToArray());
                    if (!moments.TryGetValue(key, out var vec)) { vec = new Complex[comps.Count]; moments[key] = vec; }
                    vec[cidx[Key(used.Append(aLast).ToArray())]] = t;
                }
                return;
            }
            for (int l = 0; l < n; l++)
                for (int a = 0; a <= ell - usedSum; a++)
                    Rec(depth + 1, new List<int>(ls) { l }, new List<int>(used) { a }, m * p[l][a]);
        }

        for (int l = 0; l < n; l++)
            for (int a = 0; a <= ell; a++)
                Rec(1, new List<int> { l }, new List<int> { a }, p[l][a]);
        return moments;
    }

    private static ComplexMatrix BuildHamiltonian(string[] templates, int n)
    {
        var terms = templates.Select(t => new PauliTerm(ParseTemplate(t), Complex.One)).ToList();
        return terms.ChainKBody(n);
    }

    private static PauliLetter[] ParseTemplate(string s) => s.Select(ch => ch switch
    {
        'I' => PauliLetter.I,
        'X' => PauliLetter.X,
        'Y' => PauliLetter.Y,
        'Z' => PauliLetter.Z,
        _ => throw new ArgumentException($"bad Pauli letter '{ch}' in template '{s}'"),
    }).ToArray();

    private static List<int[]> Compositions(int total, int parts)
    {
        var result = new List<int[]>();
        void Rec(int remaining, int slot, int[] acc)
        {
            if (slot == parts - 1) { acc[slot] = remaining; result.Add((int[])acc.Clone()); return; }
            for (int first = 0; first <= remaining; first++) { acc[slot] = first; Rec(remaining - first, slot + 1, acc); }
        }
        Rec(total, 0, new int[parts]);
        return result;
    }

    private static List<int[]> Tuples(int baseSize, int count)
    {
        var result = new List<int[]> { new int[count] };
        for (int slot = 0; slot < count; slot++)
        {
            var next = new List<int[]>();
            foreach (var t in result)
                for (int v = 0; v < baseSize; v++) { var copy = (int[])t.Clone(); copy[slot] = v; next.Add(copy); }
            result = next;
        }
        return result;
    }

    private static double Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0.0;
        double c = 1.0;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    private static string Key(int[] t) => string.Join(",", t);
}
