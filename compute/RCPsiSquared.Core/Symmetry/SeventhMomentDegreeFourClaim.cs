using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using Bond = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The seventh centred moment of the generator under cooling and heating at fourth order in the rates, a sum
/// over the triangles of the hopping graph (<c>docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md</c>, Theorems 1
/// and 2; <c>docs/ANALYTICAL_FORMULAS.md</c>, F137, the paragraph that closes the part of degree four).
///
/// <para>The setting of <see cref="FluxAreaFifthMomentClaim"/>: N qubits, H = Σ_bonds J(X_aX_b + Y_aY_b) +
/// D(X_aY_b − Y_aX_b) + Δ·Z_aZ_b, cooling σ⁻ at γ↓_l and heating σ⁺ at γ↑_l on site l, no fields, no dephasing,
/// c = Σ_l (γ↓_l + γ↑_l), t_ab = J + iD the amplitude of the hop a → b, s_l = (γ↓_l + γ↑_l)/2 the height of site l,
/// ω_ab = γ↓_a·γ↑_b − γ↓_b·γ↑_a. The part of Tr((L + c/2)⁷) of degree four in the rates is</para>
///
/// <para>  Σ_triangles 4^N·Im(t_uv·t_vw·t_wu)·Σ_cyc ω_{l,l+1}·(630·|s|² − 420·(s_l² + s_{l+1}²) − 210·s_{l+2}²),</para>
///
/// <para>each triangle walked u → v → w → u, the cyclic sum over its three pairs with l + 2 the opposite corner, |s|²
/// over every site of the graph; on the triangle alone, 13440·Im·Σ_cyc ω_{l,l+1}·(s_l² + s_{l+1}² + 2·s_{l+2}²).
/// With <see cref="TransportSeventhMomentClaim"/>'s part of degree two this is the whole seventh trace, the parts of
/// odd degree and of degree zero and six being zero (F137's three maps and its argument at first order in H). In the integer form the
/// witness reads, with S_l = γ↓_l + γ↑_l: the part of degree four of Tr((2L + c)⁷) is
/// 6720·4^N·Σ_triangles Im·Σ_cyc ω_{l,l+1}·(3·|S|² − 2·(S_l² + S_{l+1}²) − S_{l+2}²).</para>
///
/// <para>On the triangle, on a line of rate points γ↑ = e + m·γ↓ the whole Tr((L + c/2)⁷) is
/// −e·[3360·Im·(1 + m)²·V(γ↓) + 10752·⟨γ↓, b⟩], V(x) = (x₀ − x₁)(x₁ − x₂)(x₂ − x₀) and b_l = Σ_k β_lk the flow of
/// the parent's β; on a line γ↓ = e it is e·[3360·Im·V(γ↑) + 10752·⟨γ↑, b⟩] (Theorem 2). The shear
/// δ_l → δ_l + ν·s_l keeps every ω and every s, so it keeps this part, the part of degree two and the fifth trace, and
/// not the ninth, so not the spectrum.</para>
///
/// <para>Derived in the proof: at this order only the words whose three hops go once round a triangle survive, so the part is Im of the loop
/// times one polynomial in the rates, the same for every triangle, evaluated by an exact fit; as closed walks on the
/// Pauli letters (one arrow I → Z, three hops, the heights summed as h₃) the sites off a triangle add one constant to
/// every height, and over their letters that leaves the product of the fifth and second traces, nothing connected.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/SeventhMomentDegreeFourWitness.cs</c>,
/// <c>inspect --root seventhfour</c>: the part of degree four as the words with four rate letters,
/// 7·Σ over the five necklaces of Tr(A·B^g₁·A·B^g₂·A·B^g₃), with 2L + c = A + B, A the generator at zero rates and B the
/// part linear in the rates, in exact Gaussian integers at N = 3 and 4, beside the proof's gate
/// <c>simulations/seventh_moment_degree_four.py</c> (traces modulo primes, the degrees split by scaling the rates, a
/// fit over all 126 monomials, a walk sum, N up to 5). Parents: <see cref="FluxAreaFifthMomentClaim"/> (μ₅, the
/// product's other factor) and <see cref="TransportSeventhMomentClaim"/> (the part of degree two and the flow b).</para></summary>
public sealed class SeventhMomentDegreeFourClaim : Claim
{
    public FluxAreaFifthMomentClaim Fifth { get; }
    public TransportSeventhMomentClaim DegreeTwo { get; }

    /// <summary>The part of degree four in the rates of Tr((2L + c)⁷) = 128·Tr((L + c/2)⁷), by the closed form:
    /// 6720·4^N·Σ_triangles Im(loop)·Σ_cyc ω_{l,l+1}·(3·|S|² − 2·(S_l² + S_{l+1}²) − S_{l+2}²), S = γ↓ + γ↑.</summary>
    public static BigInteger SeventhDegreeFour(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        CheckRates(n, cool, heat);
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(n, bonds);
        var s = new BigInteger[n];
        var all = BigInteger.Zero;
        for (int x = 0; x < n; x++)
        {
            s[x] = (BigInteger)cool[x] + heat[x];
            all += s[x] * s[x];
        }
        BigInteger Omega(int a, int b) => (BigInteger)cool[a] * heat[b] - (BigInteger)cool[b] * heat[a];
        var total = BigInteger.Zero;
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++)
                for (int w = v + 1; w < n; w++)
                {
                    var im = (t[u, v] * t[v, w] * t[w, u]).Im;
                    if (im.IsZero) continue;
                    int[] c = { u, v, w };
                    var inner = BigInteger.Zero;
                    for (int i = 0; i < 3; i++)
                    {
                        int l = c[i], k = c[(i + 1) % 3], j = c[(i + 2) % 3];
                        inner += Omega(l, k) * (3 * all - 2 * (s[l] * s[l] + s[k] * s[k]) - s[j] * s[j]);
                    }
                    total += im * inner;
                }
        return 6720 * BigInteger.Pow(4, n) * total;
    }

    /// <summary>The whole Tr((2L + c)⁷) from the bonds and the rates alone: the parent's part of degree two plus the part
    /// of degree four.</summary>
    public static BigInteger SeventhCentredTrace(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat) =>
        TransportSeventhMomentClaim.SeventhDegreeTwo(n, bonds, cool, heat) + SeventhDegreeFour(n, bonds, cool, heat);

    internal static void CheckRates(int n, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        if (cool is null) throw new ArgumentNullException(nameof(cool));
        if (heat is null) throw new ArgumentNullException(nameof(heat));
        if (cool.Count != n) throw new ArgumentException($"{cool.Count} cooling rates for {n} sites", nameof(cool));
        if (heat.Count != n) throw new ArgumentException($"{heat.Count} heating rates for {n} sites", nameof(heat));
    }

    public SeventhMomentDegreeFourClaim(FluxAreaFifthMomentClaim fifth, TransportSeventhMomentClaim degreeTwo)
        : base("Under cooling and heating alone, for bonds J(XX + YY) + D(XY − YX) + Δ·ZZ on any graph without fields, the part of " +
               "Tr((L + c/2)⁷) of degree four in the rates is Σ_triangles 4^N·Im(loop)·Σ_cyc ω_{l,l+1}·(630·|s|² − 420·(s_l² + s_{l+1}²) − 210·s_{l+2}²), " +
               "13440·Im·Σ_cyc ω_{l,l+1}·(s_l² + s_{l+1}² + 2·s_{l+2}²) on the triangle; with the part of degree two the whole seventh trace",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md (Theorems 1 and 2); docs/ANALYTICAL_FORMULAS.md (F137)")
    {
        Fifth = fifth ?? throw new ArgumentNullException(nameof(fifth));
        DegreeTwo = degreeTwo ?? throw new ArgumentNullException(nameof(degreeTwo));
    }

    public override string DisplayName =>
        "Thermal seventh moment at fourth order in the rates: a sum over the triangles, the whole seventh trace closed";

    public override string Summary =>
        "cooling γ↓_l and heating γ↑_l on every site, H hopping J + iD with ZZ on any graph, no fields: the part of Tr((L + c/2)⁷) of degree " +
        "four is Σ_triangles 4^N·Im(loop)·Σ_cyc ω_{l,l+1}·(630·|s|² − 420·(s_l² + s_{l+1}²) − 210·s_{l+2}²), s the heights, ω the areas two rate " +
        "points span with the origin, twice; with the part of degree two the whole seventh trace from the bonds and the rates " +
        "(Tier 1 derived; live witness: inspect --root seventhfour).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the derivation",
                summary: "only the loop: a turn about Z per site commutes with loss and gain, so at three factors of H the part is Im(loop) times one " +
                         "polynomial in the rates, the same for every triangle (the flux-even products drop by F137's loss and gain exchanged), " +
                         "evaluated by an exact fit; as walks on the Pauli letters (the arrow I → Z once, three hops, the heights h₃ of the visited " +
                         "strings) the other sites add one constant w to every height, h₃(x + w) = h₃ + 6w·h₂ + 15w²·h₁ + 20w³, and over their " +
                         "letters only 21·(s_q²/2)·Tr₅ is left, the product, nothing connected");
            yield return new InspectableNode("the product and the touch rule",
                summary: "κ₇ = μ₇ − 21·μ₅·μ₂ with μ₂'s rate part |s|²/2 and μ₅ = 60·Im·Σω: the connected part is −26880·4^(N−3)·Im·Σ_cyc ω_{l,l+1}·(s_l² + " +
                         "s_{l+1}² + ½·s_{l+2}²); read beside Theorem B a squared height counts like minus the squared amplitude of a bond at that site " +
                         "(|2t|² or (2Δ)²), once at a pair's two sites and half at the opposite corner, except the pair's own edge, ZZ on the triangle's " +
                         "edges and, past the triangle, the borrowed triangles and the rings of four and five; a reading of two closed forms, not derived");
            yield return new InspectableNode("a line of rate points",
                summary: "on the triangle, on γ↑ = e + m·γ↓ the whole Tr((L + c/2)⁷) is −e·[3360·Im·(1 + m)²·V(γ↓) + 10752·⟨γ↓, b⟩], b the flow of β; " +
                         "on γ↓ = e it is e·[3360·Im·V(γ↑) + 10752·⟨γ↑, b⟩]; the shear δ → δ + ν·s takes a line to one arrow, the form of a temperature " +
                         "gradient as polynomials (the sheared rates need not be a bath's), and does not keep the ninth trace: " +
                         "F137's 430080 is 13440·γ⁴·V(n̄) on the unit triangle, the transport proof's 322560 is the line form at ⟨γ↓, b⟩ = 0");
            yield return new InspectableNode("what this does not say",
                summary: "why the walk classes sum so simply, and why the words with three arrows cancel (the fit shows it); the touch rule as a " +
                         "derivation; with a field, dephasing or a bond that does not conserve the excitation number the whole seventh trace is not the two " +
                         "closed forms (a Z field on one site left the part of degree four at the sum over the triangles where it was read; whether always, " +
                         "and why, is not shown). " +
                         "Whether the spectrum pairs is the equal triangle's claim (EqualTrianglePairingClaim) and open on every other triangle");
        }
    }
}
