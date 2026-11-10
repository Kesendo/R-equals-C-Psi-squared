using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The fifth centred moment of the generator under cooling and heating, for hopping bonds with ZZ on any
/// graph: flux times the area of the rate triangle (<c>docs/ANALYTICAL_FORMULAS.md</c>, F137, the paragraph "The
/// value −656·sin Φ is one point of a closed form").
///
/// <para>The setting: N qubits, H = Σ_bonds J(X_aX_b + Y_aY_b) + D(X_aY_b − Y_aX_b) + Δ·Z_aZ_b, cooling
/// σ⁻ = |0⟩⟨1| at γ↓_l and heating σ⁺ at γ↑_l on site l, no fields, no dephasing, c = Σ_l (γ↓_l + γ↑_l). With
/// t_ab = J + iD the amplitude of the hop a → b (H carries 2·t_ab·σ⁺_bσ⁻_a) and t_ba its conjugate,</para>
///
/// <para>  Tr((L + c/2)⁵) = 60·4^N · Σ_triangles Im(t_uv·t_vw·t_wu) · det[(1, 1, 1), (γ↓_u, γ↓_v, γ↓_w), (γ↑_u, γ↑_v, γ↑_w)],</para>
///
/// <para>the sum over the triangles of the hopping graph, each walked in one direction for both factors: the product
/// of the three |t| times the sine of the flux, times twice the signed area of the triangle the three rate points
/// (γ↓, γ↑) span. Δ does not enter. The trace is zero on a triangle at flux 0 or π or where the three rate points
/// lie on one line (one ratio, cooling alone, a temperature gradient at one coupling). In the integer form the
/// witness reads: Tr((2L + c)⁵) = 1920·4^N·Σ_triangles Im(loop)·det.</para>
///
/// <para>Derived in F137: with s_l = (γ↓_l + γ↑_l)/2 and δ_l = (γ↓_l − γ↑_l)/2 three maps (ρ ↦ Z⊗N·ρ with the
/// generator's transpose, the spin flip, ρ ↦ ρᵀ) leave, at the fifth power, the third order in H with one s_l and
/// one δ_k; its coefficient is a second Duhamel order, τ·∫(|f⁺_lk|² − |f⁻_lk|²), an excitation carried between the
/// two sites against a hole, and at lowest order the direct hop interferes with the detour over a third site, ket
/// on one arc of the triangle and bra on the other.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/FluxAreaFifthMomentWitness.cs</c>,
/// <c>inspect --root fifthmoment</c>, the generator in exact Gaussian integers from the matrix units, beside the
/// registry's gate <c>simulations/fifth_moment_flux_area.py</c> (floats on integer inputs, by Kronecker products).
/// F137 itself is not typed, so this claim has no typed parent. Adjacent, not parents:
/// <see cref="ThermalThirdMomentClaim"/> (the third moment in the same setting, zero for these bonds),
/// <see cref="TransportSeventhMomentClaim"/> (the order above, the seventh trace at second order in the rates),
/// <see cref="CooledBondPowerSumsClaim"/> (the parity count of the rates under cooling alone, where this trace is
/// zero), <see cref="MomentTowerPumpChannelClaim"/> (the arrow I → Z at first order; here it enters once, beside one
/// height).</para></summary>
public sealed class FluxAreaFifthMomentClaim : Claim
{
    /// <summary>One bond: J·(X_aX_b + Y_aY_b) + D·(X_aY_b − Y_aX_b) + Delta·Z_aZ_b, integers. J + iD is the amplitude
    /// of the hop a → b.</summary>
    public readonly record struct Bond(int A, int B, long J, long D, long Delta);

    /// <summary>One triangle u &lt; v &lt; w of the hopping graph: Im(t_uv·t_vw·t_wu) and
    /// det[(1, 1, 1), γ↓, γ↑] with the columns in the order u, v, w.</summary>
    public readonly record struct Triangle(int U, int V, int W, BigInteger ImLoop, BigInteger RateDeterminant)
    {
        public BigInteger Product => ImLoop * RateDeterminant;
    }

    /// <summary>t[x, y], the amplitude of the hop x → y: J + iD on a bond stored as (x, y), J − iD on one stored as
    /// (y, x); repeated bonds on a pair add up.</summary>
    public static GaussianInteger[,] HopAmplitudes(int n, IReadOnlyList<Bond> bonds)
    {
        if (bonds is null) throw new ArgumentNullException(nameof(bonds));
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "at least one site");
        var t = new GaussianInteger[n, n];
        foreach (var b in bonds)
        {
            if (b.A < 0 || b.A >= n || b.B < 0 || b.B >= n || b.A == b.B)
                throw new ArgumentException($"({b.A}, {b.B}) is not a bond on {n} sites", nameof(bonds));
            t[b.A, b.B] += new GaussianInteger(b.J, b.D);
            t[b.B, b.A] += new GaussianInteger(b.J, -(BigInteger)b.D);
        }
        return t;
    }

    /// <summary>The triangles of the hopping graph with a nonzero loop, each once, walked u → v → w → u.</summary>
    public static IReadOnlyList<Triangle> Triangles(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        if (cool is null) throw new ArgumentNullException(nameof(cool));
        if (heat is null) throw new ArgumentNullException(nameof(heat));
        if (cool.Count != n) throw new ArgumentException($"{cool.Count} cooling rates for {n} sites", nameof(cool));
        if (heat.Count != n) throw new ArgumentException($"{heat.Count} heating rates for {n} sites", nameof(heat));
        var t = HopAmplitudes(n, bonds);
        var list = new List<Triangle>();
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++)
                for (int w = v + 1; w < n; w++)
                {
                    var loop = t[u, v] * t[v, w] * t[w, u];
                    if (loop == GaussianInteger.Zero) continue;
                    BigInteger du = cool[u], dv = cool[v], dw = cool[w], hu = heat[u], hv = heat[v], hw = heat[w];
                    var det = (dv * hw - dw * hv) - (du * hw - dw * hu) + (du * hv - dv * hu);
                    list.Add(new Triangle(u, v, w, loop.Im, det));
                }
        return list;
    }

    /// <summary>Tr((2L + c)⁵) = 32·Tr((L + c/2)⁵), the closed form: 1920·4^N·Σ_triangles Im(loop)·det.</summary>
    public static BigInteger FifthCentredTrace(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        var sum = BigInteger.Zero;
        foreach (var tri in Triangles(n, bonds, cool, heat)) sum += tri.Product;
        return 1920 * BigInteger.Pow(4, n) * sum;
    }

    /// <summary>H as Pauli strings, site 0 first.</summary>
    public static IReadOnlyList<ThermalThirdMomentClaim.Term> Terms(int n, IReadOnlyList<Bond> bonds)
    {
        if (bonds is null) throw new ArgumentNullException(nameof(bonds));
        var list = new List<ThermalThirdMomentClaim.Term>();
        string Pair(int a, char pa, int b, char pb)
        {
            var s = new string('I', n).ToCharArray();
            s[a] = pa;
            s[b] = pb;
            return new string(s);
        }
        foreach (var b in bonds)
        {
            if (b.A < 0 || b.A >= n || b.B < 0 || b.B >= n || b.A == b.B)
                throw new ArgumentException($"({b.A}, {b.B}) is not a bond on {n} sites", nameof(bonds));
            if (b.J != 0)
            {
                list.Add(new(Pair(b.A, 'X', b.B, 'X'), b.J));
                list.Add(new(Pair(b.A, 'Y', b.B, 'Y'), b.J));
            }
            if (b.D != 0)
            {
                list.Add(new(Pair(b.A, 'X', b.B, 'Y'), b.D));
                list.Add(new(Pair(b.A, 'Y', b.B, 'X'), checked(-b.D)));
            }
            if (b.Delta != 0) list.Add(new(Pair(b.A, 'Z', b.B, 'Z'), b.Delta));
        }
        return list;
    }

    public FluxAreaFifthMomentClaim()
        : base("Under cooling and heating alone, for bonds J(XX + YY) + D(XY − YX) + Δ·ZZ on any graph without fields, " +
               "Tr((L + c/2)⁵) = 60·4^N·Σ_triangles Im(t_uv·t_vw·t_wu)·det[(1, 1, 1), γ↓, γ↑], t = J + iD the hop amplitude",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md (F137, the fifth centred trace as flux times the area of the rate triangle)")
    {
    }

    public override string DisplayName =>
        "Thermal fifth moment: flux times the area of the rate triangle";

    public override string Summary =>
        "cooling γ↓_l and heating γ↑_l on every site, H hopping J + iD with ZZ on any graph, no fields: " +
        "Tr((L + c/2)⁵) = 60·4^N·Σ_triangles Im(t_uv·t_vw·t_wu)·det[(1, 1, 1), γ↓, γ↑]; per triangle the sine of the flux times " +
        "twice the area the three rate points (γ↓, γ↑) span, Δ not entering (Tier 1 derived; live witness: inspect --root fifthmoment).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the derivation",
                summary: "L + c/2 = −i·ad_H + Σ_l (s_l·𝓗_l + δ_l·𝓐_l), 𝓗(ρ) = σ⁺ρσ⁻ + σ⁻ρσ⁺ the height and 𝓐(ρ) = σ⁻ρσ⁺ − σ⁺ρσ⁻ + ½{Z, ρ} the arrow; " +
                         "three maps make every trace even in the rates together and its degree in H plus its degree in δ even, one order of H with four rates " +
                         "has no trace, and the coefficient of s_l·δ_k at the third order in H is the second Duhamel order τ·∫(|f⁺_lk|² − |f⁻_lk|²): the direct " +
                         "hop against the detour over a third site, −120·4^N·Σ_j Im(t_lk·t_kj·t_jl) per ordered pair");
            yield return new InspectableNode("where it is silent",
                summary: "on a triangle at flux 0 or π, and where the three rate points lie on one line: one ratio (a line through the origin), cooling " +
                         "alone (the axis γ↑ = 0), a temperature gradient at one coupling (γ↓ − γ↑ the same on every site); the determinant does not change " +
                         "when all rate points are moved by one vector, and a bond hanging off a triangle does not enter");
            yield return new InspectableNode("the instances",
                summary: "F137's triangle with γ↓ = (½, 1, 3⁄2), γ↑ = (⅛, ⅓, ⅕): determinant −41⁄240, Tr((L + c/2)⁵) = −656·sin Φ; gate M11's control " +
                         "with hoppings 1, 2, 3i: −3936; M11's triangle at one ratio with a pendant site at another: zero, its break sitting at a " +
                         "higher trace; no ring of four has a triangle");
            yield return new InspectableNode("what this does not say",
                summary: "a vanishing fifth trace decides nothing (the temperature gradient at one coupling breaks at the seventh, in its part " +
                         "of degree four in the rates, SeventhMomentDegreeFourClaim; M11's pendant setting and its ring of four in the part of degree two, TransportSeventhMomentClaim); with a field, a " +
                         "uniform one included, or with a bond that does not conserve the excitation number the closed form fails (the witness's " +
                         "controls). No dephasing, no other jump");
        }
    }
}
