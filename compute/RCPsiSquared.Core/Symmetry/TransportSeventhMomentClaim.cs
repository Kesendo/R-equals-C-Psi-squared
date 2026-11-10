using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using Bond = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The seventh centred moment of the generator under cooling and heating at second order in the rates, for
/// hopping bonds with ZZ on any graph (<c>docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md</c>, Theorem B;
/// <c>docs/ANALYTICAL_FORMULAS.md</c>, F137, the paragraph "Step (3) and the order after it have a proof of their own").
///
/// <para>The setting of <see cref="FluxAreaFifthMomentClaim"/>: N qubits, H = Σ_bonds J(X_aX_b + Y_aY_b) +
/// D(X_aY_b − Y_aX_b) + Δ·Z_aZ_b, cooling σ⁻ at γ↓_l and heating σ⁺ at γ↑_l on site l, no fields, no dephasing,
/// c = Σ_l (γ↓_l + γ↑_l), t_ab = J + iD the amplitude of the hop a → b. The part of Tr((L + c/2)⁷) of degree two in
/// the rates is</para>
///
/// <para>  −168·4^N · Σ_{l&lt;k} β_lk · (γ↓_l·γ↑_k − γ↓_k·γ↑_l),</para>
///
/// <para>336·4^N·β_lk the coefficient of s_l·δ_k with s_l = (γ↓_l + γ↑_l)/2 and δ_l = (γ↓_l − γ↑_l)/2, and</para>
///
/// <para>  β_lk = Σ_j Im(t_lk·t_kj·t_jl)·(W_lkj + Q_lkj)
///            − 15·Σ_{i,j} |t_li|²·Im(t_ik·t_kj·t_ji) − 15·Σ_{i,j} |t_kj|²·Im(t_lj·t_ji·t_il)
///            + 15·Σ Im(t_lk·t_kx·t_xy·t_yz·t_zl) − 15·Σ Im(t_lx·t_xk·t_ky·t_yz·t_zl)
///            + 5·Σ_{i,j} Im(t_lk·t_kj·t_ji·t_il)·(Δ_lk − Δ_li − Δ_ij − Δ_jk)
///            + 5·Σ_{i,j} Im(t_li·t_ik·t_kj·t_jl)·(Δ_lj + Δ_jk),</para>
///
/// <para>the inner sites distinct and off {l, k}, with W_lkj = 30·S − 10·(d_l + d_k) − 5·d_j − 2·|t_lk|², S the sum
/// of |t|² over the bonds and d_v the sum of |t_vw|² over all sites w, and Q_lkj a quadratic form in the ZZ
/// couplings: 15·Δ² of a ZZ bond off the triangle, 5·Δ² of one with one end at l or k and the other off it, 10·Δ² of
/// one with one end at j and the other off it, 5·(Δ_lj² + Δ_jk²) + (10⁄3)·Δ_lk² + 2·Δ_lj·Δ_jk − (8⁄3)·Δ_lk·(Δ_lj + Δ_jk).
/// By lines: the triangle (l, k, j) by the imaginary part of its loop times a weight; a site joined to a corner of a
/// triangle reading that triangle's loop in the corner's place; the rings of five through l and k; the rings of
/// four, which enter only through ZZ. β_lk = −β_kl. In the integer form the witness reads, with s = γ↓ + γ↑ and
/// a = γ↓ − γ↑: the part of degree two of Tr((2L + c)⁷) is 3584·4^N·Σ_{l≠k} s_l·a_k·(3·β_lk).</para>
///
/// <para>Derived in the proof: the coefficient of s_l·δ_k in Tr e^{τ(L + c/2)} is τ·∫₀^τ (|f⁺_lk|² − |f⁻_lk|²) dx for
/// every H that conserves the excitation number (its Theorem A, step (3) of F137), an excitation carried from l to k
/// against a hole; at the seventh power the amplitudes are read on a circle of time, each site of a path excited
/// from one hop point to the next, the ZZ couplings as phases averaged over the spins that look on.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/TransportSeventhMomentWitness.cs</c>,
/// <c>inspect --root seventhmoment</c>, the part of degree two as 7·Σ_{p ≤ 2} Tr(B·A^p·B·A^(5−p)) with 2L + c = A + B, A the
/// generator at zero rates and B the part linear in the rates, in exact Gaussian integers at N = 3 and 4, beside the proof's
/// gate <c>simulations/seventh_moment_transport.py</c> (traces modulo ten primes, the degrees split by scaling the
/// rates, N up to 5, where the rings of five and a ZZ bond off a triangle are read). F137 itself is not typed, so this claim has no typed parent.
/// Adjacent, not parents: <see cref="FluxAreaFifthMomentClaim"/> (the order below, the direct hop against the detour
/// over one site, Δ not entering), <see cref="CooledBondPowerSumsClaim"/> (the seventh under cooling alone, where
/// every rate point lies on one axis and this part is zero on any graph; its seventh lives on bonds that do not
/// conserve the excitation number).</para></summary>
public sealed class TransportSeventhMomentClaim : Claim
{
    /// <summary>3·β_lk, an integer: three times the bracket of the closed form for the ordered pair (l, k).</summary>
    public static BigInteger ThriceBeta(int n, IReadOnlyList<Bond> bonds, int l, int k)
    {
        if (l < 0 || l >= n || k < 0 || k >= n || l == k) throw new ArgumentException($"({l}, {k}) is not an ordered pair of sites on {n}");
        var (t, zz) = Couplings(n, bonds);
        return ThriceBeta(n, t, zz, l, k);
    }

    /// <summary>The part of degree two in the rates of Tr((2L + c)⁷) = 128·Tr((L + c/2)⁷), by the closed form:
    /// 3584·4^N·Σ_{l≠k} (γ↓_l + γ↑_l)·(γ↓_k − γ↑_k)·(3·β_lk).</summary>
    public static BigInteger SeventhDegreeTwo(int n, IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        if (cool is null) throw new ArgumentNullException(nameof(cool));
        if (heat is null) throw new ArgumentNullException(nameof(heat));
        if (cool.Count != n) throw new ArgumentException($"{cool.Count} cooling rates for {n} sites", nameof(cool));
        if (heat.Count != n) throw new ArgumentException($"{heat.Count} heating rates for {n} sites", nameof(heat));
        var (t, zz) = Couplings(n, bonds);
        var total = BigInteger.Zero;
        for (int l = 0; l < n; l++)
            for (int k = 0; k < n; k++)
            {
                if (l == k) continue;
                BigInteger s = (BigInteger)cool[l] + heat[l], a = (BigInteger)cool[k] - heat[k];
                if (s.IsZero || a.IsZero) continue;
                total += s * a * ThriceBeta(n, t, zz, l, k);
            }
        return 3584 * BigInteger.Pow(4, n) * total;
    }

    private static (GaussianInteger[,] T, BigInteger[,] ZZ) Couplings(int n, IReadOnlyList<Bond> bonds)
    {
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(n, bonds);
        var zz = new BigInteger[n, n];
        foreach (var b in bonds)
        {
            zz[b.A, b.B] += b.Delta;
            zz[b.B, b.A] += b.Delta;
        }
        return (t, zz);
    }

    private static BigInteger Sq(GaussianInteger z) => z.Re * z.Re + z.Im * z.Im;

    private static BigInteger ThriceBeta(int n, GaussianInteger[,] t, BigInteger[,] d, int l, int k)
    {
        var bondSum = BigInteger.Zero;
        var deg = new BigInteger[n];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < n; y++)
            {
                deg[x] += Sq(t[x, y]);
                if (x < y) bondSum += Sq(t[x, y]);
            }
        var rest = Enumerable.Range(0, n).Where(v => v != l && v != k).ToArray();
        var b = BigInteger.Zero;
        foreach (int j in rest)
        {
            var loop = (t[l, k] * t[k, j] * t[j, l]).Im;
            if (loop.IsZero) continue;
            var w = 30 * bondSum - 10 * (deg[l] + deg[k]) - 5 * deg[j] - 2 * Sq(t[l, k]);
            var q = BigInteger.Zero;              // three times Q_lkj
            for (int x = 0; x < n; x++)
                for (int y = x + 1; y < n; y++)
                {
                    if (d[x, y].IsZero) continue;
                    bool xIn = x == l || x == k || x == j, yIn = y == l || y == k || y == j;
                    int weight;
                    if ((x == l && y == k) || (x == k && y == l)) weight = 10;
                    else if (xIn && yIn) weight = 15;                          // the two sides through j
                    else if (!xIn && !yIn) weight = 45;                        // off the triangle
                    else if (x == j || y == j) weight = 30;                    // one end at j, the other off it
                    else weight = 15;                                          // one end at l or k, the other off it
                    q += weight * d[x, y] * d[x, y];
                }
            q += 6 * d[l, j] * d[j, k] - 8 * d[l, k] * (d[l, j] + d[j, k]);
            b += loop * (3 * w + q);
        }
        foreach (int i in rest)
            foreach (int j in rest)
            {
                if (i == j) continue;
                b -= 45 * Sq(t[l, i]) * (t[i, k] * t[k, j] * t[j, i]).Im;                      // the copy at l
                b -= 45 * Sq(t[k, j]) * (t[l, j] * t[j, i] * t[i, l]).Im;                      // the copy at k
                var near = (t[l, k] * t[k, j] * t[j, i] * t[i, l]).Im;                         // the ring of four through the edge l–k
                b += 15 * near * (d[l, k] - d[l, i] - d[i, j] - d[j, k]);
                var far = (t[l, i] * t[i, k] * t[k, j] * t[j, l]).Im;                          // the ring of four with l and k opposite
                b += 15 * far * (d[l, j] + d[j, k]);
            }
        foreach (int x in rest)
            foreach (int y in rest)
                foreach (int z in rest)
                {
                    if (x == y || y == z || x == z) continue;
                    b += 45 * (t[l, k] * t[k, x] * t[x, y] * t[y, z] * t[z, l]).Im;            // the ring of five through the edge l–k
                    b -= 45 * (t[l, x] * t[x, k] * t[k, y] * t[y, z] * t[z, l]).Im;            // the ring of five with one site between l and k
                }
        return b;
    }

    public TransportSeventhMomentClaim()
        : base("Under cooling and heating alone, for bonds J(XX + YY) + D(XY − YX) + Δ·ZZ on any graph without fields, the part of " +
               "Tr((L + c/2)⁷) of degree two in the rates is −168·4^N·Σ_{l<k} β_lk·(γ↓_l·γ↑_k − γ↓_k·γ↑_l), β_lk in closed form on the " +
               "triangles through l and k, the triangles one bond away from l or k, the rings of five and, through ZZ alone, the rings of four",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md (Theorem B); docs/ANALYTICAL_FORMULAS.md (F137)")
    {
    }

    public override string DisplayName =>
        "Thermal seventh moment at second order in the rates: the transport formula on triangles and rings";

    public override string Summary =>
        "cooling γ↓_l and heating γ↑_l on every site, H hopping J + iD with ZZ on any graph, no fields: the part of Tr((L + c/2)⁷) of " +
        "degree two in the rates is −168·4^N·Σ_{l<k} β_lk·(γ↓_l·γ↑_k − γ↓_k·γ↑_l); each triangle through l and k by the imaginary part of " +
        "its loop times a weight (30, 25, 20, 15, 8 on the |t|², a quadratic form in Δ), a site joined to a corner of a triangle reading that triangle in the corner's place, rings of five, rings of four through ZZ alone " +
        "(Tier 1 derived; live witness: inspect --root seventhmoment).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the derivation",
                summary: "the coefficient of s_l·δ_k in Tr e^{τ(L + c/2)} is τ·∫₀^τ (|f⁺_lk|² − |f⁻_lk|²) dx for every H that conserves the excitation " +
                         "number (Theorem A of the proof, step (3) of F137); at order five in H the amplitudes are read on a circle of time with two marks, " +
                         "each site of a simple path excited from one hop point to the next, ZZ a phase averaged over the spins that look on; the fifth trace is the order " +
                         "below, the direct hop against the detour over one site");
            yield return new InspectableNode("by lines",
                summary: "a triangle (l, k, j): Im(t_lk·t_kj·t_jl) times 30·|t|² of a bond off its corners, 25 of one at j alone, 20 of one at l or k alone, " +
                         "15 on the sides through j, 8 on the edge l–k, plus the quadratic form in Δ; a site l joined to a corner i of a triangle (i, k, j): that " +
                         "triangle's loop in i's place, −15·|t_li|²; a ring of five: ±15 its loop, on a ring alone 2520·4^N·Im(loop)·(pentagram − pentagon) " +
                         "of the five rate points; a ring of four: 5·Im(loop) times a linear form in the Δ on its edges, zero without ZZ");
            yield return new InspectableNode("the instances",
                summary: "gate M11's two breaks at N = 4 with a fifth trace of zero break here, all of degree two: the unit triangle at ratio ¼ with a " +
                         "fourth site at ratio 1⁄9, Tr((L + c/2)⁷) = −44800; the ring of four at Δ = 7⁄10, flux π/2 and four ratios, −20380416⁄25, through " +
                         "its ZZ alone. On a triangle alone with three equal |t| and one Δ the three weights are equal and the part is flux times the " +
                         "area of the rate triangle again");
            yield return new InspectableNode("what this does not say",
                summary: "the part of degree four is not derived: the 430080 of the temperature gradient on the unit triangle is all of degree four, and a " +
                         "triangle with unequal weights (hoppings 1, 1 + i, 2i; γ↓ = (3, 5, 2), γ↑ = γ↓ − 2) has a part of degree two of zero and " +
                         "322560 of degree four, so a part of degree two of zero decides nothing. With a field or a bond that does not conserve the " +
                         "excitation number the closed form does not hold. No dephasing, no other jump; the rings of five and a ZZ bond off a triangle are read by the gate at " +
                         "N = 5, not by the live witness");
        }
    }
}
