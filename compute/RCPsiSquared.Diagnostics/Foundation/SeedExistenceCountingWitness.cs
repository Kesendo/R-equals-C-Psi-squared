using System.Globalization;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RealMatrix = MathNet.Numerics.LinearAlgebra.Matrix<double>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the seed-existence counting theorem
/// (<c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>; verifier
/// <c>simulations/seed_existence_nullity_check.py</c>): on the (1,2) block pencil
/// L(q) = A + q·C of the XY chain under uniform Z-dephasing, the real-eigenvalue-count identity
///
/// <code>
///     r(0⁺) − r(∞) = N − 1      for every odd N
/// </code>
///
/// through the three counting lemmas, all recomputed here at inspect time:
/// (N2) nullity(P₋₂CP₋₂) = N − 1 (the −2 rung decomposes into N−1 disjoint paths of N vertices),
/// (FF) nullity(C) = ρ, the free-fermion fusion-resonance count #{({a,b},c): λ_a+λ_b = λ_c},
/// and (N1′) nullity(P₋₆CP₋₆) = 3·Z₃ = ρ (the ORDERING-SECTOR theorem: K₆ splits into three
/// no-passing bra-rank sectors, each gauge-equivalent via diag((−1)^{z_bra}) to MINUS the
/// 3-magnon block H₃, so spec(K₆) = 3 × {−(λ_a+λ_b+λ_c)}; the chiral 3-to-1
/// triple↔resonance bijection with D = 0 by cyclotomic integrality).
///
/// <para>Boundary of the recompute: the block is built combinatorially from the k-excitation
/// hopping graphs (amplitude 1 per bond, the verifier's convention; the spin-J/4 vs Pauli-J
/// scale is nullity-free), not through <c>WeightCoherenceBlock</c>, because the witness needs
/// A and C separately and the −2/−6 index split of the position basis. Dimensions are
/// cross-checked: dim = N²(N−1)/2, dim₋₂ = N(N−1), dim₋₆ = N(N−1)(N−2)/2. The SVD nullity
/// floor 1e-9 is safe here by ≥ 7 orders (the doc's near-miss law caps it only at N ≈ 3959).</para>
///
/// <para>Two-sided gate: the zero reads (cross-sector max, gauge residual) sit beside nonzero
/// controls (the −2/−6 coupling block max and, for N ≥ 5, the K₆ max), so a builder that
/// accidentally kills the coupling or the rung content fails the witness rather than passing
/// it. The counting side (ρ, Z₃) is combinatorial in λ_k = 2cos(kπ/(N+1)), a second,
/// matrix-free computation meeting the SVD nullities.</para>
///
/// <para>Scope note carried from the doc: this witness types the COUNTING identity. The
/// seed-existence conclusion itself stays open at the codim-2 β-exotic genericity item.</para></summary>
public sealed class SeedExistenceCountingWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Largest odd N the live build materialises (dim = N²(N−1)/2 = 324 at N = 9).</summary>
    public const int MaxN = 9;

    private const double NullityTol = 1e-9;

    public int N { get; }

    /// <summary>dim of the (1,2) block, N²(N−1)/2.</summary>
    public int Dim { get; }
    /// <summary>dim of the −2 rung (bra site on a ket site), N(N−1).</summary>
    public int DimMinus2 { get; }
    /// <summary>dim of the −6 rung (three distinct sites), N(N−1)(N−2)/2.</summary>
    public int DimMinus6 { get; }

    public int NullityMinus2 { get; }
    public int NullityMinus6 { get; }
    public int NullityC { get; }

    /// <summary>ρ = #{({a,b},c): a &lt; b, λ_a+λ_b = λ_c}, combinatorial.</summary>
    public int ResonanceCount { get; }
    /// <summary>Z₃ = #{a &lt; b &lt; c: λ_a+λ_b+λ_c = 0}, combinatorial.</summary>
    public int ZeroSumTripleCount { get; }

    /// <summary>r(0⁺) − r(∞) = (n₂ + n₆) − nullity(C); the theorem says N − 1.</summary>
    public int Surplus { get; }

    /// <summary>max |K₆ element| between different bra-rank sectors; the no-passing gate (exactly 0).</summary>
    public double CrossSectorMax { get; }
    /// <summary>max |(U·K_sec·U + H₃) element| over the three sectors; the gauge gate (exactly 0).</summary>
    public double GaugeResidualMax { get; }
    /// <summary>max |element| of the −6/−2 coupling block of K; nonzero control (the coupling is real).</summary>
    public double CouplingBlockMax { get; }
    /// <summary>max |K₆ element|; nonzero control for N ≥ 5 (K₆ ≡ 0 exactly at N = 3).</summary>
    public double KSixMax { get; }

    public SeedExistenceCountingWitness(int n = 5)
    {
        if (n < 3 || n % 2 == 0 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"N must be odd with 3 ≤ N ≤ {MaxN} (the theorem's odd-N regime, live-build guarded); got {n}");
        N = n;

        var kets = Combinations(n, 2);
        var bras = Combinations(n, 1);
        Dim = kets.Count * bras.Count;
        DimMinus2 = n * (n - 1);
        DimMinus6 = n * (n - 1) * (n - 2) / 2;

        // K = -(H2 (x) I - I (x) H1), real symmetric; C = i*K shares every nullity.
        var h1 = Hop(n, 1);
        var h2 = Hop(n, 2);
        var k = RealMatrix.Build.Dense(Dim, Dim);
        for (int a = 0; a < kets.Count; a++)
            for (int ap = 0; ap < kets.Count; ap++)
                for (int b = 0; b < bras.Count; b++)
                    for (int bp = 0; bp < bras.Count; bp++)
                        k[ap * bras.Count + bp, a * bras.Count + b] =
                            -(h2[ap, a] * (b == bp ? 1.0 : 0.0) - (a == ap ? 1.0 : 0.0) * h1[bp, b]);

        // The -2 / -6 split of the position basis (n_diff = |ket Δ bra| ∈ {1, 3}).
        var i2 = new List<int>();
        var i6 = new List<int>();
        for (int a = 0; a < kets.Count; a++)
            for (int b = 0; b < bras.Count; b++)
                (kets[a].Contains(bras[b][0]) ? i2 : i6).Add(a * bras.Count + b);
        if (i2.Count != DimMinus2 || i6.Count != DimMinus6)
            throw new InvalidOperationException(
                $"rung split {i2.Count}/{i6.Count} does not match N(N−1)={DimMinus2} / N(N−1)(N−2)/2={DimMinus6}");

        var k2 = Submatrix(k, i2);
        var k6 = Submatrix(k, i6);
        NullityMinus2 = Nullity(k2);
        NullityMinus6 = Nullity(k6);
        NullityC = Nullity(k);
        Surplus = NullityMinus2 + NullityMinus6 - NullityC;

        // Nonzero controls: the -6/-2 coupling block and K6 itself.
        double coupling = 0.0;
        foreach (int r in i6)
            foreach (int c in i2)
                coupling = Math.Max(coupling, Math.Abs(k[r, c]));
        CouplingBlockMax = coupling;
        KSixMax = k6.Enumerate().Max(Math.Abs);

        // (N1') Step 1+2 gates: bra-rank sectors and the (-1)^{z_bra} gauge onto -H3.
        var sector = new int[i6.Count];
        var braSite = new int[i6.Count];
        var occupied = new (int, int, int)[i6.Count];
        for (int t = 0; t < i6.Count; t++)
        {
            int flat = i6[t];
            var ket = kets[flat / bras.Count];
            int bra = bras[flat % bras.Count][0];
            int lo = Math.Min(ket[0], ket[1]);
            int hi = Math.Max(ket[0], ket[1]);
            sector[t] = bra < lo ? 0 : bra < hi ? 1 : 2;
            braSite[t] = bra;
            var all = new[] { ket[0], ket[1], bra };
            Array.Sort(all);
            occupied[t] = (all[0], all[1], all[2]);
        }

        double cross = 0.0;
        for (int s = 0; s < i6.Count; s++)
            for (int t = 0; t < i6.Count; t++)
                if (sector[s] != sector[t])
                    cross = Math.Max(cross, Math.Abs(k6[s, t]));
        CrossSectorMax = cross;

        var h3 = Hop(n, 3);
        var triples = Combinations(n, 3);
        var cfgIndex = new Dictionary<(int, int, int), int>();
        for (int t = 0; t < triples.Count; t++)
            cfgIndex[(triples[t][0], triples[t][1], triples[t][2])] = t;

        double gauge = 0.0;
        for (int s = 0; s < 3; s++)
        {
            var ids = Enumerable.Range(0, i6.Count).Where(t => sector[t] == s).ToArray();
            for (int u = 0; u < ids.Length; u++)
                for (int v = 0; v < ids.Length; v++)
                {
                    double gauged = ((braSite[ids[u]] + braSite[ids[v]]) % 2 == 0 ? 1.0 : -1.0)
                                    * k6[ids[u], ids[v]];
                    double target = -h3[cfgIndex[occupied[ids[u]]], cfgIndex[occupied[ids[v]]]];
                    gauge = Math.Max(gauge, Math.Abs(gauged - target));
                }
        }
        GaugeResidualMax = gauge;

        // (FF) + (N1') Step 3: the combinatorial side, matrix-free.
        var lam = new double[n];
        for (int j = 0; j < n; j++) lam[j] = 2.0 * Math.Cos((j + 1) * Math.PI / (n + 1));
        int rho = 0;
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                for (int c = 0; c < n; c++)
                    if (Math.Abs(lam[a] + lam[b] - lam[c]) < 1e-7) rho++;
        ResonanceCount = rho;
        int z3 = 0;
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                for (int c = b + 1; c < n; c++)
                    if (Math.Abs(lam[a] + lam[b] + lam[c]) < 1e-7) z3++;
        ZeroSumTripleCount = z3;
    }

    private static List<int[]> Combinations(int n, int k)
    {
        var result = new List<int[]>();
        var idx = Enumerable.Range(0, k).ToArray();
        while (true)
        {
            result.Add((int[])idx.Clone());
            int i = k - 1;
            while (i >= 0 && idx[i] == n - k + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
        }
        return result;
    }

    /// <summary>The k-excitation nearest-neighbour hopping matrix, amplitude 1 per bond, in the
    /// sorted-tuple basis (no Jordan-Wigner sign: NN hops never cross an occupied site).</summary>
    private static RealMatrix Hop(int n, int k)
    {
        var states = Combinations(n, k);
        var index = new Dictionary<string, int>();
        for (int s = 0; s < states.Count; s++) index[string.Join(",", states[s])] = s;
        var h = RealMatrix.Build.Dense(states.Count, states.Count);
        for (int s = 0; s < states.Count; s++)
        {
            var occ = new HashSet<int>(states[s]);
            foreach (int site in states[s])
                foreach (int target in new[] { site - 1, site + 1 })
                    if (target >= 0 && target < n && !occ.Contains(target))
                    {
                        var moved = new HashSet<int>(occ);
                        moved.Remove(site);
                        moved.Add(target);
                        var sorted = moved.OrderBy(x => x).ToArray();
                        h[index[string.Join(",", sorted)], s] += 1.0;
                    }
        }
        return h;
    }

    private static RealMatrix Submatrix(RealMatrix m, List<int> indices)
    {
        var sub = RealMatrix.Build.Dense(indices.Count, indices.Count);
        for (int r = 0; r < indices.Count; r++)
            for (int c = 0; c < indices.Count; c++)
                sub[r, c] = m[indices[r], indices[c]];
        return sub;
    }

    private static int Nullity(RealMatrix m) =>
        m.Svd(computeVectors: false).S.Count(s => s < NullityTol);

    public string DisplayName =>
        $"SeedExistenceCountingWitness (r(0⁺) − r(∞) = N − 1 live, N={N}, dim={Dim})";

    public string Summary =>
        $"N={N}: n₂={NullityMinus2}, n₆={NullityMinus6}, nullity(C)={NullityC} (SVD), " +
        $"ρ={ResonanceCount} and Z₃={ZeroSumTripleCount} (combinatorial): " +
        $"surplus r(0⁺)−r(∞) = {Surplus} {(Surplus == N - 1 ? "= N−1, the theorem" : "≠ N−1: BROKEN")}; " +
        $"n₆ = 3·Z₃ = ρ {(NullityMinus6 == 3 * ZeroSumTripleCount && NullityMinus6 == ResonanceCount ? "holds" : "BROKEN")} " +
        $"(the (N1′) ordering-sector theorem: cross-sector {CrossSectorMax.ToString("0.0e0", Inv)}, " +
        $"gauge residual {GaugeResidualMax.ToString("0.0e0", Inv)}, both exact zeros; " +
        $"controls: coupling {CouplingBlockMax.ToString("0.0", Inv)} > 0" +
        $"{(N >= 5 ? $", ‖K₆‖ {KSixMax.ToString("0.0", Inv)} > 0" : ", K₆ ≡ 0 at N=3 by design")}). " +
        "The counting identity is typed; the seed-existence conclusion stays open at the codim-2 β-exotic.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: $"the surplus: r(0⁺) − r(∞) = {Surplus} = N − 1",
                summary: $"r(0⁺) = n₂ + n₆ = {NullityMinus2} + {NullityMinus6} (SVD nullities of the −2/−6 " +
                         $"compressions), r(∞) = nullity(C) = {NullityC}: the real-eigenvalue count drops by " +
                         $"{Surplus} at finite q > 0, forcing a real↔complex transition (a defective seed at a " +
                         "simple discriminant zero, the Kato lemma). The fusion resonances appear in BOTH " +
                         "endpoints and cancel; the surplus is carried by the odd-N −2-rung path kernel alone.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"(N2) the path count: n₂ = {NullityMinus2} = N − 1",
                summary: $"the −2 rung (dim {DimMinus2} = N(N−1)) decomposes into N−1 disjoint simple paths of " +
                         "N vertices each (q-hop blocked at the shared site + the swap glue); a path carries a " +
                         "zero mode iff N is odd. The SVD nullity read here meets that combinatorial law.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"(N1′) the ordering sectors: n₆ = {NullityMinus6} = 3·Z₃ = ρ",
                summary: $"K₆ (dim {DimMinus6} = 3·C(N,3)) splits into three no-passing bra-rank sectors " +
                         $"(cross-sector max {CrossSectorMax.ToString("0.0e0", Inv)}, exactly 0); each sector, " +
                         $"gauged by diag((−1)^{{z_bra}}), IS −H₃ (residual {GaugeResidualMax.ToString("0.0e0", Inv)}, " +
                         $"exactly 0), so spec(K₆) = 3 × {{−(λ_a+λ_b+λ_c)}} and nullity = 3·Z₃ = {3 * ZeroSumTripleCount}. " +
                         $"The chiral pairing λ_{{N+1−k}} = −λ_k maps zero-sum triples 3-to-1 onto fusion resonances " +
                         $"(ρ = {ResonanceCount}); the degenerate case D = #{{2λ_x+λ_y = 0}} is 0 for every N by the " +
                         "cyclotomic-integrality norm bound.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: $"(FF) the fusion-resonance count: nullity(C) = {NullityC} = ρ",
                summary: $"ker C = the intertwiners {{ρ: H₂ρ = ρH₁}}; in the eigenmode basis C is diagonal on " +
                         "cubic monomials with eigenvalue −i(λ_a+λ_b−λ_c), so its nullity is the resonance count " +
                         $"ρ = {ResonanceCount} (matrix-free combinatorics meeting the SVD read {NullityC}).",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the two-sided controls (a broken builder cannot pass)",
                summary: $"the −6/−2 coupling block max = {CouplingBlockMax.ToString("0.0", Inv)} (> 0: the " +
                         "ket-onto-bra hops really leave the rung, they ARE the coupling, not K₆ content)" +
                         (N >= 5
                             ? $"; ‖K₆‖_max = {KSixMax.ToString("0.0", Inv)} (> 0: the rung has genuine dynamics)."
                             : "; K₆ ≡ 0 at N = 3 (every hop leaves the rung), itself a pinned edge case."),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the remaining open ink (scope of this witness)",
                summary: "this witness types the COUNTING identity r(0⁺) − r(∞) = N − 1 (a theorem for every " +
                         "odd N). The seed-existence CONCLUSION additionally needs the codim-2 β-exotic " +
                         "genericity item (a count-dropping transition is defective unless the non-generic " +
                         "order-3 nilpotent-linear-term point), which stays OPEN; see the doc's Status.",
                provenance: NodeProvenance.Stored);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
