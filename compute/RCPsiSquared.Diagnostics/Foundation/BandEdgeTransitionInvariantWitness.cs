using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F124 computed live (typed home: <c>BandEdgeTransitionInvariantClaim</c>). For the open chain
/// with the band-edge carrier ψ_1, the FULL single-excitation bond-transition matrix
/// M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ (rows = bonds, cols = ALL N modes k=1..N) satisfies, exactly and for every N,
///
/// <para>  ‖M‖_F² + λ_min(M Mᵀ) = z = 2,   with ‖M‖_F² = 2 − E and λ_min = E,</para>
///
/// where z=2 is the chain's coordination number and E = c₀² + c_{N-1}² = (4/(N+1))·sin²(π/(N+1)) is the
/// carrier's weight on the two free ends. The non-trivial half is λ_min = E: a staggered (zone-boundary,
/// q=π) bond modulation Σ_b(−1)^b V_b couples to the band-edge carrier ONLY through the Dirichlet ends:
/// the bulk telescopes away via the conserved discrete-energy envelope Q_a = c_a²+c_{a+1}²−E₁ c_a c_{a+1} = c₀²
/// (E₁ = 2cos(π/(N+1)), the band edge / ClockHandLadder), an SSH/Peierls edge effect; the same E is the
/// deficit of the carrier's degree-weighted norm from z. One boundary quantity c₀² fixes both.
///
/// <para>This witness recomputes the identity over a chain sweep, the carrier-selection hazard (only the
/// band-edge carrier makes the staggered mode the genuine MINIMUM: an interior carrier keeps it an
/// eigenvector but not the least, so the sum drops below 2), the frame reading (λ_min = σ_min²(M), the
/// condition number κ = λ_max/λ_min, the K-partner null column ⟨ψ_N|V_b|ψ_1⟩ ≡ 0), and the object/topology
/// scope (the decoder's location dictionary k=2..N gives λ_min=0, sum≠2; the odd ring frustrates, sum>2;
/// the star breaks the trace half, ‖M‖_F²=N/2). The C# twin of simulations/_handshake_M_checksum.py and
/// simulations/_handshake_F124_adversarial.py.</para></summary>
public sealed class BandEdgeTransitionInvariantWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-9;

    /// <summary>The chain sweep the identity is read across (N=3 is the stated lower endpoint).</summary>
    private static readonly int[] ChainSweep = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    public enum Topo { Chain, Ring, Star }

    /// <summary>The full reading of the band-edge transition matrix on one (topology, N, carrier) point.</summary>
    public readonly record struct Reading(
        int N, double Fro2, double LamMin, double LamMax, double SigmaMinSq,
        double StagResid, double StagRayleigh, bool OffDiagAllPositive, double KPartnerColNorm)
    {
        public double Sum => Fro2 + LamMin;
        public double Kappa => LamMin > 0 ? LamMax / LamMin : double.PositiveInfinity;
        /// <summary>The staggered bond wave (−1)^b is the genuine minimum: an eigenvector AND its
        /// Rayleigh quotient equals λ_min (not merely some larger eigenvalue).</summary>
        public bool StaggeredIsGenuineMinimum =>
            StagResid < 1e-8 && Math.Abs(StagRayleigh - LamMin) < 1e-8;
    }

    // ---- topology ----

    private static List<(int I, int J)> Edges(Topo topo, int n)
    {
        var e = new List<(int, int)>();
        switch (topo)
        {
            case Topo.Chain: for (int i = 0; i < n - 1; i++) e.Add((i, i + 1)); break;
            case Topo.Ring: for (int i = 0; i < n; i++) e.Add((i, (i + 1) % n)); break;
            case Topo.Star: for (int i = 1; i < n; i++) e.Add((0, i)); break;
            default: throw new ArgumentOutOfRangeException(nameof(topo));
        }
        return e;
    }

    private static Matrix<double> Adjacency(Topo topo, int n)
    {
        var a = Matrix<double>.Build.Dense(n, n);
        foreach (var (i, j) in Edges(topo, n)) { a[i, j] = 1.0; a[j, i] = 1.0; }
        return a;
    }

    /// <summary>The full transition matrix M[b,k] = ψ_k(i)·ψ_1(j) + ψ_k(j)·ψ_1(i) over bonds b and ALL
    /// N adjacency eigenmodes k, ordered DESCENDING in adjacency eigenvalue (so column 0 = band edge =
    /// carrier k=1, column N−1 = bottom = the K-partner ψ_N). The carrier ψ_1 is the <paramref name="carrierRank"/>-th
    /// mode from the top (rank 0 = the nodeless band edge).</summary>
    public static Matrix<double> BuildM(Topo topo, int n, int carrierRank)
    {
        var evd = Adjacency(topo, n).Evd(Symmetricity.Symmetric);
        var vals = evd.EigenValues.Select(c => c.Real).ToArray();
        var order = Enumerable.Range(0, n).OrderByDescending(i => vals[i]).ToArray();
        var vecs = evd.EigenVectors;
        var carrier = vecs.Column(order[carrierRank]);
        var edges = Edges(topo, n);
        var m = Matrix<double>.Build.Dense(edges.Count, n);
        for (int b = 0; b < edges.Count; b++)
        {
            var (i, j) = edges[b];
            for (int k = 0; k < n; k++)
            {
                var mode = vecs.Column(order[k]);
                m[b, k] = mode[i] * carrier[j] + mode[j] * carrier[i];
            }
        }
        return m;
    }

    /// <summary>The closed form for the endpoint leakage E = c₀²+c_{N-1}² = (4/(N+1))·sin²(π/(N+1)).</summary>
    public static double EndpointClosedForm(int n) =>
        (4.0 / (n + 1)) * Math.Pow(Math.Sin(Math.PI / (n + 1)), 2);

    /// <summary>The conserved discrete-energy envelope Q_a = c_a²+c_{a+1}²−E₁ c_a c_{a+1}, E₁=2cos(π/(N+1)).
    /// Returns (isConstantAlongChain, equalsC0Squared): Part 2's bulk-cancellation heart, the reason
    /// λ_min = 2Q₀ = 2c₀² = E.</summary>
    public static (bool Constant, bool EqualsC0Sq) EnvelopeCheck(int n)
    {
        double norm = Math.Sqrt(2.0 / (n + 1));
        var c = Enumerable.Range(0, n).Select(i => norm * Math.Sin(Math.PI * (i + 1) / (n + 1))).ToArray();
        double e1 = 2.0 * Math.Cos(Math.PI / (n + 1));
        var q = Enumerable.Range(0, n - 1).Select(a => c[a] * c[a] + c[a + 1] * c[a + 1] - e1 * c[a] * c[a + 1]).ToArray();
        bool constant = q.All(qa => Math.Abs(qa - q[0]) < 1e-12);
        bool eqC0 = Math.Abs(q[0] - c[0] * c[0]) < 1e-12;
        return (constant, eqC0);
    }

    /// <summary>The full reading on one (topology, N, carrier) point: Frobenius energy, the MMᵀ spectrum
    /// floor/ceiling, σ_min²(M), the staggered-mode residual + Rayleigh quotient, the band-edge off-diagonal
    /// positivity, and the K-partner column norm.</summary>
    public static Reading Analyse(Topo topo, int n, int carrierRank = 0)
    {
        var m = BuildM(topo, n, carrierRank);
        double fro2 = m.PointwiseMultiply(m).RowSums().Sum();
        var g = m * m.Transpose();                       // bonds × bonds Gram
        var gEvals = g.Evd(Symmetricity.Symmetric).EigenValues.Select(c => c.Real).OrderBy(x => x).ToArray();
        double lamMin = gEvals[0];
        double lamMax = gEvals[^1];
        double sigmaMinSq = Math.Pow(m.Svd().S.Minimum(), 2);

        int nb = m.RowCount;
        var stag = Vector<double>.Build.Dense(nb, b => (b % 2 == 0) ? 1.0 : -1.0);
        var gStag = g * stag;
        double rayleigh = stag.DotProduct(gStag) / stag.DotProduct(stag);
        double resid = (gStag - rayleigh * stag).L2Norm();

        bool offPos = true;
        for (int a = 0; a < nb - 1 && offPos; a++) offPos &= g[a, a + 1] > 0;

        double kpartner = m.Column(n - 1).L2Norm();      // ⟨ψ_N|V_b|ψ_1⟩ over bonds (band-edge K-partner)
        return new Reading(n, fro2, lamMin, lamMax, sigmaMinSq, resid, rayleigh, offPos, kpartner);
    }

    /// <summary>λ_min of the DECODER'S LOCATION DICTIONARY (k=2..N, the strength channel k=1 dropped):
    /// the K-partner column is null there, so M_loc loses a column ⟹ M_loc M_locᵀ is rank-deficient ⟹
    /// λ_min = 0. The clean "2" needs the strength column; this is the object guard.</summary>
    public static double LocationDictionaryLamMin(int n)
    {
        var m = BuildM(Topo.Chain, n, 0);
        var loc = m.SubMatrix(0, m.RowCount, 1, n - 1);  // drop column 0 (band edge / strength channel k=1)
        var g = loc * loc.Transpose();
        return g.Evd(Symmetricity.Symmetric).EigenValues.Select(c => c.Real).Min();
    }

    // ---- IInspectable ----

    public string DisplayName =>
        "BandEdgeTransitionInvariantWitness (F124: ‖M‖_F² + λ_min(MMᵀ) = 2, λ_min = E the Dirichlet-edge coupling)";

    public string Summary
    {
        get
        {
            var r5 = Analyse(Topo.Chain, 5);
            return
                "F124 computed live (typed home: BandEdgeTransitionInvariantClaim): for the open chain, the full " +
                "bond-transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ (all N modes) has ‖M‖_F² + λ_min(MMᵀ) = z = 2 exactly, " +
                $"with ‖M‖_F²=2−E and λ_min=E=(4/(N+1))sin²(π/(N+1)) (N=5: sum={r5.Sum.ToString("0.######", Inv)}, " +
                $"λ_min={r5.LamMin.ToString("0.######", Inv)}, E={EndpointClosedForm(5).ToString("0.######", Inv)}). " +
                "The real content is λ_min=E: a staggered bond modulation couples to the band-edge carrier only " +
                "through the Dirichlet ends (an SSH/Peierls edge effect), and the same E is the carrier's degree-" +
                "weighted-norm deficit from z. Carrier-selecting: only the band edge makes the staggered mode the " +
                "MINIMUM (interior carrier → sum<2). Frame reading: λ_min=σ_min², kernel = the K-partner ψ_N. " +
                "Object guard: the location dictionary k=2..N gives λ_min=0, sum≠2.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheIdentityLadder();
            yield return TheCarrierSelection();
            yield return TheFrameReading();
            yield return TheObjectAndTopologyScope();
        }
    }

    /// <summary>The headline identity across the chain sweep: sum=2, ‖M‖_F²=2−E, λ_min=E, and the staggered
    /// bond wave is the genuine minimum. The non-trivial half is λ_min=E (the Dirichlet-edge coupling).</summary>
    private InspectableNode TheIdentityLadder()
    {
        var rows = new List<IInspectable>();
        bool allOk = true;
        foreach (int n in ChainSweep)
        {
            var r = Analyse(Topo.Chain, n);
            double e = EndpointClosedForm(n);
            bool ok = Math.Abs(r.Sum - 2.0) < 1e-9
                      && Math.Abs(r.Fro2 - (2.0 - e)) < 1e-9
                      && Math.Abs(r.LamMin - e) < 1e-9
                      && r.StaggeredIsGenuineMinimum;
            allOk &= ok;
            var (qConst, qC0) = EnvelopeCheck(n);
            rows.Add(new InspectableNode($"chain N={n}",
                summary: $"‖M‖_F²={r.Fro2.ToString("0.######", Inv)} (2−E={(2 - e).ToString("0.######", Inv)}), " +
                         $"λ_min={r.LamMin.ToString("0.######", Inv)} (E={e.ToString("0.######", Inv)}), " +
                         $"sum={r.Sum.ToString("0.######", Inv)}; staggered is the genuine minimum: {r.StaggeredIsGenuineMinimum}; " +
                         $"envelope Q constant={qConst} & =c₀²={qC0} → {(ok ? "OK" : "GATE-FIRE")}"));
        }
        return new InspectableNode("the identity ‖M‖_F² + λ_min(MMᵀ) = 2 (chain, band-edge carrier, N=3..12)",
            summary: $"‖M‖_F²=2−E (degree counting, basis-independent) + λ_min=E (the Dirichlet-edge coupling: the " +
                     $"bulk telescopes via the conserved envelope Q_a=c₀², only the boundary survives) → sum = z = 2. " +
                     $"The staggered bond wave (−1)^b is the λ_min eigenvector. All N: {(allOk ? "OK" : "A GATE FIRED")}.",
            children: rows);
    }

    /// <summary>The carrier-selecting step (the documented gate-first hazard): Part 2's conserved envelope makes
    /// the staggered vector an eigenvector for ANY carrier (eigenvalue 2Q_k), but only for the nodeless band-edge
    /// carrier are the Gram off-diagonals all positive (Perron), which is what makes staggered the MINIMUM. For
    /// an interior carrier the staggered mode is still an eigenvector but NOT the least, and the sum drops below 2.</summary>
    private InspectableNode TheCarrierSelection()
    {
        var rows = new List<IInspectable>();
        foreach (int n in new[] { 7, 9, 11 })
        {
            var edge = Analyse(Topo.Chain, n, 0);
            rows.Add(new InspectableNode($"chain N={n}, band-edge carrier (rank 0)",
                summary: $"off-diagonals all > 0: {edge.OffDiagAllPositive}; staggered is the genuine minimum: " +
                         $"{edge.StaggeredIsGenuineMinimum}; sum={edge.Sum.ToString("0.######", Inv)} (= 2)"));
            for (int rank = 1; rank <= 2; rank++)
            {
                var r = Analyse(Topo.Chain, n, rank);
                bool stagEigen = r.StagResid < 1e-8;
                bool notMin = r.StagRayleigh > r.LamMin + 1e-8;
                bool sumBelow = r.Sum < 2.0 - 1e-8;
                bool expected = stagEigen && notMin && sumBelow;  // the band-edge is load-bearing for "=2"
                rows.Add(new InspectableNode($"chain N={n}, interior carrier (rank {rank})",
                    summary: $"staggered still an eigenvector (resid={r.StagResid:0.0e+0}): {stagEigen}; but NOT the " +
                             $"minimum (Rayleigh {r.StagRayleigh.ToString("0.####", Inv)} > λ_min {r.LamMin.ToString("0.####", Inv)}): {notMin}; " +
                             $"off-diag all>0: {r.OffDiagAllPositive}; sum={r.Sum.ToString("0.####", Inv)} < 2: {sumBelow} → " +
                             $"{(expected ? "band edge is load-bearing (OK)" : "GATE-FIRE")}"));
            }
        }
        return new InspectableNode("the carrier-selecting step: only the band edge makes λ_min the genuine minimum",
            summary: "the conserved envelope makes the staggered mode an eigenvector for ANY carrier (eigenvalue 2Q_k), " +
                     "but the equality sum=2 needs it to be the LEAST eigenvalue, true only for the nodeless band-edge " +
                     "carrier, whose Gram off-diagonals are all positive (Perron). An interior carrier keeps staggered " +
                     "an eigenvector but not the minimum, and the sum drops strictly below 2. This is why λ_min must be " +
                     "read as the genuine minimum (full spectrum), not merely 'staggered is an eigenvalue'.",
            children: rows);
    }

    /// <summary>The frame reading (grounding-in-the-quantum + borrowing-a-discipline): {V_bψ_1} is a deficient
    /// non-tight Riesz basis; λ_min is the lower frame bound = σ_min²(M) = the Eckart-Young distance² to
    /// rank-collapse; the end-leakage E is the conditioner, κ = λ_max/λ_min; the exact kernel is the K-partner
    /// ψ_N (⟨ψ_N|V_b|ψ_1⟩ ≡ 0, the typed KPartnerSelectionRuleClaim). Honesty: λ_min=σ_min²(M) is a definitional
    /// Evd↔SVD identity (the eigenvalues of M Mᵀ ARE the squared singular values), so it cannot fire; the
    /// FALSIFIABLE content of the frame reading is the K-partner null column (a real selection rule) and the
    /// closed-form λ_min=E (gated in the identity ladder).</summary>
    private InspectableNode TheFrameReading()
    {
        var rows = new List<IInspectable>();
        double prevKappa = 0;
        bool kappaGrows = true;
        foreach (int n in new[] { 4, 5, 6, 7, 8 })
        {
            var r = Analyse(Topo.Chain, n);
            bool frameOk = Math.Abs(r.SigmaMinSq - r.LamMin) < 1e-9 && r.KPartnerColNorm < 1e-9;
            if (n > 4 && r.Kappa <= prevKappa) kappaGrows = false;
            prevKappa = r.Kappa;
            rows.Add(new InspectableNode($"chain N={n}",
                summary: $"λ_min={r.LamMin.ToString("0.######", Inv)} = σ_min²={r.SigmaMinSq.ToString("0.######", Inv)} " +
                         $"({(Math.Abs(r.SigmaMinSq - r.LamMin) < 1e-9 ? "match" : "MISMATCH")}); κ=λ_max/λ_min={r.Kappa.ToString("0.##", Inv)}; " +
                         $"K-partner column ‖⟨ψ_N|V_b|ψ_1⟩‖={r.KPartnerColNorm:0.0e+0} ({(r.KPartnerColNorm < 1e-9 ? "null" : "NONZERO")}) → " +
                         $"{(frameOk ? "OK" : "GATE-FIRE")}"));
        }
        return new InspectableNode("the frame reading: λ_min = σ_min² = lower frame bound; kernel = the K-partner",
            summary: $"the bond-scattered carriers {{V_bψ_1}} are a deficient, non-tight Riesz basis (rank N−1); λ_min is " +
                     $"the optimal lower frame bound = σ_min²(M) (a definitional Evd↔SVD identity, not a falsifiable check). " +
                     $"The end-leakage E is the CONDITIONER: κ = λ_max/λ_min → ∞ (a theorem, since λ_min = E → 0 like N⁻³ " +
                     $"while λ_max stays O(1/N) at most; strict step-monotonicity 3.4, 5.3, 7.6 at N=4,5,6 observed here, " +
                     $"{(kappaGrows ? "monotone on N=4..8" : "NOT monotone")}), so shorter chains are better-conditioned and " +
                     $"the long-chain / ring limit (E→0) goes singular. The FALSIFIABLE kernel content: the K-partner ψ_N " +
                     $"column ⟨ψ_N|V_b|ψ_1⟩ ≡ 0 (KPartnerSelectionRuleClaim).",
            children: rows);
    }

    /// <summary>The object guard and topology scope: the decoder's location dictionary k=2..N (strength channel
    /// dropped) gives λ_min=0 and sum≠2 (the clean "2" needs the strength column, so this is the full transition
    /// matrix, not the location reading). Topology: the even ring holds DEGENERATELY (E=0, λ_min=0, sum=2, no
    /// boundary, not a second instance of the mechanism); the odd ring frustrates the staggering (sum>2); the
    /// star breaks the trace half (‖M‖_F²=N/2).</summary>
    private InspectableNode TheObjectAndTopologyScope()
    {
        var rows = new List<IInspectable>();

        // object guard: the location dictionary k=2..N
        foreach (int n in new[] { 4, 5, 6, 7, 8 })
        {
            double locLamMin = LocationDictionaryLamMin(n);
            var full = Analyse(Topo.Chain, n);
            bool guardOk = locLamMin < 1e-9 && Math.Abs(full.Sum - 2.0) < 1e-9;
            rows.Add(new InspectableNode($"location dictionary k=2..N, chain N={n}",
                summary: $"λ_min(M_loc M_locᵀ)={locLamMin:0.0e+0} (= 0: the K-partner is a null column ⟹ rank-deficient); " +
                         $"the full k=1..N matrix has λ_min={full.LamMin.ToString("0.####", Inv)}=E and sum=2. The strength " +
                         $"column lifts the floor from 0 to E → {(guardOk ? "OK" : "GATE-FIRE")}"));
        }

        // topology scope: even ring holds degenerately (E=0, λ_min=0, sum=2, no boundary)
        foreach (int n in new[] { 4, 6 })
        {
            var r = Analyse(Topo.Ring, n);
            bool degenerate = Math.Abs(r.Sum - 2.0) < 1e-9 && r.LamMin < 1e-9;
            rows.Add(new InspectableNode($"even ring N={n}",
                summary: $"‖M‖_F²={r.Fro2.ToString("0.####", Inv)}, λ_min={r.LamMin:0.0e+0} (= 0: no boundary, E=0), " +
                         $"sum={r.Sum.ToString("0.######", Inv)} = 2 DEGENERATELY (the perfect 2-colouring closes the " +
                         $"staggering with zero boundary leakage, NOT a second instance of the mechanism) → " +
                         $"{(degenerate ? "holds degenerately (OK)" : "GATE-FIRE")}"));
        }

        // topology scope: odd ring frustrates (sum>2)
        foreach (int n in new[] { 5, 7 })
        {
            var r = Analyse(Topo.Ring, n);
            bool frustrated = r.Sum > 2.0 + 1e-6;
            rows.Add(new InspectableNode($"odd ring N={n}",
                summary: $"‖M‖_F²={r.Fro2.ToString("0.####", Inv)}, λ_min={r.LamMin.ToString("0.####", Inv)}>0, " +
                         $"sum={r.Sum.ToString("0.####", Inv)} > 2 (the odd cycle cannot 2-colour the staggering) → " +
                         $"{(frustrated ? "frustrated as expected (OK)" : "GATE-FIRE")}"));
        }

        // topology scope: star breaks the trace half (‖M‖_F² = N/2)
        foreach (int n in new[] { 5, 6 })
        {
            var r = Analyse(Topo.Star, n);
            bool starOk = Math.Abs(r.Fro2 - n / 2.0) < 1e-9;
            rows.Add(new InspectableNode($"star N={n}",
                summary: $"‖M‖_F²={r.Fro2.ToString("0.####", Inv)} = N/2 = {(n / 2.0).ToString("0.####", Inv)} " +
                         $"({(starOk ? "match" : "MISMATCH")}); the hub's high degree breaks the z=2 trace half → " +
                         $"{(starOk ? "OK" : "GATE-FIRE")}"));
        }

        return new InspectableNode("the object guard + topology scope (where the clean '2' breaks)",
            summary: "the full transition matrix k=1..N is NOT the decoder's location dictionary k=2..N: dropping the " +
                     "strength column k=1 leaves the K-partner null column, M_loc M_locᵀ is rank-deficient and λ_min=0, " +
                     "sum≠2. The clean 2 needs the strength column. Topology: chain (Dirichlet) holds; the even ring holds " +
                     "degenerately (E=0, λ_min=0, sum=2, no boundary); the odd ring frustrates the staggering (λ_min>0, " +
                     "sum>2); the star breaks the trace half at the hub (‖M‖_F²=N/2).",
            children: rows);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
