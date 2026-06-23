using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Witness = RCPsiSquared.Diagnostics.Foundation.BandEdgeTransitionInvariantWitness;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>F124, the band-edge transition invariant (the reading-grammar arc's frame-theoretic capstone):
/// for the open chain with the band-edge carrier ψ_1, the FULL single-excitation bond-transition matrix
/// M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ (all N modes k=1..N) satisfies, exactly and for every N,
///
/// <para>  ‖M‖_F² + λ_min(M Mᵀ) = z = 2,   with ‖M‖_F² = 2 − E and λ_min = E,</para>
///
/// where z=2 is the chain's coordination number and E = c₀²+c_{N-1}² = (4/(N+1))·sin²(π/(N+1)) is the
/// carrier's weight on the two free ends. The non-trivial half is λ_min = E: a staggered (zone-boundary,
/// q=π) bond modulation Σ_b(−1)^b V_b couples to the band-edge carrier ONLY through the Dirichlet ends: the
/// bulk telescopes away via the conserved discrete-energy envelope Q_a = c_a²+c_{a+1}²−E₁ c_a c_{a+1} = c₀²
/// (E₁ = 2cos(π/(N+1)), the band edge / ClockHandLadder), an SSH/Peierls edge effect; and the SAME E is the
/// deficit of the carrier's degree-weighted norm from z. One boundary quantity c₀² fixes both the spectral
/// floor and the trace deficit.
///
/// <para><b>The frame reading</b> (grounding-in-the-quantum + borrowing-a-discipline converged): the
/// bond-scattered carriers {V_bψ_1} are a deficient, non-tight Riesz basis (rank N−1); λ_min is the optimal
/// lower frame bound = σ_min²(M) = the Eckart-Young squared distance to rank-collapse (this last equality is
/// a definitional Evd↔SVD identity, not a falsifiable check; the falsifiable kernel content is the K-partner
/// null column); the end-leakage E is the conditioner, condition number κ = λ_max/λ_min → ∞ (a theorem, since
/// λ_min = E → 0 like N⁻³); and the exact kernel is the K-partner ψ_N: ⟨ψ_N|V_b|ψ_1⟩ ≡ 0, the typed
/// <see cref="KPartnerSelectionRuleClaim"/>.</para>
///
/// <para><b>Carrier-selecting (the genuine-minimum hazard).</b> Part 2's conserved envelope makes the
/// staggered mode an eigenvector for ANY carrier (eigenvalue 2Q_k); only the nodeless band-edge carrier has
/// all-positive Gram off-diagonals (Perron), which is what makes the staggered mode the MINIMUM. An interior
/// carrier keeps it an eigenvector but not the least, and the sum drops strictly below 2, so λ_min must be
/// read as the genuine minimum, not merely "staggered is an eigenvalue".</para>
///
/// <para><b>Object guard.</b> M is the FULL k=1..N matrix, NOT the decoder's location dictionary k=2..N
/// (<see cref="KPartnerSelectionRuleClaim"/>): dropping the strength channel k=1 leaves the K-partner null
/// column, so M_loc M_locᵀ is rank-deficient and λ_min=0, sum≠2. The clean 2 needs the strength column.
/// Topology: chain (Dirichlet) holds; the odd ring frustrates the staggering (sum>2); the star breaks the
/// trace half (‖M‖_F²=N/2).</para>
///
/// <para>Tier1Derived; gate-exact N=3..20 (proof verified, three Python verifiers + this live battery). Two
/// typed parents, both Tier1Derived: <see cref="KPartnerSelectionRuleClaim"/> (the kernel = the K-partner;
/// the same M, here completed with the strength column) and <see cref="ClockHandLadderClaim"/> (the band edge
/// E₁ = 2cos(π/(N+1)) the conserved envelope rides on, which selects the carrier).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md</c> + <c>docs/ANALYTICAL_FORMULAS.md</c>
/// F124 + <c>simulations/handshake_M_checksum.py</c>, <c>handshake_M_topology.py</c>,
/// <c>handshake_F124_adversarial.py</c>. Live witness: <c>inspect --root transition</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/BandEdgeTransitionInvariantWitness.cs</c>).</para></summary>
public sealed class BandEdgeTransitionInvariantClaim : Claim
{
    /// <summary>One pure-linear-algebra self-check on the small dense single-excitation transition matrix.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    private readonly KPartnerSelectionRuleClaim _kPartner;
    private readonly ClockHandLadderClaim _clock;

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    /// <summary>Parent: the K-partner selection rule. The frame's exact kernel IS the K-partner ψ_N
    /// (⟨ψ_N|V_b|ψ_1⟩ ≡ 0); F124 is the same M completed with the strength column k=1. Tier1Derived.</summary>
    public KPartnerSelectionRuleClaim KPartner => _kPartner;

    /// <summary>Parent: the clock-hand ladder. The band edge E₁ = 2cos(π/(N+1)) is the carrier the conserved
    /// envelope Q rides on; it selects the unique carrier for which the sum is exactly 2. Tier1Derived.</summary>
    public ClockHandLadderClaim Clock => _clock;

    public BandEdgeTransitionInvariantClaim(KPartnerSelectionRuleClaim kPartner, ClockHandLadderClaim clock)
        : base("F124 the band-edge transition invariant: for the open chain's band-edge carrier the full bond-" +
               "transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ (all N modes) has ‖M‖_F² + λ_min(MMᵀ) = z = 2 exactly, with " +
               "‖M‖_F²=2−E and λ_min=E=(4/(N+1))sin²(π/(N+1)); the non-trivial half λ_min=E is the Dirichlet-edge " +
               "coupling (a staggered bond modulation reaches the band-edge carrier only through the two free ends, " +
               "the bulk telescoping via the conserved envelope Q=c₀²), an SSH/Peierls edge effect, and the same E " +
               "is the carrier's degree-weighted-norm deficit from z. Frame reading: {V_bψ_1} a deficient Riesz basis, " +
               "λ_min=σ_min²=the lower frame bound, kernel = the K-partner ψ_N. Only the band-edge carrier makes the " +
               "staggered mode the genuine minimum (interior carrier → sum<2); the location dictionary k=2..N gives " +
               "λ_min=0, sum≠2 (the strength column is what lifts the floor to E).",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md (the proof) + " +
               "docs/ANALYTICAL_FORMULAS.md F124 + " +
               "simulations/handshake_M_checksum.py + simulations/handshake_M_topology.py + " +
               "simulations/handshake_F124_adversarial.py (the three gate-exact verifiers, N=3..20) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/BandEdgeTransitionInvariantWitness.cs (the live witness, inspect --root transition)")
    {
        _kPartner = kPartner ?? throw new ArgumentNullException(nameof(kPartner));
        _clock = clock ?? throw new ArgumentNullException(nameof(clock));
        Cases = BuildBattery();
    }

    public string Identity =>
        "‖M‖_F² + λ_min(MMᵀ) = z = 2 for the open chain's band-edge carrier, exactly and for every N. " +
        "‖M‖_F² = 2 − E (degree counting, basis-independent) and λ_min = E = (4/(N+1))sin²(π/(N+1)) (the " +
        "Dirichlet-edge coupling); the λ_min eigenvector is the staggered bond wave (−1)^b.";

    public string TheRealContent =>
        "λ_min = E. A staggered (zone-boundary, q=π) bond modulation couples to the band-edge carrier ONLY " +
        "through the two free ends: the bulk telescopes away via the conserved discrete-energy envelope " +
        "Q_a = c_a²+c_{a+1}²−E₁ c_a c_{a+1} = c₀² (E₁ = 2cos(π/(N+1)) the band edge), leaving only the boundary " +
        "value 2c₀² = E. An SSH/Peierls edge effect: a dimerization's grip on a band-edge state is boundary-dominated.";

    public string FrameReading =>
        "{V_bψ_1} is a deficient, non-tight Riesz basis (rank N−1). λ_min is the optimal lower frame bound = " +
        "σ_min²(M) = the Eckart-Young squared distance to rank-collapse (λ_min = σ_min² is a definitional Evd↔SVD " +
        "identity, not a falsifiable gate); the end-leakage E is the CONDITIONER (κ = λ_max/λ_min → ∞ as a theorem " +
        "since λ_min = E → 0 like N⁻³, strict step-monotonicity observed on N=4..8, shorter chains better-conditioned); " +
        "the falsifiable kernel content is the K-partner ψ_N null column (⟨ψ_N|V_b|ψ_1⟩ ≡ 0, KPartnerSelectionRuleClaim). " +
        "F124 is a completeness ⊕ conditioning identity: the Parseval deficit (z − trace) and the lower frame bound A " +
        "are both E.";

    public string ObjectGuard =>
        "M is the FULL transition matrix k=1..N, NOT the decoder's location dictionary k=2..N. Dropping the " +
        "strength channel k=1 leaves the K-partner null column, so M_loc M_locᵀ is rank-deficient and λ_min=0, " +
        "sum≠2. The clean 2 needs the strength column. The '2' itself is the coordination number z (a contract " +
        "riding on ‖V_b‖²=2), not a conservation constant: the odd ring frustrates (sum>2), the star breaks the " +
        "trace half (‖M‖_F²=N/2).";

    public override string DisplayName =>
        "F124 the band-edge transition invariant (‖M‖_F² + λ_min(MMᵀ) = 2, λ_min = E the Dirichlet-edge coupling, Tier1Derived)";

    public override string Summary =>
        "for the open chain's band-edge carrier, the full bond-transition matrix has ‖M‖_F² + λ_min(MMᵀ) = z = 2 " +
        "exactly (‖M‖_F²=2−E, λ_min=E); the real content λ_min=E is the Dirichlet-edge coupling (an SSH/Peierls " +
        "edge effect, the bulk telescoping via the conserved envelope), one boundary quantity c₀² fixing both the " +
        "spectral floor and the trace deficit; frame reading λ_min=σ_min²=the lower frame bound, kernel = the " +
        $"K-partner ψ_N; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the identity (‖M‖_F² + λ_min = z = 2)", summary: Identity);
            yield return new InspectableNode("the real content (λ_min = E, the Dirichlet-edge coupling)", summary: TheRealContent);
            yield return new InspectableNode("the frame reading (λ_min = σ_min² = lower frame bound; kernel = K-partner)", summary: FrameReading);
            yield return new InspectableNode("the object guard (full k=1..N, not the location dictionary k=2..N)", summary: ObjectGuard);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return _kPartner;   // typed parent edge (the kernel = the K-partner)
            yield return _clock;      // typed parent edge (the band edge that selects the carrier)
        }
    }

    public static BandEdgeTransitionInvariantClaim Build() =>
        new(new KPartnerSelectionRuleClaim(new ChiralMirrorTrajectoryClaim()), ClockHandLadderClaim.Build());

    public static BandEdgeTransitionInvariantClaim Shared { get; } = Build();

    /// <summary>Gate-first battery on the small dense transition matrix (no Liouvillian, no propagation): the
    /// headline identity, the genuine-minimum carrier-selection hazard, the frame identities, and the object /
    /// topology breakages. The FALSIFIABLE gates (these can fire on a wrong physics or a wrong generalization)
    /// are: the identity sum=2 and λ_min=E closed form (1-3), the genuine-minimum carrier selection (4-5), the
    /// K-partner null column (the kernel half of 6), κ growth (7), and the topology breakages (9-11: even ring
    /// degenerate, odd ring sum&gt;2, star N/2). Two sub-checks are STRUCTURAL identities that demonstrate the
    /// object distinction rather than test it: λ_min=σ_min² (the σ-half of case 6 is a definitional Evd↔SVD
    /// identity) and the location-dictionary λ_min=0 (case 8, a rank corollary of the K-partner null column).
    /// Computation delegates to the live witness.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();
        var chain = new[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

        string Tol(bool ok) => ok ? "OK" : "FAIL";

        // 1. the headline identity sum=2, all N
        bool sumOk = chain.All(n => Math.Abs(Witness.Analyse(Witness.Topo.Chain, n).Sum - 2.0) < 1e-9);
        cases.Add(new BatteryCase("identity ‖M‖_F² + λ_min(MMᵀ) = 2 (chain band-edge, N=3..12)",
            "the full transition matrix's Frobenius energy + Gram floor sum to the coordination number z=2",
            "OK", Tol(sumOk)));

        // 2. the real content λ_min = E, all N
        bool lamOk = chain.All(n =>
            Math.Abs(Witness.Analyse(Witness.Topo.Chain, n).LamMin - Witness.EndpointClosedForm(n)) < 1e-9);
        cases.Add(new BatteryCase("λ_min(MMᵀ) = E = (4/(N+1))sin²(π/(N+1)) (the Dirichlet-edge coupling, N=3..12)",
            "the spectral floor equals the carrier's weight on the two free ends, in closed form",
            "OK", Tol(lamOk)));

        // 3. the trace half ‖M‖_F² = 2 − E, all N
        bool froOk = chain.All(n =>
            Math.Abs(Witness.Analyse(Witness.Topo.Chain, n).Fro2 - (2.0 - Witness.EndpointClosedForm(n))) < 1e-9);
        cases.Add(new BatteryCase("‖M‖_F² = 2 − E (degree-weighted norm deficit, N=3..12)",
            "the same boundary quantity E is the deficit of the carrier's degree-weighted norm from z",
            "OK", Tol(froOk)));

        // 4. the staggered bond wave is the GENUINE minimum (band-edge), all N
        bool stagOk = chain.All(n => Witness.Analyse(Witness.Topo.Chain, n).StaggeredIsGenuineMinimum);
        cases.Add(new BatteryCase("staggered (−1)^b is the genuine λ_min eigenvector (band-edge, N=3..12)",
            "the λ_min eigenvector is the staggered bond wave, and its Rayleigh quotient IS the minimum (full spectrum)",
            "OK", Tol(stagOk)));

        // 5. carrier-selection hazard: interior carrier keeps staggered an eigenvector but NOT the min, sum<2
        bool selectOk = new[] { 7, 9, 11 }.All(n =>
        {
            var r = Witness.Analyse(Witness.Topo.Chain, n, 1);   // rank-1 interior carrier
            return r.StagResid < 1e-8 && r.StagRayleigh > r.LamMin + 1e-8 && r.Sum < 2.0 - 1e-8;
        });
        cases.Add(new BatteryCase("carrier-selecting: interior carrier → staggered an eigenvector but NOT the min, sum<2 (N=7,9,11)",
            "the band edge is load-bearing for '=2': an interior carrier keeps staggered an eigenvector yet not the least, so the sum drops below 2",
            "OK", Tol(selectOk)));

        // 6. frame identity: λ_min = σ_min²(M) AND the K-partner column is null, all N
        bool frameOk = chain.All(n =>
        {
            var r = Witness.Analyse(Witness.Topo.Chain, n);
            return Math.Abs(r.SigmaMinSq - r.LamMin) < 1e-9 && r.KPartnerColNorm < 1e-9;
        });
        cases.Add(new BatteryCase("frame reading: λ_min = σ_min²(M) [definitional] and ⟨ψ_N|V_b|ψ_1⟩ ≡ 0 [falsifiable] (the K-partner kernel, N=3..12)",
            "λ_min = σ_min² is a definitional Evd↔SVD identity (cannot fire); the falsifiable half is the K-partner null column, the exact kernel of the frame",
            "OK", Tol(frameOk)));

        // 7. the conditioner: κ = λ_max/λ_min grows with N (theorem: λ_min ~ N⁻³ → 0; step-monotonicity empirical here)
        var kappas = new[] { 4, 5, 6, 7, 8 }.Select(n => Witness.Analyse(Witness.Topo.Chain, n).Kappa).ToArray();
        bool kappaGrows = Enumerable.Range(1, kappas.Length - 1).All(i => kappas[i] > kappas[i - 1]);
        cases.Add(new BatteryCase("the conditioner: κ = λ_max/λ_min grows with N (N=4..8)",
            "κ → ∞ is a theorem (λ_min = E → 0 like N⁻³); strict step-monotonicity on N=4..8 is the observed instance; the long-chain limit goes singular",
            "OK", Tol(kappaGrows)));

        // 8. object guard: the location dictionary k=2..N gives λ_min = 0 (a structural corollary of the K-partner rule)
        bool guardOk = new[] { 4, 5, 6, 7, 8 }.All(n => Witness.LocationDictionaryLamMin(n) < 1e-9);
        cases.Add(new BatteryCase("object guard: location dictionary k=2..N has λ_min = 0 (N=4..8)",
            "a structural corollary of the K-partner null column (dropping the strength channel k=1 leaves it ⟹ M_loc M_locᵀ rank-deficient); demonstrates that the clean 2 needs the strength column",
            "OK", Tol(guardOk)));

        // 9. topology: the even ring holds degenerately (sum=2 AND λ_min=0: E=0, no boundary)
        bool evenRingOk = new[] { 4, 6 }.All(n =>
        {
            var r = Witness.Analyse(Witness.Topo.Ring, n);
            return Math.Abs(r.Sum - 2.0) < 1e-9 && r.LamMin < 1e-9;
        });
        cases.Add(new BatteryCase("topology: even ring holds degenerately, sum=2 AND λ_min=0 (N=4,6)",
            "the even ring 2-colours perfectly with zero boundary leakage (E=0): sum=2 holds but DEGENERATELY (λ_min=0, no boundary), not a second instance of the Dirichlet-edge mechanism",
            "OK", Tol(evenRingOk)));

        // 10. topology: the odd ring frustrates the staggering (sum > 2)
        bool oddRingOk = new[] { 5, 7 }.All(n => Witness.Analyse(Witness.Topo.Ring, n).Sum > 2.0 + 1e-6);
        cases.Add(new BatteryCase("topology: odd ring frustrates the staggering, sum > 2 (N=5,7)",
            "the '2' is the coordination number, not a universal constant; the odd cycle cannot 2-colour the staggered mode",
            "OK", Tol(oddRingOk)));

        // 11. topology: the star breaks the trace half (‖M‖_F² = N/2)
        bool starOk = new[] { 5, 6 }.All(n =>
            Math.Abs(Witness.Analyse(Witness.Topo.Star, n).Fro2 - n / 2.0) < 1e-9);
        cases.Add(new BatteryCase("topology: star breaks the trace half, ‖M‖_F² = N/2 (N=5,6)",
            "the hub's high degree changes the degree-weighted norm; the star is outside the z=2 chain/ring scope",
            "OK", Tol(starOk)));

        return cases;
    }
}
