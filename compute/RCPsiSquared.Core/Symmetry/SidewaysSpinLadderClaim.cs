using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The sideways spin ladder (the <c>sideways_spin_ladder</c> open arc, measured 2026-08-07,
/// typed 2026-08-09): the staggered second ladder
///
/// <code>
///     S⁺(ρ) = Σ_l (−1)^l c_l† ρ c_l†      (p, q̃) → (p+1, q̃−1)
/// </code>
///
/// is §6's V of <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> with the stagger restored, which is
/// F142's spin raising (<c>docs/proofs/PROOF_FROZEN_BAND_SO4.md</c> Lemma 2.3): the (−1)^l is the whole
/// difference, and it is what makes S⁺ commute with the Hamiltonian part. Two halves with two standings:
///
/// <para><b>The intertwining half is DERIVED</b> (F142 Lemma 2.1 + Lemma 2.3): [D, S⁺] = 0 for every rate
/// profile γ_l (the ladder never moves the disagreement set), and, for REAL SYMMETRIC h (Lemma 2.3's
/// hypothesis, load-bearing: its proof uses h_{la} = h_{al}), [K, S⁺] = 0 iff Σ h Σ = −h with
/// Σ = diag((−1)^l), i.e. iff the hopping matrix connects only opposite-parity sites. SCOPE, stated once
/// and not softened: S⁺ does NOT cover F125's H class. On-site terms and any same-parity hopping are
/// excluded; arbitrary nearest-neighbour bond profiles and arbitrary γ_l are admitted (adjacent sites
/// always have opposite parity, and the 220-rung sweep ran at non-uniform J_b and γ_l); the break-inputs
/// below are the measured face of that boundary. S⁺ preserves p + q̃, the parity fact behind §6 item
/// (ii)'s route-unavailability remark: band blocks have p + q̃ odd, diagonal cores p + q̃ even, so no
/// staggered ladder connects the two at any N (the exclusion itself is closed by the rate window and the
/// R4 certificates, not by any ladder property).</para>
///
/// <para><b>The multiplet half: the S⁺ chains are MEASURED and GATED, N=5, 7 and 9; the η chains are
/// INFERRED.</b> The F125 fold family is exactly the interiors of the two S⁺ chains at p + q̃ = N∓1
/// (SHAPE gate), each an sl(2) chain of interior length N−2 carrying spin ℓ = (N−3)/2, and the S⁺
/// transport norms along a chain are the Clebsch-Gordan coefficients √(ℓ(ℓ+1) − m(m+1)) with no free
/// parameter, derived from the closed form and measured at N=7 as 2, √6, √6, 2 (the arc entry records
/// the prediction preceding that run). The two η chains carry the same spin by the weight bookkeeping,
/// since fold + band = 4N−8 = 4(N−2) is the multiplet-dimension accounting (COUNT gate); their transport
/// norms are UNTESTED, and σ_min is not them (next paragraph). The chain terminates at the eigensolver
/// floor at the step into a boundary block (highest weight S⁺|ℓ,ℓ⟩ = 0, a RATIO gate, not an exact zero:
/// the vector comes from an eigensolver while the intertwining residual is an operator identity and IS
/// exactly 0.0). The falsifiable next case N=9, ℓ=3 was CONFIRMED 2026-08-09, both chains, to six
/// decimals: norms √6, √10, √12, √12, √10, √6, walked past the dense wall (middle blocks 10584² and
/// 15876²) by one dense LU shift-invert per block plus inverse iteration
/// (<c>SidewaysSpinLadderSparse</c>, gate Category SLOW_SIDEWAYS9), CONTROL still compared to 0.0
/// exactly and exactly 0.0. So the multiplet half is measured at THREE odd N: 5, 7, 9.</para>
///
/// <para><b>The σ_min fence.</b> σ_min of a rung map is NOT the multiplet's CG coefficient: it reads the
/// weakest direction of the WHOLE rung map, several multiplets per block, and DISAGREES with the CG value
/// on the two middle N=7 rungs (1.414214 against 2.449490), rungs F125 never pinned. F125's pinned values
/// (1 at N=4, √2 at N=5) numerically EQUAL the CG value there, coinciding by smallness of the block, so
/// they must not be quoted as a confirmation of the multiplet reading; the gate holds the rejection.</para>
///
/// <para><b>Break-inputs</b> (measured, five geometries at N=5 and N=6): a next-nearest bond breaks S⁺,
/// an on-site potential breaks it, the ring breaks it at odd N and holds at even N, the star breaks it at
/// the labelling tried (F142 Corollary 2.4: bipartite, but centre-against-leaves, not even-against-odd);
/// Φ held on every row (it needs no condition on h). Read as hopping matrices, NOT as spin chains: the
/// distinction caught a real build error (σ⁺_iσ⁻_j is the fermionic hop only for adjacent sites).</para>
///
/// <para><b>Typed parent.</b> <see cref="SpectatorIntertwinerClaim"/> (F125, whose W is the η sibling Φ
/// and whose fold family is what the chains carry). Evidence:
/// <c>simulations/eta_ladder_chain.py</c> (the Python gate, exits non-zero) and
/// <c>compute/RCPsiSquared.Diagnostics.Tests/Foundation/SidewaysSpinLadderGateTests.cs</c> and the N=9
/// walk <c>SidewaysSpinLadderSparseTests.cs</c> (Category SLOW_SIDEWAYS9) (the C# gates,
/// Categories SIDEWAYS + SLOW_SIDEWAYS, CONTROL compared to 0.0 exactly); break-inputs:
/// <c>simulations/eta_ladder_breakinput.py</c>; rung sweep: <c>simulations/eta_ladder_blocks.py</c>;
/// arc entry: <c>compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs</c>.</para></summary>
public sealed class SidewaysSpinLadderClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: F125, whose W = Φ is the η sibling and whose fold family
    // the S⁺ chains explain as a multiplet.
    public SpectatorIntertwinerClaim Intertwiner { get; }

    /// <summary>The spin carried by each of the four transport chains at odd N: ℓ = (N−3)/2.</summary>
    public static double ChainSpin(int n) => (n - 3) / 2.0;

    /// <summary>The Clebsch-Gordan transport norm √(ℓ(ℓ+1) − m(m+1)) at weight m along a chain.</summary>
    public static double CgTransportNorm(double spin, double m) => Math.Sqrt(spin * (spin + 1) - m * (m + 1));

    /// <summary>F125's orbit size at ODD N, derived combinatorially in PROOF_CODIM1_BY_ADDITIVITY; that it
    /// equals the multiplet dimension 4(N−2) is the accounting the multiplet reading rests on. Odd N only:
    /// the proof gives 4N−12 at even N, where the chain reading is untested.</summary>
    public static int OrbitSizeOddN(int n) => (n & 1) == 1
        ? 4 * n - 8
        : throw new ArgumentOutOfRangeException(nameof(n),
            $"4N−8 is the ODD-N orbit size (even N is 4N−12, chain reading untested); got N={n}");

    public SidewaysSpinLadderClaim(SpectatorIntertwinerClaim intertwiner)
        : base("The sideways spin ladder: S⁺(ρ) = Σ_l (−1)^l c_l†ρc_l† (§6's V of PROOF_CODIM1_BY_ADDITIVITY " +
               "with the stagger restored = F142's spin raising) intertwines the Liouvillian exactly, " +
               "(p,q̃)→(p+1,q̃−1), for any rate profile γ_l and any REAL SYMMETRIC pure-hopping h with " +
               "ΣhΣ = −h (derived, F142 Lemmas 2.1 + 2.3; NOT F125's H class: on-site terms and same-parity " +
               "hopping excluded, nearest-neighbour bond disorder and arbitrary γ_l admitted, and the " +
               "measured break-inputs are the boundary: next-nearest bond, on-site potential, odd ring, star " +
               "all break it, even ring and the disordered open chain hold); S⁺ preserves p+q̃, so no " +
               "staggered ladder connects the odd-sum band blocks to the even-sum diagonal cores at any N " +
               "(the parity fact behind §6 item (ii); the core exclusion itself is closed by the rate window " +
               "and the R4 certificates); the F125 fold family is the interiors of the two S⁺ chains at " +
               "p+q̃ = N∓1, each carrying spin ℓ = (N−3)/2 with Clebsch-Gordan transport norms " +
               "√(ℓ(ℓ+1)−m(m+1)), no free parameter, measured and gated at N=5 and N=7 (2, √6, √6, 2 at " +
               "N=7); the two η chains carry the same spin by the 4N−8 = 4(N−2) odd-N multiplet accounting, " +
               "their transport norms untested; chain death = highest-weight annihilation at the eigensolver " +
               "floor, a ratio gate; σ_min is NOT the CG coefficient (differs on the middle N=7 rungs, √2 vs " +
               "√6; F125's pinned 1 and √2 coincide with CG by smallness and confirm nothing); the N=9, ℓ=3 " +
               "prediction √6, √10, √12, √12, √10, √6 CONFIRMED 2026-08-09 on both chains to six decimals " +
               "(LU shift-invert past the dense wall, gate SLOW_SIDEWAYS9), so the multiplet half is measured " +
               "at N = 5, 7, 9",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_FROZEN_BAND_SO4.md + docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md + " +
               "simulations/eta_ladder_chain.py + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/SidewaysSpinLadderGateTests.cs + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/SidewaysSpinLadderSparseTests.cs")
    {
        Intertwiner = intertwiner ?? throw new ArgumentNullException(nameof(intertwiner));
    }

    public override string DisplayName =>
        "Sideways spin ladder: S⁺ = Σ (−1)^l c_l†(·)c_l† carries the F125 fold family as one spin-(N−3)/2 multiplet";

    public override string Summary =>
        "S⁺(ρ) = Σ_l (−1)^l c_l†ρc_l† intertwines L exactly for Σ-odd real-symmetric pure hopping and any γ_l " +
        "(derived, F142 L2.1+2.3; scope NARROWER than F125's H class, break-inputs measured); preserves p+q̃ " +
        "(the parity fact behind §6 item (ii)); the fold family = the two S⁺ chain interiors at p+q̃ = N∓1, " +
        "spin ℓ = (N−3)/2 per chain, odd-N orbit 4N−8 = 4(N−2); S⁺ transport norms = CG coefficients, " +
        "measured N=5,7,9 (√2,√2; 2,√6,√6,2; √6,√10,√12,√12,√10,√6, the N=9 case predicted before it was " +
        "walked), the η chains inferred from the count, untested; σ_min is NOT a CG confirmation " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("chain spin ℓ = (N−3)/2 at N=7 (= 2)", ChainSpin(7));
            yield return InspectableNode.RealScalar("CG norm at N=7, edge rung (ℓ=2, m=−2) (= 2)",
                CgTransportNorm(2.0, -2.0));
            yield return InspectableNode.RealScalar("CG norm at N=7, middle rung (ℓ=2, m=−1) (= √6)",
                CgTransportNorm(2.0, -1.0));
            yield return new InspectableNode("the derived half: [D,S⁺] = 0 always; [K,S⁺] = 0 iff ΣhΣ = −h (real symmetric h)",
                summary: "F142 Lemma 2.1: neither ladder moves the disagreement set, at any profile γ_l. " +
                         "Lemma 2.3, stated for REAL SYMMETRIC h (the hypothesis is load-bearing, the proof uses " +
                         "h_{la} = h_{al}): [K,S⁺] = Σ h_{al}((−1)^a + (−1)^l)c†_{a↑}c_{l↓}, zero iff " +
                         "h connects only opposite-parity sites. The star is the sufficient-not-necessary fence " +
                         "(Corollary 2.4): symmetric spectrum, wrong two-colouring, S⁺ does not commute there.");
            yield return new InspectableNode("the parity fact behind §6 item (ii): S⁺ preserves p+q̃",
                summary: "any d-raising ladder of this shape maps (p,q̃)→(p+1,q̃−1): band blocks have p+q̃ odd, " +
                         "diagonal cores p+q̃ even, so no staggered ladder carries the band's cap onto the core at " +
                         "any N. This repairs item (ii)'s route-unavailability remark, not the exclusion: the core " +
                         "exclusion is closed by the rate window and the R4 certificates and never depended on any " +
                         "ladder property.");
            yield return new InspectableNode("the measured half (S⁺ chains) and the inferred half (η chains)",
                summary: "MEASURED: fold family = interiors of the two S⁺ chains at p+q̃ = N∓1 (SHAPE gate), " +
                         "fold + band = 4N−8 (COUNT gate), S⁺ norms = √(ℓ(ℓ+1)−m(m+1)) at 1e-6 (LADDER gate), " +
                         "terminal step at the eigensolver floor with a survivors negative control; N=5 and N=7 " +
                         "in C# and Python, N=9 in C# (LU shift-invert, SLOW_SIDEWAYS9), CONTROL residual " +
                         "compared to 0.0 EXACTLY everywhere. INFERRED: the two η " +
                         "chains carry the same spin by the 4N−8 = 4(N−2) accounting; their transport norms are " +
                         "untested (the η evidence in hand is rank + spectral containment, print-only).");
            yield return new InspectableNode("the σ_min fence (a rejected confirmation, kept refutable)",
                summary: "σ_min of a rung map reads the weakest direction present, not one multiplet's CG value: " +
                         "equal at N=4, N=5 and the outer N=7 rungs, DIFFERENT on the middle N=7 rungs " +
                         "(1.414214 vs 2.449490). Gated so the rejection stays measurable.");
            yield return new InspectableNode("N=9 confirmed; even N is what stays open",
                summary: "N=9, ℓ=3 CONFIRMED 2026-08-09: norms √6, √10, √12, √12, √10, √6 on both chains to " +
                         "six decimals, seven sectors per chain, 28 = 4·9−8; the 10584²/15876² middle blocks " +
                         "walked by one dense LU shift-invert each (SidewaysSpinLadderSparse, gate " +
                         "SLOW_SIDEWAYS9, ~2 min). The measurement is robust to the fold pair's closeness " +
                         "(members −12.878060/−12.880829, split 2.77e-3, eigenvector angle ~1.8e-3): both " +
                         "members AND every in-plane mix transport at the same CG norm, verified to nine " +
                         "decimals in review. OPEN: even N. No even-N defective seed is RECORDED " +
                         "(RealDefectiveSeeds lists odd N, lower bounds), and at even N the four members " +
                         "around half filling are simultaneously band and fold image " +
                         "(PROOF_CODIM1_BY_ADDITIVITY), so the even-N chain reading is untested, not predicted.");
        }
    }
}
