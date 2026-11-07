using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F171, the grey rooms inside the coloured blocks: under Z dephasing at one rate γ on every site,
/// the η ladder Φ and the spin ladder S⁻ of <c>docs/proofs/PROOF_FROZEN_BAND_SO4.md</c> §2 carry the boundary
/// block (k, 0) of the XY chain, at any real couplings, one-to-one into every joint-popcount block (p, q)
/// whose Hamming window holds k, onto a subspace of dimension C(N, k) that the Hamiltonian part keeps and on
/// which the dissipator is −2γk: a grey room, where the motion comes through whole, e^{tL} = e^{−2γkt}·e^{tL_H},
/// and the modes do not move with γ. In the block (1, 1) of the uniform open chain the two images, C(N, 2) at
/// the distance 2 and the identity at 0, are all there is. On the even ring the image is kept at every even k,
/// and at odd k, with the wrap bond nonzero, exactly when one coupling-free map keeps it, the bra's hop across
/// the wrap bond alone, which at N = 4, 6, 8 it does only where the image fills a corner (the losses exact
/// there, and for (5, 2) at N = 10). On the odd ring the image is kept when it reaches its block from a corner
/// by the ladder its parity of k allows; that the rest is lost is measured at N = 3, 5, 7, not proved.
/// Recomputed live by <see cref="GreyRoomLadderImagesWitness"/> (<c>inspect --root greyroom</c>); the sober
/// twin, the count alone, is <c>MirrorWorld.GreyRoom</c>.
///
/// <para><b>Words, fenced at the door.</b> A ROOM here is a subspace the Hamiltonian part never leaves on
/// which the dissipator is one number, not <see cref="FrozenDivisorClaim"/>'s even and odd rooms of an
/// involution. A CORNER is a block with an empty or a full side, F153's pinned block, where a ladder image fills
/// its block; NOT <see cref="FrozenDivisorClaim"/>'s four corners p, q ∈ {1, N − 1}, nor
/// <see cref="RingFrozenBandClaim"/>'s corner, the block (1, 1), which this claim calls the block (1, 1). The
/// LADDERS are PROOF_FROZEN_BAND_SO4's two, and a HEIGHT is the
/// letter cube's k_Z, the Hamming distance of a cell. The WINDOW is <see cref="PinnedBlockFloorClaim"/>'s
/// disagreement window, listed.</para>
///
/// <para><b>The chain of reasons, in four steps.</b> Under equal light a cell |a⟩⟨b| is charged −2γ times
/// its Hamming distance and moved nowhere (the Absorption cell cost, <see cref="AbsorptionTheoremClaim"/>
/// through <see cref="PinnedBlockFloorClaim"/>), so the dissipator is one number on each Hamming plane of a
/// block, and a subspace the Hamiltonian part keeps inside one plane is grey (PROOF_GREY_GLASS_LADDER_IMAGES
/// Lemma 2a and 2b). A map that commutes with both parts carries grey rooms onto grey rooms (its Lemma 2c).
/// Neither ladder moves a disagreement set, so both commute with the dissipator (PROOF_FROZEN_BAND_SO4
/// Lemma 2.1); Φ commutes with the Hamiltonian part for every single-particle hopping
/// (<see cref="SpectatorIntertwinerClaim"/>, F125: W is Φ) and S⁺, hence its adjoint S⁻, for a real symmetric
/// hopping exactly when ΣhΣ = −h, which nearest-neighbour hopping on a chain always meets
/// (PROOF_FROZEN_BAND_SO4 Lemma 2.3). The boundary block (k, 0) is grey as a whole (F153, the pinning
/// criterion) and is a highest weight of the spin ladder and a lowest weight of the η ladder; the two ladders
/// commute (the two SU(2)s of the doubled picture), so S⁻^s Φ^r in either order is one-to-one on it for
/// s ≤ k and r ≤ N − k, and with s = (k − p + q)/2 and r = q − s it lands in (p, q) exactly when k lies in the
/// block's window, |p − q| ≤ k ≤ min(p + q, 2N − p − q) in steps of two
/// (<see cref="PinnedBlockFloorClaim.DisagreementWindow"/>). The image is then a grey room of dimension
/// C(N, k) at −2γk.</para>
///
/// <para><b>The block (1, 1), and why nothing more is grey there</b> (PROOF_GREY_GLASS_LADDER_IMAGES §4).
/// The room at −4γ is the sum over the frequencies of ad_h of the kernel of the diagonal map on each
/// eigenspace. The involution U v_(a,b) = (−1)^(a+b) v_(b,a) of PROOF_R90_FROZEN_DIVISOR §9.1, typed in
/// <see cref="FrozenDivisorClaim"/>, commutes with ad_h and keeps every diagonal entry, so it keeps each
/// eigenspace, and the diagonal map kills U's odd part. On the even part the uniform chain's spectrum inversion
/// Σψ_a = ±ψ_{N+1−a} (<see cref="ChiralKClaim"/>) turns each eigenspace into a forest of orbits whose diagonals
/// are differences of type-I cosine vectors, and §4 shows those differences independent, so the diagonal map
/// is one-to-one on the even part. The room at 0 is the identity alone. So C(N, 2) + 1, at every N, on the
/// uniform open chain.</para>
///
/// <para><b>The rings.</b> On p particles the spin ring with wrap w is the fermion ring with wrap sign
/// w·(−1)^(p−1) (PROOF_RING_FROZEN_BAND Lemma 1, typed in <see cref="RingFrozenBandClaim"/> as
/// <see cref="RingFrozenBandClaim.SectorWrapSign"/>), so a block's two sides see one flux when p ≡ q (even k)
/// and opposite fluxes at odd k. On the even ring Σh_fΣ = −h_f for either flux, so both ladders need one flux
/// on both sides, which the target has at even k; at odd k the target's two sides disagree unless one is
/// empty or full, and off those corners the Hamiltonian part differs from the hopping that keeps the image by
/// J_w times the bra's hop across the wrap bond alone, a coupling-free map: the image is kept exactly when
/// that map keeps it, which at N = 4, 6, 8 it does not (the proof's §6, G12). On the odd ring Σh_fΣ = −h_{−f}:
/// Φ needs one flux on both sides and S⁻ opposite ones, so no flux serves both and the argument joins the two
/// ladders only at a corner a single ladder fills: kept when (k even and s ∈ {0, k}) or (k odd and
/// r ∈ {0, N − k}); that every other image is lost is measured at N = 3, 5, 7, at the uniform ring and at
/// couplings of both signs (the proof's §7, G13), not proved. The η ladder's half of the ring rule, W carrying
/// a ring block (p, q) with p ≡ q onto the π-flux ring's (p + 1, q + 1), is gated apart in
/// <c>SpectatorRingParityGateTests</c>. <see cref="ImageProvedKept"/> is this rule, and it says where F171
/// PROVES the image kept: false is "not proved kept", a loss where the proof measures one and an open row
/// beyond.</para>
///
/// <para><b>Typed parents, and what each gives.</b> <see cref="SpectatorIntertwinerClaim"/> supplies the η
/// ladder Φ with its exact intertwining of both parts for any quadratic particle-conserving H.
/// <see cref="PinnedBlockFloorClaim"/> supplies the starts, the boundary blocks (k, 0) grey as wholes, and the
/// window combinatorics the landing rule is read in; its parent the Absorption Theorem is the cell cost behind
/// the lemma. <see cref="ChiralKClaim"/> supplies the spectrum inversion the (1, 1) count rests on.
/// <see cref="FrozenDivisorClaim"/> supplies the involution U whose odd half is the distance-2 image in (1, 1),
/// and its ⌊N/2⌋ floor is the count the Heisenberg (1, 1) room meets, adjacent and not docked.
/// <see cref="RingFrozenBandClaim"/> supplies Lemma 1, the flux the particle number sets, on which the whole
/// ring rule rests; it is typed as an odd-N statement about the frozen band at uniform bonds, and its Lemma 1 holds
/// at every N and every coupling,
/// an argument this claim carries as that claim carried F140's floor. None alone gives the statement: the first
/// has the ladder without the rooms, the second the rooms without the ladder, the others one block or one
/// ring each.</para>
///
/// <para><b>What is not typed here.</b> The spin ladder's intertwining half, PROOF_FROZEN_BAND_SO4 Lemma 2.3,
/// carries no typed Tier-1-derived claim: <see cref="SidewaysSpinLadderClaim"/> holds it, derived, beside a
/// measured multiplet half and is tiered Candidate for that half, so it is adjacent and not a parent (a
/// parent's tier caps the child's), though every image with s &gt; 0 rests on it. The mirror products at the
/// distance N/2, the half-chain reading and its uniqueness (the proof's §5), the antipodal map and the
/// sum-of-full-pairs criterion (§8) are proved there and not typed; what they add to the ladder images is
/// measured (8 against 6 at N = 4, 22 against 20 at N = 6) and completeness beyond them is measured, not
/// proved. The typed Krylov room of F157's blind seat (<see cref="SeatCutBlindnessClaim"/>) and the untyped
/// node-pair room m² + 1 (PROOF_NODE_PAIR_RESOLVENT §8; <c>NodePairResolventClaim</c> types that proof's bond
/// theorems) are rooms at μ = 0 under one lit seat; this claim reads rooms at every height under equal light.
/// The ladder images of the one-sided sectors were named before their count by <c>OneSidedSectorClosedForm</c>,
/// <c>BlockSpectrumOpenQuestions</c> and experiments/N10_BLOCK_SPECTRUM_PUSH.md, which point here. The
/// Heisenberg chain is outside: on a connected chain its ZZ term loses both ladders' images from N = 3 on (the
/// one-particle ZZ diagonal differs at the ends; a read, and with zero bonds making equal dimers the diagonal is
/// constant and the spin image survives), and its (1, 1) room is the ⌊N/2⌋ frozen modes at N = 3 to 8, measured,
/// and certified at N = 3 to 6 by the witness, whose mod-p room meets F140's floor there: under uniform light a
/// frozen mode at exactly −4γ lies on the plane of distance 2 in the kernel of the Hamiltonian part, so F140's
/// ⌊N/2⌋ at every reversal-invariant coupling is a lower bound on that room.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F171 + <c>docs/proofs/PROOF_GREY_GLASS_LADDER_IMAGES.md</c>
/// + <c>docs/THE_GREY_GLASS.md</c> (the page, the measured rooms, the open questions) +
/// <c>simulations/grey_glass_gate.py</c> (G7 to G13) + <c>simulations/grey_glass_ring_rooms.py</c> (R1 to R4).
/// </para></summary>
public sealed class GreyRoomLadderImagesClaim : Claim
{
    /// <summary>The η ladder Φ, F125: W(ρ) = Σ_l c_l†ρc_l intertwines both parts of L for any quadratic
    /// particle-conserving H.</summary>
    public SpectatorIntertwinerClaim Spectator { get; }

    /// <summary>The starts: the 4N boundary blocks are grey as wholes, F153, and its window combinatorics.</summary>
    public PinnedBlockFloorClaim Pinned { get; }

    /// <summary>The spectrum inversion Σψ_a = ±ψ_{N+1−a} of the uniform chain, behind the (1, 1) count.</summary>
    public ChiralKClaim ChiralK { get; }

    /// <summary>The involution U of PROOF_R90_FROZEN_DIVISOR §9.1, F140, whose odd half is the distance-2 image
    /// in (1, 1).</summary>
    public FrozenDivisorClaim FrozenDivisor { get; }

    /// <summary>The flux the particle number sets on the spin ring, PROOF_RING_FROZEN_BAND Lemma 1, F170.</summary>
    public RingFrozenBandClaim RingBand { get; }

    /// <summary>C(N, k), the dimension of the boundary block (k, 0) and of each of its ladder images.</summary>
    public static int LadderDimension(int n, int k)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), n, "N ≥ 0");
        if (k < 0 || k > n) return 0;
        long c = 1;
        for (int i = 1; i <= k; i++) c = c * (n - k + i) / i;
        return checked((int)c);
    }

    /// <summary>The Hamming distances the block (p, q) holds, |p − q| ≤ k ≤ min(p + q, 2N − p − q) in steps of
    /// two: F153's disagreement window, listed.</summary>
    public static IReadOnlyList<int> Window(int n, int p, int q)
    {
        var (min, max) = PinnedBlockFloorClaim.DisagreementWindow(n, p, q);
        var w = new List<int>();
        for (int k = min; k <= max; k += 2) w.Add(k);
        return w;
    }

    /// <summary>A block with an empty or a full side, F153's pinned block, at one Hamming distance: where a
    /// ladder image fills its block.</summary>
    public static bool IsCorner(int n, int p, int q) => PinnedBlockFloorClaim.IsPinned(n, p, q);

    /// <summary>The steps from (k, 0) to (p, q): s = (k − p + q)/2 down the spin ladder and r = q − s up the η
    /// ladder; both lie in [0, k] and [0, N − k] exactly when k is in the window, and the call throws off it.
    /// Symmetric under the transpose in what the rule below reads (s ↔ k − s leaves {0, k} alone, r is
    /// unchanged).</summary>
    public static (int Spin, int Eta) LadderSteps(int n, int p, int q, int k)
    {
        if (!Window(n, p, q).Contains(k))
            throw new ArgumentOutOfRangeException(nameof(k), k, $"the distance {k} is not in the window of ({p}, {q}) at N = {n}");
        int s = (k - p + q) / 2;
        return (s, q - s);
    }

    /// <summary>Where F171 PROVES the ladder image of (k, 0) kept by the Hamiltonian part in the block (p, q):
    /// always on the chain, at any real couplings; a ring whose wrap bond is zero is the chain; on the even
    /// ring at even k, and at odd k only where the image fills a corner; on the odd ring when (k even and
    /// s ∈ {0, k}) or (k odd and r ∈ {0, N − k}). False is NOT PROVED KEPT: a loss where the proof measures one
    /// (the even ring at N = 4, 6, 8 and (5, 2) at N = 10, exact; the odd ring at N = 3, 5, 7, at the uniform
    /// ring and couplings of both signs, measured), an open row beyond. An XY statement: on the Heisenberg book
    /// nothing is claimed, and this rule is not asked.</summary>
    public static bool ImageProvedKept(int n, int p, int q, int k, bool ring, bool wrapBondNonzero)
    {
        if (!Window(n, p, q).Contains(k)) return false;
        if (!ring || !wrapBondNonzero) return true;
        var (s, r) = LadderSteps(n, p, q, k);
        if (n % 2 == 0) return k % 2 == 0 || IsCorner(n, p, q);
        return k % 2 == 0 ? (s == 0 || s == k) : (r == 0 || r == n - k);
    }

    /// <summary>The proved lower bound on the grey room at −2γk in (p, q): C(N, k) where the image is proved kept.</summary>
    public static int LadderBound(int n, int p, int q, int k, bool ring, bool wrapBondNonzero) =>
        ImageProvedKept(n, p, q, k, ring, wrapBondNonzero) ? LadderDimension(n, k) : 0;

    /// <summary>The whole grey content of the block (1, 1) of the uniform open chain: C(N, 2) at the distance 2
    /// and the identity at 0, nothing else, at every N ≥ 2.</summary>
    public static int UniformOpenChainBlock11Room(int n) =>
        n >= 2 ? LadderDimension(n, 2) + 1
            : throw new ArgumentOutOfRangeException(nameof(n), n, "the block (1, 1) needs N ≥ 2");

    public GreyRoomLadderImagesClaim(SpectatorIntertwinerClaim spectator, PinnedBlockFloorClaim pinned, ChiralKClaim chiralK,
        FrozenDivisorClaim frozenDivisor, RingFrozenBandClaim ringBand)
        : base(
            "F171: under Z dephasing at one rate γ on every site, the η ladder Φ and the spin ladder S⁻ carry the " +
            "boundary block (k, 0) of the XY chain, at any real couplings, one-to-one into every joint-popcount " +
            "block (p, q) whose Hamming window |p − q| ≤ k ≤ min(p + q, 2N − p − q) holds k, onto a subspace of " +
            "dimension C(N, k) kept by the Hamiltonian part on which the dissipator is −2γk: a grey room, " +
            "e^{tL} = e^{−2γkt}·e^{tL_H} on it, its modes not moving with γ. In the block (1, 1) of the uniform " +
            "open chain the two images are all there is, C(N, 2) + 1. On the even ring the image is kept at every " +
            "even k, and at odd k, with the wrap bond nonzero, exactly when the bra's hop across the wrap bond alone " +
            "keeps it, a coupling-free map, which at N = 4, 6, 8 it does only where the image fills a corner, and for " +
            "(5, 2) at N = 10 not at all (the losses exact there). On the odd ring the image is kept when it reaches " +
            "its block from a corner by the ladder of its parity of k, s ∈ {0, k} at even k and r ∈ {0, N − k} at odd " +
            "k; that the rest is lost is measured at N = 3, 5, 7, at the uniform ring and at couplings of both signs, " +
            "not proved. SCOPE: XY couplings, real; equal light (with positive rates not all equal only the four single " +
            "cells stay grey among whole blocks); the even ring's losses are proved at N = 4, 6, 8 and for (5, 2) at " +
            "N = 10 and claimed at no other N, the odd ring's are measured and proved at none; completeness beyond " +
            "the ladder images and the reflection's mirror products at the distance N/2 is measured, not proved.",
            Tier.Tier1Derived,
            "docs/ANALYTICAL_FORMULAS.md F171 + docs/proofs/PROOF_GREY_GLASS_LADDER_IMAGES.md + docs/THE_GREY_GLASS.md + " +
            "simulations/grey_glass_gate.py + simulations/grey_glass_ring_rooms.py + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/GreyRoomLadderImagesWitness.cs (inspect --root greyroom)")
    {
        Spectator = spectator ?? throw new ArgumentNullException(nameof(spectator));
        Pinned = pinned ?? throw new ArgumentNullException(nameof(pinned));
        ChiralK = chiralK ?? throw new ArgumentNullException(nameof(chiralK));
        FrozenDivisor = frozenDivisor ?? throw new ArgumentNullException(nameof(frozenDivisor));
        RingBand = ringBand ?? throw new ArgumentNullException(nameof(ringBand));
    }

    public override string DisplayName =>
        "F171: the ladder images of the boundary blocks are grey rooms in every coloured block";

    public override string Summary =>
        "under equal light the dissipator is one number per Hamming plane, both ladders commute with both parts " +
        "on the XY chain, and (k, 0) is a highest weight of one ladder and a lowest weight of the other: every " +
        "block holds, at each distance k of its window, a grey room of dimension at least C(N, k), exactly " +
        "C(N, 2) + 1 in all in the (1, 1) of the uniform open chain; on the rings the flux the particle number " +
        "sets decides which images are proved kept, the rest lost where the proof computes or measures it; live as " +
        $"exact intertwining identities and GF(p) ranks ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the four steps",
                summary: "the light charges a cell its Hamming distance and moves nothing, so it is one number per " +
                         "plane and a kept subspace inside a plane is grey; a map commuting with both parts carries " +
                         "grey rooms onto grey rooms; Φ and S⁻ commute with the dissipator (no disagreement set moves) " +
                         "and, on the XY chain, with the Hamiltonian part (Φ for any hopping, S⁻ when ΣhΣ = −h), and " +
                         "with each other; (k, 0) is grey as a whole and a highest weight of the spin ladder, a lowest " +
                         "weight of the η ladder, so S⁻^s Φ^r is one-to-one on it exactly across the window");
            yield return new InspectableNode("the rule on the rings, written out",
                summary: "the spin ring on p particles is the fermion ring with wrap sign w·(−1)^(p−1), so a block's " +
                         "sides share a flux at even k and differ at odd k; even ring: Σh_fΣ = −h_f, both ladders want " +
                         "one flux, kept at even k, at odd k exactly when the bra's hop across the wrap bond keeps it, " +
                         "which at N = 4, 6, 8 it does only in a corner (exact; (5, 2) at N = 10 lost); odd ring: " +
                         "Σh_fΣ = −h_{−f}, Φ wants one flux and S⁻ opposite ones, kept when it comes from a corner by " +
                         "the ladder of k's parity, the rest measured lost at N = 3, 5, 7; e.g. at N = 5 the image of " +
                         $"(2, 0) in (1, 1) is {(ImageProvedKept(5, 1, 1, 2, ring: true, wrapBondNonzero: true) ? "proved kept" : "not proved kept, and measured lost")}");
            yield return new InspectableNode($"the block (1, 1) of the uniform open chain: C(N, 2) + 1 = {UniformOpenChainBlock11Room(5)} at N = 5",
                summary: "the room at −4γ is, per frequency of ad_h, the kernel of the diagonal map; the diagonal map " +
                         "kills U's odd part, and on the even part the spectrum inversion makes each eigenspace a forest " +
                         "of orbits whose diagonals are differences of type-I cosine vectors, independent (the proof's " +
                         "§4); the room at 0 is the identity alone; so the two ladder images and nothing else, at every N");
            yield return new InspectableNode("what is typed and what stays outside",
                summary: "typed here: the images, their count, their keeping on the chain and the rings' rule; the " +
                         "spin ladder's Lemma 2.3 carries no typed Tier-1-derived claim (SidewaysSpinLadderClaim holds " +
                         "it beside a measured multiplet half, tiered Candidate, adjacent); the mirror products at N/2, " +
                         "the half-chain reading, the antipodal pairs and the sum-of-full-pairs criterion are proved in " +
                         "the source and not typed; the rooms above the images (8 against 6, 22 against 20) and the " +
                         "rings' losses beyond the N named are reads; the Heisenberg chain is outside (its (1, 1) reads " +
                         "⌊N/2⌋ at N = 3 to 8, the frozen modes). The sober twin: MirrorWorld.GreyRoom, the count " +
                         "alone, run mode greyroom N [ring]");
            yield return Spectator;
            yield return Pinned;
            yield return ChiralK;
            yield return FrozenDivisor;
            yield return RingBand;
        }
    }
}
