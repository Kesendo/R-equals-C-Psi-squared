using System;
using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The empirical Q=20 instance refinement of <see cref="DefectReadingEquivarianceClaim"/> (the
/// structural parent). The decoder projects an observed reading onto the PAINTED per-site α-dictionary
/// {F_b} (painted = the per-site profile read off a PROPAGATED run at Q=20, as opposed to the BARE algebraic
/// dictionary; least-squares plain Euclidean cosine, distinct from the parent's mode-basis Gram metric).
///
/// <para>The chain mirror. Under R (site i ↦ N − 1 − i, bond b ↦ N − 2 − b) the chain, the Z watching and
/// the carrier ψ₁ are symmetric, so for the protocol the letter of bond N − 2 − b is R applied to the letter
/// of bond b, and the painters' reliability flags are mirror-consistent. Split a letter into its mirror-even
/// part and its mirror-odd part (seesaw-led: to first order in δJ the response to (V_b − V_{N−2−b})/2, which
/// collects every R-odd channel, the parent's k = 2 seesaw first), w(b) the mirror-odd share of ‖F_b‖²; then
/// cos(F_b, F_{N−2−b}) = 1 − 2·w(b) for every chain-mirror pair, and at even N the centre bond's letter is
/// mirror-even, w = 0. The identity is exact for the protocol; the code's α fit realises its premise on the
/// reliable sites to two stopping widths (|F_{N−2−b} − R·F_b| ≤ 2·1e-7/δJ = 1e-5). In the bare mode basis the mirror-odd share is
/// the even-k mode weight, so 1 − 2w is the parent's Σ(−1)^{k−1}·w_k read on the painted dictionary.</para>
///
/// <para>The reading, over the sites reliable for every bond (the painters' rule: |f| ≤ 10 and the guard at
/// δJ/2 agreeing; a mirror-symmetric mask, so w and the identity are taken over it): the worst signed cos is
/// −0.976, −0.541, −0.965, −0.962, −0.940, −0.927, −0.992, −0.950, −0.914, −0.961 at N = 3..12. The worst pair
/// is a chain-mirror pair at N = 3, 6, 9, 11, 12 and a pair across the centre (a &lt; (N − 2)/2 &lt; b, not a
/// mirror pair) at N = 5, 7, 8, 10; at N = 5 that is (0,2) and (1,3), the pair the decoder demo plants
/// ("bond 3 weakened reads like bond 1 strengthened"), while N = 5's painted chain-mirror pairs read −0.783 and
/// −0.894. At N = 4 the three bonds hold one mirror pair, (0,2) with w = 0.707 and cos = −0.414, and two pairs
/// with the mirror-even centre; the worst is (0,1) at −0.541. Anti-collinear letters are the rule at Q = 20, not
/// a property of odd N. Over all sites N = 6 reads −0.378 and N = 10 −0.701; those values are carried by sites
/// the painters' rule rejects (|f| up to 216 at N = 6, 42 at N = 10), and what the rule catches there is a far
/// minimum of the α fit at δJ = 0.02, not an unreadable site: the exactly stretched clean purity has the same far
/// minima (α = 5.33 and 0.155 on N = 6 bond 1, sites 4 and 1) and loses them at δJ = 0.01. The linear letter (the δJ → 0 limit of the
/// exactly stretched fit, a projection with no fit at all; the painters' own fit, which clamps beyond the grid,
/// meets it to about 1 % at small δJ, and at the protocol's δJ = 0.02 only to 6 % at N = 4 and 17 % at N = 5) is ordinary on those sites and anti-collinear over all sites, −0.965, −0.536, −0.983,
/// −0.973, −0.928, −0.960, −0.968, −0.931, −0.905, −0.915 at N = 3..12. The bare single-excitation dictionary
/// reads −0.52..−0.74 at N = 4..9 (N = 3 is the −1.0 single-pair case).</para>
///
/// <para>The carrier channels, on the linear letter. By linearity the linear letter is a sum over mode pairs; the
/// pairs that touch the carrier give Σ_k c_k(b)·u_k, with c_k(b) = ⟨ψ_k|V_b|ψ₁⟩ in closed form (k = 2..N F124's
/// location dictionary M[b, k], k = 1 the strength) and u_k the letter of the channel (1, k), fixed per N and
/// mirror-odd for even k. The u_k are the first-principles counterpart of the R_k the handshake_decoder arc fitted;
/// their sizes stay without a closed form (the handshake_Rk facet, IsDeadEnd, stands). The carrier part leaves a
/// remainder of 8 % to 17 % in norm over the dictionary at N = 4..12 (single edge letters up to 42 % at N = 12). A
/// carrier energy shift, the channel (1, 1), commutes with the carrier and acts only on what the watching has
/// already moved out of it (|u₁| 0.2..0.8 against |u₂| 4..39). On the carrier part the seesaw is the even-k sum and
/// the mirror-even half the odd-k sum. The linear letter's w is the balance of its two halves, and no single
/// criterion sets it: ‖even‖/‖odd‖ runs on a continuum, 0.12 at N = 6 bond 1 and 0.21 at N = 12 bond 3, both with c₃
/// near its node (at (b + 3/2)/(N + 1) ≈ 1/3), and 0.31 at N = 8 bond 2, where |c₃| = 1/6 and bond 1, with the same
/// |c₃|, reads 0.39. What separates N = 8's two bonds is spread over the channels and the remainder and is not
/// read. The carrier model reproduces the linear letter's w within 0.06 at N = 5..10 and within 0.1 at N = 4 and
/// N = 11..13; the painted letters at δJ = 0.02 it does not reproduce (misses up to 0.42), the painting departing from
/// the linear letter there.</para>
///
/// <para>Anchors: hypotheses/HANDSHAKE_GEOMETRY.md + SingleExcitationPainterDictionary.cs (the dictionary in
/// the single-excitation block; its Python predecessor is simulations/handshake_rk_block.py) + DefectDecoder.cs
/// (DictionaryWorstAntiCos, the live quantity, read over all sites, which is the reliable-site reading at the
/// decoder's N ≤ 5 since every site is reliable there) + ReadingPowerWitness.cs (the live node at the inspect
/// root decoder). Verifiers: SingleExcitationPainterDictionaryTests (the handshake with the full pipeline at
/// N = 3..5, the premise, the identity, the table to N = 12) and DictionaryAntiCollinearityInvestigationTests
/// (the full pipeline at N = 3..6).</para></summary>
public sealed class DecoderAntiCollinearitySeesawClaim : Claim
{
    private readonly DefectReadingEquivarianceClaim _equivariance;

    public DecoderAntiCollinearitySeesawClaim(DefectReadingEquivarianceClaim equivariance)
        : base("Decoder α-dictionary anti-collinearity is the seesaw's (a Q=20 painting effect): for every chain-mirror " +
               "pair cos(F_b, F_{N−2−b}) = 1 − 2·w(b), w the letter's mirror-odd weight, and the centre bond at even N is " +
               "mirror-even; over the sites reliable for every bond the worst signed cos is ≤ −0.91 at N = 3 and 5..12 and " +
               "−0.541 at N = 4, the worst pair a chain-mirror pair or a pair across the centre: the empirical Q=20 instance " +
               "of the structural DefectReadingEquivarianceClaim",
               Tier.Tier2Empirical,
               "hypotheses/HANDSHAKE_GEOMETRY.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SingleExcitationPainterDictionary.cs (the block dictionary) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (DictionaryWorstAntiCos, the live quantity) + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/SingleExcitationPainterDictionaryTests.cs (the verifier)")
    {
        _equivariance = equivariance ?? throw new ArgumentNullException(nameof(equivariance));
    }

    public string MirrorPairIdentity =>
        "The chain, the Z watching and ψ₁ are symmetric under R (i ↦ N − 1 − i), so the letter of bond N − 2 − b is R " +
        "applied to the letter of bond b (exact for the protocol; the α fit realises it to 1e-5 on the reliable sites). Splitting a letter into " +
        "its mirror-even and mirror-odd (seesaw-led) parts, cos(F_b, F_{N−2−b}) = (‖even‖² − ‖odd‖²)/‖F_b‖² = 1 − 2·w(b). " +
        "At even N the centre bond is its own mirror: w = 0. In the bare mode basis 1 − 2w is the parent's Σ(−1)^{k−1}w_k.";

    public string Reading =>
        "Worst signed cos over the sites reliable for every bond (γ = 0.05, Q = 20): N=3 −0.976, N=4 −0.541, N=5 −0.965, " +
        "N=6 −0.962, N=7 −0.940, N=8 −0.927, N=9 −0.992, N=10 −0.950, N=11 −0.914, N=12 −0.961; a chain-mirror pair at " +
        "N = 3, 6, 9, 11, 12 (e.g. N=12 bond 3, w = 0.98, 1 − 2w = −0.96), a pair across the centre at N = 5, 7, 8, 10 " +
        "(N = 5: (0,2) and (1,3), while its chain-mirror pairs read −0.783 and −0.894); at N = 4 the only mirror pair " +
        "reads −0.414 and the worst is (0,1) with the mirror-even centre.";

    public string RejectedSites =>
        "Over all sites N = 6 reads −0.378 and N = 10 −0.701; those values are carried by sites the painters' reliability " +
        "rule rejects (|f| up to 216 at N = 6, 42 at N = 10), where the α fit at δJ = 0.02 has a far minimum, the exactly " +
        "stretched purity included. The linear letter, the δJ → 0 limit with no fit, is ordinary there (1.05 and −0.97 on N = 6 bond 1) and anti-collinear " +
        "over all sites at every N = 3, 5..12 (−0.905 to −0.983; −0.536 at N = 4). The bare single-excitation dictionary " +
        "reads −0.52..−0.74 at N = 4..9 (N = 3 is the −1.0 single-pair case).";

    public string CarrierChannels =>
        "By linearity the linear letter is a sum over mode pairs; the pairs that touch ψ₁ give Σ_k c_k(b)·u_k, c_k(b) = " +
        "⟨ψ_k|V_b|ψ₁⟩ in closed form (k = 2..N F124's M[b, k], k = 1 the strength) and u_k the fixed letter of the " +
        "channel (1, k), mirror-odd for even k, the first-principles counterpart of the handshake_decoder R_k; the carrier " +
        "part leaves 8–17 % in norm over the dictionary at N = 4..12. The linear letter's w is the balance of its halves " +
        "with no single criterion: ‖even‖/‖odd‖ is 0.12 at N = 6 bond 1 and 0.21 at N = 12 bond 3, c₃ near its node in " +
        "both, and 0.31 against 0.39 at N = 8's bonds 2 and 1, which share |c₃| = 1/6; what separates those two is not " +
        "read. The model reproduces the linear letter's w within 0.06 at N = 5..10, within 0.1 at N = 4 and N = 11..13; " +
        "the sizes of the u_k are computed, with no closed form.";

    public override string DisplayName =>
        "Decoder α-dictionary anti-collinearity: the seesaw's, cos = 1 − 2w on chain-mirror pairs (Q=20 painting, Tier2Empirical)";

    public override string Summary =>
        "The decoder's painted α-dictionary is anti-collinear at every N read (worst signed cos ≤ −0.91 at N = 3, 5..12 over " +
        "the sites reliable for every bond, −0.541 at N = 4); for every chain-mirror pair cos(F_b, F_{N−2−b}) = 1 − 2w(b), w the " +
        "letter's mirror-odd weight, and the centre bond at even N is mirror-even; the worst pair is a chain-mirror pair or a " +
        "pair across the centre; which bonds carry a nearly pure seesaw has no single criterion in the carrier channels " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The chain-mirror identity, cos = 1 − 2w", summary: MirrorPairIdentity);
            yield return new InspectableNode("The reading over reliable sites, N = 3..12", summary: Reading);
            yield return new InspectableNode("The rejected sites, the linear letter and the bare dictionary", summary: RejectedSites);
            yield return new InspectableNode("The carrier channels, F_b ≈ Σ c_k(b)·u_k", summary: CarrierChannels);
            yield return _equivariance; // the typed parent edge (the structural cousin)
        }
    }
}
