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
/// the painters' rule rejects (|f| up to 216 at N = 6, 42 at N = 10). The bare single-excitation dictionary
/// reads −0.52..−0.74 at N = 4..9 (N = 3 is the −1.0 single-pair case). Which bonds carry a nearly pure seesaw is
/// not derived: a mode picture (the defect couples ψ₁ to ψ₃ with a strength that samples ψ₃ at the bond; w ≈ 1
/// where that coupling is small) places it at N = 6, 9 and 12, and N = 8's bond 2, w = 0.94 off any node, keeps
/// it a sighting.</para>
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
        "rule rejects (|f| up to 216 at N = 6, 42 at N = 10). The bare single-excitation dictionary reads −0.52..−0.74 at " +
        "N = 4..9 (N = 3 is the −1.0 single-pair case).";

    public override string DisplayName =>
        "Decoder α-dictionary anti-collinearity: the seesaw's, cos = 1 − 2w on chain-mirror pairs (Q=20 painting, Tier2Empirical)";

    public override string Summary =>
        "The decoder's painted α-dictionary is anti-collinear at every N read (worst signed cos ≤ −0.91 at N = 3, 5..12 over " +
        "the sites reliable for every bond, −0.541 at N = 4); for every chain-mirror pair cos(F_b, F_{N−2−b}) = 1 − 2w(b), w the " +
        "letter's mirror-odd weight, and the centre bond at even N is mirror-even; the worst pair is a chain-mirror pair or a " +
        "pair across the centre; which bonds carry a nearly pure seesaw is open " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The chain-mirror identity, cos = 1 − 2w", summary: MirrorPairIdentity);
            yield return new InspectableNode("The reading over reliable sites, N = 3..12", summary: Reading);
            yield return new InspectableNode("The rejected sites and the bare dictionary", summary: RejectedSites);
            yield return _equivariance; // the typed parent edge (the structural cousin)
        }
    }
}
