using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>A field conditioned on a cooled bit, whose spectrum pairs where none of F137's carriers applies
/// (<c>docs/ANALYTICAL_FORMULAS.md</c>, F137, the paragraph "At unequal rates the arrow can cancel a field").
///
/// <para>The setting: two qubits, H = h·(X₀ + X₀Z₁) + w·Z₁ = 2h·X₀·|0⟩⟨0|₁ + w·Z₁, site 0 cooled at a and heated at
/// b, s₀ = (a + b)/2, site 1 cooled alone at g. H conserves Z₁ and site 1 moves only from |1⟩ to |0⟩, so L is block
/// triangular in site 1's ket-bra sectors: only site 1's jump leaves a sector, and only from (1, 1) to (0, 0).
/// The centre is −(s₀ + g/2). The field-free block (1, 1) holds −g, −s₀ − g twice and −2s₀ − g; the field's block
/// (0, 0) holds 0, −s₀ and a pair with mean −3s₀/2. Its 0 and −s₀ pair with −2s₀ − g and −s₀ − g at every rate; the
/// remaining −g and −s₀ − g pair with each other, and the field's pair sits about the centre, both exactly when
/// s₀ = g; the one-sided blocks (0, 1) and (1, 0) pair with each other at every rate, h and w. The third moment of
/// <see cref="ThermalThirdMomentClaim"/> reads the converse: Tr(L + c/2) = 0, so by Newton's identities the coefficient of x¹³ in the
/// characteristic polynomial of L + c/2 is −Tr((L + c/2)³)/3 = −8h²·(a + b − 2g), so the whole spectrum pairs exactly
/// when s₀ = g, at every h ≠ 0 and any w. The field term 48·s₀·h² and the arrow term −48·g·h² cancel there; the
/// arrow is the second rung t₂(1) = Tr(Z₁·H²) = 8h² of the pump channel.</para>
///
/// <para>H does not commute with Z⊗N, no map of the cube's one-site closure applied alike on every site
/// reflects it at h = 1, w = 0, site 0 at 1/1 and site 1 cooled at 1 (the sextant's search there, cited from MirrorWorld and not recomputed by
/// this claim's witness), and none of F137's carriers applies. Mirrored by conjugation with X₁ (X₀ − X₀Z₁, site 1 heated
/// alone), turned to Y by a quarter turn of site 0 about Z, which the bath keeps. The setting is narrow: at the rates where the controls try them,
/// with the third moment still cancelled, heating on site 1, the field left on in the other sector, a ZZ term, or a
/// second neighbour conditioned on site 1 break at the fifth sum.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/ConditionedFieldPairingWitness.cs</c>,
/// <c>inspect --root conditionedfield</c>, which proves the block identities on a grid of parameter points one wider
/// than their degree (a polynomial of degree at most four in each of five variables that vanishes on such a grid is
/// zero), beside the registry's gate <c>simulations/thermal_third_moment.py</c>, which proves them with h, w and the
/// rates as symbols. Parent: <see cref="ThermalThirdMomentClaim"/>, whose law at this H is the x¹³ coefficient.
/// Adjacent: <see cref="MomentTowerPumpChannelClaim"/> (whose battery takes X₀ + X₀Z₁ as its rung-2 case),
/// <see cref="PalindromeComplementConnectionClaim"/> (Theorem 4 conditions on the bits of undephased sites, the same
/// blocks; here the conditioning site is damped).</para></summary>
public sealed class ConditionedFieldPairingClaim : Claim
{
    public ThermalThirdMomentClaim Moment { get; }

    /// <summary>The coefficient of y¹³ in det(y − (2L + c)), −Tr((2L + c)³)/3 by Newton's identities (the trace of 2L + c
    /// is zero), read off the parent's law at this H: −64h²·(a + b − 2g), eight times the coefficient of x¹³ in
    /// det(x − (L + c/2)), which is −8h²·(a + b − 2g). The frequency w does not enter.</summary>
    public static BigInteger ThirteenthCoefficientOfTwoLPlusC(long h, long a, long b, long g)
    {
        var cubed = ThermalThirdMomentClaim.CubedCentredTrace(
            new[] { new ThermalThirdMomentClaim.Term("XI", h), new ThermalThirdMomentClaim.Term("XZ", h) },
            new[] { a, g }, new[] { b, 0L });
        return -cubed / 3;   // exact: the law carries the factor 12
    }

    public ConditionedFieldPairingClaim(ThermalThirdMomentClaim moment)
        : base("H = h(X₀ + X₀Z₁) + w·Z₁, site 0 cooled at a and heated at b, site 1 cooled alone at g: L is block triangular in " +
               "site 1's ket-bra sectors, and the spectrum pairs about −c/2 exactly when (a + b)/2 = g, at every h ≠ 0 and any w, " +
               "none of F137's carriers applying",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md (F137, the conditioned field, at unequal rates the arrow can cancel a field)")
    {
        Moment = moment ?? throw new ArgumentNullException(nameof(moment));
    }

    public override string DisplayName =>
        "Conditioned field: X₀(1 + Z₁) pairs exactly when s₀ = γ↓₁ at h ≠ 0, outside F137's carriers";

    public override string Summary =>
        "H = h(X₀ + X₀Z₁) + w·Z₁, site 0 at a/b, site 1 cooled alone at g: block triangular in site 1's ket-bra sectors, the field-free " +
        "block −g, −s₀ − g (twice), −2s₀ − g, the field's block 0, −s₀ and a pair of mean −3s₀/2, the one-sided blocks paired with each " +
        "other, so the spectrum pairs at s₀ = g; the coefficient of x¹³ is −8h²(a + b − 2g), so at h ≠ 0 only there (Tier 1 derived; parent " +
        "ThermalThirdMomentClaim; live witness: inspect --root conditionedfield).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the triangularity",
                summary: "H keeps site 1's bit on both sides and only site 1's jump changes it, from (1, 1) to (0, 0); every other entry of L " +
                         "keeps the sector (ket bit, bra bit) of site 1, so the characteristic polynomial is the product of the four blocks'");
            yield return new InspectableNode("the blocks",
                summary: "field-free (1, 1): −g, −s₀ − g twice, −2s₀ − g; the field's (0, 0): 0, −s₀ and a pair of mean −3s₀/2; about the " +
                         "centre −s₀ − g/2, 0 and −s₀ pair with −2s₀ − g and −s₀ − g at every rate, −g and the other −s₀ − g pair with each other " +
                         "and the field's pair sits about the centre, both exactly when s₀ = g; the one-sided (0, 1) and (1, 0) pair with each " +
                         "other at every a, b, g, h, w");
            yield return new InspectableNode("the converse",
                summary: "the x¹³ coefficient is −Tr((L + c/2)³)/3, and the parent's law gives Tr((L + c/2)³) = 48·s₀·h² − 48·g·h² at this H, " +
                         "so off s₀ = g with h ≠ 0 the third moment is not zero and the spectrum does not pair");
            yield return new InspectableNode("outside the known carriers",
                summary: "H does not commute with Z⊗N; no map of the one-site closure applied alike on every site reflects it at h = 1, w = 0, site 0 " +
                         "at 1/1 and site 1 cooled at 1 (cited: sextant decide 2 " +
                         "XI:1,XZ:1 1/1,1/0 loss+gain), and none of F137's carriers applies; the triangularity with the four blocks proves the " +
                         "pairing, the third moment its converse");
            yield return new InspectableNode("narrow",
                summary: "at the rates where the controls try them, with the third moment cancelled, heating on site 1, a ZZ term, the field left on in the other sector (2X₀ + X₀Z₁) or a " +
                         "second neighbour conditioned on site 1 break at the fifth sum; mirrored by X₁ (X₀ − X₀Z₁, site 1 heated alone) and turned " +
                         "to Y₀(1 + Z₁) it holds");
        }
    }
}
