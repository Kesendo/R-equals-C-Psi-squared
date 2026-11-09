using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The third centred moment of the generator under cooling and heating, closed for every H
/// (<c>docs/ANALYTICAL_FORMULAS.md</c>, F137, the paragraph "The third moment, for every H, under cooling and
/// heating").
///
/// <para>The setting: N qubits, H a sum of Pauli strings with real coefficients, cooling σ⁻ = |0⟩⟨1| at γ↓_l and
/// heating σ⁺ at γ↑_l on site l, nothing else; s_l = (γ↓_l + γ↑_l)/2 and c = Σ_l (γ↓_l + γ↑_l), so −c/2 is the
/// centre Tr L / 4^N. Then</para>
///
/// <para>  Tr((L + c/2)³) = 12·4^(N−1)·Σ_l s_l·(h_X,l² + h_Y,l²) − 3·2^(N−1)·Σ_l (γ↓_l − γ↑_l)·Tr(Z_l·H²),</para>
///
/// <para>h_X,l and h_Y,l the coefficients of X_l and Y_l in H, the trace read on the traceless part of H (a
/// constant does not enter L but does enter Tr(Z_l·H²)). L is read as written: a uniform Z field, which the frame
/// rotating with the qubits removes where H conserves the excitation number, still enters the arrow term through ZZ
/// terms (the −768 of the instances node). Write L + c/2 = A + B with A = −i[H, ·] and B the centred dissipator, the dissipator plus
/// c/2; Tr A³ = 0, Tr B³ = 0 and Tr(A·B²) = 0, so the moment is 3·Tr(A²·B). B's diagonal gives the field term (on
/// one site the bath gives the letters I, X, Y, Z the diagonal entries 0, −s, −s, −2s, so B's sit at +s, 0, 0, −s, and a
/// transverse field turns a pair with offset s/2); B's arrows, D[σ∓](I) = ±Z, give the second, A² undoing an arrow on site l only
/// through two terms of H whose product is a multiple of Z_l. An anticommuting pair, product ±iZ_l, cancels
/// against its reversed product; the commuting pairs, a string with I at l and the same string with Z at l, are
/// what Tr(Z_l·H²) collects.</para>
///
/// <para>In the integer form the witness reads: Tr((2L + c)³) = 12·4^N·(F − 2·R), with F = Σ_l (γ↓_l + γ↑_l)·
/// (h_X,l² + h_Y,l²) and R = Σ_l (γ↓_l − γ↑_l)·P_l, P_l the sum of c_T·c_T' over the strings T with I at l and
/// T' = T with Z at l, so that Tr(Z_l·H²) = 2^(N+1)·P_l.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/ThermalThirdMomentWitness.cs</c>,
/// <c>inspect --root thirdmoment</c>, the generator built in exact Gaussian integers on its own road beside the
/// registry's gate <c>simulations/thermal_third_moment.py</c>. F137 itself is not typed, so this claim has no typed
/// parent. Adjacent, not parents: <see cref="MomentTowerPumpChannelClaim"/> (the arrow term is
/// −(3/2)·4^N·slope₂, the second rung of the pump-slope law), <see cref="CorrelatedZDephasingPriceClaim"/>
/// (F166's third moment under correlated Z dephasing, the twin with Tr(H²·Z_k·Z_l) where this one has
/// Tr(Z_l·H²)), <see cref="CooledBondPowerSumsClaim"/> (p₃ = 0 at every N for its five-letter class by the phase
/// count; for the nine-letter bonds of WHICH_BONDS_PAIR_UNDER_COOLING, which that page reads off clusters fitting in
/// three sites, this law is a second route to p₃ = 0 at every N), <see cref="ThreeLadderHingeClaim"/> (under Z dephasing its 3·Tr(Q·A²) holds the cell-free
/// face as the first moment squared, the t₀·t₂ term dying with Tr(Z_l) = 0; under the bath the arrow I → Z revives
/// it, and 3·Tr(A²·B) holds the transverse field term and the rung t₂(l)), <see cref="PalindromeTwoEndCountClaim"/> (at equal rates
/// the bath is Pauli dephasing along X and Y, and a transverse field's break is F158's shortest odd word
/// Tr(H·X_l)). The field term is F87's cell-free face 6·4^N·Σ c²·γ read on the two transverse axes with the
/// offset s/2.</para></summary>
public sealed class ThermalThirdMomentClaim : Claim
{
    /// <summary>One Pauli string of H, site 0 first, with an integer coefficient.</summary>
    public readonly record struct Term(string Letters, long Coefficient)
    {
        public override string ToString() =>
            Coefficient == 1 ? Letters : Coefficient == -1 ? "−" + Letters : $"{(Coefficient < 0 ? "−" : "")}{BigInteger.Abs(Coefficient)}·{Letters}";
    }

    /// <summary>F and R of the integer form, Tr((2L + c)³) = 12·4^N·(F − 2·R): F the field part
    /// Σ_l (γ↓_l + γ↑_l)·(h_X,l² + h_Y,l²), R the arrow part Σ_l (γ↓_l − γ↑_l)·P_l. Repeated strings are summed and
    /// the identity string is dropped.</summary>
    public static (BigInteger Field, BigInteger Arrow) Parts(IReadOnlyList<Term> terms, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        if (cool is null) throw new ArgumentNullException(nameof(cool));
        if (heat is null) throw new ArgumentNullException(nameof(heat));
        int n = cool.Count;
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(cool), n, "at least one site");
        if (heat.Count != n) throw new ArgumentException($"{heat.Count} heating rates for {n} sites", nameof(heat));
        for (int l = 0; l < n; l++)
        {
            if (cool[l] < 0) throw new ArgumentOutOfRangeException(nameof(cool), "rates are nonnegative");
            if (heat[l] < 0) throw new ArgumentOutOfRangeException(nameof(heat), "rates are nonnegative");
        }
        var h = new Dictionary<string, BigInteger>();
        string identity = new('I', n);
        foreach (var t in terms)
        {
            if (t.Letters is null || t.Letters.Length != n || t.Letters.Any(ch => "IXYZ".IndexOf(ch) < 0))
                throw new ArgumentException($"'{t.Letters}' is not a Pauli string on {n} sites", nameof(terms));
            if (t.Letters == identity) continue;
            h[t.Letters] = (h.TryGetValue(t.Letters, out var o) ? o : BigInteger.Zero) + t.Coefficient;
        }
        BigInteger field = BigInteger.Zero, arrow = BigInteger.Zero;
        foreach (var (s, c) in h)
        {
            int support = s.Count(ch => ch != 'I');
            if (support == 1)
            {
                int l = s.IndexOfAny(new[] { 'X', 'Y', 'Z' });
                if (s[l] != 'Z') field += (new BigInteger(cool[l]) + heat[l]) * c * c;
            }
            for (int m = 0; m < n; m++)
            {
                if (s[m] != 'I') continue;
                var partner = s.Substring(0, m) + "Z" + s.Substring(m + 1);
                if (h.TryGetValue(partner, out var cp)) arrow += (new BigInteger(cool[m]) - heat[m]) * c * cp;
            }
        }
        return (field, arrow);
    }

    /// <summary>Tr((2L + c)³) = 8·Tr((L + c/2)³), the closed form: 12·4^N·(F − 2·R).</summary>
    public static BigInteger CubedCentredTrace(IReadOnlyList<Term> terms, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        var (f, r) = Parts(terms, cool, heat);
        return 12 * BigInteger.Pow(4, cool.Count) * (f - 2 * r);
    }

    public ThermalThirdMomentClaim()
        : base("Under cooling and heating alone, for every H a sum of Pauli strings with real coefficients, " +
               "Tr((L + c/2)³) = 12·4^(N−1)·Σ_l s_l(h_X,l² + h_Y,l²) − 3·2^(N−1)·Σ_l (γ↓_l − γ↑_l)·Tr(Z_l·H²), H taken traceless",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md (F137, the third moment, for every H, under cooling and heating)")
    {
    }

    public override string DisplayName =>
        "Thermal third moment: a field term and an arrow term, for every H";

    public override string Summary =>
        "cooling γ↓_l and heating γ↑_l on every site, s_l = (γ↓_l + γ↑_l)/2, H any real sum of Pauli strings taken traceless: " +
        "Tr((L + c/2)³) = 12·4^(N−1)·Σ_l s_l(h_X,l² + h_Y,l²) − 3·2^(N−1)·Σ_l (γ↓_l − γ↑_l)·Tr(Z_l·H²); a transverse field on a damped " +
        "site enters with the sum of its rates, the arrow through the commuting term pairs with product ±Z_l (Tier 1 derived; " +
        "live witness: inspect --root thirdmoment).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the derivation",
                summary: "L + c/2 = A + B, A = −i[H, ·], B the centred dissipator; Tr A³ = 0, Tr B³ = 0 and Tr(A·B²) = 0 (B moves letters only from I to Z, and A could undo that only through a Z-string, which commutes there), " +
                         "so the moment is 3·Tr(A²·B); B's diagonal (centred +s, 0, 0, −s on I, X, Y, Z) gives the field term, a transverse field " +
                         "turning a pair with offset s/2, and B's arrows D[σ∓](I) = ±Z give Tr(Z_l·H²), the commuting term pairs whose product is ±Z_l");
            yield return new InspectableNode("where the two terms already live",
                summary: "the field term is F87's cell-free face 6·4^N·Σ c²·γ on the two transverse axes (ThreeLadderHingeClaim holds that face " +
                         "under Z dephasing, where the t₀·t₂ term dies with Tr(Z_l) = 0 and the bath's arrow revives it); the arrow term is −(3/2)·4^N·slope₂, the " +
                         "second rung of the pump-slope law (MomentTowerPumpChannelClaim, PROOF_MOMENT_TOWER_PUMP_CHANNEL §1); F166's third moment " +
                         "under correlated Z dephasing is the twin, Tr(H²·Z_k·Z_l) in the place of Tr(Z_l·H²)");
            yield return new InspectableNode("the instances",
                summary: "Heisenberg N = 3 under cooling 1: a uniform X field 1 gives 288 (a formal setting at zero splitting), a uniform Z field 1 " +
                         "gives −768 as written (the lab frame of the field-free chain, which pairs in F137's frame); every bond of the nine " +
                         "letters of WHICH_BONDS_PAIR_UNDER_COOLING has a vanishing third moment at every N, no bond term being a field and no two " +
                         "having a product ±Z_l");
            yield return new InspectableNode("equal rates",
                summary: "γ·(D[σ⁻] + D[σ⁺]) = (γ/2)·(D[X] + D[Y]): at equal rates on every site the bath is Pauli dephasing along X and Y, the arrow " +
                         "term is gone, and a transverse field on a damped site breaks the pairing, F158's shortest odd word Tr(H·X_l) (PROOF_PALINDROME_TWO_END_COUNT §(f5))");
            yield return new InspectableNode("what this does not say",
                summary: "a vanishing third moment decides nothing; where the field term and the arrow term cancel, the fifth sum can still break " +
                         "(ConditionedFieldPairingClaim names the controls). The law is about cooling and heating alone: no dephasing, no other jump");
        }
    }
}
