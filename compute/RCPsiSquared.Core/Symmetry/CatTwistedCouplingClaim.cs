using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F168 under a twisted coupling (Tier 1 derived, <c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>,
/// section "Under a twisted coupling"): a hopping with a twist δ_e on each bond e = (a, b), oriented
/// from a to b, J_e cos δ_e·(X_aX_b + Y_aY_b) + J_e sin δ_e·(X_aY_b − Y_aX_b), ZZ terms beside it, a
/// field h_l·(cos φ_l·X + sin φ_l·Y) at an angle of its own on every bit (h_l ≠ 0) and Z watching at
/// rates γ_l. Turning every bit by its angle leaves the residual twist β_e = δ_e − (φ_b − φ_a) on each
/// bond, the twisted current j_e = h_a·h_b·J_e·sin β_e, and the cat's phase ϑ = θ − Σφ_l.
///
/// <list type="bullet">
///   <item>Where every β_e is a multiple of π, F168 holds with Nφ replaced by Σφ_l; on a ring that
///   needs the flux Σδ_e to be a multiple of π.</item>
///   <item>Otherwise <c>m_7 = 96·Σ_l γ_l·(div j)_l²</c>, (div j)_l the current leaving bit l less the
///   current arriving: proven from N = 3 on where sin ϑ = 0 and at every phase from N = 6 on; at
///   N = 4 and 5 at other phases it rests on the gates and on the live witness.</item>
/// </list>
///
/// <para>The proof's four steps: the first move turns the twist's commutator into its
/// anticommutator, so only a twist acting once on the ket and once on the bra counts; each side needs a
/// field step before and after it and one rate falls between (order 7); the watching charges only where
/// ket and bra disagree, and with the loop sums zero only the bits where both hold their excitation
/// survive (the 16); the orderings weigh 12 + 12 − 9 − 9 = 6.</para>
///
/// <para>Live witness: <c>inspect --root cattwist</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/CatTwistedCouplingWitness.cs</c>), which reads the
/// moments word by word in the original frame (twists and own angles at once), a different road from
/// the page's ℚ(i) matrix gates in <c>simulations/the_wave_between.py</c> section 11.</para>
///
/// <para>Typed parent: <see cref="CatTurnedFieldSelfEncounterClaim"/> (F168, of which this is the
/// widening: where the twist is followed it is F168 in the turned frame).</para></summary>
public sealed class CatTwistedCouplingClaim : Claim
{
    public CatTurnedFieldSelfEncounterClaim CatTurnedField { get; }

    /// <summary>The weight of the order-7 law: 6 orderings times 16 from the amplitudes.</summary>
    public const int OrderSevenWeight = 96;

    /// <summary>sin β_e for the residual twist β_e = δ_e − (φ_b − φ_a), from rational points of the
    /// unit circle (cos, sin) of δ_e, φ_a and φ_b, exact.</summary>
    public static BigRational SinResidualTwist((BigRational C, BigRational S) delta,
                                               (BigRational C, BigRational S) phiA,
                                               (BigRational C, BigRational S) phiB)
    {
        foreach (var (c, s) in new[] { delta, phiA, phiB })
            if (c * c + s * s != BigRational.One)
                throw new ArgumentException("every angle must be a rational point of the unit circle");
        // φ_b − φ_a, then δ − (φ_b − φ_a)
        var dc = phiB.C * phiA.C + phiB.S * phiA.S;
        var ds = phiB.S * phiA.C - phiB.C * phiA.S;
        return delta.S * dc - delta.C * ds;
    }

    /// <summary>(div j)_l: the current leaving bit l less the current arriving, for bonds (A, B, j)
    /// oriented from A to B.</summary>
    public static BigRational[] Divergence(int n, IEnumerable<(int A, int B, BigRational Current)> bonds)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1");
        var div = new BigRational[n];
        for (int l = 0; l < n; l++) div[l] = BigRational.Zero;
        foreach (var (a, b, j) in bonds)
        {
            div[a] += j;
            div[b] -= j;
        }
        return div;
    }

    /// <summary>m_7 = 96·Σ_l γ_l·(div j)_l², exact.</summary>
    public static BigRational OrderSevenMoment(IReadOnlyList<BigRational> rates, IReadOnlyList<BigRational> divergence)
    {
        if (rates.Count != divergence.Count) throw new ArgumentException("one rate per bit");
        BigRational sum = BigRational.Zero;
        for (int l = 0; l < rates.Count; l++) sum += rates[l] * divergence[l] * divergence[l];
        return new BigRational(OrderSevenWeight) * sum;
    }

    public CatTwistedCouplingClaim(CatTurnedFieldSelfEncounterClaim catTurnedField)
        : base(
            "The cat under a twisted coupling: where the field angles follow the twist F168 holds with Nφ → Σφ_l " +
            "(on a ring only at flux 0 or π); otherwise m_7 = 96·Σ_l γ_l·(div j)_l², the look seeing where the " +
            "twisted current j_e = h_a·h_b·J_e·sin β_e piles up at a watched bit",
            Tier.Tier1Derived,
            "experiments/THE_CAT_UNDER_A_TURNED_FIELD.md + docs/ANALYTICAL_FORMULAS.md (F168, the index entry)")
    {
        CatTurnedField = catTurnedField ?? throw new ArgumentNullException(nameof(catTurnedField));
    }

    public override string DisplayName =>
        "F168 under a twisted coupling: m_7 = 96·Σγ_l(div j)_l²";

    public override string Summary =>
        "A hopping with a twist (DM part) and a field angle of its own on every bit: with the residual twist " +
        "β_e = δ_e − (φ_b − φ_a) and the current j_e = h_a h_b J_e sin β_e, F168 holds with Σφ_l where every β_e is a " +
        "multiple of π (a ring needs flux 0 or π), and otherwise m_7 = 96·Σγ_l(div j)_l², proven where sin(Σφ_l − θ) = 0 " +
        "from N = 3 on and at every phase from N = 6 on, gated at N = 4, 5 off it; a ring's source-free circulation " +
        "escapes order 7 (Tier 1 derived; live witness: inspect --root cattwist).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("both sides",
                summary: "the first move R keeps the even steps up to the adjoint and turns the twist's commutator Λ − Ω into its anticommutator Λ + Ω (X^⊗N·T·X^⊗N = −T); at ϑ = 0 or π, averaging the two readings of m_k cancels every term with an even number of bra steps, and parity kills the odd ones: only a twist once on the ket and once on the bra counts");
            yield return new InspectableNode("the order",
                summary: "the twist moves nothing on |0…0⟩ or |1…1⟩, so each side needs a field step before it and after it, and without a rate the look is even in t: 2 twists + 4 field steps + 1 rate, order 7 at every N ≥ 3");
            yield return new InspectableNode("where the rate falls",
                summary: "the rate charges σ − 2γ_u − 2γ_v + 4γ_u·[u = v] (the cube's height, weighted); the ket's loop sums Σκ = ⟨f, T₁f⟩ = 0, so σ and the single charges drop and only the bits where ket and bra hold their excitation together survive: 4Σγκ² = −16Σγ(div j)²");
            yield return new InspectableNode("the count",
                summary: "with the rate after both twist steps 12 orderings, before both 12, mixed 9 each with the sign flipped: 12 + 12 − 9 − 9 = 6, and with step 1's minus m_7 = 96·Σγ(div j)²");
            yield return new InspectableNode("the torque",
                summary: "with an untwisted hopping and the fields at their own angles, (div j)_l = Σ_m J_lm h_l h_m sin(φ_m − φ_l) up to one sign: the look asks whether every watched bit's field lies along the field its neighbours' fields add up to");
        }
    }
}
