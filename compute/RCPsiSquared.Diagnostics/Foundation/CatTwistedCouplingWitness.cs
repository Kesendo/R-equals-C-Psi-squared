using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F168 under a twisted coupling (<c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>,
/// section "Under a twisted coupling"): the cat under Z dephasing, a field at an angle of its own on
/// every bit and a hopping with a twist on every bond, its centred odd moments read exactly word by word
/// (<see cref="CatPauliMoments"/>).
///
/// <para>The page's gates work in the frame where every field lies along X, or with one of the two
/// (twist, own angles) at a time, as density matrices over ℚ(i). This witness stays in the original
/// frame with both at once, as Pauli words over the rationals, and computes the residual twist
/// β_e = δ_e − (φ_b − φ_a) only for the claim's closed form. Bit 1 is unwatched (γ_1 = 0), so a source
/// there is invisible to the law, and the strengths carry both signs.</para>
///
/// <para>Children: the law (m_1, m_3, m_5 = 0 and m_7 = 96·Σγ_l(div j)_l² at the phase θ = Σφ_l, nonzero);
/// the other phase (the same m_7 at θ = Σφ_l + arccos 3/5 from N = 4 on, and at N = 3 a different
/// number, where F168's own term and a term linear in the twist join); following the twist (the angles
/// turned along every bond, the ring closed at flux π: every odd moment through 2N + 3 vanishes at
/// θ = Σφ_l, and at θ = 0 the first is F168's with Nφ replaced by Σφ_l); the control (the same weights
/// on each bond's own current squared differ); a constant current (fields along X, j_e the same on every
/// bond: on the chain its sources sit at the two ends, m_7 = 96·(γ_0 + γ_{N−1})·j², on the ring it has
/// none, m_7 = 0 and m_9 ≠ 0).</para>
///
/// <para>Guard: N ∈ {3, 4, 5}. Typed home <see cref="CatTwistedCouplingClaim"/>.</para></summary>
public sealed class CatTwistedCouplingWitness : IInspectable
{
    public const int MinN = 3;
    public const int MaxN = 5;

    private static readonly BigRational[] Strengths = { R(3, 2), R(-1), R(4, 5), R(-5, 3), R(2, 3) };
    private static readonly BigRational[] Rates = { R(2, 5), R(0), R(1, 3), R(3, 4), R(1, 6) };
    private static readonly BigRational[] Couplings = { R(1), R(-3, 2), R(2, 3), R(5, 4), R(1, 2) };
    private static readonly BigRational[] Anisotropies = { R(-1, 3), R(2), R(1, 2), R(0), R(-3) };
    private static readonly (BigRational C, BigRational S)[] Angles =
        { (R(5, 13), R(-12, 13)), (R(3, 5), R(4, 5)), (R(-8, 17), R(15, 17)), (R(24, 25), R(7, 25)), (R(-4, 5), R(-3, 5)) };
    private static readonly (BigRational C, BigRational S)[] Twists =
        { (R(20, 29), R(21, 29)), (R(-3, 5), R(4, 5)), (R(12, 13), R(-5, 13)), (R(7, 25), R(-24, 25)), (R(15, 17), R(8, 17)) };
    private static readonly (BigRational C, BigRational S) OffPhase = (R(3, 5), R(4, 5));
    private static readonly BigRational ConstantCurrent = R(2, 5);

    public int N { get; }
    public bool Ring { get; }
    private int BondCount => Ring ? N : N - 1;
    private (int A, int B) Bond(int e) => (e, (e + 1) % N);

    /// <summary>m_1 … m_7 at θ = Σφ_l, and the claim's 96·Σγ(div j)².</summary>
    public IReadOnlyDictionary<int, BigRational> MatchedMoments { get; }
    public BigRational LawValue { get; }
    /// <summary>m_1 … m_7 at θ = Σφ_l + arccos 3/5.</summary>
    public IReadOnlyDictionary<int, BigRational> OffPhaseMoments { get; }
    /// <summary>The same weights on each bond's own current squared, Σ_e 96·(γ_a + γ_b)·j_e², the control.</summary>
    public BigRational PerBondValue { get; }
    /// <summary>The angles following the twist (ring at flux π): odd moments through 2N + 3 at θ = Σφ_l.</summary>
    public IReadOnlyDictionary<int, BigRational> FollowingMoments { get; }
    /// <summary>The same following setting at θ = 0: odd moments through 2N + 1, and F168's law with Σφ_l.</summary>
    public IReadOnlyDictionary<int, BigRational> FollowingOffMoments { get; }
    public BigRational FollowingLawValue { get; }
    /// <summary>Fields along X, the same current on every bond: m_1 … m_9, and the end-source value on the chain.</summary>
    public IReadOnlyDictionary<int, BigRational> ConstantCurrentMoments { get; }
    public BigRational EndSourceValue { get; }

    public CatTwistedCouplingWitness(int n = 4, bool ring = false)
    {
        if (n < MinN || n > MaxN) throw new ArgumentOutOfRangeException(nameof(n), n, $"N must lie in {MinN}..{MaxN}");
        N = n;
        Ring = ring;
        var angles = Angles.Take(n).ToArray();
        var twists = Twists.Take(BondCount).ToArray();
        var sumPhi = Sum(angles);

        // the law, in the original frame
        var terms = Terms(angles, twists);
        var a = CatPauliMoments.CatAmounts(n, sumPhi.C, sumPhi.S);
        MatchedMoments = CatPauliMoments.OddMoments(n, a, terms, Rates, 7);
        var currents = Enumerable.Range(0, BondCount).Select(e =>
        {
            var (p, q) = Bond(e);
            var sinBeta = CatTwistedCouplingClaim.SinResidualTwist(twists[e], angles[p], angles[q]);
            return (p, q, Strengths[p] * Strengths[q] * Couplings[e] * sinBeta);
        }).ToArray();
        var div = CatTwistedCouplingClaim.Divergence(n, currents);
        LawValue = CatTwistedCouplingClaim.OrderSevenMoment(Rates.Take(n).ToArray(), div);
        PerBondValue = currents.Aggregate(BigRational.Zero, (acc, c) =>
            acc + new BigRational(CatTwistedCouplingClaim.OrderSevenWeight) * (Rates[c.p] + Rates[c.q]) * c.Item3 * c.Item3);
        var off = Rotate(sumPhi, OffPhase);
        OffPhaseMoments = CatPauliMoments.OddMoments(n, CatPauliMoments.CatAmounts(n, off.C, off.S), terms, Rates, 7);

        // following the twist: φ_{l+1} = φ_l + δ_l; on the ring the last twist closes the flux at π
        var followTwists = twists.ToArray();
        if (ring)
        {
            var partial = Sum(followTwists.Take(BondCount - 1));
            followTwists[BondCount - 1] = Rotate((R(-1), R(0)), (partial.C, -partial.S));
        }
        var followAngles = new (BigRational C, BigRational S)[n];
        followAngles[0] = angles[0];
        for (int l = 1; l < n; l++) followAngles[l] = Rotate(followAngles[l - 1], followTwists[l - 1]);
        var followTerms = Terms(followAngles, followTwists);
        var followSum = Sum(followAngles);
        FollowingMoments = CatPauliMoments.OddMoments(n, CatPauliMoments.CatAmounts(n, followSum.C, followSum.S),
                                                      followTerms, Rates, 2 * n + 3);
        FollowingOffMoments = CatPauliMoments.OddMoments(n, CatPauliMoments.CatAmounts(n, BigRational.One, BigRational.Zero),
                                                         followTerms, Rates, 2 * n + 1);
        BigRational sigma = BigRational.Zero, prod = BigRational.One;
        for (int l = 0; l < n; l++) { sigma += Rates[l]; prod *= Strengths[l] * Strengths[l]; }
        FollowingLawValue = CatTurnedFieldSelfEncounterClaim.LeadingOddMoment(n, sigma, prod, followSum.S * followSum.S);

        // a constant current: fields along X, plain hopping, DM coefficient j / (h_a h_b) on every bond
        var cTerms = new List<(PauliMask, BigRational)>();
        for (int l = 0; l < n; l++) cTerms.Add((new PauliMask(1UL << l, 0), Strengths[l]));
        for (int e = 0; e < BondCount; e++)
        {
            var (p, q) = Bond(e);
            AddBond(cTerms, p, q, Couplings[e], ConstantCurrent / (Strengths[p] * Strengths[q]), Anisotropies[e]);
        }
        ConstantCurrentMoments = CatPauliMoments.OddMoments(n, CatPauliMoments.CatAmounts(n, BigRational.One, BigRational.Zero),
                                                            cTerms, Rates, 9);
        EndSourceValue = ring ? BigRational.Zero
            : new BigRational(CatTwistedCouplingClaim.OrderSevenWeight) * (Rates[0] + Rates[n - 1]) * ConstantCurrent * ConstantCurrent;
    }

    private static BigRational R(long p, long q = 1) => new(new BigInteger(p), new BigInteger(q));

    private static (BigRational C, BigRational S) Rotate((BigRational C, BigRational S) a, (BigRational C, BigRational S) b) =>
        (a.C * b.C - a.S * b.S, a.S * b.C + a.C * b.S);

    private static (BigRational C, BigRational S) Sum(IEnumerable<(BigRational C, BigRational S)> angles) =>
        angles.Aggregate((BigRational.One, BigRational.Zero), (acc, x) => Rotate(acc, x));

    private static void AddBond(List<(PauliMask, BigRational)> terms, int p, int q,
                                BigRational hop, BigRational dm, BigRational zz)
    {
        ulong a = 1UL << p, b = 1UL << q;
        terms.Add((new PauliMask(a | b, 0), hop));              // X_p X_q
        terms.Add((new PauliMask(a | b, a | b), hop));          // Y_p Y_q
        terms.Add((new PauliMask(a | b, b), dm));               // X_p Y_q
        terms.Add((new PauliMask(a | b, a), -dm));              // −Y_p X_q
        terms.Add((new PauliMask(0, a | b), zz));               // Z_p Z_q
    }

    /// <summary>H = Σ_l h_l(cos φ_l X + sin φ_l Y) + Σ_e J_e[cos δ_e (XX + YY) + sin δ_e (XY − YX)] + Δ_e ZZ.</summary>
    private List<(PauliMask T, BigRational Coefficient)> Terms((BigRational C, BigRational S)[] angles,
                                                               (BigRational C, BigRational S)[] twists)
    {
        var terms = new List<(PauliMask, BigRational)>();
        for (int l = 0; l < N; l++)
        {
            ulong bit = 1UL << l;
            terms.Add((new PauliMask(bit, 0), Strengths[l] * angles[l].C));
            terms.Add((new PauliMask(bit, bit), Strengths[l] * angles[l].S));
        }
        for (int e = 0; e < BondCount; e++)
        {
            var (p, q) = Bond(e);
            AddBond(terms, p, q, Couplings[e] * twists[e].C, Couplings[e] * twists[e].S, Anisotropies[e]);
        }
        return terms;
    }

    private static bool BelowSevenVanish(IReadOnlyDictionary<int, BigRational> m) => m[1].IsZero && m[3].IsZero && m[5].IsZero;

    public bool LawHolds => BelowSevenVanish(MatchedMoments) && MatchedMoments[7] == LawValue && !LawValue.IsZero;
    /// <summary>From N = 4 on the same m_7 at the other phase; at N = 3 a different one.</summary>
    public bool OffPhaseReadingHolds => N >= 4
        ? BelowSevenVanish(OffPhaseMoments) && OffPhaseMoments[7] == LawValue
        : OffPhaseMoments[7] != LawValue;
    public bool ControlDiffers => PerBondValue != LawValue;
    public bool FollowingVanishes => FollowingMoments.Values.All(m => m.IsZero);
    public bool FollowingLawHolds =>
        FollowingOffMoments.Where(kv => kv.Key < 2 * N + 1).All(kv => kv.Value.IsZero)
        && FollowingOffMoments[2 * N + 1] == FollowingLawValue && !FollowingLawValue.IsZero;
    public bool ConstantCurrentReadingHolds => BelowSevenVanish(ConstantCurrentMoments)
        && ConstantCurrentMoments[7] == EndSourceValue
        && (Ring ? !ConstantCurrentMoments[9].IsZero : !EndSourceValue.IsZero);

    public string DisplayName => $"CatTwistedCouplingWitness (N={N}, {(Ring ? "ring" : "chain")})";

    public string Summary =>
        $"N={N} {(Ring ? "ring" : "chain")}, every bond twisted, every field at its own angle, bit 1 unwatched: " +
        $"m_7 = {MatchedMoments[7]} {(LawHolds ? "is" : "is NOT")} 96·Σγ(div j)² at θ = Σφ; " +
        $"at another phase {(OffPhaseReadingHolds ? (N >= 4 ? "the same" : "a different number, as at N = 3 it must be") : "NOT as expected")}; " +
        $"the angles following the twist {(FollowingVanishes ? "leave every odd moment zero" : "do NOT leave every odd moment zero")} " +
        $"and at θ = 0 give F168's law with Σφ {(FollowingLawHolds ? "exactly" : "NOT")}; a constant current " +
        $"{(ConstantCurrentReadingHolds ? (Ring ? "is invisible at order 7 and breaks at m_9" : "shows its two end sources") : "does NOT read as predicted")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the law",
                summary: $"θ = Σφ_l: m_1, m_3, m_5 = {MatchedMoments[1]}, {MatchedMoments[3]}, {MatchedMoments[5]}; " +
                         $"m_7 = {MatchedMoments[7]} against 96·Σγ_l(div j)_l² = {LawValue}: {(LawHolds ? "equal, and not zero" : "DIFFERENT or zero")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the other phase",
                summary: $"θ = Σφ_l + arccos 3/5: m_7 = {OffPhaseMoments[7]}" +
                         (N >= 4
                             ? $" {(OffPhaseReadingHolds ? "equals the law, free of the phase" : "does NOT equal the law")}"
                             : $" {(OffPhaseReadingHolds ? "differs from the law: at N = 3 order 7 is also F168's 2N + 1, and a term linear in the twist joins" : "equals the law, unexpectedly")}") + ".",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "following the twist",
                summary: $"φ_{{l+1}} = φ_l + δ_l{(Ring ? ", the ring closed at flux π" : "")}: at θ = Σφ_l every odd moment through " +
                         $"m_{2 * N + 3} is {(FollowingVanishes ? "exactly 0" : "NOT 0")}; at θ = 0, m_{2 * N + 1} = {FollowingOffMoments[2 * N + 1]} " +
                         $"{(FollowingLawHolds ? "is" : "is NOT")} F168's (σ/3)(2N+1)!·∏h_l²·sin²(Σφ_l) = {FollowingLawValue}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the control",
                summary: $"the same weights on each bond's own current squared would give {PerBondValue}: " +
                         $"{(ControlDiffers ? "a different number, so the gate sees the divergence" : "the same number, the control is void")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "a constant current",
                summary: $"fields along X, h_a h_b·ξ_e = {ConstantCurrent} on every bond: m_7 = {ConstantCurrentMoments[7]}" +
                         (Ring
                             ? $" (no source on the ring) and m_9 = {ConstantCurrentMoments[9]}: {(ConstantCurrentReadingHolds ? "invisible at order 7, broken at order 9" : "NOT as predicted")}."
                             : $" against 96·(γ_0 + γ_{N - 1})·j² = {EndSourceValue}: {(ConstantCurrentReadingHolds ? "the chain's current has its sources at the two ends" : "NOT as predicted")}."),
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
