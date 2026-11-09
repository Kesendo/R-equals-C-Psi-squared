using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The palindrome with the jump on one seat alone
/// (<c>docs/proofs/PROOF_ONE_DEPHASED_SEAT.md</c>): the Heisenberg triangle at N = 3, with letter fields and with
/// planar fields, and at N = 4 and N = 7 palindromes with neither a colouring nor a site symmetry fixing the dephased
/// seat. Theorem numbers
/// are those of the complement connection (<see cref="PalindromeComplementConnectionClaim"/>), which reads every
/// site dephased and then undephased sites under its own hypotheses.
///
/// <para><b>Common triangle carrier.</b> With jump Z on site 0, Heisenberg couplings s=J01, t=J02, r=J12,
/// all real and nonzero, and field h(X1+Y2), the common Z0-odd commutant of the bond Hamiltonian and X1+Y2
/// is a line iff (s-t)(st+sr+tr)=0, zero otherwise. Its generator squares to (s²+t²)³ I. This classifies
/// one invertible carrier independent of h, not every fixed-h palindrome. Exact symbolic gate:
/// <c>simulations/triangle_common_carrier_gate.py</c>.</para>
///
/// <para><b>The triangle's palindromic locus at N = 3.</b> Heisenberg triangle, jump Z on site 0 only, a letter field
/// of free magnitude on every site, at most one coupling zero: a row pairs exactly when it is coloured, carries the
/// 1↔2 swap with a π rotation about an axis perpendicular to z (J01 = J02), or lies on e₂ = J01J02 + J01J12 + J02J12 = 0
/// with no field on 0 and fields of equal magnitude perpendicular to z on 1 and 2.
/// Necessity by an exact computer-assisted elimination: <c>simulations/n3_palindrome_covering.py</c> (certificate) and
/// <c>simulations/n3_palindrome_covering_gate.py</c> (independent gate); live at <c>inspect --root trianglelocus</c>
/// (<c>TrianglePalindromeLocusWitness</c>), which recomputes both ends exactly over Q at points.</para>
///
/// <para><b>The triangle with planar fields.</b> Same triangle and jump, no field on 0, fields f₁, f₂ (the section's f_u,
/// f_v) of any direction perpendicular to z, at most one coupling zero: a row pairs exactly when f₁ ∥ f₂, when J01 = J02 and |f₁| = |f₂|,
/// when e₂ = 0 and |f₁| = |f₂| at any angle, or on R1 = J12(J01 − J02)² − J01J02(J01 + J02) = 0, R2 = f₁·f₂ + 2e₂ = 0,
/// R3 = (J01 + J12)|f₁|² + (J02 + J12)|f₂|² + 2(J01 + J02 + J12)(f₁·f₂) = 0; with site 0 at an end only the first.
/// Both directions exact and computer-assisted: <c>simulations/n3_planar_field_locus_gate.py</c>; live at
/// <c>inspect --root trianglelocus</c> with <c>--fu</c>, <c>--fv</c>, which recomputes both ends exactly over Q at
/// points and on a planar grid.</para>
///
/// <para><b>Blind dephased seats, and palindromes beyond colourings and site symmetries.</b> On the triangle the Krylov determinant of the
/// dephased seat under the signed one-magnon Laplacian is −(s − t)·e₂: s = t a site symmetry fixing it, e₂ = 0 a rank
/// drop that blinds every seat (F157, <see cref="SeatCutBlindnessClaim"/>). At N = 4, H = e·σ1·σ2 + a·(σ0 − σ3)·(σ1 − σ2) with fields 0, −h·X1, h·X2, g·Y3
/// pairs for every real a, g ≠ 0 and every real e, h, for h ≠ 0 with neither a colouring nor a site symmetry fixing site 0, its carrier
/// in closed form and depending on h (exact, computer-assisted: <c>simulations/n4_pair_channel_family.py</c>, elimination and dense gate);
/// at N = 7 a construction, a dephased centre joined to two scalene weighted triangles with different weight sets, equivalent as spin representations, pairs with no rank drop and no site symmetry
/// (<c>simulations/blind_seat_carriers_gate.py</c>). The sheet and blindness: on the triangle with Heisenberg bonds of
/// which at most one is zero, the jump Z on site 0 alone and nonzero planar fields f₁, f₂ on sites 1 and 2 neither
/// parallel nor antiparallel, a far element of odd height (every string with an odd number of letters X or Y) exists
/// with no field on site 0 exactly when site 0 is blind and |f₁| = |f₂|, and with a nonzero planar field f₀ on site 0 exactly
/// when J01 = J02, |f₁| = |f₂| and f₀ lies on the line of f₁ + f₂ (<c>simulations/blind_seat_sheet_gate.py</c>; live in
/// <c>TrianglePalindromeLocusWitness</c>, node "the sheet and blindness").</para></summary>
public sealed class OneDephasedSeatClaim : Claim
{
    /// <summary>F158, whose Lemma 1's second corollary puts both ends on spans of Pauli strings and whose Lemma 3 with its
    /// theorem lets an invertible far element carry the palindrome: every carrier here is such an element.</summary>
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>The complement connection, whose far equations, written before Theorem 6's spanning-tree hypothesis, the
    /// section on planar fields solves on a component without an invertible spanning tree, outside that hypothesis; its
    /// Theorems 2, 4, 5 and 7 are read beside the sections, adjacent rather than built on.</summary>
    public PalindromeComplementConnectionClaim ComplementConnection { get; }

    /// <summary>F157, the blind seat: on the triangle the dephased seat's Krylov determinant −(s − t)·e₂ vanishes
    /// exactly where the seat is blind, and the conjecture asks, under Heisenberg bonds with letter fields, that every
    /// palindrome without a colouring at one dephased seat sit on a blind one.</summary>
    public SeatCutBlindnessClaim BlindSeat { get; }

    public OneDephasedSeatClaim(
        PalindromeTwoEndCountClaim twoEndCount,
        PalindromeComplementConnectionClaim complementConnection,
        SeatCutBlindnessClaim blindSeat)
        : base("With the jump Z on one dephased seat alone, under Heisenberg bonds. On the triangle on (d, u, v) with " +
               "couplings s = J_du, t = J_dv, r = J_uv and e2 = st + sr + tr: at real nonzero couplings the common space " +
               "of carriers that serve the whole field line h(X_u + Y_v) is a line when (s - t)e2 = 0 and zero otherwise " +
               "(by hand, with an exact gate). With letter fields and at most one coupling zero, a row pairs exactly when it " +
               "is coloured, carries the u<->v swap composed with a pi rotation about an axis perpendicular to z (s = t), or " +
               "lies on e2 = 0 with no field on d and fields of equal magnitude perpendicular to z on u and v (an exact, " +
               "computer-assisted elimination with a stored certificate and an independent gate). With no field on d and " +
               "fields of any direction perpendicular to z, at most one coupling zero, it pairs exactly with parallel or " +
               "antiparallel fields (one of them zero included), at s = t with fields of equal length, at e2 = 0 with fields of equal length at any angle, and on " +
               "one further locus of three polynomial relations (both directions exact and computer-assisted). The dephased " +
               "seat's Krylov determinant under the signed one-magnon Laplacian is -(s - t)e2, so the swap and reciprocal " +
               "classes sit on a blind seat (F157). Beyond the triangle, with one dephased seat, Heisenberg bonds and " +
               "letter fields: a family at N = 4 pairs for every real a, g nonzero and, at a nonzero field h, with neither a " +
               "colouring nor a site symmetry fixing d, through a closed-form carrier (exact, computer-assisted); at N = 7 a " +
               "construction, a dephased centre joined to two scalene weighted triangles with different weight sets, " +
               "equivalent as spin representations, " +
               "pairs with neither a colouring nor any site symmetry and with no rank drop of the Laplacian; at N = 5, for a, h nonzero, the N = 4 carrier survives no attachment of a fifth spin by " +
               "isotropic bonds as a product W4 x Q. That every palindrome there without a colouring needs a blind " +
               "dephased seat is a conjecture. In the widened N = 4 family (the pair seen from 0 and 3 through its " +
               "difference only, the family's fields, J03 = 0, h, g and the two difference couplings nonzero) a far element even under the pair's swap with " +
               "a nonzero singlet block exists exactly where the pair carries a mode that vanishes at 0 and at 3. " +
               "Where a row pairs through one far element, for Pauli jumps, a mirror that " +
               "factorises over the sites, alone or followed by a site permutation, exists only where that element is a " +
               "product times a site permutation fixing the watched seat; a mirror conditioned on the bits of an " +
               "undephased site needs a letter there that H conserves, a second element of the near end, so none " +
               "arises. On the triangle, with at most one coupling zero, no field on d and nonzero planar fields on u " +
               "and v neither parallel nor antiparallel, a far element of odd height exists exactly when d is blind and " +
               "the two fields have one length; with a nonzero planar field on d added, exactly when s = t, the two fields have " +
               "one length and the field on d lies on the line of their sum.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ONE_DEPHASED_SEAT.md + " +
               "docs/ANALYTICAL_FORMULAS.md (F165, the index entry; F157, the blind seat) + " +
               "docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md (Theorem 6's far equations) + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (F158, Lemma 1's second corollary and Lemma 3) + " +
               "experiments/THE_PALINDROME_AS_A_COLOURING.md (the colouring rule)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
        ComplementConnection = complementConnection ?? throw new ArgumentNullException(nameof(complementConnection));
        BlindSeat = blindSeat ?? throw new ArgumentNullException(nameof(blindSeat));
    }

    public override string DisplayName =>
        "F165: with one dephased seat, palindromes beyond colourings and site symmetries (the one dephased seat)";

    public override string Summary =>
        "the triangle's loci exact at N = 3 (letter fields: colourings, the swap class, the reciprocal class on e2 = 0; " +
        "planar fields with none on d: one locus more), the dephased seat blind on (s - t)e2 = 0; at N = 4 and N = 7 " +
        $"palindromes with neither a colouring nor a site symmetry fixing the dephased seat ({Tier.Label()}); that under " +
        "Heisenberg bonds with letter fields a palindrome without a colouring needs its dephased seat blind is a " +
        "conjecture";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("a common carrier for the triangle's whole field line",
                summary: "Heisenberg triangle, jump Z0 only, field h(X1+Y2), real nonzero s=J01, t=J02, r=J12: " +
                         "the common odd commutant of the bond Hamiltonian and X1+Y2 is a line exactly when " +
                         "(s-t)(st+sr+tr)=0, including the branch intersection, and zero otherwise. Its generator " +
                         "squares to (s^2+t^2)^3 I. This classifies one invertible carrier independent of h, not " +
                         "every fixed-h palindrome. Exact symbolic gate: simulations/triangle_common_carrier_gate.py.");

            yield return new InspectableNode("the triangle's palindromic locus at N = 3",
                summary: "Same triangle and jump, a letter field of free real magnitude on every site, at most one " +
                         "coupling zero: a row pairs exactly when it is coloured (A), carries the 1<->2 swap with a pi " +
                         "rotation about an axis perpendicular to z (B, J01 = J02), or lies on e2 = J01J02 + J01J12 + J02J12 = 0 " +
                         "with no field on 0 and fields of equal magnitude perpendicular to z on 1 and 2 (C). Sufficiency from F138, the symmetry element " +
                         "and the common carrier; necessity by exact elimination over the ten word orbits in two charts, " +
                         "certificate and independent gate: simulations/n3_palindrome_covering.py and " +
                         "simulations/n3_palindrome_covering_gate.py; live: inspect --root trianglelocus.");

            yield return new InspectableNode("the triangle with planar fields",
                summary: "Same triangle and jump, no field on 0, fields f1, f2 of any direction perpendicular to z, at " +
                         "most one coupling zero: a row pairs exactly when f1 || f2 (A), J01 = J02 and |f1| = |f2| (B), " +
                         "e2 = 0 and |f1| = |f2| at any angle (C), or on R1 = J12(J01-J02)^2 - J01J02(J01+J02) = 0, " +
                         "R2 = f1.f2 + 2e2 = 0, R3 = (J01+J12)|f1|^2 + (J02+J12)|f2|^2 + 2(J01+J02+J12)(f1.f2) = 0 (E); " +
                         "with site 0 at an end only A. Both directions exact and computer-assisted: " +
                         "simulations/n3_planar_field_locus_gate.py; live: inspect --root trianglelocus --J 1,2,6 " +
                         "--fu 8,0 --fv -5,3.");

            yield return new InspectableNode("blind dephased seats, and palindromes beyond colourings and site symmetries",
                summary: "On the triangle the Krylov determinant of the dephased seat under the signed one-magnon " +
                         "Laplacian is -(s-t)e2, the common-carrier factor: s = t is a site symmetry fixing d (B), e2 = 0 a " +
                         "rank drop of the Laplacian that blinds every seat (C) (F157). At N = 4, H = e s1.s2 + " +
                         "a (s0 - s3).(s1 - s2) with fields 0, -hX1, hX2, gY3 pairs for every real a, g != 0 and every real e, h, " +
                         "near = far = 1, for h != 0 with neither a colouring nor a site symmetry fixing site 0 (exact, computer-assisted); its " +
                         "carrier, in closed form, depends on h. At N = 7 a construction (a dephased centre joined by one coupling c > 0 to all six " +
                         "sites of two scalene weighted triangles with different weight sets, equivalent as spin representations, " +
                         "no bonds between them, fields " +
                         "h X on one and h Y on the other, h != 0) pairs with no rank drop and no site symmetry, the centre " +
                         "still blind. Gates: " +
                         "simulations/n4_pair_channel_family.py, simulations/blind_seat_carriers_gate.py. That under Heisenberg " +
                         "bonds with letter fields a palindrome without a colouring at one dephased seat needs that seat " +
                         "blind is a conjecture; " +
                         "off the letters at N = 3 blindness is not necessary, a Heisenberg row pairing with the seat not blind " +
                         "(simulations/n3_offletter_locus_gate.py, PROOF_ONE_DEPHASED_SEAT's Open list; those rows, E and D, pair through far " +
                         "elements that are not of odd height, of no definite height at the rows read). On the triangle with Heisenberg bonds of which at most one is zero, the jump Z " +
                         "on site 0 alone and nonzero planar fields f1, f2 on sites 1 and 2 neither parallel nor antiparallel, a " +
                         "far element of odd height (every string with an odd number of letters X or Y) exists, with no field on " +
                         "site 0, exactly when site 0 is blind and |f1| = |f2|, and with a nonzero planar field f0 on site 0 exactly when " +
                         "J01 = J02, |f1| = |f2| and f0 lies on the line of f1 + f2 (the paragraph on the sheet and blindness; gate " +
                         "simulations/blind_seat_sheet_gate.py; live: inspect --root trianglelocus --J 3,6,-2 --fu 3,4 --fv 5,0, with --fd x,y " +
                         "for a field on site 0, node 'the sheet and blindness'). In the widened N = 4 family (the pair (1, 2) seen " +
                         "from 0 and 3 through its difference only, b0 = (J01 - J02)/2 and b3 = (J13 - J23)/2, the family's " +
                         "fields, J03 = 0, h, g, b0, b3 nonzero) a far element even under the pair's swap has the sheet IIIY, " +
                         "and one with a nonzero singlet block exists exactly where the pair carries a mode that vanishes " +
                         "at 0 and at 3, b3 = -b0.");

            yield return new InspectableNode("the carriers' locality",
                summary: "Where a row pairs through one far element W, for Pauli jumps, a mirror that factorises over " +
                         "the sites, alone or followed by a site permutation, sends the identity to W, so it exists only " +
                         "where W is a product times a site permutation fixing the watched seat, and then W has rank one " +
                         "across that seat and the rest. A mirror conditioned on the bits of an undephased site needs a " +
                         "letter there that H conserves, a second element of the near end, so where near = 1 none arises. " +
                         "Read at the rows the proof names (C, D, E and the N = 4 family have rank 2 at the watched seat, " +
                         "so neither kind of mirror carries them); live: the operator Schmidt ranks of W across the three " +
                         "cuts at inspect --root trianglelocus. The same lemma in rank form, for a thermal steady state, is " +
                         "MirrorLocalityBoundClaim's (finding 1 of experiments/THE_CUBE_SITE_BY_SITE.md), adjacent, not a parent.");

            yield return new InspectableNode("scope, and what would falsify it",
                summary: "The triangle's loci need Heisenberg bonds with at most one coupling zero and the jump Z on d " +
                         "alone, letter fields for the first and, for the second, no field on d and fields perpendicular to " +
                         "z; a counterexample is a row of the setting that pairs outside the classes or a row of a class " +
                         "that does not pair. The N = 4 proposition falls to a real row of the family with a, g != 0 that " +
                         "does not pair; the conjecture to a palindrome under Heisenberg bonds with letter fields and one " +
                         "dephased seat, without a colouring, whose dephased seat is not blind. The carriers' locality needs " +
                         "Pauli jumps and a row that pairs through one far element; the sheet and blindness, at most one " +
                         "coupling zero, nonzero fields on u and v perpendicular to z, neither parallel nor antiparallel, and " +
                         "no field on d or a planar one; the widened N = 4 rule (h, g, b0, b3 nonzero) falls to a far element even " +
                         "under the pair's swap with a nonzero singlet block at b3 != -b0. " +
                         "Fields off the letters at N = 3 beyond the planar fields with none on d, other bonds at N = 3, a " +
                         "carrier at N = 5 that entangles the fifth spin and every N beyond are open (the proof's Open list).");
        }
    }
}
