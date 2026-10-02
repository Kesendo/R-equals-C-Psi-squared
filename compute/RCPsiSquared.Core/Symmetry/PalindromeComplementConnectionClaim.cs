using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>With every site dephased, the palindrome is a flat connection
/// (<c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>). With exactly one single-site Pauli jump
/// on every site and H a real combination of Pauli strings, turn each jump to Z by a proper rotation of its
/// site's letters; then F158's near end is the diagonal matrices commuting with H and its far end is
/// X^⊗N times the diagonals D with
///
/// <code>
///     H_x̄ȳ · d_y = d_x · H_xy     for all bitstrings x, y     (x̄ the bitwise complement)
/// </code>
///
/// <para><b>Theorem 2.</b> dim ker L is the number of components of the hopping graph Γ_H (an edge x–y
/// where H_xy ≠ 0), dim ker(L + 2σ) the number of GOOD components of Γ = Γ_H joined with its image under
/// x ↦ x̄ (equal moduli on every edge, H_xx = H_x̄x̄, trivial holonomy of d_y/d_x = H_xy/H_x̄ȳ), and the
/// palindrome holds exactly when every component of Γ is good. Colourings are the flat sections that are
/// characters.</para>
///
/// <para><b>Theorem 1.</b> Hopping bonds J(XX + YY), every J ≠ 0, with any complement-invariant diagonal
/// part (Heisenberg, XXZ, XY), on a connected graph, Z on every site, letter fields of any magnitudes:
/// the palindrome holds exactly when no field lies along Z and the fields do not use both X and Y.</para>
///
/// <para><b>Theorem 3.</b> Heisenberg bonds, every J ≠ 0, on a connected graph, every site dephased along
/// one letter with at least two letters among the axes, fields of any direction: the near end is
/// one-dimensional, and the palindrome holds exactly when two axes occur and every field lies along the
/// third letter c, carried by c^⊗N alone; on every broken row the far end is zero. Theorems 1, 3 and 5 are the
/// three classes in which F138's converse holds exactly.</para>
///
/// <para><b>Theorem 4.</b> At least one site dephased, every dephased site with one single-site Pauli jump,
/// and every undephased site keeping a letter H conserves: turned to Z, H is block
/// diagonal in their bits, H = ⊕_σ H_σ, and near = Σ_{σ,τ} hom(H_τ, H_σ), far = Σ_{σ,τ} hom(H_τ, H̄_σ),
/// each hom(A, B) = dim{g : B_xy·g_y = g_x·A_xy} a count of good components.</para>
///
/// <para><b>Theorem 5.</b> Heisenberg bonds on a connected graph, letter fields, exactly one undephased site:
/// the palindrome holds exactly when a colouring exists, F138's converse a third time; with two undephased
/// sites it fails (F138 (b), the SWAP rows).</para>
///
/// <para><b>Theorem 6.</b> Undephased sites that keep no letter, wherever each component has a spanning tree of
/// invertible blocks H_xy of H over the dephased bits (for H and for its complement image, on the same support):
/// the connection is a local system; transported along a spanning tree to one base point per component, both ends are the commutant and
/// the intertwiners of finitely many generators of size 2^|F|, and the palindrome holds exactly when the
/// generator tuples of H and of its complement image are simultaneously similar, decided by traces of words,
/// of length at most 3 with one undephased site.</para>
///
/// <para><b>Theorem 7.</b> H conserving a nonconstant diagonal charge Q = Σ w_l Z_l of the frame, dephased and
/// undephased bits mixed or not (the total magnetization of XXZ bonds with Z jumps and Z fields): both ends are sums over pairs
/// of charge sectors, near = Σ dim Hom(ρ_q', ρ_q), far = Σ dim Hom(ρ_q', ρ̄_q), a count of good components where
/// every sector space is at most one-dimensional; Theorem 4 is a charge on the undephased sites alone.</para>
///
/// <para><b>Common triangle carrier.</b> With jump Z on site 0, Heisenberg couplings s=J01, t=J02, r=J12,
/// all real and nonzero, and field h(X1+Y2), the common Z0-odd commutant of the bond Hamiltonian and X1+Y2
/// is a line iff (s-t)(st+sr+tr)=0, zero otherwise. Its generator squares to (s²+t²)³ I. This classifies
/// one invertible carrier independent of h, not every fixed-h palindrome. Exact symbolic gate:
/// <c>simulations/triangle_common_carrier_gate.py</c>.</para>
///
/// <para>Gate: <c>compute/MirrorWorld/EndCount.cs</c>'s <c>ComplementConnection</c> (Theorems 1 to 3) and
/// <c>SectorConnection</c> (Theorem 4) are held against exact end-count verdicts and certified counts;
/// <c>LocalSystem</c> (Theorem 6) and <c>ChargeSectors</c> (Theorem 7) compare modular upper counts.
/// These gates live in <c>compute/MirrorWorld.Tests/EndCountTests.cs</c>; Theorem 5 compares the exact verdict directly. Live lab: <c>inspect --root complement</c>
/// (<c>ComplementConnectionWitness</c>), which shares no code with the gate, reads Theorems 2 and 4 on its own
/// graph (with <c>--bonds ZZ</c> for undephased sites that keep a letter), prints Theorem 1's, 3's or 5's rule
/// beside the reading, and reads <c>twoendstrings</c> beside itself. The exact routine is
/// <see cref="ComplementConnectionGraph"/>, which <c>PalindromeSoftCertifier.DecideAtN</c> uses too, deciding
/// any real template set under Z dephasing at one N where the N-free strategies leave a row undetermined.</para></summary>
public sealed class PalindromeComplementConnectionClaim : Claim
{
    /// <summary>The typed parent: this claim is F158 read where every jump is a single-site letter, turned to
    /// Z: both ends become problems on the bitstrings of the dephased sites, scalar where every site is
    /// dephased or every undephased site keeps a letter.</summary>
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>Theorem 3's rule for Heisenberg bonds on a connected graph, every site dephased along
    /// mixed axes (at least two distinct letters among <paramref name="axes"/>): the palindrome holds
    /// exactly when two axes occur and every field letter is the third. Field letters are given as a
    /// set; a field of general direction enters as the letters of its nonzero components.</summary>
    public static bool MixedAxesPalindrome(IEnumerable<char> axes, IEnumerable<char> fieldLetters)
    {
        var a = axes.Select(char.ToUpperInvariant).Distinct().ToList();
        if (a.Count < 2 || a.Any(c => c is not ('X' or 'Y' or 'Z')))
            throw new ArgumentException("Theorem 3 needs at least two distinct axes among X, Y, Z.", nameof(axes));
        var free = "XYZ".Where(c => !a.Contains(c)).ToList();
        return free.Count == 1 && fieldLetters.Select(char.ToUpperInvariant).All(f => f == free[0]);
    }

    /// <summary>Theorem 5's rule for Heisenberg bonds on a connected graph with exactly one undephased site
    /// and letter fields: the palindrome holds exactly when some letter p is the axis of no dephased site
    /// and every field letter, the undephased site's included, is p. With every site dephased under one
    /// common axis it is Theorem 1's rule, and under mixed axes Theorem 3's.</summary>
    public static bool OneUndephasedSitePalindrome(IEnumerable<char> axes, IEnumerable<char> fieldLetters)
    {
        var a = axes.Select(char.ToUpperInvariant).Distinct().ToList();
        if (a.Count == 0 || a.Any(c => c is not ('X' or 'Y' or 'Z')))
            throw new ArgumentException("Theorem 5 needs at least one dephased site with a letter axis.", nameof(axes));
        var f = fieldLetters.Select(char.ToUpperInvariant).Distinct().ToList();
        return "XYZ".Any(p => !a.Contains(p) && f.All(l => l == p));
    }

    public PalindromeComplementConnectionClaim(PalindromeTwoEndCountClaim twoEndCount)
        : base("With every site dephased the palindrome is a flat connection: with one single-site Pauli jump on " +
               "every site, turned to Z, the far end of F158 is X^N times the flat sections of d_y/d_x = " +
               "H_xy/H_x̄ȳ on the hopping graph over bitstrings; near count = components of the hopping graph, " +
               "far count = good components of its union with the complement image, palindrome iff every " +
               "component is good (Theorem 2). Corollaries, two of the three classes where F138's converse holds exactly: " +
               "hopping bonds J(XX + YY), every J nonzero, with a complement-invariant diagonal part, under uniform Z with letter " +
               "fields, palindrome iff no Z field and not both X and Y " +
               "(Theorem 1); Heisenberg bonds, every J nonzero, under mixed axes with fields of any direction, near end " +
               "one-dimensional, palindrome iff two axes occur and every field lies along the third letter " +
               "(Theorem 3). Both need a connected graph. With the dephased sites (at least one) each under one single-site " +
               "jump and every undephased site keeping a letter H conserves: H splits into " +
               "sectors of their bits, near = sum over sector pairs of hom(H_tau, H_sigma), far = the same against the complement " +
               "image, each hom a count of good components (Theorem 4). Heisenberg bonds on a connected graph with exactly one " +
               "undephased site and letter fields: palindrome iff a colouring exists (Theorem 5); with two it fails (F138 (b)). " +
               "On the Heisenberg triangle with Z0, field h(X1+Y2), and real nonzero couplings s=J01, t=J02, r=J12, " +
               "a common invertible carrier independent of h exists iff (s-t)(st+sr+tr)=0; its space is a line.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md (Theorems 1-7 and common triangle carrier) + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (F158, section (f10)) + " +
               "docs/ANALYTICAL_FORMULAS.md (F138, whose converse holds exactly in three classes) + " +
               "experiments/THE_PALINDROME_AS_A_COLOURING.md (the colouring rule and the frame theorem)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "With every site dephased, the palindrome is a flat connection (the complement connection)";

    public override string Summary =>
        "near = components of the hopping graph, far = its good components joined with the complement image; " +
        "F138's converse exact under one axis (Theorem 1), under mixed axes (Theorem 3) and with one undephased site " +
        $"under Heisenberg bonds (Theorem 5); undephased sites that keep a letter split into sectors (Theorem 4) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("why every site dephased makes the ends diagonal",
                summary: "F158's second corollary to Lemma 1 puts the near end on the strings commuting with every " +
                         "jump and the far end on those anticommuting. With a Z on every site the first are the Z " +
                         "strings, which span the diagonal matrices, and the second are X^N times them. So both " +
                         "ends are equations on functions over the 2^N bitstrings: constancy along the edges of " +
                         "the hopping graph at the near end, the complement connection at the far end.");

            yield return new InspectableNode("Lemma A, and why mixed axes leave one near element",
                summary: "Under mixed axes the turn to Z makes a bond between sites of different axes flip both " +
                         "bits on every basis state (its third-letter term alone) and flip single bits somewhere; " +
                         "same-axis bonds swap bits. Swaps and pair flips connect each popcount-parity class " +
                         "(a leaf induction over the clusters of same-axis bonds), a single flip joins the two, " +
                         "and the hopping graph is connected: dim ker L = 1, so dim ker(L + 2σ) is at most 1.");

            yield return new InspectableNode("from a flat section to a colouring (Theorem 3)",
                summary: "Pair flips and swaps make the section a character times one constant per parity; a " +
                         "single flip seen from both parities makes the constant a sign (at N = 2 the constant " +
                         "d_x·d_x̄ does it, and on the three-axis triangle the pair flips close a cycle of " +
                         "holonomy −1, the sign X^⊗3 gives A_1·A_2·A_3). A character is one lit string commuting " +
                         "with H, a colouring, and a colouring of a connected Heisenberg graph is c^⊗N with c no " +
                         "axis and every field along c.");

            yield return new InspectableNode("undephased sites that keep a letter (Theorem 4)",
                summary: "When every undephased site keeps a letter that H conserves, turning those letters to Z makes H " +
                         "block diagonal in their bits, H = sum over sectors of H_sigma. Both ends split into scalar " +
                         "problems read between two sectors: near = sum of hom(H_tau, H_sigma), far = sum of " +
                         "hom(H_tau, bar H_sigma), each hom the good components of the graph of the two blocks. The " +
                         "cross-sector terms are what an undephased site adds. A Heisenberg or XX + YY bond touching the " +
                         "site keeps no letter, and there the connection is matrix valued: Theorem 6 reads it where a " +
                         "spanning tree carries invertible blocks. Gate: " +
                         "EndCount.SectorConnection in MirrorWorld; live: inspect --root complement --bonds ZZ with " +
                         "undephased sites, which prints the cross-sector part apart.");

            yield return new InspectableNode("a common carrier for the triangle's whole field line",
                summary: "Heisenberg triangle, jump Z0 only, field h(X1+Y2), real nonzero s=J01, t=J02, r=J12: " +
                         "the common odd commutant of the bond Hamiltonian and X1+Y2 is a line exactly when " +
                         "(s-t)(st+sr+tr)=0, including the branch intersection, and zero otherwise. Its generator " +
                         "squares to (s^2+t^2)^3 I. This classifies one invertible carrier independent of h, not " +
                         "every fixed-h palindrome. Exact symbolic gate: simulations/triangle_common_carrier_gate.py.");

            yield return new InspectableNode("scope, and the counterexamples at its edges",
                summary: "Theorems 1 to 3: every site dephased by one single-site Pauli jump. Theorem 4: undephased " +
                         "sites where each keeps a letter H conserves; Theorem 5: one undephased site under Heisenberg bonds with " +
                         "letter fields, where no letter is kept; Theorem 6: any undephased sites, wherever every component has " +
                         "a spanning tree whose blocks of H and of its complement image are invertible. A component without " +
                         "such a tree lies outside Theorem 6; Theorem 7 counts every row whose conserved diagonal charge has " +
                         "sectors of dimension at most one, and only splits the others. Theorems 1, 3 and 5 need a connected graph: two " +
                         "disjoint Heisenberg bonds with axes X, Y and Y, Z and fields along each bond's third " +
                         "letter use three axes and pair, by ZZXX. Theorem 3 needs isotropic bonds; XXZ and XY " +
                         "under mixed axes are outside it.");

            yield return new InspectableNode("live lab (the witness)",
                summary: "ComplementConnectionWitness turns the jumps (and the letters undephased sites keep) on " +
                         "PauliMask, walks the 2^N bitstrings over the Gaussian integers, counts the hopping and good " +
                         "components (Theorem 2) or the sector sums of hom counts (Theorem 4), prints Theorem 1's, 3's " +
                         "or 5's rule beside the verdict (Theorem 5's beside the string route, since no graph reads a " +
                         "site that keeps no letter), and reads twoendstrings on the same row: inspect --root complement.");
        }
    }
}
