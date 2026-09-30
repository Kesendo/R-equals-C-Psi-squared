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
/// third letter c, carried by c^⊗N alone; on every broken row the far end is zero. Theorems 1 and 3 are the
/// two classes in which F138's converse holds exactly.</para>
///
/// <para>Gate: <c>compute/MirrorWorld/EndCount.cs</c>'s <c>ComplementConnection</c>, held in
/// <c>compute/MirrorWorld.Tests/EndCountTests.cs</c> against the end count's exact verdict. Live lab:
/// <c>inspect --root complement</c> (<c>ComplementConnectionWitness</c>), which shares no code with the
/// gate and reads <c>twoendstrings</c> beside itself.</para></summary>
public sealed class PalindromeComplementConnectionClaim : Claim
{
    /// <summary>The typed parent: this claim is F158 read on one class of jumps, both ends of which
    /// become diagonal problems on the bitstrings once every jump is a Z.</summary>
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

    public PalindromeComplementConnectionClaim(PalindromeTwoEndCountClaim twoEndCount)
        : base("With every site dephased the palindrome is a flat connection: with one single-site Pauli jump on " +
               "every site, turned to Z, the far end of F158 is X^N times the flat sections of d_y/d_x = " +
               "H_xy/H_x̄ȳ on the hopping graph over bitstrings; near count = components of the hopping graph, " +
               "far count = good components of its union with the complement image, palindrome iff every " +
               "component is good (Theorem 2). Corollaries, the two classes where F138's converse holds exactly: " +
               "hopping bonds J(XX + YY), every J nonzero, with a complement-invariant diagonal part, under uniform Z with letter " +
               "fields, palindrome iff no Z field and not both X and Y " +
               "(Theorem 1); Heisenberg bonds, every J nonzero, under mixed axes with fields of any direction, near end " +
               "one-dimensional, palindrome iff two axes occur and every field lies along the third letter " +
               "(Theorem 3). Both corollaries need a connected graph",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md (primary: Theorems 1-3) + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (F158, section (f10)) + " +
               "docs/ANALYTICAL_FORMULAS.md (F138, whose converse holds exactly in the two classes) + " +
               "experiments/THE_PALINDROME_AS_A_COLOURING.md (the colouring rule and the frame theorem)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "With every site dephased, the palindrome is a flat connection (the complement connection)";

    public override string Summary =>
        "near = components of the hopping graph, far = its good components joined with the complement image; " +
        $"F138's converse exact under one axis (Theorem 1) and under mixed axes (Theorem 3) ({Tier.Label()})";

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

            yield return new InspectableNode("scope, and the counterexamples at its edges",
                summary: "Every site dephased by one single-site Pauli jump; an undephased site lies outside " +
                         "(the near span stops being abelian). Theorems 1 and 3 need a connected graph: two " +
                         "disjoint Heisenberg bonds with axes X, Y and Y, Z and fields along each bond's third " +
                         "letter use three axes and pair, by ZZXX. Theorem 3 needs isotropic bonds; XXZ and XY " +
                         "under mixed axes are outside it.");

            yield return new InspectableNode("live lab (the witness)",
                summary: "ComplementConnectionWitness turns the jumps on PauliMask, walks the 2^N bitstrings over " +
                         "the Gaussian integers, counts the hopping and good components, prints Theorem 1's or " +
                         "Theorem 3's prediction beside the graph's verdict, and reads twoendstrings on the same " +
                         "row: inspect --root complement.");
        }
    }
}
