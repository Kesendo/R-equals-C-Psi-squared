using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>F1PalindromeOrbitPairing (Tier 1 derived; 2026-05-22): the F1 palindrome
/// conjugation Π is an order-4 operator that, on the joint-popcount sector labels
/// (p_c, p_r), acts as the whole-sector permutation (p_c, p_r) ↦ (N − p_r, p_c). It
/// strictly subsumes the X⊗N sector pairing of <see cref="XGlobalChargeConjugationPairing"/>:
/// Π² = X⊗N sends (p_c, p_r) ↦ (N − p_c, N − p_r), the X⊗N rule.
///
/// <para>F1 (<c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>) is the identity
/// <c>Π·L·Π⁻¹ = −L − 2Σγ·I</c>: a sector and its Π-image have spectra related by the
/// reflection <c>λ ↦ −2Σγ − λ</c>. The whole-sector permutation organises the (N+1)²
/// joint-popcount sectors into Π-orbits, each of size 4 except the single Π-fixed sector
/// (N/2, N/2) present at even N. One eigendecomposition per orbit therefore feeds three
/// follower sectors: the Π²-image (X⊗N partner) by a verbatim spectrum copy, and the Π-
/// and Π³-images by the F1 reflection.</para>
///
/// <para>Algorithmic gain over <see cref="XGlobalChargeConjugationPairing"/>: where X⊗N
/// alone halves the number of distinct eigendecompositions, the Π-orbit quarters it (one
/// per orbit-of-4 instead of one per X⊗N-pair-of-2). Block sizes are unchanged; this is
/// sector-orbiting, not sector-splitting. The Π-orbit subsumes the X⊗N pair, so a builder
/// wiring this primitive replaces the X⊗N partition rather than layering on top of it;
/// the X⊗N <see cref="XGlobalChargeConjugationPairing.PairSector"/> rule is reused here to
/// express and check the Π² sub-map internally.</para>
///
/// <para>Π-fixedness: (p_c, p_r) = (N − p_r, p_c) forces p_c = N − p_r and p_r = p_c, hence
/// p_c = p_r = N/2. Only the central sector at even N is Π-fixed; it is also X⊗N-self-
/// paired, so a builder gets no F1 gain there, matching the prior X⊗N behaviour exactly.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> (F1; Π² = X⊗N is F1²);
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> (the Π·L·Π⁻¹ = −L − 2Σγ·I proof);
/// <c>compute/RCPsiSquared.Core/Symmetry/PiOperator.cs</c> (Π in the 4^N Pauli-string
/// basis, order-4 with Π² = X⊗N); <c>compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs</c>
/// (the X⊗N pairing this orbit subsumes).</para></summary>
public sealed class F1PalindromeOrbitPairing : Claim
{
    private readonly SymmetryFamilyInventory _inventory;

    public F1PalindromeOrbitPairing(SymmetryFamilyInventory inventory)
        : base("F1PalindromeOrbitPairing: the F1 palindrome Π permutes joint-popcount sectors as the whole-sector cycle (p_c, p_r) ↦ (N − p_r, p_c) in orbits of size 4; one eigendecomposition per orbit feeds the X⊗N follower (verbatim) and the two Π/Π³ followers (reflected through λ ↦ −2Σγ − λ).",
               Tier.Tier1Derived,
               "docs/proofs/MIRROR_SYMMETRY_PROOF.md (the F1 identity Π·L·Π⁻¹ = −L − 2Σγ·I); docs/ANALYTICAL_FORMULAS.md (F1 and the F1² corollary Π² = X⊗N); Π is order-4 with Π² = X⊗N so the whole-sector cycle subsumes the XGlobalChargeConjugationPairing X⊗N pair; verified bit-exact vs dense and vs un-halved per-block in BlockSpectrumF1OrbitPairingTests")
    {
        _inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
    }

    /// <summary>How a follower sector's spectrum is obtained from its orbit primary.</summary>
    public enum F1FollowerKind
    {
        /// <summary>The Π²-image (the X⊗N partner) of the primary: spectrum copied verbatim,
        /// since X⊗N is a genuine symmetry (X⊗N·L·X⊗N⁻¹ = L).</summary>
        XnCopy = 0,
        /// <summary>The Π- or Π³-image of the primary: spectrum reflected through
        /// <c>λ ↦ −2Σγ − λ</c>, the F1 palindrome map.</summary>
        F1Reflect = 1,
    }

    /// <summary>A follower sector: an index into the primary list together with the rule
    /// (<see cref="F1FollowerKind"/>) for deriving its spectrum from that primary's.</summary>
    public readonly record struct F1Follower(int PrimaryIndex, F1FollowerKind Kind);

    /// <summary>Π-image of joint-popcount sector (p_c, p_r): (N − p_r, p_c). Order 4;
    /// applying it four times is the identity. Π² is the X⊗N rule (N − p_c, N − p_r).</summary>
    public static (int PCol, int PRow) PiImage(int N, int pCol, int pRow)
    {
        if (N < 0 || pCol < 0 || pCol > N || pRow < 0 || pRow > N)
            throw new ArgumentOutOfRangeException(
                $"N={N}, pCol={pCol}, pRow={pRow}: require 0 ≤ pCol, pRow ≤ N");
        return (N - pRow, pCol);
    }

    /// <summary>True if sector (p_c, p_r) is Π-fixed (its own singleton orbit). The fixedness
    /// condition (p_c, p_r) = (N − p_r, p_c) reduces to 2·p_c == N AND 2·p_r == N, so only the
    /// central (N/2, N/2) sector at even N is Π-fixed; at odd N no sector is.</summary>
    public static bool IsPiFixed(int N, int pCol, int pRow)
    {
        if (N < 0 || pCol < 0 || pCol > N || pRow < 0 || pRow > N)
            throw new ArgumentOutOfRangeException(
                $"N={N}, pCol={pCol}, pRow={pRow}: require 0 ≤ pCol, pRow ≤ N");
        return 2 * pCol == N && 2 * pRow == N;
    }

    /// <summary>Number of distinct spectral classes after Π-orbit grouping at given N: one
    /// eigendecomposition per orbit. At odd N every orbit has size 4, so the count is
    /// (N+1)²/4. At even N the (N+1)² sectors split into the single Π-fixed (N/2, N/2)
    /// sector plus ((N+1)² − 1)/4 orbits of size 4, giving ((N+1)² − 1)/4 + 1.
    ///
    /// <para>This is at most the <see cref="XGlobalChargeConjugationPairing.DistinctSpectralClasses"/>
    /// count: the Π-orbit groups four sectors where X⊗N groups two.</para></summary>
    public static int DistinctSpectralClasses(int N)
    {
        int total = (N + 1) * (N + 1);
        // Odd N: every orbit has size 4, total is divisible by 4.
        // Even N: exactly one Π-fixed sector; the remaining total - 1 sectors form orbits of 4.
        return (N % 2 == 0) ? (total - 1) / 4 + 1 : total / 4;
    }

    /// <summary>Partition a list of joint-popcount sectors into "primary" sectors (compute
    /// eig) and "follower" sectors (derive spectrum from the orbit primary). Π-fixed sectors
    /// are always primary; in each size-4 orbit the lex-smallest (PCol, PRow) sector is the
    /// primary, the other three are followers.
    ///
    /// <para>Returns indices into the input sectors list. The dictionary maps follower-index
    /// → <see cref="F1Follower"/>: the primary's index plus the derivation rule. The Π²-image
    /// of the primary is an <see cref="F1FollowerKind.XnCopy"/> follower (k = 2 applications
    /// of <see cref="PiImage"/> from the primary); the Π- and Π³-images are
    /// <see cref="F1FollowerKind.F1Reflect"/> followers (k ∈ {1, 3}).</para>
    ///
    /// <para>Optional <paramref name="sectorSize"/> projection: when supplied, primary
    /// sectors are sorted descending by size so the most expensive eig starts first under
    /// <see cref="Parallel.ForEach"/>, overlapping the largest block's wall-time with smaller
    /// blocks running concurrently. Mirrors
    /// <see cref="XGlobalChargeConjugationPairing.PartitionByXNPairing{TSector}"/>.</para></summary>
    public static (List<int> Primaries, Dictionary<int, F1Follower> FollowerToPrimary) PartitionByPiOrbit<TSector>(
        int N,
        IReadOnlyList<TSector> sectors,
        Func<TSector, (int PCol, int PRow)> getLabel,
        Func<TSector, int>? sectorSize = null)
    {
        if (sectors is null) throw new ArgumentNullException(nameof(sectors));
        if (getLabel is null) throw new ArgumentNullException(nameof(getLabel));

        var sectorIndexByLabel = new Dictionary<(int, int), int>(sectors.Count);
        for (int i = 0; i < sectors.Count; i++)
            sectorIndexByLabel[getLabel(sectors[i])] = i;

        var primaries = new List<int>(sectors.Count);
        var followerToPrimary = new Dictionary<int, F1Follower>();
        for (int i = 0; i < sectors.Count; i++)
        {
            var (pCol, pRow) = getLabel(sectors[i]);
            if (IsPiFixed(N, pCol, pRow))
            {
                primaries.Add(i);
                continue;
            }

            // Walk the size-4 Π-orbit {label, Π·label, Π²·label, Π³·label} and pick the
            // lex-smallest label as the orbit primary.
            var orbit = new (int PCol, int PRow)[4];
            orbit[0] = (pCol, pRow);
            for (int k = 1; k < 4; k++)
                orbit[k] = PiImage(N, orbit[k - 1].PCol, orbit[k - 1].PRow);

            int primaryK = 0;
            for (int k = 1; k < 4; k++)
            {
                bool kSmaller = orbit[k].PCol < orbit[primaryK].PCol
                    || (orbit[k].PCol == orbit[primaryK].PCol && orbit[k].PRow < orbit[primaryK].PRow);
                if (kSmaller) primaryK = k;
            }

            if (primaryK == 0)
            {
                // This sector is the orbit primary.
                primaries.Add(i);
            }
            else
            {
                // This sector is a follower; its offset from the primary in the 4-cycle is
                // (0 - primaryK) mod 4. k = 2 ⇒ the primary's Π²-image (X⊗N partner) ⇒
                // XnCopy; k ∈ {1, 3} ⇒ Π/Π³-image ⇒ F1Reflect.
                int offsetFromPrimary = (4 - primaryK) % 4;
                var primaryLabel = orbit[primaryK];
                int primaryIndex = sectorIndexByLabel[primaryLabel];
                var kind = offsetFromPrimary == 2
                    ? F1FollowerKind.XnCopy
                    : F1FollowerKind.F1Reflect;
                followerToPrimary[i] = new F1Follower(primaryIndex, kind);
            }
        }

        if (sectorSize is not null)
            primaries.Sort((a, b) => sectorSize(sectors[b]).CompareTo(sectorSize(sectors[a])));

        return (primaries, followerToPrimary);
    }

    public override string DisplayName =>
        "F1PalindromeOrbitPairing: Π orbits joint-popcount sectors (p_c, p_r) ↦ (N − p_r, p_c) in 4-cycles; quarters the eig-call count";

    public override string Summary =>
        $"F1 Π-orbit sector-grouping under chain XY+Z-deph; (N+1)² sectors collapse to ≈ (N+1)²/4 distinct spectral classes ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Π-image formula",
                summary: "(p_c, p_r) ↦ (N − p_r, p_c); order 4");
            yield return new InspectableNode("Π² = X⊗N",
                summary: "(p_c, p_r) ↦ (N − p_c, N − p_r): the XGlobalChargeConjugationPairing rule");
            yield return new InspectableNode("follower derivation",
                summary: "Π²-image: spectrum copied verbatim; Π/Π³-image: reflected through λ ↦ −2Σγ − λ");
            yield return new InspectableNode("N=8 distinct classes",
                summary: $"{DistinctSpectralClasses(8)} (vs {XGlobalChargeConjugationPairing.DistinctSpectralClasses(8)} for X⊗N alone, 81 unpaired)");
            yield return new InspectableNode("N=10 distinct classes",
                summary: $"{DistinctSpectralClasses(10)} (vs {XGlobalChargeConjugationPairing.DistinctSpectralClasses(10)} for X⊗N alone, 121 unpaired)");
        }
    }
}
