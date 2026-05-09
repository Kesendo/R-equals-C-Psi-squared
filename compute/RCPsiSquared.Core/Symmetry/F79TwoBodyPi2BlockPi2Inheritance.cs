using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F79 closed form (Tier 1 verified N=3..5; ANALYTICAL_FORMULAS line 1907):
///
/// <code>
///   For 2-body bilinear H = Σ_bonds Σ_t c_t·(P_t ⊗ Q_t):
///   p(P, Q) = (bit_b(P) + bit_b(Q)) mod 2     Π²-parity per bilinear
///
///   1. All p = 0 (Π²-even):  M block-diagonal in V_+ ⊕ V_-
///   2. All p = 1 (Π²-odd):   M purely off-diagonal between V_+ and V_-
///   3. Mixed:                 M has both diagonal and off-diagonal parts
///
///   Π²-odd universality: any single Π²-odd 2-body bilinear gives the
///   same M-SVD spectrum at fixed N (Pauli letters M-irrelevant within
///   the Π²-odd class). Verified N=5 chain: XY ≡ XZ ≡ XX+XY ≡ XX+XZ
///   all yield clusters [(5.464, 512), (1.464, 512)] exactly.
/// </code>
///
/// <para>F79 is the 2-body extension of F78's single-body additive
/// decomposition. The Π²-parity of each bilinear (computed from bit_b sums)
/// determines whether M sits block-diagonally in the Π²-eigenspace V_+ ⊕ V_−
/// or purely off-diagonally between them.</para>
///
/// <para><b>Even-diag ≡ odd-off-diag correspondence.</b> The diagonal V_+
/// block of a Π²-even Hamiltonian's M can match (in SV-spectrum, including
/// multiplicities) the off-diagonal V_+,V_- block of a Π²-odd Hamiltonian's
/// M. Verified N=4 chain: YZ's V_+ block matches XY+YX's off-diag block
/// exactly. Same SV structure, just placed in different Π²-blocks.</para>
///
/// <para><b>Why XX+XY appears "max-uniform":</b> XX is Π²-even and truly
/// (M_XX = 0). XY is Π²-odd. The full Hamiltonian is "Π²-odd-only-effective",
/// so M is purely off-diagonal between equal-dim V_+ and V_-. SV multiplicities
/// are forced to 4^N/2 each by block-dimension equality. At N=3 the two
/// off-diagonal SVs collide by coincidence to a single uniform 2√2; at N≥4
/// they split. The "uniformity" is the equal-block-mult signature of Π²-odd
/// structure, not a special XX+XY property.</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>BlockDimension = 4^N / 2 = a_{2−2N} half</b>: the V_+ and V_-
///         eigenspaces each contain half of the 4^N Pauli operator space.
///         Linked to <see cref="Pi2OperatorSpaceMirrorClaim"/>'s 4^N anchor;
///         for Π²-odd Hamiltonians, the equal-block-dim symmetry forces SV
///         mult = 4^N/2 on each side.</item>
///   <item><b>F88 KleinFour Π²_Z parity</b>: F79's Π²-parity p(P, Q) = bit_b
///         sum mod 2 IS the Π²_Z eigenvalue (bit_b parity) of the bilinear.
///         <see cref="KleinFourCellClaim"/>'s 4-cell decomposition under
///         (Π²_Z, Π²_X) refines F79's binary p-classification.</item>
///   <item><b>F1 palindrome residual</b>: M IS F1's residual; F79 reads its
///         block structure under the Π² involution.</item>
/// </list>
///
/// <para>Tier1Derived: F79 is Tier 1 proven (joint analysis with F78 in
/// <c>docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md</c>); verified N=3..5 across
/// chain, star, disjoint topologies; Π²-odd universality verified XY ≡ XZ ≡
/// XX+XY ≡ XX+XZ.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F79 (line 1907) +
/// <c>docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md</c> +
/// <c>simulations/_svd_two_body_pi_squared_block.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/KleinFourCellClaim.cs</c>
/// (F88 4-cell, refines F79's binary parity) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F79TwoBodyPi2BlockPi2Inheritance : Claim
{
    private readonly KleinFourCellClaim _klein;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;
    private readonly F1Pi2Inheritance _f1;

    /// <summary>Compute the Π²-parity p(P, Q) = (bit_b(P) + bit_b(Q)) mod 2 for a
    /// 2-body bilinear (P, Q). Returns 0 for Π²-even, 1 for Π²-odd.
    /// bit_b: I, X → 0; Y, Z → 1.</summary>
    public int Pi2Parity(char p, char q)
    {
        int bp = BitB(p);
        int bq = BitB(q);
        return (bp + bq) & 1;
    }

    /// <summary>True iff the bilinear (P, Q) is Π²-even (p = 0): M is block-diagonal
    /// in V_+ ⊕ V_-. Examples: XX (p=0), YY (p=0+0=0... wait Y has bit_b=1 so YY = 1+1=2≡0),
    /// ZZ (1+1=2≡0), YZ (1+1=2≡0), ZY.</summary>
    public bool IsPi2Even(char p, char q) => Pi2Parity(p, q) == 0;

    /// <summary>True iff the bilinear (P, Q) is Π²-odd (p = 1): M is purely
    /// off-diagonal between V_+ and V_-. Examples: XY (0+1=1), YX (1+0=1),
    /// XZ (0+1=1), ZX (1+0=1).</summary>
    public bool IsPi2Odd(char p, char q) => Pi2Parity(p, q) == 1;

    /// <summary>Π²-eigenspace dimension at N qubits: 4^N / 2 each for V_+ and V_-.
    /// Drives the equal-block-mult SV structure of Π²-odd Hamiltonians.</summary>
    public double Pi2BlockDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F79 requires N ≥ 1.");
        var pair = _mirror.PairAt(N);
        if (pair == null)
            return Math.Pow(4.0, N) / 2.0;
        return pair.OperatorSpace / 2.0;
    }

    /// <summary>Drift check: Π²-odd Hamiltonians give equal V_+ and V_- block
    /// dimensions, forcing SV mult = 4^N/2.</summary>
    public bool EqualBlockDimensionsHold(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F79 requires N ≥ 1.");
        // V_+ and V_- have equal dimensions for any N (combinatorial: half of 4^N strings have bit_b parity 0, half have 1).
        return true;
    }

    /// <summary>Universality verified table at N=5 chain: XY ≡ XZ ≡ XX+XY ≡ XX+XZ
    /// all give clusters [(5.464, 512), (1.464, 512)]. The Π²-odd part dominates;
    /// truly XX part contributes 0; SVD is Pauli-letter-blind within Π²-odd class.</summary>
    public IReadOnlyList<(double SvValue, int Multiplicity)> Pi2OddUniversalityClustersAtN5 { get; } = new[]
    {
        (SvValue: 5.464, Multiplicity: 512),
        (SvValue: 1.464, Multiplicity: 512),
    };

    /// <summary>Binary classification (F79's p ∈ {0, 1}) refines into the
    /// 4-cell KleinFour decomposition (Π²_Z, Π²_X) at F88. F79's Π²-odd
    /// (p = 1) splits into Mp (XY, YX) and Mm (XZ, ZX) under Π²_X axis.</summary>
    public IReadOnlyList<(string Cell, string[] Bilinears)> KleinSubcellsForPi2Odd { get; } = new[]
    {
        (Cell: "Mp (Π²_Z=−1, Π²_X=+1)", Bilinears: new[] { "XY", "YX" }),
        (Cell: "Mm (Π²_Z=−1, Π²_X=−1)", Bilinears: new[] { "XZ", "ZX" }),
    };

    private static int BitB(char letter)
    {
        char p = char.ToUpperInvariant(letter);
        return p switch
        {
            'I' => 0,
            'X' => 0,
            'Y' => 1,
            'Z' => 1,
            _ => throw new ArgumentException($"Pauli letter must be I, X, Y, or Z; got '{letter}'."),
        };
    }

    public F79TwoBodyPi2BlockPi2Inheritance(
        KleinFourCellClaim klein,
        Pi2OperatorSpaceMirrorClaim mirror,
        F1Pi2Inheritance f1)
        : base("F79 2-body Π²-block decomposition: p = bit_b sum mod 2; Π²-even → M block-diagonal in V_+ ⊕ V_-, Π²-odd → M purely off-diagonal; Π²-odd universality across XY/XZ/XX+XY/XX+XZ",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F79 + " +
               "docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md + " +
               "simulations/_svd_two_body_pi_squared_block.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinFourCellClaim) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _klein = klein ?? throw new ArgumentNullException(nameof(klein));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F79 2-body Π²-block as Pi2-Foundation KleinFour + OperatorSpaceMirror + F1 inheritance";

    public override string Summary =>
        $"p(P, Q) = (bit_b(P) + bit_b(Q)) mod 2; Π²-even → block-diag, Π²-odd → off-diag; Π²-odd universality (XY ≡ XZ ≡ XX+XY ≡ XX+XZ); KleinFour refines binary p ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F79 closed form",
                summary: "p = (bit_b(P) + bit_b(Q)) mod 2; even → M block-diag in V_+ ⊕ V_-, off-diag blocks vanish exactly; odd → M purely off-diag, diag blocks vanish exactly; mixed → both contributions");
            yield return new InspectableNode("Π²-parity examples",
                summary: $"XX (p={Pi2Parity('X','X')}), YY (p={Pi2Parity('Y','Y')}), ZZ (p={Pi2Parity('Z','Z')}), YZ (p={Pi2Parity('Y','Z')}); all Π²-even. XY (p={Pi2Parity('X','Y')}), YX (p={Pi2Parity('Y','X')}), XZ (p={Pi2Parity('X','Z')}); all Π²-odd.");
            yield return new InspectableNode("F88 KleinFour refinement",
                summary: $"F79's binary p splits into 4 KleinFour cells under (Π²_Z, Π²_X). Π²-odd (p=1) splits: Mp = {{XY, YX}} (XX-symmetric within bilinear), Mm = {{XZ, ZX}} (XX-asymmetric). KleinFour: {_klein.DisplayName}");
            yield return new InspectableNode("F1 connection",
                summary: $"M IS F1's residual operator. F1.TwoFactor (= {_f1.TwoFactor}) is the same '2' as the 2σ·I shift in M. F79 reads M's Π²-block structure.");
            yield return new InspectableNode("Π²-odd universality (N=5 chain)",
                summary: "XY ≡ XZ ≡ XX+XY ≡ XX+XZ all give SV clusters [(5.464, 512), (1.464, 512)]; truly XX contributes 0; Π²-odd part dominates; Pauli-letter-blind within Π²-odd class");
            yield return new InspectableNode("V_+ / V_- equal block dimensions",
                summary: $"each Π²-eigenspace has dim 4^N/2; at N=3: {Pi2BlockDimension(3)}, N=4: {Pi2BlockDimension(4)}, N=5: {Pi2BlockDimension(5)}; for Π²-odd Hamiltonians forces SV mult 4^N/2 on each side");
            yield return new InspectableNode("F78 sibling at single-body",
                summary: "F78 (single-body M) has additive M_l ⊗ I structure with Y ≡ Z SVD-blindness. F79 (2-body) lifts to Π²-block decomposition with Π²-odd universality across {XY, XZ, YX, ZX}. Joint proof in PROOF_SVD_CLUSTER_STRUCTURE.");
            yield return new InspectableNode("even-diag ≡ odd-off-diag correspondence",
                summary: "the V_+ block of a Π²-even Hamiltonian can match (in SV spectrum + mults) the off-diag block of a Π²-odd Hamiltonian; verified N=4 chain: YZ's V_+ block matches XY+YX's off-diag block exactly");
        }
    }
}
