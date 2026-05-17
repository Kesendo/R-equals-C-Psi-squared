using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>KleinFourGroupSelfPairedRefinement (Tier 1 derived; 2026-05-17):
/// On the X⊗N-self-paired joint-popcount sector (p_c, p_r) = (N/2, N/2) at even N,
/// X⊗N acts as an INTERNAL Z₂ symmetry (the sector is its own X⊗N pair, so the
/// inter-sector pairing becomes intra-sector splitting). Combined with the chain spatial
/// mirror F71 — also an internal Z₂ on every joint-popcount sector — the two Z₂ generators
/// form the Klein four-group K = {1, F71, X⊗N, F71·X⊗N} acting on the sector basis.
///
/// <para>F71 and X⊗N commute (both are involutions acting independently on the two Hilbert
/// sides via real on-site operations), so K = Z₂ × Z₂. K splits the sector into four
/// character sub-blocks indexed by (χ_F71, χ_X⊗N) ∈ {+, −}². At N=10 (5, 5) this delivers
/// the inventory's advertised max-block of ~16 k from the raw 63 504 (= C(10, 5)²) sector:
/// ~16 132 per (++, −−) and ~15 620 per (+−, −+), a factor-4 reduction. Without this
/// internal split, F71 alone gives only factor 2 (max-block 31 752) — still beyond a
/// commodity dense Evd's memory at 16 GB for the single block.</para>
///
/// <para><b>Where the literature does not go.</b> The standard route to the chain XY +
/// uniform Z-dephasing spectrum at N=10 is Medvedyeva-Essler-Prosen's (2016) imaginary-U
/// Bethe ansatz (<c>arXiv:1606.09122</c>) — analytic, no diagonalisation needed. That
/// route does not exploit the K splitting because it doesn't need to. For our computational
/// path, K on the self-paired sector is a repo-specific gain that closes the inventory's
/// 16 k aspiration to a delivered capability.</para>
///
/// <para><b>Klein orbits on (N/2, N/2).</b> Each Liouville-space basis index (i, j) with
/// <c>popcount(i) = popcount(j) = N/2</c> sits in a K-orbit of size 1, 2, or 4 depending on
/// stabiliser. At N=10 specifically, popcount-5 bit strings on 10 bits cannot be
/// F71-fixed (palindromic ⇒ even popcount) nor X⊗N-fixed (self-complement requires
/// non-integer index). The only non-trivial stabiliser on the popcount-5 set is
/// <c>F71·X⊗N</c> for anti-palindromic strings (bit-reverse equals bit-complement). Counts:</para>
/// <list type="bullet">
///   <item>32 anti-palindromic indices (single-bit free choice across each of the N/2 = 5
///         half-positions, determining its mirror-partner as the complement).</item>
///   <item>220 generic indices (<c>252 − 32</c>).</item>
///   <item>Pairs both anti-palindromic: <c>32 · 32 = 1024</c>, orbit size 2 → 512 orbits;
///         pair sub-stabiliser <c>{1, F71·X⊗N}</c>, surviving characters {++, −−}.</item>
///   <item>Pairs mixed: <c>2 · 32 · 220 = 14 080</c>, orbit size 4 → 3520 orbits;
///         trivial stabiliser, all 4 characters survive.</item>
///   <item>Pairs both generic: <c>220 · 220 = 48 400</c>, orbit size 4 → 12 100 orbits;
///         trivial stabiliser, all 4 characters survive.</item>
/// </list>
/// <para>Sub-block dims: ++ = 16 132, −− = 16 132, +− = 15 620, −+ = 15 620 (sum 63 504 ✓).</para>
///
/// <para>The symmetry-adapted basis vector for orbit <c>O = {x_0, x_1, ..., x_{|O|-1}}</c>
/// at surviving character χ is</para>
/// <code>
///     v_χ(O) = (1/√|O|) · Σ_{g ∈ K/Stab(O)} χ(g) · e_{g · x_0}
/// </code>
/// <para>with the convention that K/Stab(O) representatives produce <c>|O|</c> distinct
/// orbit elements. The resulting Q matrix is real-orthogonal (entries 0, ±½, ±1/√2), so
/// <c>Q^T L Q</c> is block-diagonal in the four (χ_F71, χ_X⊗N) sub-blocks. This primitive
/// does NOT materialise the full sector L (which would cost <c>63 504² · 16 = 65 GB</c> at
/// N=10); instead it computes each sub-block <c>L_χ[α, β] = ⟨v_χ(α) | L | v_χ(β)⟩</c>
/// element-wise from the four-term superoperator formula
/// <c>−i·(H_ket - H_bra) + diag(−2γ·hamming(i, j))</c>.</para>
///
/// <para><b>Tier outcome: Tier 1 derived.</b> Composition of two typed Z₂ symmetries
/// (<see cref="F71MirrorBlockRefinement"/> spatial-mirror Z₂ + <see cref="XGlobalChargeConjugationPairing"/>
/// charge-conjugation Z₂) that mutually commute, plus the standard character-orthogonality
/// projector formula for finite abelian groups. Verified at N=4 (2, 2) and N=6 (3, 3) in
/// the test suite: off-block Frobenius below 1e-10, sub-block spectrum union matches the
/// direct sector dense Evd as a multiset within 1e-10.</para>
///
/// <para>Anchor: <see cref="F71MirrorBlockRefinement"/> + <see cref="XGlobalChargeConjugationPairing"/> +
/// <see cref="PerBlockLiouvillianBuilder"/>; standard finite-group representation theory
/// (Klein four-group character table, isotypic decomposition).</para>
/// </summary>
public sealed class KleinFourGroupSelfPairedRefinement : Claim
{
    public int N { get; }
    public int M { get; }   // = N / 2, the self-paired popcount

    /// <summary>K-orbits on the sector basis. Each orbit is a sorted list of
    /// <c>(rowIndex, colIndex)</c> Hilbert-side index pairs (so the operator is
    /// <c>|row⟩⟨col|</c>; popcount(row) = popcount(col) = M). The first entry is the
    /// canonical orbit representative (smallest pair under lex order); subsequent entries
    /// are the other K-orbit members.</summary>
    public IReadOnlyList<KleinOrbit> Orbits { get; }

    /// <summary>Sub-block dimensions indexed by Klein character. Sum equals the sector
    /// dimension <c>C(N, N/2)²</c>.</summary>
    public IReadOnlyDictionary<KleinCharacter, int> SubBlockDims { get; }

    public int SectorDim => SubBlockDims.Values.Sum();

    public static KleinFourGroupSelfPairedRefinement Build(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        if ((N & 1) != 0) throw new ArgumentException(
            $"N must be even (self-paired sector requires N/2 ∈ ℤ); got N={N}.", nameof(N));

        int m = N / 2;
        var orbits = EnumerateOrbits(N, m);

        var subBlockDims = new Dictionary<KleinCharacter, int>
        {
            [KleinCharacter.PlusPlus] = 0,
            [KleinCharacter.PlusMinus] = 0,
            [KleinCharacter.MinusPlus] = 0,
            [KleinCharacter.MinusMinus] = 0,
        };
        foreach (var orbit in orbits)
            foreach (var chi in orbit.SurvivingCharacters)
                subBlockDims[chi]++;

        return new KleinFourGroupSelfPairedRefinement(N, m, orbits, subBlockDims);
    }

    /// <summary>Build the dense L sub-block for one Klein character. Element-wise
    /// construction: <c>L_χ[α, β] = (1/√(|O_α|·|O_β|)) · Σ_{g_α, g_β} χ(g_α)*·χ(g_β) ·
    /// ⟨e_{g_α·x_α} | L | e_{g_β·x_β}⟩</c>, where the inner matrix element follows from
    /// the chain-XY + Z-dephasing Liouvillian formula.</summary>
    public ComplexMatrix BuildSubBlockL(KleinCharacter chi, IReadOnlyList<double> gammaPerSite)
    {
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != N)
            throw new ArgumentException(
                $"gammaPerSite length {gammaPerSite.Count} != N {N}", nameof(gammaPerSite));

        var ordered = Orbits.Where(o => o.SurvivingCharacters.Contains(chi)).ToList();
        int dim = ordered.Count;
        var M_chi = Matrix<Complex>.Build.Dense(dim, dim);

        // Per (i, j) → set of (i', j') with nonzero L matrix element. Matching
        // PauliHamiltonian.XYChain(N, J) convention H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}):
        // each bond's (XX + YY) acts as |01⟩↔|10⟩ swap with amplitude 2, so per-bond H
        // matrix element is (J/2)·2 = J.
        //   - i' = i ⊕ bondFlip[b] for each bond b with bit_b(i) ≠ bit_{b+1}(i), j' = j
        //     contributes −i·J at L[(i', j), (i, j)].
        //   - i' = i, j' = j ⊕ bondFlip[b] similarly contributes +i·J.
        //   - i' = i, j' = j (diagonal): −2·Σ_l γ_l · (n_l^i XOR n_l^j).
        // (Hamiltonian uses J = 1 implicitly — caller scales if needed via gamma; future
        // extension: pass J as parameter. For now this primitive targets the canonical
        // J = 1 setup matching XyJordanWignerModes.)
        const double J = 1.0;
        int[] bondFlip = Enumerable.Range(0, N - 1).Select(b => (1 << b) | (1 << (b + 1))).ToArray();
        int dHilbert = 1 << N;

        for (int colOrbit = 0; colOrbit < dim; colOrbit++)
        {
            var oCol = ordered[colOrbit];
            double invSqrtSizeCol = 1.0 / Math.Sqrt(oCol.Size);

            for (int rowOrbit = 0; rowOrbit < dim; rowOrbit++)
            {
                var oRow = ordered[rowOrbit];
                double invSqrtSizeRow = 1.0 / Math.Sqrt(oRow.Size);

                Complex sum = Complex.Zero;
                for (int gColIdx = 0; gColIdx < oCol.Size; gColIdx++)
                {
                    double chiCol = ChiSign(chi, gColIdx, oCol.GroupElementIndices);
                    var (iCol, jCol) = oCol.Members[gColIdx];

                    for (int gRowIdx = 0; gRowIdx < oRow.Size; gRowIdx++)
                    {
                        double chiRow = ChiSign(chi, gRowIdx, oRow.GroupElementIndices);
                        var (iRow, jRow) = oRow.Members[gRowIdx];

                        // L[(iRow, jRow), (iCol, jCol)] for chain XY + Z-deph:
                        //   −i · H[iRow, iCol] · δ_{jRow, jCol}
                        // + i · H[jCol, jRow] · δ_{iRow, iCol}
                        // − 2 · Σ_l γ_l · (n_l^iCol XOR n_l^jCol) · δ_{iRow, iCol} · δ_{jRow, jCol}
                        Complex elem = Complex.Zero;
                        if (jRow == jCol)
                        {
                            // Hamiltonian on ket (row) side.
                            for (int b = 0; b < N - 1; b++)
                            {
                                int bitB = (iCol >> b) & 1;
                                int bitBp1 = (iCol >> (b + 1)) & 1;
                                if (bitB == bitBp1) continue;
                                int iColSwapped = iCol ^ bondFlip[b];
                                if (iColSwapped == iRow)
                                {
                                    elem += new Complex(0, -J);
                                }
                            }
                        }
                        if (iRow == iCol)
                        {
                            // Hamiltonian on bra (col) side: contributes +i·H[jCol, jRow] (note transposition for the bra commutator).
                            for (int b = 0; b < N - 1; b++)
                            {
                                int bitB = (jRow >> b) & 1;
                                int bitBp1 = (jRow >> (b + 1)) & 1;
                                if (bitB == bitBp1) continue;
                                int jRowSwapped = jRow ^ bondFlip[b];
                                if (jRowSwapped == jCol)
                                {
                                    elem += new Complex(0, +J);
                                }
                            }
                            if (jRow == jCol)
                            {
                                // Diagonal dissipator contribution.
                                int xorIJ = iCol ^ jCol;
                                double diss = 0.0;
                                for (int l = 0; l < N; l++)
                                    if (((xorIJ >> l) & 1) != 0)
                                        diss -= 2 * gammaPerSite[l];
                                elem += new Complex(diss, 0);
                            }
                        }

                        sum += chiRow * chiCol * elem;
                    }
                }
                M_chi[rowOrbit, colOrbit] = sum * invSqrtSizeRow * invSqrtSizeCol;
            }
        }

        return M_chi;
    }

    /// <summary>Sign χ(g) where g is the gIndex-th element of the orbit's group-element
    /// enumeration. The enumeration is ordered consistently for size-2 and size-4 orbits;
    /// see <see cref="KleinOrbit"/>.</summary>
    private static double ChiSign(KleinCharacter chi, int gIndex, IReadOnlyList<KleinGroupElement> groupElements)
    {
        var g = groupElements[gIndex];
        bool flipF71 = (g & KleinGroupElement.F71) == KleinGroupElement.F71;
        bool flipXN = (g & KleinGroupElement.XN) == KleinGroupElement.XN;
        int sign = 1;
        if (flipF71 && (chi == KleinCharacter.MinusPlus || chi == KleinCharacter.MinusMinus)) sign = -sign;
        if (flipXN && (chi == KleinCharacter.PlusMinus || chi == KleinCharacter.MinusMinus)) sign = -sign;
        return sign;
    }

    /// <summary>Enumerate K-orbits on the (m, m) joint-popcount sector. Each orbit is
    /// labelled by its canonical (smallest lex) member; orbit members are listed in
    /// canonical group-element order so character signs apply consistently.</summary>
    private static List<KleinOrbit> EnumerateOrbits(int N, int m)
    {
        var sectorIndices = EnumeratePopcount(N, m).ToArray();
        var seen = new HashSet<(int, int)>();
        var orbits = new List<KleinOrbit>();

        foreach (var i in sectorIndices)
        {
            foreach (var j in sectorIndices)
            {
                if (seen.Contains((i, j))) continue;

                int iF71 = ReverseBits(i, N);
                int jF71 = ReverseBits(j, N);
                int iXN = ((1 << N) - 1) ^ i;
                int jXN = ((1 << N) - 1) ^ j;
                int iBoth = ReverseBits(iXN, N);
                int jBoth = ReverseBits(jXN, N);

                var members = new List<((int row, int col) Pair, KleinGroupElement G)>
                {
                    ((i, j), KleinGroupElement.Identity),
                    ((iF71, jF71), KleinGroupElement.F71),
                    ((iXN, jXN), KleinGroupElement.XN),
                    ((iBoth, jBoth), KleinGroupElement.F71 | KleinGroupElement.XN),
                };

                // Distinct orbit members in canonical order: identity first, then
                // distinct group elements encountered as new (row, col) pairs.
                var orbitPairs = new List<(int row, int col)>();
                var orbitGroupElements = new List<KleinGroupElement>();
                var membersSeen = new HashSet<(int, int)>();
                foreach (var (pair, g) in members)
                {
                    if (membersSeen.Add(pair))
                    {
                        orbitPairs.Add(pair);
                        orbitGroupElements.Add(g);
                    }
                }

                foreach (var pair in orbitPairs) seen.Add(pair);

                // Stabiliser characters: those χ with χ(g) = +1 for every g such that
                // g · (i, j) = (i, j). The non-identity stabiliser elements are exactly
                // those Klein elements that DON'T appear as orbit-member-generators.
                var stabilizerSet = new HashSet<KleinGroupElement> { KleinGroupElement.Identity };
                foreach (var (pair, g) in members)
                    if (pair == (i, j))
                        stabilizerSet.Add(g);

                var surviving = new HashSet<KleinCharacter>();
                foreach (var chi in new[] { KleinCharacter.PlusPlus, KleinCharacter.PlusMinus, KleinCharacter.MinusPlus, KleinCharacter.MinusMinus })
                {
                    bool ok = true;
                    foreach (var g in stabilizerSet)
                    {
                        bool flipF71 = (g & KleinGroupElement.F71) == KleinGroupElement.F71;
                        bool flipXN = (g & KleinGroupElement.XN) == KleinGroupElement.XN;
                        int sign = 1;
                        if (flipF71 && (chi == KleinCharacter.MinusPlus || chi == KleinCharacter.MinusMinus)) sign = -sign;
                        if (flipXN && (chi == KleinCharacter.PlusMinus || chi == KleinCharacter.MinusMinus)) sign = -sign;
                        if (sign != 1) { ok = false; break; }
                    }
                    if (ok) surviving.Add(chi);
                }

                orbits.Add(new KleinOrbit(orbitPairs, orbitGroupElements, surviving));
            }
        }

        return orbits;
    }

    private static IEnumerable<int> EnumeratePopcount(int N, int m)
    {
        int max = 1 << N;
        for (int i = 0; i < max; i++)
            if (System.Numerics.BitOperations.PopCount((uint)i) == m)
                yield return i;
    }

    private static int ReverseBits(int x, int N)
    {
        int y = 0;
        for (int b = 0; b < N; b++) y |= ((x >> b) & 1) << (N - 1 - b);
        return y;
    }

    private KleinFourGroupSelfPairedRefinement(int n, int m, IReadOnlyList<KleinOrbit> orbits,
        IReadOnlyDictionary<KleinCharacter, int> subBlockDims)
        : base($"Klein 4-group internal symmetry on self-paired (N/2, N/2) sector at N={n}: " +
               $"K = {{1, F71, X⊗N, F71·X⊗N}} splits {subBlockDims.Values.Sum()} basis elements into 4 sub-blocks " +
               $"of dims [++={subBlockDims[KleinCharacter.PlusPlus]}, +-={subBlockDims[KleinCharacter.PlusMinus]}, " +
               $"-+={subBlockDims[KleinCharacter.MinusPlus]}, --={subBlockDims[KleinCharacter.MinusMinus]}].",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs (spatial-mirror Z₂) + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs (charge-conjugation Z₂) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/PerBlockLiouvillianBuilder.cs (per-block L source); " +
               "finite abelian group character table for Z₂ × Z₂ isotypic decomposition.")
    {
        N = n; M = m;
        Orbits = orbits;
        SubBlockDims = subBlockDims;
    }

    public override string DisplayName =>
        $"Klein 4-group refinement on (N/2, N/2)={N}/2 sector at N={N}";

    public override string Summary =>
        $"sector dim {SectorDim} → 4 sub-blocks " +
        $"{{++:{SubBlockDims[KleinCharacter.PlusPlus]}, +-:{SubBlockDims[KleinCharacter.PlusMinus]}, " +
        $"-+:{SubBlockDims[KleinCharacter.MinusPlus]}, --:{SubBlockDims[KleinCharacter.MinusMinus]}}} " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("m (= N/2)", M);
            yield return InspectableNode.RealScalar("sector dim C(N, N/2)²", SectorDim);
            yield return InspectableNode.RealScalar("orbit count", Orbits.Count);
            foreach (var (chi, dim) in SubBlockDims)
                yield return InspectableNode.RealScalar($"sub-block dim ({chi})", dim);
        }
    }
}

/// <summary>One K-orbit on a (m, m) self-paired sector. <see cref="Members"/> lists the
/// distinct (row, col) Hilbert-index pairs in the orbit, in the canonical order matching
/// <see cref="GroupElementIndices"/> (so <c>Members[k] = GroupElementIndices[k] · Members[0]</c>).
/// Orbit size is <c>Members.Count</c>, one of {1, 2, 4}; <see cref="SurvivingCharacters"/>
/// lists the Klein characters χ with non-vanishing symmetry-adapted basis vector on
/// this orbit (those trivial on the orbit's stabiliser).</summary>
public sealed record KleinOrbit(
    IReadOnlyList<(int Row, int Col)> Members,
    IReadOnlyList<KleinGroupElement> GroupElementIndices,
    IReadOnlySet<KleinCharacter> SurvivingCharacters)
{
    public int Size => Members.Count;
}

/// <summary>The four elements of the Klein four-group K = Z₂ × Z₂ generated by F71
/// (spatial mirror) and X⊗N (global charge conjugation). Encoded as bit flags so
/// composition is XOR.</summary>
[Flags]
public enum KleinGroupElement
{
    Identity = 0,
    F71 = 1,
    XN = 2,
}

/// <summary>The four 1-dimensional characters of the Klein four-group on the (χ_F71, χ_X⊗N)
/// basis: ++ trivial, +- flips X⊗N sign, -+ flips F71 sign, -- flips both.</summary>
public enum KleinCharacter
{
    PlusPlus,
    PlusMinus,
    MinusPlus,
    MinusMinus,
}
