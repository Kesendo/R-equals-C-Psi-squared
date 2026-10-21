using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Mat = RCPsiSquared.Core.Numerics.GaussianRational[,];

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the chain part of F167 (<c>experiments/THE_SODIUM_CHAIN.md</c>) at two
/// nuclei: each spin 3/2 a qubit pair (|3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩, |10⟩, |11⟩),
/// per nucleus H_k = q_k·ZZ_k + ω_k·I_z,k (q_k = ω_Q,k·√6/2, the coefficient of ZZ, as in
/// <see cref="SodiumTwoMirrorsWitness"/>), jumps ZZ at 3J₀/2, ZX and ZY at 3J₁/2, XI and YI at
/// 3J₂/2, a real two-body coupling between the nuclei; Hilbert 16, Liouville 256, cell basis
/// 16i + j of |i⟩⟨j|, Q = ZZ₁ZZ₂, S(ρ) = Q·ρᵀ.
///
/// <para>Exact over ℚ(i), no double in any gate. The coupling is read from the SPIN MATRICES:
/// I_± held as A + √3·B (the page's verifier reads the strings by sympy), the two-nucleus secular
/// dipolar coupling 3I_zI_z − I·I decomposed into strings with a rational and a √3 part, and the
/// truncations taken as the projection onto the energy shells of q₁ZZ₁ + q₂ZZ₂. The √3 of the full
/// coupling is not rounded and not replaced: the diagonal similarity Φ(ρ) = DρD⁻¹,
/// D = diag(3^(−b/2)) with b the number of inner levels, multiplies the cell entry (r, c) of L by
/// 3^((β_r − β_c)/2), β = b(bra) − b(ket) on the cell |i⟩⟨j|, and the witness checks that every
/// rational entry sits on an even β difference and every √3 entry on an odd one, so ΦLΦ⁻¹ is a
/// rational matrix with the spectrum, the kernels and the invariant hulls (Φ keeps every cell) of
/// the page's generator. Where a rational stand-in c for √3 is used (the lemma's control and a
/// second kernel count), it is named so: the page proves both for any non-zero coefficient.</para>
///
/// <para>The parameters are the page verifier's, in its units: d = 4 (the coupling d·I_z⊗I_z or
/// d·(3I_zI_z − I·I)), q = (12, 20) and ω = (20, 28) for the Ising coupling, q = (12, 20) and
/// ω = 20 for the unequal truncation, q = 12 and ω = 20 for the equal one (its coupling at 2d, as
/// the verifier builds it, twice the secular part), q = 0 and ω = 20 for
/// the full coupling; identical rates (J₀, J₁, J₂) = (6, 4, 2) and, for the Ising coupling with
/// rates per nucleus, (J₀, J₁, J₂) and (10, 2, 4). User inputs: the identical rates and d, and c.</para>
///
/// <para>Gates, each with a control that can fail: (1) the dictionary I_z, I_x, I_y, T₂₀ against
/// the string forms, the 22 strings of the dipolar coupling, its eight √3 strings exactly those
/// anticommuting with Q, IXIX and IYIY the rational ones anticommuting with each ZZ_k, and the two
/// truncations as shell projections; (2) the lemma S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ_k(J₁,k + J₂,k)
/// entry by entry for the Ising coupling (rates per nucleus and identical) and both truncations,
/// the √3 strings at coefficient c breaking it; (3) the ladder: [L, P_k] = 0 for the Ising
/// coupling and the unequal truncation and not for the equal one, every bright-set block's χ
/// palindromic about −s_B and about no neighbouring rung; (4) the full coupling at J₀ = 0: F158's
/// far span (the strings anticommuting with every jump are Q alone, none once the J₀ jumps are
/// on), its two ends as exact kernel dimensions at √3 and at c, and the characteristic polynomial
/// at √3 not palindromic about −σ, the equal truncation pairing by both routes as the control;
/// (5) the explicit sharp spaces: E₁ (twelve operators, σσσ on three qubits) and E₂ (the
/// signed satellite flip-flop) built from strings, each invariant under every parameter piece of L
/// separately and annihilated by the √3 coupling from both sides, so E_m ⊆ U_m ⊆ K_m at every
/// parameter with identical rates, one offset and an equal splitting; the controls: rates per
/// nucleus, a splitting or an offset on one nucleus break both, the exchange-symmetric order ±4
/// partners leave the kernel, a bare column breaks the rational coupling. Readings, not gates:
/// the exact dimensions of K_m (equal truncation) and of K_m and U_m (full coupling at √3) by
/// Krylov stabilisation, and whether each is palindromic about −s_m.</para>
///
/// <para>Typed home <see cref="SodiumChainMirrorLadderClaim"/>.</para></summary>
public sealed class SodiumChainWitness : IInspectable
{
    /// <summary>The largest |input| accepted; every comparison is exact at any size, this only
    /// keeps the inspect-time arithmetic small.</summary>
    public const long MaxMagnitude = 1000;

    private const int Nuclei = 2, Levels = 4, HDim = 16, Dim = 256;

    /// <summary>The second nucleus's rates for the Ising coupling with rates per nucleus.</summary>
    public static readonly (long J0, long J1, long J2) SecondNucleusRates = (10, 2, 4);

    /// <summary>(q₁, q₂, ω₁, ω₂) for the Ising coupling: two isotopes or two sites, every frequency apart.</summary>
    public static readonly (long Q1, long Q2, long Omega1, long Omega2) IsingFrequencies = (12, 20, 20, 28);

    /// <summary>(q₁, q₂, ω) for the unequal truncation, |q₁| ≠ |q₂|, one offset.</summary>
    public static readonly (long Q1, long Q2, long Omega) UnequalFrequencies = (12, 20, 20);

    /// <summary>(q, ω) for the equal truncation, identical nuclei.</summary>
    public static readonly (long Q, long Omega) EqualFrequencies = (12, 20);

    /// <summary>The offset of the full coupling, computed at zero splitting.</summary>
    public const long FullOmega = 20;

    /// <summary>The page's identity of what K₀ holds beyond U₀ at its verifier's point (J₀, J₁, J₂) =
    /// (6, 4, 2), d = 4: χ(K₀) = χ(U₀)·x(x + 48)²((x + 36)² + 1600).</summary>
    public static readonly (long J0, long J1, long J2, long D) PagePoint = (6, 4, 2, 4);

    /// <summary>One bright-set block of the ladder: its size, its rung 2s_B, whether χ pairs about
    /// −s_B, and the neighbouring rungs 2s_{B△{k}} with whether χ pairs about each.</summary>
    public sealed record LadderBlock(int BrightMask, int Size, long TwoS, bool Pairs,
        IReadOnlyList<(long TwoS, bool Pairs)> Neighbours)
    {
        public bool Holds => Pairs && Neighbours.All(n => !n.Pairs);
        public string Label => BrightMask switch { 0 => "∅", 1 => "{1}", 2 => "{2}", _ => "{1, 2}" };
    }

    /// <summary>One invariant hull reading: m bright nuclei, the coordinate size, the exact
    /// dimension, the rung 2s_m, and whether χ of L on it pairs about −s_m (null when empty).</summary>
    public sealed record HullReading(int M, int CellCount, int Dimension, long TwoS, bool? Pairs);

    public long J0 { get; }
    public long J1 { get; }
    public long J2 { get; }
    public long D { get; }
    public long C { get; }

    // gate 1: the dictionary
    public bool IzMatches { get; }
    public bool IxMatches { get; }
    public bool IyMatches { get; }
    public bool T20Matches { get; }
    public bool SqrtThreePartNeeded { get; }
    public int DipolarStringCount { get; }
    public IReadOnlyList<string> SqrtThreeStrings { get; }
    public bool PartsDisjoint { get; }
    public bool SqrtThreeAreTheQAnticommuting { get; }
    public IReadOnlyList<string> SatelliteSatelliteStrings { get; }
    public bool SatelliteSatelliteAnticommuteWithEach { get; }
    public bool CouplingsReal { get; }
    public bool UnequalTruncationIsShellProjection { get; }
    public bool EqualTruncationIsShellProjection { get; }
    public bool TruncationsDiffer { get; }
    public bool DictionaryHolds => IzMatches && IxMatches && IyMatches && T20Matches && SqrtThreePartNeeded
        && DipolarStringCount == 22 && SqrtThreeStrings.Count == 8 && PartsDisjoint && SqrtThreeAreTheQAnticommuting
        && SatelliteSatelliteStrings.SequenceEqual(new[] { "IXIX", "IYIY" }) && SatelliteSatelliteAnticommuteWithEach
        && CouplingsReal && UnequalTruncationIsShellProjection && EqualTruncationIsShellProjection && TruncationsDiffer;

    // gate 2: the lemma
    public bool LemmaIsingPerNucleus { get; }
    public bool LemmaIsingIdentical { get; }
    public bool LemmaUnequal { get; }
    public bool LemmaEqual { get; }
    public bool SqrtThreeStringsBreakLemma { get; }
    public bool LemmaHolds => LemmaIsingPerNucleus && LemmaIsingIdentical && LemmaUnequal && LemmaEqual;

    // gate 3: the ladder
    public bool IsingKeepsBrightness { get; }
    public bool UnequalKeepsBrightness { get; }
    public bool EqualKeepsBrightness { get; }
    public IReadOnlyList<LadderBlock> IsingBlocks { get; }
    public IReadOnlyList<LadderBlock> UnequalBlocks { get; }
    public bool LadderHolds => IsingKeepsBrightness && UnequalKeepsBrightness && !EqualKeepsBrightness
        && IsingBlocks.All(b => b.Holds) && UnequalBlocks.All(b => b.Holds);

    // gate 4: the full coupling at J₀ = 0
    public IReadOnlyList<string> FarSpanAtJ0Zero { get; }
    public int FarSpanWithJ0Jumps { get; }
    public bool RestCommutesWithQ { get; }
    public bool SqrtThreeAnticommutesWithQ { get; }
    public bool SimilarityIsRational { get; }
    public long TwoSigma { get; }
    public int NearKernel { get; }
    public int FarKernel { get; }
    public int NearKernelAtC { get; }
    public int FarKernelAtC { get; }
    public bool CharpolyPairsAboutSigma { get; }
    public int ControlNearKernel { get; }
    public int ControlFarKernel { get; }
    public bool ControlCharpolyPairs { get; }
    public bool GlobalBreakHolds => FarSpanAtJ0Zero.SequenceEqual(new[] { "ZZZZ" }) && FarSpanWithJ0Jumps == 0
        && RestCommutesWithQ && SqrtThreeAnticommutesWithQ && SimilarityIsRational
        // the near ends always hold the identity; the far ends' zero carries the gate
        && NearKernel >= 1 && FarKernel == 0 && NearKernelAtC >= 1 && FarKernelAtC == 0 && !CharpolyPairsAboutSigma
        && ControlNearKernel == ControlFarKernel && ControlCharpolyPairs;

    // readings
    public IReadOnlyList<HullReading> EqualSharp { get; }
    public IReadOnlyList<HullReading> FullSharp { get; }
    public IReadOnlyList<HullReading> FullKernelHull { get; }
    public int SqrtThreeRank { get; }
    public bool HullsInvariant { get; }

    /// <summary>χ(K₀)/χ(U₀) = x(x + 48)²((x + 36)² + 1600) exactly, read at <see cref="PagePoint"/>
    /// only (null elsewhere, where the page states nothing).</summary>
    public bool? K0BeyondU0MatchesPage { get; }

    // gate 5: the explicit sharp spaces E₁ ⊆ U₁ ⊆ K₁ and E₂ ⊆ U₂ ⊆ K₂ at every parameter with identical rates
    public int E1Dimension { get; }
    public int E2Dimension { get; }
    public bool E1InOneBrightCells { get; }
    public bool E2InTwoBrightCells { get; }
    public bool E1InvariantUnderEveryPiece { get; }
    public bool E2InvariantUnderEveryPiece { get; }
    public bool E1InSqrtThreeKernel { get; }
    public bool E2InSqrtThreeKernel { get; }
    public bool E2OnItsRung { get; }
    public bool E1BreaksWithRatesPerNucleus { get; }
    public bool E2BreaksWithRatesPerNucleus { get; }
    public bool E1TopBlockKeepsRatesPerNucleus { get; }
    public bool SymmetricPartnersKeepRationalCoupling { get; }
    public bool SymmetricPartnersLeaveSqrtThreeKernel { get; }
    public bool BareColumnBreaksRationalCoupling { get; }
    public bool UnequalSplittingBreaksBoth { get; }
    public bool UnequalOffsetBreaksBoth { get; }
    /// <summary>At these parameters dim K_m = dim U_m = dim E_m for m = 1, 2, so E_m = U_m = K_m here.</summary>
    public bool ExplicitAreTheSharpSpacesHere { get; }
    public bool ExplicitHolds => E1Dimension == 12 && E2Dimension == 1 && E1InOneBrightCells && E2InTwoBrightCells
        && E1InvariantUnderEveryPiece && E2InvariantUnderEveryPiece && E1InSqrtThreeKernel && E2InSqrtThreeKernel && E2OnItsRung
        && E1BreaksWithRatesPerNucleus && E2BreaksWithRatesPerNucleus && E1TopBlockKeepsRatesPerNucleus
        && SymmetricPartnersKeepRationalCoupling && SymmetricPartnersLeaveSqrtThreeKernel && BareColumnBreaksRationalCoupling
        && UnequalSplittingBreaksBoth && UnequalOffsetBreaksBoth;

    /// <summary>I_z,tot is an eigenoperator of L (identical rates) exactly when J₁ = J₂; it lies in the
    /// dark cells, so then in K₀, and outside the kernel of the √3 coupling.</summary>
    public bool IzTotalIsEigenoperator { get; }
    public bool IzTotalInSqrtThreeKernel { get; }

    private readonly Mat _rat;     // the dipolar coupling 3I_zI_z − I·I, rational part
    private readonly Mat _root;    // its √3 part, the coefficient of √3

    public SodiumChainWitness(long j0 = 6, long j1 = 4, long j2 = 2, long d = 4, long c = 1)
    {
        foreach (var (v, name) in new[] { (j0, nameof(j0)), (j1, nameof(j1)), (j2, nameof(j2)), (d, nameof(d)), (c, nameof(c)) })
            if (v > MaxMagnitude || v < -MaxMagnitude)
                throw new ArgumentOutOfRangeException(name, v, $"|input| ≤ {MaxMagnitude}");
        if (j0 < 1) throw new ArgumentOutOfRangeException(nameof(j0), j0, "J₀ ≥ 1 (the ladder's rungs are 3J₀ apart; J₀ = 0 is gate 4's own point)");
        if (j1 < 1 || j2 < 1) throw new ArgumentOutOfRangeException(nameof(j1), (j1, j2), "J₁, J₂ ≥ 1 (F158's break needs them on)");
        if (d == 0) throw new ArgumentOutOfRangeException(nameof(d), d, "a coupling d ≠ 0");
        if (c == 0) throw new ArgumentOutOfRangeException(nameof(c), c, "a stand-in c ≠ 0 for √3");
        (J0, J1, J2, D, C) = (j0, j1, j2, d, c);

        // ---- gate 1: spin 3/2 from the spin formula, entries in ℚ(√3)
        var x = new BigRational(15, 4);
        var iz = Zero(Levels);
        var ipA = Zero(Levels);
        var ipB = Zero(Levels);
        for (int k = 0; k < Levels; k++)
        {
            var m = new BigRational(3 - 2 * k, 2);
            iz[k, k] = m;
            if (k == 0) continue;
            var (a, b) = SqrtInQ3(x - m * (m + 1));       // I_+|m_k⟩ = √(X − m(m + 1))·|m_{k−1}⟩
            ipA[k - 1, k] = a;
            ipB[k - 1, k] = b;
        }
        var ip = new Root3(ipA, ipB);
        var im = new Root3(Transpose(ipA), Transpose(ipB));
        var half = new BigRational(1, 2);
        var ix = (ip + im).Scale(half);
        var iy = (ip + im.Scale(-1)).Scale(-GaussianRational.I * half);  // (I_+ − I_−)/(2i)
        IzMatches = Equal(iz, Add(Str("ZI"), Scale(Str("IZ"), half)));
        IxMatches = Equal(ix.A, Scale(Add(Str("XX"), Str("YY")), half)) && Equal(ix.B, Scale(Str("IX"), half));
        IyMatches = Equal(iy.A, Scale(Sub(Str("YX"), Str("XY")), half)) && Equal(iy.B, Scale(Str("IY"), half));
        T20Matches = Equal(Sub(Scale(Mul(iz, iz), 3), Scale(Identity(Levels), x)), Scale(Str("ZZ"), 3)); // T₂₀ = D/√6, D = 3·ZZ
        SqrtThreePartNeeded = !IsZero(ix.B) && !IsZero(iy.B);

        var izR = new Root3(iz, Zero(Levels));
        var dip = Root3.Kron(izR, izR).Scale(2) + Root3.Kron(ix, ix).Scale(-1) + Root3.Kron(iy, iy).Scale(-1); // 3IzIz − I·I
        _rat = dip.A;
        _root = dip.B;
        var ratStrings = Decompose(_rat);
        var rootStrings = Decompose(_root);
        PartsDisjoint = !ratStrings.Keys.Intersect(rootStrings.Keys).Any()
                        && Equal(Recompose(ratStrings), _rat) && Equal(Recompose(rootStrings), _root);
        var all = ratStrings.Keys.Union(rootStrings.Keys).ToHashSet();
        DipolarStringCount = all.Count;
        SqrtThreeStrings = rootStrings.Keys.OrderBy(s => s, StringComparer.Ordinal).ToArray();
        SqrtThreeAreTheQAnticommuting = all.Where(w => Anticommutes(w, "ZZZZ")).ToHashSet().SetEquals(rootStrings.Keys)
                                        && !all.Where(w => Anticommutes(w, "ZZII")).ToHashSet().SetEquals(rootStrings.Keys);
        SatelliteSatelliteStrings = ratStrings.Keys.Where(w => Anticommutes(w, "ZZII") || Anticommutes(w, "IIZZ"))
            .OrderBy(s => s, StringComparer.Ordinal).ToArray();
        SatelliteSatelliteAnticommuteWithEach = SatelliteSatelliteStrings.All(w => Anticommutes(w, "ZZII") && Anticommutes(w, "IIZZ"));
        var each = Recompose(ratStrings.Where(kv => !SatelliteSatelliteStrings.Contains(kv.Key)).ToDictionary(kv => kv.Key, kv => kv.Value));
        CouplingsReal = IsRealSymmetric(_rat) && IsRealSymmetric(_root);

        var (uq1, uq2, uOmega) = UnequalFrequencies;
        var (eq, eOmega) = EqualFrequencies;
        UnequalTruncationIsShellProjection = Equal(Secular(_rat, uq1, uq2), each) && IsZero(Secular(_root, uq1, uq2));
        var kappa = ratStrings["IXIX"];
        var flip = Scale(Add(Str("IXIX"), Str("IYIY")), kappa);
        var equalForm = Add(each, Mul(flip, Scale(Sub(Identity(HDim), Str("ZZZZ")), half)));   // (IXIX + IYIY)·(1 − ZZ₁ZZ₂)/2
        EqualTruncationIsShellProjection = ratStrings["IYIY"] == kappa && Equal(Secular(_rat, eq, eq), equalForm)
                                           && IsZero(Secular(_root, eq, eq));
        TruncationsDiffer = !Equal(each, equalForm);

        // ---- the couplings and Hamiltonians
        var dd = (GaussianRational)(BigRational)d;
        var ising = Scale(Kron(iz, iz), dd);
        var unequal = Scale(each, dd);
        var equal = Scale(equalForm, 2 * dd);       // the verifier's equal2: twice the secular part
        var fullRat = Scale(_rat, dd);
        var fullRoot = Scale(_root, dd);           // H₃ = √3·fullRoot

        var (iq1, iq2, iw1, iw2) = IsingFrequencies;
        var hIsing = Add(Local(iz, iq1, iq2, iw1, iw2), ising);
        var hUnequal = Add(Local(iz, uq1, uq2, uOmega, uOmega), unequal);
        var hEqual = Add(Local(iz, eq, eq, eOmega, eOmega), equal);
        var hFullRest = Add(Local(iz, 0, 0, FullOmega, FullOmega), fullRat);
        var hFullAtC = Add(hFullRest, Scale(fullRoot, (GaussianRational)(BigRational)c));

        var identical = new[] { (j0, j1, j2), (j0, j1, j2) };
        var perNucleus = new[] { (j0, j1, j2), SecondNucleusRates };

        // ---- gate 2: the lemma
        LemmaIsingPerNucleus = Lemma(hIsing, perNucleus);
        LemmaIsingIdentical = Lemma(hIsing, identical);
        LemmaUnequal = Lemma(hUnequal, identical);
        LemmaEqual = Lemma(hEqual, identical);
        SqrtThreeStringsBreakLemma = !Lemma(hFullAtC, identical);

        // ---- gate 3: the ladder
        var lIsing = Generator(hIsing, perNucleus);
        var lUnequal = Generator(hUnequal, identical);
        var lEqual = Generator(hEqual, identical);
        IsingKeepsBrightness = KeepsBrightness(lIsing);
        UnequalKeepsBrightness = KeepsBrightness(lUnequal);
        EqualKeepsBrightness = KeepsBrightness(lEqual);
        IsingBlocks = IsingKeepsBrightness ? Ladder(lIsing, perNucleus) : Array.Empty<LadderBlock>();
        UnequalBlocks = UnequalKeepsBrightness ? Ladder(lUnequal, identical) : Array.Empty<LadderBlock>();

        // ---- gate 4: the full coupling at J₀ = 0
        var zeroJ0 = new[] { (0L, j1, j2), (0L, j1, j2) };
        var jumps = new List<string>();
        for (int k = 0; k < Nuclei; k++)
            foreach (var j in new[] { "ZX", "ZY", "XI", "YI" }) jumps.Add(Place(k, j));
        var allStrings = AllStrings(2 * Nuclei);
        FarSpanAtJ0Zero = allStrings.Where(w => jumps.All(j => Anticommutes(w, j))).ToArray();
        var withJ0 = jumps.Concat(Enumerable.Range(0, Nuclei).Select(k => Place(k, "ZZ"))).ToArray();
        FarSpanWithJ0Jumps = allStrings.Count(w => withJ0.All(j => Anticommutes(w, j)));
        var q = Str("ZZZZ");
        RestCommutesWithQ = Equal(Mul(hFullRest, q), Mul(q, hFullRest));
        SqrtThreeAnticommutesWithQ = Equal(Mul(fullRoot, q), Scale(Mul(q, fullRoot), -1)) && !IsZero(fullRoot);

        TwoSigma = SodiumChainMirrorLadderClaim.LemmaConstant(zeroJ0);   // 2σ = 6Σ(J₁ + J₂) at J₀ = 0
        var (lFull0, rational0) = Rationalised(hFullRest, fullRoot, zeroJ0);
        SimilarityIsRational = rational0;
        (NearKernel, FarKernel) = TwoEnds(lFull0, TwoSigma);
        (NearKernelAtC, FarKernelAtC) = TwoEnds(Generator(hFullAtC, zeroJ0), TwoSigma);
        CharpolyPairsAboutSigma = SpectrumPairs(lFull0, Enumerable.Range(0, Dim).ToArray(), TwoSigma);
        var lEqual0 = Generator(hEqual, zeroJ0);
        (ControlNearKernel, ControlFarKernel) = TwoEnds(lEqual0, TwoSigma);
        ControlCharpolyPairs = SpectrumPairs(lEqual0, Enumerable.Range(0, Dim).ToArray(), TwoSigma);

        // ---- readings: the invariant hulls
        bool invariant = true;
        EqualSharp = Hulls(lEqual, identical, Array.Empty<Mat>(), ref invariant, out _);
        var (lFull, rational) = Rationalised(hFullRest, fullRoot, identical);
        SimilarityIsRational &= rational;
        FullSharp = Hulls(lFull, identical, Array.Empty<Mat>(), ref invariant, out var restrictionsK);
        var h3 = SimilarH3(fullRoot);
        SqrtThreeRank = Rank(fullRoot);
        FullKernelHull = Hulls(lFull, identical, new[] { LeftMultiplication(h3), RightMultiplication(h3) }, ref invariant, out var restrictionsU);
        HullsInvariant = invariant;
        if ((j0, j1, j2, d) == PagePoint)
        {
            var scale = CommonScale(restrictionsK[0].Concat(restrictionsU[0]));
            K0BeyondU0MatchesPage = PageQuotientMatches(Chi(restrictionsK[0], scale), Chi(restrictionsU[0], scale), scale);
        }

        // ---- gate 5: the explicit sharp spaces, for every parameter
        var e1 = ExplicitE1();
        var e2 = Columns(new[] { Strings(("IXZX", 1), ("IYZY", 1), ("ZXIX", 1), ("ZYIY", 1)) });
        var e1Top = TakeColumns(e1, 8);
        E1Dimension = Rank(e1);
        E2Dimension = Rank(e2);
        E1InOneBrightCells = InCells(e1, 1);
        E2InTwoBrightCells = InCells(e2, 2);
        var izTot = Add(Kron(iz, Identity(Levels)), Kron(Identity(Levels), iz));
        var zeroH = Zero(HDim);
        var ratePieces = new[] { (1L, 0L, 0L), (0L, 1L, 0L), (0L, 0L, 1L) }
            .Select(r => Generator(zeroH, new[] { r, r })).ToArray();
        // one offset, an equal splitting on both nuclei (zero included), the coupling's two parts
        var hamiltonianPieces = new[] { HamiltonianPart(izTot), HamiltonianPart(_rat), HamiltonianPart(_root),
                                        HamiltonianPart(Add(Str("ZZII"), Str("IIZZ"))) };
        var pieces = ratePieces.Concat(hamiltonianPieces).ToArray();
        E1InvariantUnderEveryPiece = pieces.All(p => Invariant(e1, p));
        E2InvariantUnderEveryPiece = pieces.All(p => Invariant(e2, p));
        var left3 = LeftMultiplication(_root);
        var right3 = RightMultiplication(_root);
        E1InSqrtThreeKernel = IsZero(Mul(left3, e1)) && IsZero(Mul(right3, e1));
        E2InSqrtThreeKernel = IsZero(Mul(left3, e2)) && IsZero(Mul(right3, e2));
        var minusSix = (GaussianRational)(BigRational)(-6);
        E2OnItsRung = ratePieces.All(p => Equal(Mul(p, e2), Scale(e2, minusSix)))
                      && hamiltonianPieces.All(p => IsZero(Mul(p, e2)))
                      && SodiumChainMirrorLadderClaim.TwiceRung(Nuclei, 2, j0, j1, j2) == 12 * (j0 + j1 + j2);   // the cross-dock to the claim's rung
        var perNucleusPieces = Enumerable.Range(0, Nuclei)
            .SelectMany(k => new[] { (1L, 0L, 0L), (0L, 1L, 0L), (0L, 0L, 1L) }
                .Select(r => Generator(zeroH, k == 0 ? new[] { r, (0L, 0L, 0L) } : new[] { (0L, 0L, 0L), r })))
            .ToArray();
        E1BreaksWithRatesPerNucleus = !perNucleusPieces.All(p => Invariant(e1, p));
        E2BreaksWithRatesPerNucleus = !perNucleusPieces.All(p => Invariant(e2, p));
        E1TopBlockKeepsRatesPerNucleus = perNucleusPieces.Concat(hamiltonianPieces).All(p => Invariant(e1Top, p));
        // controls on the coupling legs: the exchange-symmetric order ±4 partners keep the rational coupling and leave V;
        // one bare TQ ⊗ satellite column breaks the rational coupling
        var symmetric = Columns("IZ".SelectMany(letter => new[] { false, true }
            .Select(imag => Add(SigmaTriple(new[] { 1, 2, 3 }, 0, letter, imag), SigmaTriple(new[] { 0, 1, 3 }, 2, letter, imag)))).ToArray());
        SymmetricPartnersKeepRationalCoupling = Invariant(symmetric, hamiltonianPieces[1]);
        SymmetricPartnersLeaveSqrtThreeKernel = !(IsZero(Mul(left3, symmetric)) && IsZero(Mul(right3, symmetric)));
        BareColumnBreaksRationalCoupling = !Invariant(Columns(new[] { SigmaTriple(new[] { 1, 2, 3 }, 0, 'I', false) }), hamiltonianPieces[1]);
        // the scope: a splitting or an offset on one nucleus alone breaks both spaces
        var oneSplitting = HamiltonianPart(Str("ZZII"));
        var oneOffset = HamiltonianPart(Kron(iz, Identity(Levels)));
        UnequalSplittingBreaksBoth = !Invariant(e1, oneSplitting) && !Invariant(e2, oneSplitting);
        UnequalOffsetBreaksBoth = !Invariant(e1, oneOffset) && !Invariant(e2, oneOffset);
        // E_m ⊆ U_m needs the cells, the invariance and the kernel; equal dimension then makes E_m = U_m = K_m here
        ExplicitAreTheSharpSpacesHere = E1InOneBrightCells && E2InTwoBrightCells && E1InvariantUnderEveryPiece && E2InvariantUnderEveryPiece
                                        && E1InSqrtThreeKernel && E2InSqrtThreeKernel && HullsInvariant
                                        && FullSharp[1].Dimension == E1Dimension && FullKernelHull[1].Dimension == E1Dimension
                                        && FullSharp[2].Dimension == E2Dimension && FullKernelHull[2].Dimension == E2Dimension;

        // the total magnetisation: an eigenoperator of L exactly when J₁ = J₂ (the longitudinal rates 6J₁, 6J₂ meet)
        var izVec = Columns(new[] { izTot });
        var lRates = Generator(zeroH, identical);
        IzTotalIsEigenoperator = Rank(HStack(izVec, Mul(lRates, izVec))) == 1 && hamiltonianPieces.All(p => IsZero(Mul(p, izVec)));
        IzTotalInSqrtThreeKernel = IsZero(Mul(left3, izVec)) && IsZero(Mul(right3, izVec));
    }

    // ===== the explicit sharp spaces =====

    /// <summary>The twelve operators of E₁, written in σ = X + iY on three qubits and I or Z on the
    /// fourth: Re and Im of σσσ on the qubits (1, 2, 3) and on (1, 3, 4), the remaining qubit I or Z
    /// (the whole block of coherence order ±5: TQ on one nucleus, DQ on the other); then Re and Im
    /// of σσσ on (2, 3, 4) with I or Z on qubit 1 minus the same on (1, 2, 4) with that letter on
    /// qubit 3 (TQ ⊗ satellite − satellite ⊗ TQ, order ±4).</summary>
    private static Mat ExplicitE1()
    {
        var cols = new List<Mat>();
        foreach (var (qs, free) in new[] { (new[] { 0, 1, 2 }, 3), (new[] { 0, 2, 3 }, 1) })
            foreach (var letter in "IZ")
                foreach (var imag in new[] { false, true })
                    cols.Add(SigmaTriple(qs, free, letter, imag));
        foreach (var letter in "IZ")
            foreach (var imag in new[] { false, true })
                cols.Add(Sub(SigmaTriple(new[] { 1, 2, 3 }, 0, letter, imag), SigmaTriple(new[] { 0, 1, 3 }, 2, letter, imag)));
        return Columns(cols);
    }

    /// <summary>Re (or Im) of σ ⊗ σ ⊗ σ on the qubits <paramref name="qs"/>, σ = X + iY, the letter on
    /// the qubit <paramref name="free"/>: the strings with n_Y letters Y carry the real (imaginary)
    /// part of iⁿ.</summary>
    private static Mat SigmaTriple(int[] qs, int free, char letter, bool imag)
    {
        var m = Zero(HDim);
        for (int bits = 0; bits < 1 << qs.Length; bits++)
        {
            int ny = BitOperations.PopCount((uint)bits);
            int coef = imag ? (ny % 4 == 1 ? 1 : ny % 4 == 3 ? -1 : 0) : (ny % 4 == 0 ? 1 : ny % 4 == 2 ? -1 : 0);
            if (coef == 0) continue;
            var word = new char[4];
            word[free] = letter;
            for (int t = 0; t < qs.Length; t++) word[qs[t]] = ((bits >> t) & 1) == 1 ? 'Y' : 'X';
            m = Add(m, Scale(Str(new string(word)), (BigRational)coef));
        }
        return m;
    }

    private static Mat Strings(params (string Word, long Coefficient)[] terms)
    {
        var m = Zero(HDim);
        foreach (var (w, c) in terms) m = Add(m, Scale(Str(w), (BigRational)c));
        return m;
    }

    /// <summary>The operators as columns of a 256 × k matrix, cell 16i + j holding ρ[i, j].</summary>
    private static Mat Columns(IReadOnlyList<Mat> ops)
    {
        var r = Zero(Dim, ops.Count);
        for (int k = 0; k < ops.Count; k++)
            for (int i = 0; i < HDim; i++)
                for (int j = 0; j < HDim; j++) r[HDim * i + j, k] = ops[k][i, j];
        return r;
    }

    private static Mat TakeColumns(Mat a, int count)
    {
        var r = Zero(a.GetLength(0), count);
        for (int i = 0; i < a.GetLength(0); i++)
            for (int k = 0; k < count; k++) r[i, k] = a[i, k];
        return r;
    }

    private static Mat HStack(Mat a, Mat b)
    {
        int n = a.GetLength(0), p = a.GetLength(1), q = b.GetLength(1);
        var r = Zero(n, p + q);
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < p; k++) r[i, k] = a[i, k];
            for (int k = 0; k < q; k++) r[i, p + k] = b[i, k];
        }
        return r;
    }

    /// <summary>M maps the column span of E into itself: rank [E | M·E] = rank E, exactly.</summary>
    private static bool Invariant(Mat e, Mat m) => Rank(HStack(e, Mul(m, e))) == Rank(e);

    private static bool InCells(Mat e, int brightCount)
    {
        for (int r = 0; r < Dim; r++)
            for (int k = 0; k < e.GetLength(1); k++)
                if (!e[r, k].IsZero && BrightCount(r) != brightCount) return false;
        return true;
    }

    // ===== the generator =====

    private static string Place(int nucleus, string pair) => nucleus == 0 ? pair + "II" : "II" + pair;

    /// <summary>Σ_k (q_k·ZZ_k + ω_k·I_z,k) on the two nuclei.</summary>
    private static Mat Local(Mat iz, long q1, long q2, long w1, long w2) =>
        Add(Add(Scale(Str("ZZII"), (BigRational)q1), Scale(Str("IIZZ"), (BigRational)q2)),
            Add(Scale(Kron(iz, Identity(Levels)), (BigRational)w1), Scale(Kron(Identity(Levels), iz), (BigRational)w2)));

    /// <summary>L = −i[H, ·] + Σ_k Σ_P γ_P (P·ρ·P − ρ), the jumps of nucleus k ZZ at 3J₀/2, ZX and ZY
    /// at 3J₁/2, XI and YI at 3J₂/2; row-major cells 16i + j.</summary>
    private static Mat Generator(Mat h, IReadOnlyList<(long J0, long J1, long J2)> rates, bool flipJ0 = false)
    {
        var l = HamiltonianPart(h);
        for (int k = 0; k < Nuclei; k++)
        {
            var (a0, a1, a2) = rates[k];
            if (flipJ0) a0 = -a0;
            foreach (var (s, rate) in new[] { ("ZZ", a0), ("ZX", a1), ("ZY", a1), ("XI", a2), ("YI", a2) })
            {
                if (rate == 0) continue;
                var g = (GaussianRational)new BigRational(3 * rate, 2);
                var p = Str(Place(k, s));
                for (int i = 0; i < HDim; i++)
                    for (int a = 0; a < HDim; a++)
                    {
                        if (p[i, a].IsZero) continue;
                        for (int b = 0; b < HDim; b++)
                            for (int j = 0; j < HDim; j++)
                                if (!p[b, j].IsZero) l[HDim * i + j, HDim * a + b] += g * p[i, a] * p[b, j];
                    }
                for (int r = 0; r < Dim; r++) l[r, r] -= g;
            }
        }
        return l;
    }

    private static Mat HamiltonianPart(Mat h)
    {
        var l = Zero(Dim);
        var minusI = -GaussianRational.I;
        for (int a = 0; a < HDim; a++)
            for (int b = 0; b < HDim; b++)
            {
                if (h[a, b].IsZero) continue;
                for (int j = 0; j < HDim; j++)
                {
                    l[HDim * a + j, HDim * b + j] += minusI * h[a, b];             // −i·H·ρ: (Hρ)_{a,j} ∋ H_{a,b}·ρ_{b,j}
                    l[HDim * j + b, HDim * j + a] += GaussianRational.I * h[a, b]; // +i·ρ·H: (ρH)_{j,b} ∋ ρ_{j,a}·H_{a,b}
                }
            }
        return l;
    }

    private static Mat LeftMultiplication(Mat h)
    {
        var l = Zero(Dim);
        for (int i = 0; i < HDim; i++)
            for (int k = 0; k < HDim; k++)
                if (!h[i, k].IsZero)
                    for (int j = 0; j < HDim; j++) l[HDim * i + j, HDim * k + j] = h[i, k];
        return l;
    }

    private static Mat RightMultiplication(Mat h)
    {
        var l = Zero(Dim);
        for (int k = 0; k < HDim; k++)
            for (int j = 0; j < HDim; j++)
                if (!h[k, j].IsZero)
                    for (int i = 0; i < HDim; i++) l[HDim * i + j, HDim * i + k] = h[k, j];
        return l;
    }

    private static int Inner(int level) => level == 1 || level == 2 ? 1 : 0;
    private static int InnerCount(int state) => Inner(state >> 2) + Inner(state & 3);

    /// <summary>The bright set of the cell |i⟩⟨j|: bit k set when nucleus k is outer on one side and inner on the other.</summary>
    private static int BrightMask(int cell)
    {
        int i = cell / HDim, j = cell % HDim;
        return (Inner(i >> 2) != Inner(j >> 2) ? 1 : 0) | (Inner(i & 3) != Inner(j & 3) ? 2 : 0);
    }

    private static int BrightCount(int cell) => BitOperations.PopCount((uint)BrightMask(cell));

    /// <summary>S·L·S⁻¹ entry by entry: S = D_Q∘T, (S L S⁻¹)[a, b] = Q(a)·L[Ta, Tb]·Q(b), Q(|i⟩⟨j|) = (−1)^b(i).</summary>
    private static Mat ConjugateByS(Mat l)
    {
        var r = Zero(Dim);
        for (int a = 0; a < Dim; a++)
            for (int b = 0; b < Dim; b++)
            {
                var v = l[Transposed(a), Transposed(b)];
                if (v.IsZero) continue;
                r[a, b] = (InnerCount(a / HDim) + InnerCount(b / HDim)) % 2 == 0 ? v : -v;
            }
        return r;
    }

    private static int Transposed(int cell) => (cell % HDim) * HDim + cell / HDim;

    private static bool Lemma(Mat h, IReadOnlyList<(long J0, long J1, long J2)> rates)
    {
        var left = ConjugateByS(Generator(h, rates));
        var right = Scale(Generator(h, rates, flipJ0: true), -1);
        long shift = SodiumChainMirrorLadderClaim.LemmaConstant(rates);
        for (int r = 0; r < Dim; r++) right[r, r] -= (GaussianRational)(BigRational)shift;
        return Equal(left, right);
    }

    private static bool KeepsBrightness(Mat l)
    {
        for (int r = 0; r < Dim; r++)
            for (int c = 0; c < Dim; c++)
                if (!l[r, c].IsZero && BrightMask(r) != BrightMask(c)) return false;
        return true;
    }

    private static IReadOnlyList<LadderBlock> Ladder(Mat l, IReadOnlyList<(long J0, long J1, long J2)> rates)
    {
        var blocks = new List<LadderBlock>();
        for (int mask = 0; mask < 1 << Nuclei; mask++)
        {
            var idx = Enumerable.Range(0, Dim).Where(r => BrightMask(r) == mask).ToArray();
            long twoS = SodiumChainMirrorLadderClaim.TwiceCentre(rates, mask);
            var chi = BlockCharpoly(l, idx, out var scale);
            var neighbours = Enumerable.Range(0, Nuclei)
                .Select(k => SodiumChainMirrorLadderClaim.TwiceCentre(rates, mask ^ (1 << k)))
                .Select(t => (t, IsPalindromic(chi, scale * t))).ToArray();
            blocks.Add(new LadderBlock(mask, idx.Length, twoS, IsPalindromic(chi, scale * twoS), neighbours));
        }
        return blocks;
    }

    /// <summary>The full coupling in the frame Φ: L_rest + √3·L₃ with every rational entry multiplied
    /// by 3^((β_r − β_c)/2) and every √3 entry by 3^((β_r − β_c + 1)/2); rational exactly when the
    /// parities fall as stated (returned).</summary>
    private static (Mat L, bool Rational) Rationalised(Mat hRest, Mat root, IReadOnlyList<(long J0, long J1, long J2)> rates)
    {
        var rest = Generator(hRest, rates);
        var l3 = HamiltonianPart(root);
        var r = Zero(Dim);
        bool ok = true;
        for (int a = 0; a < Dim; a++)
            for (int b = 0; b < Dim; b++)
            {
                int t = Beta(a) - Beta(b);
                bool odd = (t & 1) != 0;
                if (!rest[a, b].IsZero)
                {
                    ok &= !odd;
                    if (!odd) r[a, b] += rest[a, b] * PowerOfThree(t / 2);
                }
                if (!l3[a, b].IsZero)
                {
                    ok &= odd;
                    if (odd) r[a, b] += l3[a, b] * PowerOfThree((t + 1) / 2);
                }
            }
        return (r, ok);
    }

    /// <summary>√3·D·H₃·D⁻¹ with D = diag(3^(−b/2)): the √3 coupling in the frame Φ, rational since
    /// every entry of H₃ joins states of opposite inner-count parity.</summary>
    private static Mat SimilarH3(Mat root)
    {
        var r = Zero(HDim);
        for (int a = 0; a < HDim; a++)
            for (int c = 0; c < HDim; c++)
                if (!root[a, c].IsZero) r[a, c] = root[a, c] * PowerOfThree((InnerCount(c) - InnerCount(a) + 1) / 2);
        return r;
    }

    private static int Beta(int cell) => InnerCount(cell % HDim) - InnerCount(cell / HDim);

    private static GaussianRational PowerOfThree(int e)
    {
        var p = BigInteger.Pow(3, Math.Abs(e));
        return e >= 0 ? (GaussianRational)new BigRational(p) : (GaussianRational)new BigRational(BigInteger.One, p);
    }

    /// <summary>F158's two ends, dim ker L and dim ker(L + 2σ), as exact ranks per connected component.</summary>
    private static (int Near, int Far) TwoEnds(Mat l, long twoSigma)
    {
        int near = 0, far = 0;
        foreach (var comp in Components(l))
        {
            var sub = Submatrix(l, comp);
            near += comp.Length - Rank(sub);
            for (int i = 0; i < comp.Length; i++) sub[i, i] += (GaussianRational)(BigRational)twoSigma;
            far += comp.Length - Rank(sub);
        }
        return (near, far);
    }

    private static bool SpectrumPairs(Mat l, int[] idx, long twoS)
    {
        var chi = BlockCharpoly(l, idx, out var scale);
        return IsPalindromic(chi, scale * twoS);
    }

    // ===== invariant hulls =====

    /// <summary>For m = 0 … n: the largest L-invariant subspace inside the cells with m bright nuclei
    /// and the common kernel of <paramref name="annihilators"/>, by Krylov stabilisation of the
    /// constraint rows (a row space closed under v ↦ v·L), per connected component of the union graph;
    /// the χ of L restricted to it and its scale are returned in <paramref name="chis"/>.</summary>
    private static IReadOnlyList<HullReading> Hulls(Mat l, IReadOnlyList<(long J0, long J1, long J2)> rates, Mat[] annihilators,
        ref bool invariant, out List<Mat>[] restrictionsByM)
    {
        var comps = Components(new[] { l }.Concat(annihilators).ToArray());
        var readings = new List<HullReading>();
        restrictionsByM = new List<Mat>[Nuclei + 1];
        for (int m = 0; m <= Nuclei; m++)
        {
            var restrictions = new List<Mat>();
            int dim = 0;
            foreach (var comp in comps)
            {
                var (k, restriction, inv) = Hull(l, comp, r => BrightCount(r) == m, annihilators);
                invariant &= inv;
                dim += k;
                if (k > 0) restrictions.Add(restriction!);
            }
            var scale = CommonScale(restrictions);
            long twoS = SodiumChainMirrorLadderClaim.TwiceRung(Nuclei, m, rates[0].J0, rates[0].J1, rates[0].J2);
            int cells = Enumerable.Range(0, Dim).Count(r => BrightCount(r) == m);
            readings.Add(new HullReading(m, cells, dim, twoS, dim == 0 ? null : IsPalindromic(Chi(restrictions, scale), scale * twoS)));
            restrictionsByM[m] = restrictions;
        }
        return readings;
    }

    private static BigInteger CommonScale(IEnumerable<Mat> ms)
    {
        BigInteger scale = BigInteger.One;
        foreach (var mm in ms) scale = BigRational.Lcm(scale, DenominatorLcm(mm));
        return scale;
    }

    /// <summary>det(x − scale·M) for the direct sum of the blocks <paramref name="ms"/>.</summary>
    private static GaussianInteger[] Chi(IEnumerable<Mat> ms, BigInteger scale)
    {
        var chi = new[] { GaussianInteger.One };
        foreach (var mm in ms) chi = GaussianPolynomial.Multiply(chi, ScaledCharpoly(mm, scale));
        return chi;
    }

    private static (int Dim, Mat? Restriction, bool Invariant) Hull(Mat l, int[] comp, Func<int, bool> inside, Mat[] annihilators)
    {
        int s = comp.Length;
        var lc = Submatrix(l, comp);
        var basis = new List<GaussianRational[]>();
        var pivots = new List<int>();
        var queue = new Queue<GaussianRational[]>();

        void Offer(GaussianRational[] row)
        {
            var v = (GaussianRational[])row.Clone();
            for (int t = 0; t < basis.Count; t++)
            {
                var f = v[pivots[t]];
                if (f.IsZero) continue;
                for (int col = 0; col < s; col++)
                    if (!basis[t][col].IsZero) v[col] -= f * basis[t][col];
            }
            int piv = Array.FindIndex(v, z => !z.IsZero);
            if (piv < 0) return;
            var inv = GaussianRational.One / v[piv];
            for (int col = 0; col < s; col++) v[col] *= inv;
            foreach (var b in basis)                              // keep the basis fully reduced
            {
                var f = b[piv];
                if (f.IsZero) continue;
                for (int col = 0; col < s; col++)
                    if (!v[col].IsZero) b[col] -= f * v[col];
            }
            basis.Add(v);
            pivots.Add(piv);
            queue.Enqueue(v);
        }

        for (int p = 0; p < s; p++)
            if (!inside(comp[p]))
            {
                var e = ZeroRow(s);
                e[p] = GaussianRational.One;
                Offer(e);
            }
        foreach (var x in annihilators)
            for (int g = 0; g < Dim; g++)
            {
                var row = ZeroRow(s);
                bool any = false;
                for (int p = 0; p < s; p++)
                    if (!x[g, comp[p]].IsZero) { row[p] = x[g, comp[p]]; any = true; }
                if (any) Offer(row);
            }
        while (queue.Count > 0)
        {
            var v = queue.Dequeue();
            var w = ZeroRow(s);
            for (int i = 0; i < s; i++)
            {
                if (v[i].IsZero) continue;
                for (int j = 0; j < s; j++)
                    if (!lc[i, j].IsZero) w[j] += v[i] * lc[i, j];
            }
            Offer(w);
        }

        var free = Enumerable.Range(0, s).Where(c => !pivots.Contains(c)).ToArray();
        int k = free.Length;
        if (k == 0) return (0, null, true);
        var basisMat = Zero(s, k);
        for (int f = 0; f < k; f++)
        {
            basisMat[free[f], f] = GaussianRational.One;
            for (int t = 0; t < basis.Count; t++) basisMat[pivots[t], f] = -basis[t][free[f]];
        }
        var lb = Mul(lc, basisMat);
        var restriction = Zero(k, k);
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++) restriction[a, b] = lb[free[a], b];  // the basis is the identity on the free rows
        return (k, restriction, Equal(lb, Mul(basisMat, restriction)));
    }

    private static bool PageQuotientMatches(GaussianInteger[] chiK, GaussianInteger[] chiU, BigInteger scale)
    {
        // det(x − σM) is the χ of σM, whose eigenvalues are σ times L's: the factor is x(x + 48σ)²((x + 36σ)² + 1600σ²)
        var lin = (BigInteger v) => new[] { new GaussianInteger(v * scale, BigInteger.Zero), GaussianInteger.One };
        var x = new[] { GaussianInteger.Zero, GaussianInteger.One };
        var q = GaussianPolynomial.Multiply(x, GaussianPolynomial.Multiply(lin(48), lin(48)));
        var sq = GaussianPolynomial.Multiply(lin(36), lin(36));
        sq[0] += new GaussianInteger(1600 * scale * scale, BigInteger.Zero);
        var expected = GaussianPolynomial.Trim(GaussianPolynomial.Multiply(chiU, GaussianPolynomial.Multiply(q, sq)));
        var k = GaussianPolynomial.Trim(chiK);
        return k.Length == expected.Length && k.Zip(expected).All(z => z.First.Equals(z.Second));
    }

    // ===== characteristic polynomials =====

    /// <summary>χ of the principal submatrix on <paramref name="idx"/>, scaled to Gaussian integers by
    /// the lcm of its denominators (returned): the product over the connected components inside it.</summary>
    private static GaussianInteger[] BlockCharpoly(Mat l, int[] idx, out BigInteger scale)
    {
        var sub = Submatrix(l, idx);
        scale = DenominatorLcm(sub);
        var chi = new[] { GaussianInteger.One };
        foreach (var comp in Components(sub))
            chi = GaussianPolynomial.Multiply(chi, ScaledCharpoly(Submatrix(sub, comp), scale));
        return chi;
    }

    private static BigInteger DenominatorLcm(Mat m)
    {
        BigInteger c = BigInteger.One;
        foreach (var z in m) c = BigRational.Lcm(BigRational.Lcm(c, z.Re.Denominator), z.Im.Denominator);
        return c;
    }

    /// <summary>det(x − scale·M) over ℤ[i] (Berkowitz, division free).</summary>
    private static GaussianInteger[] ScaledCharpoly(Mat m, BigInteger scale)
    {
        int n = m.GetLength(0);
        var g = new GaussianInteger[n, n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
            {
                var re = m[a, b].Re * scale;
                var im = m[a, b].Im * scale;
                if (!re.IsInteger || !im.IsInteger) throw new InvalidOperationException("scale does not clear the denominators");
                g[a, b] = new GaussianInteger(re.Numerator, im.Numerator);
            }
        return GaussianPolynomial.Trim(GaussianMatrixCharpoly.Characteristic(g));
    }

    /// <summary>χ(−x − shift) = (−1)^deg·χ(x) identically: the spectrum is symmetric about −shift/2.</summary>
    private static bool IsPalindromic(GaussianInteger[] chi, BigInteger shift)
    {
        chi = GaussianPolynomial.Trim(chi);
        var mirrored = GaussianPolynomial.Trim(GaussianPolynomial.ComposeLinear(chi, -GaussianInteger.One, new GaussianInteger(-shift, BigInteger.Zero)));
        if (mirrored.Length != chi.Length) return false;
        bool odd = (chi.Length - 1) % 2 == 1;
        for (int i = 0; i < chi.Length; i++)
            if (!mirrored[i].Equals(odd ? -chi[i] : chi[i])) return false;
        return true;
    }

    // ===== components =====

    /// <summary>The connected components of the union of the nonzero patterns (as undirected graphs).</summary>
    private static List<int[]> Components(params Mat[] ms)
    {
        int n = ms[0].GetLength(0);
        var parent = Enumerable.Range(0, n).ToArray();
        int Find(int a) { while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; }
        foreach (var m in ms)
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                    if (a != b && !m[a, b].IsZero) parent[Find(a)] = Find(b);
        return Enumerable.Range(0, n).GroupBy(Find).Select(g => g.ToArray()).ToList();
    }

    private static Mat Submatrix(Mat l, int[] idx)
    {
        var r = Zero(idx.Length);
        for (int a = 0; a < idx.Length; a++)
            for (int b = 0; b < idx.Length; b++) r[a, b] = l[idx[a], idx[b]];
        return r;
    }

    // ===== strings =====

    private static readonly Dictionary<char, Mat> Letters = new()
    {
        ['I'] = new GaussianRational[,] { { 1, 0 }, { 0, 1 } },
        ['X'] = new GaussianRational[,] { { 0, 1 }, { 1, 0 } },
        ['Y'] = new GaussianRational[,] { { 0, -GaussianRational.I }, { GaussianRational.I, 0 } },
        ['Z'] = new GaussianRational[,] { { 1, 0 }, { 0, -1 } },
    };

    private static Mat Str(string word)
    {
        Mat m = new GaussianRational[,] { { 1 } };
        foreach (var ch in word) m = Kron(m, Letters[ch]);
        return m;
    }

    private static string[] AllStrings(int qubits)
    {
        IEnumerable<string> words = new[] { "" };
        for (int q = 0; q < qubits; q++) words = words.SelectMany(w => "IXYZ".Select(ch => w + ch));
        return words.ToArray();
    }

    private static bool Anticommutes(string a, string b) =>
        a.Zip(b).Count(z => z.First != 'I' && z.Second != 'I' && z.First != z.Second) % 2 == 1;

    /// <summary>The Pauli-string expansion of a 16 × 16 matrix, c_w = Tr(P_w·M)/16, nonzero terms only.</summary>
    private static Dictionary<string, GaussianRational> Decompose(Mat m)
    {
        var r = new Dictionary<string, GaussianRational>();
        var sixteenth = (GaussianRational)new BigRational(1, HDim);
        foreach (var w in AllStrings(4))
        {
            var c = Trace(Mul(Str(w), m)) * sixteenth;
            if (!c.IsZero) r[w] = c;
        }
        return r;
    }

    private static Mat Recompose(IReadOnlyDictionary<string, GaussianRational> terms)
    {
        var m = Zero(HDim);
        foreach (var (w, c) in terms) m = Add(m, Scale(Str(w), c));
        return m;
    }

    /// <summary>The secular part with respect to q₁ZZ₁ + q₂ZZ₂: the entries joining states of equal energy.</summary>
    private static Mat Secular(Mat m, long q1, long q2)
    {
        long Energy(int s) => q1 * (Inner(s >> 2) == 0 ? 1 : -1) + q2 * (Inner(s & 3) == 0 ? 1 : -1);
        var r = Zero(HDim);
        for (int a = 0; a < HDim; a++)
            for (int b = 0; b < HDim; b++)
                if (Energy(a) == Energy(b)) r[a, b] = m[a, b];
        return r;
    }

    private static bool IsRealSymmetric(Mat m)
    {
        for (int a = 0; a < m.GetLength(0); a++)
            for (int b = 0; b < m.GetLength(1); b++)
                if (!m[a, b].Im.IsZero || m[a, b] != m[b, a]) return false;
        return true;
    }

    // ===== ℚ(√3) for the spin matrices =====

    /// <summary>A + √3·B with A, B matrices over ℚ(i).</summary>
    private readonly record struct Root3(Mat A, Mat B)
    {
        public static Root3 operator *(Root3 u, Root3 v) =>
            new(Add(Mul(u.A, v.A), SodiumChainWitness.Scale(Mul(u.B, v.B), 3)), Add(Mul(u.A, v.B), Mul(u.B, v.A)));
        public static Root3 operator +(Root3 u, Root3 v) => new(Add(u.A, v.A), Add(u.B, v.B));
        public Root3 Scale(GaussianRational s) => new(SodiumChainWitness.Scale(A, s), SodiumChainWitness.Scale(B, s));
        public static Root3 Kron(Root3 u, Root3 v) =>
            new(Add(SodiumChainWitness.Kron(u.A, v.A), SodiumChainWitness.Scale(SodiumChainWitness.Kron(u.B, v.B), 3)),
                Add(SodiumChainWitness.Kron(u.A, v.B), SodiumChainWitness.Kron(u.B, v.A)));
    }

    private static (GaussianRational A, GaussianRational B) SqrtInQ3(BigRational r)
    {
        if (RationalSqrt(r) is { } a) return (a, GaussianRational.Zero);
        if (RationalSqrt(r / 3) is { } b) return (GaussianRational.Zero, b);
        throw new InvalidOperationException($"√{r} is not in ℚ(√3)");
    }

    private static BigRational? RationalSqrt(BigRational r)
    {
        if (r.Sign < 0) return null;
        var n = ISqrt(r.Numerator);
        var d = ISqrt(r.Denominator);
        return n * n == r.Numerator && d * d == r.Denominator ? new BigRational(n, d) : null;
    }

    private static BigInteger ISqrt(BigInteger v)
    {
        if (v.IsZero) return v;
        var x = v;
        var y = (x + 1) / 2;
        while (y < x) { x = y; y = (x + v / x) / 2; }
        return x;
    }

    // ===== exact matrices over ℚ(i) =====

    private static GaussianRational[] ZeroRow(int n)
    {
        var r = new GaussianRational[n];
        for (int i = 0; i < n; i++) r[i] = GaussianRational.Zero;
        return r;
    }

    private static Mat Zero(int n) => Zero(n, n);

    private static Mat Zero(int n, int p)
    {
        var m = new GaussianRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) m[i, j] = GaussianRational.Zero;
        return m;
    }

    private static Mat Identity(int n)
    {
        var m = Zero(n);
        for (int i = 0; i < n; i++) m[i, i] = GaussianRational.One;
        return m;
    }

    private static Mat Mul(Mat a, Mat b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), p = b.GetLength(1);
        var r = Zero(n, p);
        for (int i = 0; i < n; i++)
            for (int t = 0; t < k; t++)
            {
                if (a[i, t].IsZero) continue;
                for (int j = 0; j < p; j++)
                    if (!b[t, j].IsZero) r[i, j] += a[i, t] * b[t, j];
            }
        return r;
    }

    private static Mat Add(Mat a, Mat b) => Map2(a, b, (u, v) => u + v);
    private static Mat Sub(Mat a, Mat b) => Map2(a, b, (u, v) => u - v);
    private static Mat Scale(Mat a, GaussianRational s) => Map2(a, a, (u, _) => s * u);
    private static Mat Scale(Mat a, BigRational s) => Scale(a, (GaussianRational)s);

    private static Mat Map2(Mat a, Mat b, Func<GaussianRational, GaussianRational, GaussianRational> f)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = new GaussianRational[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[i, j] = f(a[i, j], b[i, j]);
        return r;
    }

    private static Mat Transpose(Mat a)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var r = new GaussianRational[p, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) r[j, i] = a[i, j];
        return r;
    }

    private static Mat Kron(Mat a, Mat b)
    {
        int n = a.GetLength(0), m = b.GetLength(0);
        var r = new GaussianRational[n * m, n * m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < m; k++)
                    for (int l = 0; l < m; l++) r[i * m + k, j * m + l] = a[i, j] * b[k, l];
        return r;
    }

    private static GaussianRational Trace(Mat a)
    {
        var s = GaussianRational.Zero;
        for (int i = 0; i < a.GetLength(0); i++) s += a[i, i];
        return s;
    }

    private static bool Equal(Mat a, Mat b)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1)) return false;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    private static bool IsZero(Mat a)
    {
        foreach (var z in a) if (!z.IsZero) return false;
        return true;
    }

    /// <summary>Exact rank over ℚ(i) by Gaussian elimination.</summary>
    private static int Rank(Mat a)
    {
        var m = (Mat)a.Clone();
        int rows = m.GetLength(0), cols = m.GetLength(1), rank = 0;
        for (int col = 0; col < cols && rank < rows; col++)
        {
            int piv = -1;
            for (int r = rank; r < rows; r++) if (!m[r, col].IsZero) { piv = r; break; }
            if (piv < 0) continue;
            for (int c = 0; c < cols; c++) (m[rank, c], m[piv, c]) = (m[piv, c], m[rank, c]);
            for (int r = rank + 1; r < rows; r++)
            {
                if (m[r, col].IsZero) continue;
                var f = m[r, col] / m[rank, col];
                for (int c = col; c < cols; c++) m[r, c] -= f * m[rank, c];
            }
            rank++;
        }
        return rank;
    }

    // ===== IInspectable =====

    public string DisplayName => $"SodiumChainWitness (two nuclei; J₀={J0}, J₁={J1}, J₂={J2}, d={D}, c={C})";

    private static string Verdict(bool ok, string yes, string no) => ok ? yes : no;

    private static string ShowHulls(IReadOnlyList<HullReading> hs) =>
        string.Join(", ", hs.Select(h => $"m = {h.M}: {h.Dimension} of {h.CellCount}" +
            (h.Pairs is null ? " (empty)" : h.Pairs.Value ? $", pairs about −s_{h.M} = −{h.TwoS / 2}" : $", NOT about −s_{h.M} = −{h.TwoS / 2}")));

    private static string ShowBlocks(IReadOnlyList<LadderBlock> bs) =>
        string.Join("; ", bs.Select(b => $"B = {b.Label}: {(b.Pairs ? "pairs" : "does NOT pair")} about −s_B, 2s_B = {b.TwoS}, " +
            $"and {(b.Neighbours.All(n => !n.Pairs) ? "about no neighbour" : "about a NEIGHBOUR")} ({string.Join(", ", b.Neighbours.Select(n => n.TwoS))})"));

    public string Summary =>
        $"two spin-3/2 nuclei, the coupling read from the spin matrices, exact over ℚ(i): the dictionary " +
        $"{Verdict(DictionaryHolds, "holds", "FAILS")} (22 dipolar strings, the eight √3 ones those anticommuting with Q); " +
        $"S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ(J₁ + J₂) {Verdict(LemmaHolds, "holds", "FAILS")} for the Ising coupling and both truncations, " +
        $"the √3 strings at c = {C} {Verdict(SqrtThreeStringsBreakLemma, "breaking it", "FAILING to break it")}; the ladder " +
        $"{Verdict(LadderHolds, "holds", "FAILS")} (every bright-set block about its own rung and no neighbour, the equal truncation moving brightness); " +
        $"the full coupling at J₀ = 0 {Verdict(GlobalBreakHolds, "breaks the global mirror", "does NOT break as stated")}: " +
        $"dim ker L = {NearKernel} against dim ker(L + 2σ) = {FarKernel} at √3 ({NearKernelAtC} against {FarKernelAtC} at c), " +
        $"χ {(CharpolyPairsAboutSigma ? "=" : "≠")} its mirror about −σ, 2σ = {TwoSigma}. Readings: K_m {string.Join(" + ", EqualSharp.Select(h => h.Dimension))} " +
        $"(equal truncation), K_m {string.Join(" + ", FullSharp.Select(h => h.Dimension))} and U_m {string.Join(" + ", FullKernelHull.Select(h => h.Dimension))} (full coupling), " +
        $"the invariance L·B = B·M {Verdict(HullsInvariant, "checked on every hull", "FAILING")}; the explicit E₁ ({E1Dimension}) and E₂ ({E2Dimension}) " +
        $"{Verdict(ExplicitHolds, "lie in U₁ and U₂ at every parameter with identical rates, one offset and an equal splitting", "do NOT hold as stated")}" +
        $"{(ExplicitHolds && ExplicitAreTheSharpSpacesHere ? " and are K₁ and K₂ here" : "")}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the gate: the dictionary and the couplings",
                summary: $"from the spin-3/2 matrices with I_± in ℚ(√3): I_z = ZI + ½·IZ {Verdict(IzMatches, "✓", "FAILS")}, " +
                         $"I_x = (√3/2)·IX + ½(XX + YY) {Verdict(IxMatches, "✓", "FAILS")}, I_y = (√3/2)·IY + ½(YX − XY) {Verdict(IyMatches, "✓", "FAILS")}, " +
                         $"T₂₀ = (3/√6)·ZZ {Verdict(T20Matches, "✓", "FAILS")}, the √3 parts {Verdict(SqrtThreePartNeeded, "nonzero (dropping them is wrong)", "ZERO")}; " +
                         $"3I_zI_z − I·I has {DipolarStringCount} strings, {SqrtThreeStrings.Count} with √3 ({string.Join(", ", SqrtThreeStrings)}), " +
                         $"{Verdict(PartsDisjoint, "the two parts disjoint and reproducing the matrix", "the parts NOT disjoint or NOT reproducing it")}; " +
                         $"the √3 strings {Verdict(SqrtThreeAreTheQAnticommuting, "are exactly those anticommuting with ZZ₁ZZ₂ (not those anticommuting with ZZ₁ alone)", "are NOT the Q-anticommuting ones")}; " +
                         $"the rational strings anticommuting with a ZZ_k are {string.Join(", ", SatelliteSatelliteStrings)}, " +
                         $"{Verdict(SatelliteSatelliteAnticommuteWithEach, "anticommuting with each", "NOT with each")}; both parts {Verdict(CouplingsReal, "real symmetric", "NOT real symmetric")}; " +
                         $"the shell projection {Verdict(UnequalTruncationIsShellProjection, "is", "is NOT")} the commuting part at q = ({UnequalFrequencies.Q1}, {UnequalFrequencies.Q2}) " +
                         $"and {Verdict(EqualTruncationIsShellProjection, "is", "is NOT")} that plus (IXIX + IYIY)·(1 − ZZ₁ZZ₂)/2 at q = ({EqualFrequencies.Q}, {EqualFrequencies.Q}), " +
                         $"the √3 part secular in neither; the two truncations {Verdict(TruncationsDiffer, "differ", "do NOT differ")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the lemma: the mirror image is the J₀-reversed chain",
                summary: $"S(ρ) = Q·ρᵀ, S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ_k(J₁,k + J₂,k), entry by entry on the 256 × 256 generator: " +
                         $"Ising with rates per nucleus ({J0}, {J1}, {J2}) and {SecondNucleusRates} {Verdict(LemmaIsingPerNucleus, "holds", "FAILS")}, " +
                         $"Ising with identical rates {Verdict(LemmaIsingIdentical, "holds", "FAILS")}, the unequal truncation {Verdict(LemmaUnequal, "holds", "FAILS")}, " +
                         $"the equal truncation {Verdict(LemmaEqual, "holds", "FAILS")}; control: the full coupling with its √3 strings at the rational " +
                         $"coefficient c = {C} (a stand-in for √3; the page proves the break for any non-zero coefficient) {Verdict(SqrtThreeStringsBreakLemma, "breaks it", "does NOT break it")}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the ladder: one centre per bright set",
                summary: $"[L, P_k] = 0 for every k: Ising {Verdict(IsingKeepsBrightness, "yes", "NO")}, unequal truncation {Verdict(UnequalKeepsBrightness, "yes", "NO")}, " +
                         $"equal truncation {Verdict(!EqualKeepsBrightness, "no (the satellite flip-flop moves brightness; the control)", "YES (the control FAILS)")}. " +
                         $"χ of each 64-dimensional bright-set block, Ising with rates per nucleus: {ShowBlocks(IsingBlocks)}. " +
                         $"Unequal truncation, identical rates: {ShowBlocks(UnequalBlocks)}.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the full coupling at J₀ = 0: the global mirror breaks",
                summary: $"F158's far span: the strings anticommuting with every jump are {string.Join(", ", FarSpanAtJ0Zero)}, " +
                         $"{FarSpanWithJ0Jumps} once the J₀ jumps are on; the rest of H {Verdict(RestCommutesWithQ, "commutes", "does NOT commute")} with Q, " +
                         $"the √3 part {Verdict(SqrtThreeAnticommutesWithQ, "anticommutes", "does NOT anticommute")}. In the frame Φ the generator " +
                         $"{Verdict(SimilarityIsRational, "is rational (every √3 entry on an odd β difference)", "is NOT rational")}, so at √3 itself: " +
                         $"dim ker L = {NearKernel}, dim ker(L + 2σ) = {FarKernel}, 2σ = {TwoSigma}; at the stand-in c = {C}: {NearKernelAtC} against {FarKernelAtC}; " +
                         $"route B, χ(−x − 2σ) {(CharpolyPairsAboutSigma ? "=" : "≠")} χ(x) at √3. Control, the equal truncation at J₀ = 0: " +
                         $"{ControlNearKernel} against {ControlFarKernel}, χ {(ControlCharpolyPairs ? "palindromic" : "NOT palindromic")} " +
                         $"({Verdict(GlobalBreakHolds, "both routes as the page states", "a route DISAGREES with the page")}).",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "the explicit sharp spaces, at every parameter",
                summary: $"E₁ (σ = X + iY: Re and Im of σσσ on three qubits, I or Z on the fourth; the whole order ±5 block, TQ ⊗ DQ, " +
                         $"and TQ ⊗ satellite − satellite ⊗ TQ at order ±4) has dimension {E1Dimension}, " +
                         $"{Verdict(E1InOneBrightCells, "inside the cells with one bright nucleus", "OUTSIDE the one-bright cells")}; " +
                         $"E₂ = (Z₁ + Z₃)·(X₂X₄ + Y₂Y₄), the satellite flip-flop signed by Z₁ + Z₃, dimension {E2Dimension}, " +
                         $"{Verdict(E2InTwoBrightCells, "inside the two-bright cells", "OUTSIDE the two-bright cells")}. Each is invariant under " +
                         $"every piece of L separately (J₀, J₁, J₂ with identical rates, ω·I_z,tot, an equal splitting q(ZZ₁ + ZZ₂), the rational and the √3 coupling): " +
                         $"E₁ {Verdict(E1InvariantUnderEveryPiece, "yes", "NO")}, E₂ {Verdict(E2InvariantUnderEveryPiece, "yes", "NO")}; " +
                         $"both annihilated by the √3 coupling from both sides: E₁ {Verdict(E1InSqrtThreeKernel, "yes", "NO")}, E₂ {Verdict(E2InSqrtThreeKernel, "yes", "NO")}; " +
                         $"E₂ {Verdict(E2OnItsRung, "is an eigenoperator at −s₂ = −6(J₀ + J₁ + J₂), on its rung", "is NOT on its rung")}. So E_m ⊆ U_m ⊆ K_m " +
                         $"at every parameter with identical rates, one offset and an equal splitting. Controls: with rates per nucleus E₁ {Verdict(E1BreaksWithRatesPerNucleus, "breaks", "does NOT break")} " +
                         $"and E₂ {Verdict(E2BreaksWithRatesPerNucleus, "breaks", "does NOT break")}, while the order ±5 block of E₁ " +
                         $"{Verdict(E1TopBlockKeepsRatesPerNucleus, "keeps", "does NOT keep")} its invariance; the exchange-symmetric order ±4 partners " +
                         $"{Verdict(SymmetricPartnersKeepRationalCoupling, "keep", "do NOT keep")} the rational coupling and " +
                         $"{Verdict(SymmetricPartnersLeaveSqrtThreeKernel, "leave", "do NOT leave")} the kernel of the √3 one, and one bare TQ ⊗ satellite column " +
                         $"{Verdict(BareColumnBreaksRationalCoupling, "breaks", "does NOT break")} the rational coupling; a splitting on one nucleus alone " +
                         $"{Verdict(UnequalSplittingBreaksBoth, "breaks", "does NOT break")} both spaces, and so {Verdict(UnequalOffsetBreaksBoth, "does", "does NOT")} " +
                         $"an offset on one nucleus alone. At these parameters dim K_m = dim U_m = dim E_m " +
                         $"for m = 1, 2: {Verdict(ExplicitAreTheSharpSpacesHere, "yes, so E_m = U_m = K_m here", "NO")}. " +
                         $"I_z,tot {Verdict(IzTotalIsEigenoperator, "is", "is not")} an eigenoperator of L (J₁ {(J1 == J2 ? "=" : "≠")} J₂), " +
                         $"{Verdict(IzTotalInSqrtThreeKernel, "inside", "outside")} the kernel of the √3 coupling.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "readings: the sharp subspaces",
                summary: $"exact dimensions by Krylov stabilisation over ℚ(i), at these parameters, and χ of L restricted (readings, not gates): " +
                         $"equal truncation K_m: {ShowHulls(EqualSharp)}; full coupling at √3, K_m: {ShowHulls(FullSharp)}; " +
                         $"U_m inside the two-sided kernel of the √3 coupling (rank {SqrtThreeRank} of 16): {ShowHulls(FullKernelHull)}" +
                         $"{(FullKernelHull.Zip(FullSharp).Where(z => z.First.Dimension == z.Second.Dimension && z.First.Dimension > 0).Select(z => z.First.M).ToArray() is { Length: > 0 } same ? $"; U_m = K_m (U_m ⊆ K_m, equal dimension) at m = {string.Join(", ", same)}" : "")}" +
                         $"{(K0BeyondU0MatchesPage is { } kq ? $"; χ(K₀)/χ(U₀) {(kq ? "=" : "≠")} x(x + 48)²((x + 36)² + 1600)" : "")}; " +
                         $"invariance L·B = B·M {Verdict(HullsInvariant, "checked on every hull", "FAILS")}.",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
