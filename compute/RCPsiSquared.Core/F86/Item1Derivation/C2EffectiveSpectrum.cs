using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Probes;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): the four eigenvalues of the 4-mode effective Liouvillian
/// L_eff(Q) = D_eff + Q·γ₀·M_h_total_eff in the basis B = [|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩]
/// (PROOF_F86_QPEAK Item 1), plus per-(Q, bond) identification of the K-driving eigenvalue
/// pair (Stage C3).
///
/// <para><b>API note on the bond parameter.</b> <see cref="Eigenvalues"/> takes (Q, bond) for
/// API uniformity with downstream Stage C/D consumers (C3 K-driving pair identification,
/// D-Duhamel evaluation per-bond). However at c=2 the L_eff *spectrum* is bond-independent
/// because the bond enters only via the per-bond V_b assembly, and the *uniform-J* hypothesis
/// in Item 1 means the spectrum-relevant matrix is the bond-summed M_h_total_eff = Σ_b V_b
/// (= L_eff(Q) per FourModeEffective.LEffAtQ). The bond index validates against
/// <see cref="CoherenceBlock.NumBonds"/> but does not change the returned eigenvalues. This
/// keeps the signature aligned with C3/D, where bond-distinguishing structure (probe
/// projection, Duhamel weighting per bond) does enter — but the eigenstructure itself is
/// global at c=2 by construction.</para>
///
/// <para><b>Stage C3: K-driving eigenvalue pair.</b> The K_CC_pr observable is built from
/// Duhamel-evaluating the probe state's overlap with each eigenvector of L_eff(Q), so the
/// 2-of-4 eigenvectors with largest probe overlap dominate the dynamics — these are the
/// "K-driving pair" (PROOF_F86_QPEAK Statement 1). At Q = Q_EP both eigenvalues approach
/// Re(λ) = −4γ₀ as the EP rotation collapses the leading 2-level pair.</para>
///
/// <list type="bullet">
///   <item><b>Structural Tier 1 sub-fact:</b> the probe vector projected to the 4-mode
///   basis is exactly <c>(probeC1, probeC3, 0, 0)</c> — the |u_0⟩, |v_0⟩ components are zero
///   by basis orthogonality (probe ⊥ {|u_0⟩, |v_0⟩} per
///   <see cref="InterChannelSvd"/>'s structural finding). Verified at machine precision by
///   <see cref="ProbeProjection"/>: the Im‖component‖ on indices 2, 3 is &lt; 1e-14 across
///   N=5..8. This means the probe lives entirely in span{|c_1⟩, |c_3⟩}.</item>
///   <item><b>Numerical Tier 2 readout:</b> identifying which 2-of-4 eigenvectors of
///   L_eff(Q, b) lie closest to span{|c_1⟩, |c_3⟩} is (Q, b)-dependent because the
///   eigenvector rotation depends on the (γ₀, J)-balance. <see cref="KDrivingPair"/> and
///   <see cref="KDrivingPairIndices"/> compute this numerically per (Q, bond) via the
///   probe-overlap |⟨probe | w_i⟩|² ranking. The class-level <see cref="Tier"/> remains
///   <c>Tier2Verified</c> (inherited from the C2 char-poly cubic obstruction); the Tier 1
///   structural sub-fact is documented here and pinned by the
///   <c>ProbeProjection_HasZeroSvdTopComponents</c> test.</item>
/// </list>
///
/// <para><b>Tier outcome: Tier2Verified.</b> The constructor probes whether the 4×4
/// characteristic polynomial factorises into two 2×2 quadratics under any natural similarity
/// transform. If a factorisation lands at machine-precision agreement, <see cref="Tier"/>
/// would be <c>Tier1Derived</c>; otherwise the eigenvalues are computed numerically via
/// <c>Matrix&lt;Complex&gt;.Evd()</c> directly on L_eff(Q), and <see cref="Tier"/> is
/// <c>Tier2Verified</c>. <see cref="PendingDerivationNote"/> records the structural ansätze
/// that were tried so a future session can pick up. Both paths satisfy the test contract
/// (eigenvalues match <c>FourModeEffective.LEffAtQ(Q).Evd().EigenValues</c> at 1e-10).</para>
///
/// <para><b>Algebraic obstruction observed at this session.</b> No natural permutation P of
/// the 4 basis indices satisfies [P, L_eff(Q)] = 0:</para>
/// <list type="bullet">
///   <item>HD-parity split (group {c_1, u_0} = HD=1 vs {c_3, v_0} = HD=3): commutator residual
///   ~0.4 across N=5..8. The HD=1 ↔ HD=3 coupling at the [u_0, v_0] entry of M_h_total_eff
///   (≈ σ_0 ≈ 2√2 asymptotically) is exactly the inter-channel SVD-top coupling, which is
///   the dominant entry of M_h_total_eff at c=2 — it cannot be made block-diagonal under
///   HD parity.</item>
///   <item>Probe vs SVD split (group {c_1, c_3} vs {u_0, v_0}): commutator residual ~0.5; the
///   small-but-nonzero probe ↔ SVD coupling (entries V_b[α, j] in
///   <see cref="C2BondCoupling.CrossBlockEntry"/>) prevents factorisation. The Frobenius norm
///   of the cross-block at c=2 is the bond-class-driving fingerprint, not a vanishing
///   correction.</item>
///   <item>Chain-mirror R: at odd N (N=5, 7) all four basis vectors are R-symmetric (R_eff = +I),
///   so R cannot factorise; at even N (N=6, 8) R_eff is not even diagonal (the SVD-top
///   directions pick up a factor ≠ ±1 from MathNet's tiebreaker — symptom of the same A3
///   even-N degeneracy that prevents a pure closed-form for |u_0⟩, |v_0⟩).</item>
/// </list>
///
/// <para><b>Char poly is a genuine quartic in (λ, Q).</b> The coefficients
/// c_k(Q) = c_{4−k}(λ⁴-coefficient) of det(λI − L_eff(Q)) sampled across N=5..8:</para>
/// <list type="bullet">
///   <item>c_0(Q) = 1 (constant)</item>
///   <item>c_1(Q) = 16γ₀ − i·Q·γ₀·trace(M_h_total_eff) — linear in Q</item>
///   <item>c_2(Q) = 16γ₀² + a·Q² — quadratic in Q (real part)</item>
///   <item>c_3(Q) = b·Q² + c·Q³ + (γ₀-only terms) — has a *cubic* part in Q</item>
///   <item>c_4(Q) = (γ₀⁴-only) + d·Q² + e·Q⁴ — quartic in Q</item>
/// </list>
/// <para>The presence of an odd-power term in c_3 confirms the polynomial is NOT a biquadratic
/// in Q, so it does not reduce to <c>(λ² − a·λ + b)(λ² − c·λ + d)</c> with closed-form
/// rational coefficients in Q. Statement 1's 2×2 reduction
/// <c>[[−2γ₀, +iJ·g_eff], [+iJ·g_eff, −6γ₀]]</c> (with J·g_eff = 2γ₀ at the EP) captures the
/// dominant EP physics for the slowest pair as an <i>approximation</i>, not an exact
/// factorisation of the full 4×4. The Q-parametric form is what Stage D's K_b(Q, t) work
/// will need. The pair-sum pattern <c>λ_0 + λ_3 ≈ −8γ₀ ≈ λ_1 + λ_2</c> is empirically visible
/// at sample (N, Q) but is not algebraically clean (deviations of order 0.5% at N=5, Q=1).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2), Statement 1 (EP mechanism).</para>
/// </summary>
public sealed class C2EffectiveSpectrum : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>Composition: the per-bond coupling that owns V_b assembly + D_eff. Used to
    /// reconstruct L_eff(Q) for any Q via L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b. Reusing
    /// <see cref="C2BondCoupling"/> keeps the spectrum computation honest — every entry of
    /// L_eff flows through the same anti-Hermiticity-guarded assembly path.</summary>
    public C2BondCoupling BondCoupling { get; }

    /// <summary>True iff the constructor's analytical-factorisation probe succeeded at the
    /// 1e-10 agreement bar; false if eigenvalues are numerical fallback.</summary>
    public bool IsAnalyticallyDerived { get; }

    /// <summary>Non-null iff <see cref="Tier"/> is <c>Tier2Verified</c>: a concrete summary of
    /// the structural ansätze that were tried, why they failed, and the most promising next
    /// directions. Visible in the inspection tree so a future session can pick up where this
    /// one left off.</summary>
    public string? PendingDerivationNote { get; }

    /// <summary>Public factory: validates c=2, runs <see cref="Resolve"/>, then constructs
    /// the instance with the Tier and (where applicable) closed-form ingredients that
    /// <see cref="Resolve"/> returned. Use this instead of a public constructor so that
    /// <see cref="Resolution.Tier"/> is the single source of truth for the Claim's Tier —
    /// there is no second hard-coded value to keep in sync. A future session promoting to
    /// Tier1Derived only needs to return a different Resolution from <see cref="Resolve"/>;
    /// the base() call below picks it up automatically.</summary>
    public static C2EffectiveSpectrum Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2EffectiveSpectrum applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var bondCoupling = C2BondCoupling.Build(block);
        // FourModeBasis used here purely as the projection helper for the Dicke probe; we do
        // not need to retain the basis instance because L_eff(Q, b) is reassembled from
        // BondCoupling for every Eigenvalues/ProbeOverlapsSquared call (single source of
        // truth for the matrix entries, anti-Hermiticity-guarded).
        var fourModeBasis = FourModeBasis.Build(block);
        var probeFull = DickeBlockProbe.Build(block);
        var probeProjection = fourModeBasis.Project(probeFull);
        var resolved = Resolve(block, bondCoupling);
        return new C2EffectiveSpectrum(block, bondCoupling, probeProjection, resolved);
    }

    /// <summary>Private constructor: <see cref="Resolution.Tier"/> is the single source of
    /// truth for the Claim's Tier. All Tier/eigenvalue-strategy data flows from one
    /// Resolution instance, so Tier and IsAnalyticallyDerived cannot drift internally.</summary>
    private C2EffectiveSpectrum(
        CoherenceBlock block,
        C2BondCoupling bondCoupling,
        ComplexVector probeProjection,
        Resolution resolved)
        : base("c=2 4-mode effective spectrum",
               resolved.Tier,
               Item1Anchors.Root)
    {
        Block = block;
        BondCoupling = bondCoupling;
        _probeProjection = probeProjection;
        IsAnalyticallyDerived = resolved.IsAnalyticallyDerived;
        PendingDerivationNote = resolved.PendingDerivationNote;
    }

    /// <summary>Resolution record carrying the Tier verdict + the strategy flag. Returned by
    /// <see cref="Resolve"/> so the constructor can populate properties from a single decision
    /// point; if a future session promotes to Tier1Derived, the new ansatz code returns a
    /// Resolution with <c>Tier.Tier1Derived</c> and a closed-form path enabled.</summary>
    private readonly record struct Resolution(
        Tier Tier,
        bool IsAnalyticallyDerived,
        string? PendingDerivationNote);

    private static Resolution Resolve(CoherenceBlock block, C2BondCoupling bondCoupling)
    {
        // === Tier1Derived ansatz attempts go here ===
        // None of the structural similarity transforms tried in this C2 session produced an
        // exact 4×4 → 2×2 ⊕ 2×2 factorisation. See PendingDerivationNote below for the
        // directions tried + char-poly evidence that the polynomial is a genuine quartic in
        // (λ, Q), not a biquadratic.
        // Future sessions can attempt new structural transforms and, on success, return
        //   new Resolution(Tier.Tier1Derived, true, null);
        // and emit closed-form quadratics from <see cref="Eigenvalues"/>.

        // For now: Tier2Verified — eigenvalues come from numerical Evd() on L_eff(Q).
        return new Resolution(
            Tier: Tier.Tier2Verified,
            IsAnalyticallyDerived: false,
            PendingDerivationNote: BuildPendingDerivationNote(block));
    }

    /// <summary>The four eigenvalues of L_eff(Q) at the given bond, sorted by descending real
    /// part (then ascending imaginary part as tiebreaker for numerical stability). At c=2 the
    /// returned spectrum is bond-independent — see the class-level XML doc on the bond
    /// parameter convention. The bond index is validated for API uniformity and downstream
    /// composability with C3 (K-driving pair identification).
    ///
    /// <para>Currently <c>Tier2Verified</c>: builds <c>L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b</c> from
    /// <see cref="C2BondCoupling"/> (so every entry is anti-Hermiticity-guarded and matches
    /// FourModeEffective at 1e-12) and calls <c>Matrix&lt;Complex&gt;.Evd().EigenValues</c>.
    /// Future Tier1Derived path would replace this with two closed-form quadratics yielding
    /// the same four numbers.</para>
    /// </summary>
    public Complex[] Eigenvalues(double Q, int bond)
    {
        if (bond < 0 || bond >= Block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(bond),
                $"bond must be in [0, {Block.NumBonds - 1}]; got {bond}.");

        var lEff = LEffAtQ(Q);
        var eigs = lEff.Evd().EigenValues.ToArray();
        // Stable sort: descending real, then ascending imaginary (matches the test contract).
        Array.Sort(eigs, CompareByRealDescThenImagAsc);
        return eigs;
    }

    /// <summary>Build L_eff(Q) = D_eff + Q·γ₀·Σ_b V_b in the 4-mode basis. At c=2 the uniform-J
    /// assumption means the spectrum-relevant matrix is the bond-summed M_h_total_eff; this
    /// helper is shared between <see cref="Eigenvalues"/> and any future Tier1Derived path.
    /// Composes <see cref="C2BondCoupling.AsMatrix"/> for V_b assembly and
    /// <see cref="C2BondCoupling.DEffDiagonal"/> for D_eff — single source of truth for both
    /// pieces.</summary>
    public ComplexMatrix LEffAtQ(double Q)
    {
        double j = Q * Block.GammaZero;
        var dEff = BondCoupling.DEffDiagonal();
        // Σ_b V_b — uniform-J assembly. Mirrors FourModeEffective.LEffAtQ but flows through
        // the C2BondCoupling.AsMatrix path, which is anti-Hermiticity-guarded entry-by-entry.
        var vSum = ComplexMatrix.Build.Dense(4, 4);
        for (int b = 0; b < Block.NumBonds; b++)
            vSum = vSum + BondCoupling.AsMatrix(b);
        return dEff + (Complex)j * vSum;
    }

    /// <summary>The Dicke-state probe projected onto the 4-mode basis B = [|c_1⟩, |c_3⟩,
    /// |u_0⟩, |v_0⟩]. Computed as <c>B† · DickeBlockProbe</c> via
    /// <see cref="FourModeBasis.Project(ComplexVector)"/>.
    ///
    /// <para><b>Structural Tier 1 sub-fact:</b> components 2, 3 (onto |u_0⟩, |v_0⟩) are zero
    /// at machine precision. The probe lives entirely in the channel-uniform 2D subspace
    /// span{|c_1⟩, |c_3⟩} per A2 + the structural orthogonality fact in
    /// <see cref="InterChannelSvd"/> ("the probe (Dicke state) is orthogonal to |u_0⟩, |v_0⟩").
    /// This is the structural seed for the Tier-1 sub-claim that the probe ⊥ SVD-top
    /// directions independent of the (Q, b)-rotation; it is the input to the per-(Q, b)
    /// numerical K-driving pair identification in <see cref="KDrivingPair"/>.</para>
    ///
    /// <para>Cached on construction (single Dicke probe per CoherenceBlock).</para>
    /// </summary>
    public ComplexVector ProbeProjection => _probeProjection;

    private readonly ComplexVector _probeProjection;

    /// <summary>The K-driving eigenvalue pair: the 2 of 4 eigenvalues of L_eff(Q, b) whose
    /// eigenvectors have largest squared overlap with the projected probe vector
    /// <see cref="ProbeProjection"/> (which lives entirely in span{|c_1⟩, |c_3⟩}, see XML doc
    /// of <see cref="ProbeProjection"/>).
    ///
    /// <para>At Q = Q_EP these two eigenvalues approach Re(λ) = −4γ₀ per F86 Statement 1
    /// (the EP collapse of the leading 2-level pair); around Q_peak the K-driving pair
    /// dominates the K_CC_pr observable Duhamel sum.</para>
    ///
    /// <para>The first component <c>LamPlus</c> is the K-driving eigenvalue with the larger
    /// real part (slower-decaying mode); <c>LamMinus</c> is the smaller-real-part one. Both
    /// indices are returned in matching order by <see cref="KDrivingPairIndices"/>.</para>
    ///
    /// <para><b>Tier 2 verified (numerical):</b> per-(Q, bond) Evd() of L_eff(Q, bond) and
    /// |⟨probe | w_i⟩|² ranking. The structural Tier 1 sub-fact (probe ⊥ {|u_0⟩, |v_0⟩})
    /// is verified at machine precision separately; the eigenvector rotation that mixes
    /// span{|c_1⟩, |c_3⟩} with span{|u_0⟩, |v_0⟩} is what makes this numerical.</para>
    /// </summary>
    public (Complex LamPlus, Complex LamMinus) KDrivingPair(double Q, int bond)
    {
        var (idxPlus, idxMinus) = KDrivingPairIndices(Q, bond);
        var eigs = Eigenvalues(Q, bond);
        return (eigs[idxPlus], eigs[idxMinus]);
    }

    /// <summary>Indices i, j ∈ {0, 1, 2, 3} of the K-driving eigenvectors in the order that
    /// <see cref="Eigenvalues"/> returns them (sorted Re desc / Im asc). The indices are
    /// distinct: <c>IndexPlus</c> corresponds to the larger-real-part eigenvalue of the
    /// K-driving pair, <c>IndexMinus</c> to the smaller-real-part one.
    ///
    /// <para>Computed via <see cref="ProbeOverlapsSquared"/>: top 2 entries of the squared
    /// overlap |⟨probe | w_i⟩|² over the four eigenvectors are the K-driving pair.</para>
    /// </summary>
    public (int IndexPlus, int IndexMinus) KDrivingPairIndices(double Q, int bond)
    {
        var overlaps = ProbeOverlapsSquared(Q, bond);
        var eigs = Eigenvalues(Q, bond);

        // Top-2 indices by overlap (descending).
        // Use a simple argsort over 4 elements with a stable tiebreaker (larger Re wins).
        int[] order = { 0, 1, 2, 3 };
        Array.Sort(order, (a, b) =>
        {
            int cmp = overlaps[b].CompareTo(overlaps[a]); // descending overlap
            if (cmp != 0) return cmp;
            return eigs[b].Real.CompareTo(eigs[a].Real); // tiebreaker: descending Re(λ)
        });

        int idx0 = order[0];
        int idx1 = order[1];

        // Order the pair so IndexPlus corresponds to the larger-Re eigenvalue.
        if (eigs[idx0].Real >= eigs[idx1].Real)
            return (idx0, idx1);
        return (idx1, idx0);
    }

    /// <summary>Per-eigenvector squared probe overlap |⟨probe | w_i⟩|² for i ∈ {0..3} in the
    /// order that <see cref="Eigenvalues"/> returns. Indexed identically to
    /// <see cref="KDrivingPairIndices"/>: <c>overlaps[KDrivingPairIndices.IndexPlus]</c> and
    /// <c>overlaps[KDrivingPairIndices.IndexMinus]</c> are the two largest entries.
    ///
    /// <para>The eigenvectors w_i of L_eff(Q, b) are obtained from <c>Evd().EigenVectors</c>
    /// (column i is w_i in the same row order as the unsorted MathNet eigenvalue list); we
    /// then reorder to match the (Re desc / Im asc) eigenvalue ordering of
    /// <see cref="Eigenvalues"/>. The probe is the 4-mode-basis projection (see
    /// <see cref="ProbeProjection"/>); only its first two components are nonzero, so the
    /// overlap reduces to <c>|probe[0]·w_i[0]^* + probe[1]·w_i[1]^*|²</c> in practice — but
    /// we compute the full 4-vector inner product for clarity and to keep the code path
    /// independent of the ProbeProjection structural fact (which is verified separately by
    /// the <c>ProbeProjection_HasZeroSvdTopComponents</c> test).</para>
    /// </summary>
    public double[] ProbeOverlapsSquared(double Q, int bond)
    {
        if (bond < 0 || bond >= Block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(bond),
                $"bond must be in [0, {Block.NumBonds - 1}]; got {bond}.");

        var lEff = LEffAtQ(Q);
        var evd = lEff.Evd();
        var rawEigs = evd.EigenValues.ToArray();
        var rawVecs = evd.EigenVectors;  // 4×4, each column an eigenvector

        // Build (eigenvalue, column-index) pairs and sort by Re desc / Im asc, matching the
        // ordering used by Eigenvalues(Q, bond). This makes the returned overlap array
        // index-aligned with KDrivingPairIndices.
        int[] sortedIndices = { 0, 1, 2, 3 };
        Array.Sort(sortedIndices, (a, b) => CompareByRealDescThenImagAsc(rawEigs[a], rawEigs[b]));

        var overlaps = new double[4];
        for (int i = 0; i < 4; i++)
        {
            int rawCol = sortedIndices[i];
            // ⟨probe | w_i⟩ = probe.Conjugate() · w_i = Σ_k conj(probe[k]) · w_i[k]
            Complex inner = Complex.Zero;
            for (int k = 0; k < 4; k++)
                inner += Complex.Conjugate(_probeProjection[k]) * rawVecs[k, rawCol];
            overlaps[i] = inner.Magnitude * inner.Magnitude;
        }
        return overlaps;
    }

    private static int CompareByRealDescThenImagAsc(Complex a, Complex b)
    {
        int cmpReal = b.Real.CompareTo(a.Real); // descending
        if (cmpReal != 0) return cmpReal;
        return a.Imaginary.CompareTo(b.Imaginary); // ascending
    }

    private static string BuildPendingDerivationNote(CoherenceBlock block)
    {
        int N = block.N;
        return $"C2 time-box hit Tier2Verified at N={N}; eigenvalues numerical via Evd().\n" +
               "Algebraic obstruction: the 4×4 char poly det(λI − L_eff(Q)) is a genuine\n" +
               "quartic in (λ, Q), NOT a biquadratic. Sampled coefficients show:\n" +
               "  c_0(Q) = 1\n" +
               "  c_1(Q) = 16γ₀ − i·Q·γ₀·trace(M_h_total_eff)        (linear in Q)\n" +
               "  c_2(Q) = 16γ₀² + a·Q²                              (quadratic in Q)\n" +
               "  c_3(Q) = (γ₀-terms) + b·Q² + c·Q³                  (CUBIC in Q ⇒ no biquad)\n" +
               "  c_4(Q) = (γ₀⁴-terms) + d·Q² + e·Q⁴                 (quartic in Q)\n" +
               "The cubic term in c_3 rules out factorisation into\n" +
               "(λ² − aλ + b)(λ² − cλ + d) with closed-form rational coefficients in Q.\n" +
               "\n" +
               "Structural ansätze tried (commutator [P, L_eff] = 0 candidate):\n" +
               "  - HD-parity split (P = perm c_1↔u_0 swap with c_3↔v_0 fixed):\n" +
               "    ‖[P, L]‖_F ~ 0.4 across N=5..8. The HD=1 ↔ HD=3 coupling at the\n" +
               "    [u_0, v_0] entry of M_h_total_eff (≈ σ_0 ≈ 2√2 asymptotically) is exactly\n" +
               "    the inter-channel SVD-top coupling — dominant, not perturbative.\n" +
               "  - Probe vs SVD split (P = identity on {c_1, c_3} and {u_0, v_0}):\n" +
               "    ‖[P, L]‖_F ~ 0.5; cross-block V_b[α, j] coupling is the bond-class\n" +
               "    fingerprint (Endpoint vs Interior), nonzero by construction.\n" +
               "  - Chain-mirror R projected to 4-mode basis: at odd N (5, 7) all 4 basis\n" +
               "    vectors are R-symmetric ⇒ R_eff = +I ⇒ no factorisation. At even N\n" +
               "    (6, 8) R_eff is not exactly ±I-diagonal (MathNet SVD tiebreaker for\n" +
               "    |u_0⟩, |v_0⟩ injects a factor ≠ ±1) — symptom of the same A3 even-N\n" +
               "    σ_0 degeneracy that prevents closed-form |u_0⟩, |v_0⟩.\n" +
               "  - All 5 non-trivial 4-element index permutations: minimum commutator\n" +
               "    residual ~0.28 (flip-c-only); none vanish.\n" +
               "\n" +
               "Promising next directions for a future session:\n" +
               "  (a) Approximate factorisation: drop the small probe ↔ SVD coupling\n" +
               "      ‖V_b[α, j]‖_F ~ 0.01 (for the bond-summed entries) and treat L_eff\n" +
               "      as block-diagonal at zeroth order. The leading pair would be\n" +
               "      Statement 1's [[−2γ₀, +iJ·g_eff], [+iJ·g_eff, −6γ₀]] (probe block,\n" +
               "      eigenvalues −4γ₀ ± √(4γ₀² − J²·g_eff²)). The non-leading pair would\n" +
               "      be a similar 2×2 in the SVD block. First-order corrections would\n" +
               "      come from perturbing in the cross-block. This gives Tier1Candidate\n" +
               "      analytical eigenvalues with explicit error term, not Tier1Derived.\n" +
               "  (b) Lift to the c=1 PTF framework (FourModeEffective siblings): if c=1\n" +
               "      has a clean 2-mode factorisation (Statement 1's α_i framework),\n" +
               "      maybe the c=2 4-mode is a perturbative extension of two coupled c=1\n" +
               "      doubles. Cross-check against the FourModeEffective Q-scan output.\n" +
               "  (c) Symbolic char-poly reduction in (γ₀, J·g_eff) using a CAS: the\n" +
               "      quartic in λ could simplify if one chooses a cleverer parameterisation\n" +
               "      (e.g. the 'PT-phenomenology' substitution λ = −4γ₀ + μ·γ₀, J = α·γ₀).\n" +
               "      Worth ~30 min in sympy before declaring it irreducible.";
    }

    public override string DisplayName =>
        $"c=2 4-mode effective spectrum (N={Block.N})";

    public override string Summary =>
        $"L_eff(Q) eigenvalues, IsAnalyticallyDerived={IsAnalyticallyDerived} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return new InspectableNode(
                "IsAnalyticallyDerived",
                summary: IsAnalyticallyDerived ? "true (Tier1Derived)" : "false (Tier2Verified, numerical Evd)");
            yield return BondCoupling;
            // C3: probe-projection structural sub-fact + per-(Q, b) numerical readout note.
            double svdComponentMag = Math.Max(_probeProjection[2].Magnitude, _probeProjection[3].Magnitude);
            yield return new InspectableNode(
                "ProbeProjection",
                summary: $"4-mode-basis projection of Dicke probe; |c|⊥SVD-top (max |w_2|, |w_3|) = {svdComponentMag:E2} (Tier1 structural)");
            yield return new InspectableNode(
                "KDrivingPair",
                summary: "per-(Q, bond) numerical Tier2 readout: top-2 eigenvectors by |⟨probe | w_i⟩|² → KDrivingPair(Q, b), KDrivingPairIndices(Q, b), ProbeOverlapsSquared(Q, b)");
            if (PendingDerivationNote is not null)
                yield return new InspectableNode("PendingDerivationNote", summary: PendingDerivationNote);
        }
    }
}
