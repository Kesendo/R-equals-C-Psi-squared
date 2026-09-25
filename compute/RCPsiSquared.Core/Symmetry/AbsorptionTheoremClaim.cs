using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Absorption Theorem (Tier 1 derived, <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c>),
/// with its two mathematical objects kept distinct; stated first at uniform γ₀, then per site. Under
/// uniform Z-dephasing γ₀ the dissipator is diagonal on a computational-basis pair |A⟩⟨B| and
/// contributes the cell cost 2γ₀·n_diff(A,B); for an initially isolated cell this is its initial
/// decay slope, and it is a decay rate only where |A⟩⟨B| is itself an eigenmode (|A⟩ and |B⟩
/// eigenstates of H). A Liouvillian eigenmode is generally a superposition of such cells. Its decay
/// rate follows from the Hermitian part of L and is −Re λ = 2γ₀·⟨n_XY⟩_v, where the expectation may
/// be non-integer. Thus the coefficient 2 is exact, but the interacting eigenvalue spectrum is not
/// thereby quantized in steps of 2γ₀.
///
/// <para><b>The carrier is a vector (Theorem 2, the per-site law).</b> Under site-dependent
/// Z-dephasing {γ_l} the rate of every eigenmode v_k is the γ profile paired with the mode's
/// per-site light profile:</para>
///
/// <code>
///   −Re(λ_k) = 2 · Σ_l γ_l · ⟨Δ_l⟩_k,   ⟨Δ_l⟩_k = ⟨v_k|N_l|v_k⟩ / ‖v_k‖² ∈ [0, 1],
///   N_l = (I − Z_l⊗Z_l)/2
/// </code>
///
/// <para>N_l is the site-l disagreement projector: on a coherence |A⟩⟨B| it reads the sharp bit
/// Δ_l = [A_l ≠ B_l], whether bra and ket disagree at site l, and the Hamiltonian, by rotating the
/// eigenmodes, turns that bit into the expectation ⟨Δ_l⟩. The diagonal (bra = ket at every site)
/// never pays. The uniform statement is the all-sites-equal projection, Σ_l ⟨Δ_l⟩ = ⟨n_XY⟩; a
/// uniform γ₀ leaves a site-permutation degeneracy (the N + 1 rungs of the bare dissipator), a
/// profile lifts it where the site rates differ, fully when its subset sums are distinct, and a
/// single-site coherence then reads its own γ_l (its cell cost 2γ_l, a decay rate wherever H leaves
/// both of its labels eigenstates). Read with each site as an independent dephasing channel
/// it is the per-channel decoherence-rate law of a heterogeneous substrate
/// (<c>CarrierVectorPortfolio</c> in Diagnostics). Computed here by <see cref="PerSiteLightProfile"/>
/// and <see cref="PerSiteEigenmodeDecayRate"/>; the same per-site dissipator diagonal is Step 1 of
/// <c>docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md</c>, where it feeds the F1 palindrome's independence
/// of the γ pattern.</para>
///
/// <para><b>H-independence, for any Hermitian H.</b> In the row-major vec representation
/// L_H = −i(H⊗I − I⊗Hᵀ) is anti-Hermitian for every Hermitian H, since Hermiticity alone gives
/// Hᵀ = H* (no reality assumption). So v†L_H v is purely imaginary, the Hamiltonian drops out of
/// Re λ through the Rayleigh quotient, and Herm(L) = (L + L†)/2 is the dephasing dissipator alone:
/// complex Hermitian H, Dzyaloshinskii-Moriya, transverse or Y terms, any graph (proof Step 1;
/// <c>simulations/popcount_identity_h_independence.py</c> reads the difference as exactly 0 against
/// a random complex Hermitian H). The Hamiltonian decides where a mode's light sits, never what a
/// unit of it costs.</para>
///
/// <para><b>Sitting on a rung is not a weight certificate.</b> At uniform γ₀, mixing across weight
/// sectors usually moves a mode off the integer rungs (F33's 8γ₀/3 and 10γ₀/3 at N = 3,
/// strong-coupling limits of ⟨n_XY⟩ = 4/3 and 5/3), but a mixture whose light averages to an
/// integer lands exactly on one. The {0,2}-coherence is the standing example: on the uniform XY
/// chain at N = 2 and 3 and, with the ZZ term (Heisenberg), at N = 2 only, above its own Q*, its
/// n_diff histogram is {0: ½, 2: ½}, so ⟨n_diff⟩ = 1 and it sits on −2γ₀ beside the pure
/// distance-1 modes (<c>docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md</c> §4,
/// <see cref="CoherenceHorizonClaim"/>). Where the block no longer closes (on the XY chain at N = 4
/// and 5, as <see cref="CoherenceHorizonClaim"/> reads it), its share w2 exceeds ½ and it sits below
/// that line by 2γ₀(2w2 − 1).
/// <see cref="F86.LEffMirrorAxisClaim"/> carries the same shape one rung up,
/// ½·1 + ½·3 = 2 at −4γ₀. A rung's occupants are therefore not all of pure weight.</para>
///
/// <para><b>The two ends.</b> 0 ≤ ⟨Δ_l⟩ ≤ 1 bounds every rate by 2·Σ_l γ_l = 2σ, 2γ₀N at uniform
/// γ₀. The bottom end is reached by every Hamiltonian, because the identity is a fixed point of every
/// such L (−i[H, I] = 0, and Z-dephasing is unital). The top end is reached within the
/// number-conserving XY/Heisenberg family (the XOR drain, proof §4.1), so there the ladder is spanned
/// end to end; for a generic Hermitian H it need not be: at N = 3, γ₀ = 0.05 a generic real symmetric
/// H reaches no mode at 2Nγ₀, its fastest sitting at 4.20γ₀ on the draw the proof's regime script
/// fixes (proof §4 and §6, <c>simulations/absorption_ladder_regimes.py</c>).</para>
///
/// <para><b>Four more readings of the same Rayleigh quotient</b> (proof §2 Extensions and §4.7).
/// (i) Two-sided: a left eigenvector w (w†L = λw†) carries the same weighted light as the right
/// one, since w†Lw = λ‖w‖² decomposes the same way; biorthogonal bookkeeping cannot disagree with
/// itself about absorption. (ii) The projector form: for a degenerate cluster the basis-free
/// reading is light_l(V) = Tr(Π_V·Δ_l)/dim V through the orthogonal projector Π_V onto the
/// cluster's invariant subspace (not the biorthogonal projector's own diagonal, whose entries can
/// leave [0, 1]), exact by a block-triangular trace argument; <c>SlowLightDistribution</c> in
/// Diagnostics computes it. (iii) The dephase-letter rotation: X-dephasing reads n_YZ and
/// Y-dephasing n_XZ, the letters that anticommute with the dephasing letter. (iv) The recentred
/// face: at uniform γ₀, L_D = γ₀·(Q − N·I) with Q = Σ_l Z_l⊗Z_l, so M = L + γ₀N·I = L_H + γ₀Q moves
/// the ladder's midpoint ⟨n_XY⟩ = N/2 to zero, and one diagonal enters three arguments: the
/// absorption ladder (a Rayleigh quotient), the palindrome (F1 partners carry complementary light,
/// ⟨n_XY⟩_s + ⟨n_XY⟩_f = N, gated per mode by <c>F8PartnerLightComplementarityTests</c>) and the F87
/// windowed converse (a power-sum expansion). The Pauli reading is an identity for every H exactly
/// when Herm(L_D) is diagonal in the Pauli basis, because an anti-Hermitian part of L_D drops out of Re λ as L_H
/// does. Jumps each proportional to a Pauli string are one sufficient way to get that, not a necessary
/// one: the correlated jumps √a(Z₀ + iZ₁), √b(Z₀ − iZ₁), a ≠ b, come from no Pauli-string jump set, yet
/// their Hermitian part is local Z-dephasing at rate a + b. Where Herm(L_D) is not diagonal in the Pauli
/// basis (collective dephasing through the one jump Σ_k Z_k, amplitude damping) the Rayleigh identity
/// still holds and reads the rate in the basis that diagonalizes Herm(L_D) (proof §2, "Where the
/// boundary actually runs"; each case checked in <c>simulations/absorption_ladder_regimes.py</c>).</para>
///
/// <para><b>The coefficient 2.</b> It comes straight from Step 2: Z anticommutes with X and Y, so
/// ZXZ − X = −2X, and no dimension enters. The repository also carries it as the Pi2 dyadic
/// ladder's a_0 (<see cref="Pi2DyadicLadderClaim.Term"/>(0), the root d of d² − 2d = 0), shared
/// with F1's TwoFactor, F50's DecayRateFactor and F66's UpperPoleCoefficient. That anchor records
/// that these coefficients all equal 2; whether they equal 2 for the same reason the proof leaves
/// open (§6). 2γ₀ is the rung of the bare dissipator at uniform γ₀, the proof's "absorption
/// quantum": a cell cost, not a spacing of the interacting rate spectrum.</para>
///
/// <para><b>Direct typed children.</b> In the default registry
/// (<c>KnowledgeRegistryFactory.BuildDefault</c>) thirty-four claims take this claim as a direct
/// typed parent. Fourteen of them name themselves by an F-number of
/// <c>docs/ANALYTICAL_FORMULAS.md</c>, and those numbers are <see cref="FNumberedDirectChildren"/>:
/// F25, F33, F50, F55, F64–F68, F74 and F89 (the Pi2 inheritances), F135
/// (<see cref="RecordParityLawClaim"/>), F153 (<c>PinnedBlockFloorClaim</c>) and F154
/// (<c>CompressedDensityLocusClaim</c>). The rule is mechanical and reads the child alone: a child
/// counts when its class name or its Name opens with the number, or when the title of its Name, the
/// text before the first colon, closes with the number in parentheses. A number cited in passing does
/// not count (<see cref="ClockHandLadderClaim"/>'s Name mentions the F2b band edge after its title),
/// and neither does an entry that names a child as its typed claim without the child carrying the
/// number: the F122 entry names <see cref="StructuralCeilingClaim"/>, the F86 table
/// <c>TPeakLaw</c> and <see cref="F86.LEffMirrorAxisClaim"/>, and none of them is in the list.
/// <c>AbsorptionTheoremClaimRegistrationTests</c> applies the rule to every direct child the default
/// registry returns, checks each number's entry heading in the formula registry, and asserts equality
/// with the list and the two counts. Transitive descendants
/// (F3 through F50, the F2b inheritance claim through F65) are not direct children;
/// <c>knowledge descendants AbsorptionTheoremClaim</c> walks every child, numbered or not, and
/// <c>compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md</c> carries the per-formula one-liners. The
/// proof keeps its own list, of the results the theorem unifies, and the two lists differ in both
/// directions.</para>
///
/// <para><b>Verification.</b> Per mode at N = 2..5 on the Heisenberg chain at J = 1, γ₀ = 0.05:
/// 13 + 60 + 251 + 1,018 = 1,342 active modes, the ratio α/(2γ₀⟨n_XY⟩) at mean 1.000000 with CV
/// 0.0000 to the printed digits (<c>simulations/absorption_theorem_discovery.py</c> Step 6,
/// <c>simulations/results/absorption_theorem_discovery.txt</c>; the 18 kernel modes, where the
/// ratio is 0/0, are the ⟨n_XY⟩ = 0 end of the law rather than exceptions to it). Over all 1,024
/// modes of the N = 5 Heisenberg chain the largest deviation is 2.4·10⁻¹⁴
/// (<c>simulations/path_d_bell_pair_absorption.py</c>, recorded in
/// <c>experiments/COCKPIT_SCALING.md</c> §11). The per-site law on the N = 3 Heisenberg chain with
/// a non-uniform γ: <c>simulations/absorption_gamma_vector.py</c>. In C#,
/// <c>AbsorptionTheoremClaimTests</c> checks Herm(L) = diag(−2·Σ_l γ_l·Δ_l) entry by entry and
/// exactly for a random complex Hermitian H, and the per-site law on every right and every left
/// eigenvector at N = 2 and 3 with a non-uniform γ, the deviation held in one band around
/// eps·‖L‖_F, and within a factor 16 of itself, while the coupling runs over six decades, with
/// mutations that must fail (the uniform-γ formula on a non-uniform profile, the reversed site
/// order).</para>
///
/// <para><b>Hardware, and what it does not reach.</b> IBM Torino Q52 single-qubit tomography
/// (Confirmations <c>absorption_theorem_ratio_torino</c>, <c>experiments/IBM_ABSORPTION_THEOREM.md</c>)
/// reads excess/(2γ*) = 1.03 on the free-evolution T2* baseline. The proof's fence (§3): at N = 1
/// with n_XY = 1, γ* is extracted from the same coherence envelope the ratio reads, so the ratio is
/// 1 by construction up to the T1 each side subtracts (241.4 μs from the population fit in the
/// excess, 221.2 μs from the calibration in 2γ*), and the 3 % is that difference: the consistency of
/// two fits to one decay, not the ladder. What the theorem predicts beyond one site is additivity: a
/// coherence that differs at sites a and b decays at the sum of its two sites' rates on the same
/// device, 2(γ_a + γ_b), which is twice a single-site rate only when the site rates are equal
/// (uniform γ₀). That needs N ≥ 2, and no registered rate ratio with an error bar reads it. The
/// price-pair campaign (ibm_marrakesh, 2026-07-04, four runs) pre-registered the additivity as its P1.
/// Its registered verdict (<c>price_pair_locality_marrakesh_july2026</c>, in both Confirmations
/// registries) is that the dephasing covariances were local in the clean run-3 session, the premise the
/// additivity follows from under the local-channel model; the direct rate test P1 was masked on the
/// product-state singles by a coherent nearest-neighbour ZZ
/// (<c>experiments/PRICE_PAIR_HARDWARE_PREDICTION.md</c>). So the premise has a registered reading and
/// the rate law itself does not. The staircase null test (ibm_marrakesh, 2026-07-26, two flights)
/// tested the other half of the per-site law, that a site where bra and ket agree adds no dephasing
/// cost, together with Δ(11) = Δ(10) + Δ(01), the additivity of two spectators' relaxation
/// contributions to a rate difference, not a coherence that differs at two sites. Both flights
/// returned B-BLOCK-INVALID (its in-situ T1 reference failed bracket consistency in the first; in the
/// second the residual condition on the far spectator failed, which the audit traces to a mid-batch
/// T1 excursion), so by its frozen truth table that claim is neither confirmed
/// nor falsified, the spectator additivity holding both times
/// (<c>experiments/STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md</c>). A draft pre-registration for a
/// clean reading of the rate law exists, unfrozen and with nothing in it registered
/// (<c>experiments/ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md</c>).</para>
///
/// <para>Anchors: the proof (Theorems 1 and 2, the Extensions, §4.7),
/// <c>docs/ANALYTICAL_FORMULAS.md</c> (entry AT and the entries in <see cref="FNumberedDirectChildren"/>),
/// <c>experiments/ABSORPTION_THEOREM_DISCOVERY.md</c>, <c>experiments/IBM_ABSORPTION_THEOREM.md</c>,
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> (the cell costs and the Hamming-complement sum)
/// and <see cref="Pi2DyadicLadderClaim"/>.</para></summary>
public sealed class AbsorptionTheoremClaim : Claim
{
    /// <summary>The F-numbers of <c>docs/ANALYTICAL_FORMULAS.md</c> by which the direct typed children
    /// of this claim in the default registry name themselves: a child's class name or Name opens with
    /// the number, or the title of its Name (the text before the first colon) closes with it in
    /// parentheses. Fourteen children, fourteen numbers. Pinned by equality in
    /// <c>AbsorptionTheoremClaimRegistrationTests</c>, which applies this rule to every direct child
    /// <c>KnowledgeRegistryFactory.BuildDefault</c> returns.</summary>
    public static IReadOnlyList<string> FNumberedDirectChildren { get; } = new[]
    {
        "F25", "F33", "F50", "F55", "F64", "F65", "F66", "F67", "F68", "F74", "F89",
        "F135", "F153", "F154",
    };

    public Pi2DyadicLadderClaim Ladder { get; }

    /// <summary>The exact coefficient 2 multiplying both the basis-pair dissipator
    /// cost and the eigenmode expectation value.</summary>
    public double DissipatorCoefficient => Ladder.Term(0);

    /// <summary>Dissipator cost of one bra-ket disagreement at uniform γ₀: <c>2γ₀</c>. This is a
    /// basis-cell quantity, not the smallest nonzero eigenmode decay-rate step.</summary>
    public double SingleDisagreementCellCost(double gammaZero)
    {
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero;
    }

    /// <summary>Diagonal dissipator cost / isolated-cell initial decay slope of |A⟩⟨B| at uniform
    /// γ₀: <c>2γ₀·n_diff(A,B)</c> (per site, <see cref="PerSiteEigenmodeDecayRate"/> on the sharp
    /// profile of |A⟩⟨B|).</summary>
    public double BasisPairDissipatorCost(int nDiff, double gammaZero)
    {
        if (nDiff < 0)
            throw new ArgumentOutOfRangeException(nameof(nDiff), nDiff, "n_diff must be >= 0.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * nDiff;
    }

    /// <summary>Decay rate of a Liouvillian eigenmode under uniform γ₀ from its normalized light
    /// expectation: <c>−Re λ = 2γ₀·⟨n_XY⟩_v</c>. The expectation is continuous. This is the
    /// all-sites-equal projection of <see cref="PerSiteEigenmodeDecayRate"/>.</summary>
    public double EigenmodeDecayRate(double averageNXy, double gammaZero)
    {
        if (!double.IsFinite(averageNXy) || averageNXy < 0)
            throw new ArgumentOutOfRangeException(nameof(averageNXy), averageNXy,
                "average n_XY must be finite and >= 0.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * averageNXy;
    }

    /// <summary>Theorem 2, the per-site law: the decay rate of a Liouvillian eigenmode under
    /// site-dependent Z-dephasing, <c>−Re λ = 2·Σ_l γ_l·⟨Δ_l⟩_v</c>, from the γ profile and the
    /// mode's per-site light profile ⟨Δ_l⟩_v ∈ [0, 1] (<see cref="PerSiteLightProfile"/>). Exact for
    /// every right and every left eigenvector and any Hermitian H, because Herm(L) is the dephasing
    /// dissipator alone. On a sharp profile (a single coherence |A⟩⟨B|, bits Δ_l ∈ {0, 1}) it is the
    /// cell cost 2·Σ_{l: A_l ≠ B_l} γ_l, a decay rate where that coherence is an eigenmode; at
    /// uniform γ₀ it reduces to <see cref="EigenmodeDecayRate"/> with ⟨n_XY⟩ = Σ_l ⟨Δ_l⟩. The sites
    /// are summed in ascending order l = 0..N−1, the order the Liouvillian builder uses, so on a sharp
    /// profile this reproduces the builder's diagonal bit for bit.</summary>
    public double PerSiteEigenmodeDecayRate(IReadOnlyList<double> gammaPerSite, IReadOnlyList<double> lightPerSite)
    {
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (lightPerSite is null) throw new ArgumentNullException(nameof(lightPerSite));
        if (gammaPerSite.Count == 0)
            throw new ArgumentException("the γ profile needs at least one site.", nameof(gammaPerSite));
        if (lightPerSite.Count != gammaPerSite.Count)
            throw new ArgumentException(
                $"light profile has {lightPerSite.Count} sites, γ profile has {gammaPerSite.Count}.",
                nameof(lightPerSite));
        double weighted = 0.0;
        for (int l = 0; l < gammaPerSite.Count; l++)
        {
            double gamma = gammaPerSite[l];
            if (!double.IsFinite(gamma) || gamma < 0.0)
                throw new ArgumentOutOfRangeException(nameof(gammaPerSite), gamma,
                    $"γ_{l} must be finite and >= 0.");
            double light = lightPerSite[l];
            if (!double.IsFinite(light) || light < 0.0 || light > 1.0)
                throw new ArgumentOutOfRangeException(nameof(lightPerSite), light,
                    $"⟨Δ_{l}⟩ must be finite and lie in [0, 1].");
            weighted += gamma * light;
        }
        return DissipatorCoefficient * weighted;
    }

    /// <summary>The per-site light profile ⟨Δ_l⟩_v = ⟨v|N_l|v⟩/‖v‖², N_l = (I − Z_l⊗Z_l)/2, of a
    /// Liouville-space vector of length 4^N. Layout: the row-major vec of
    /// <see cref="Lindblad.PauliDephasingDissipator"/>, index x = a·d + b for the coherence |a⟩⟨b|,
    /// d = 2^N; site l is the l-th Kronecker factor from the left (<see cref="Pauli.PauliString.SiteOp"/>),
    /// i.e. bit N−1−l of a basis index. Only the bits of a ⊕ b enter, so the column-stacked vec of the
    /// same operator gives the same profile. N_l is diagonal, so the reading is the weighted count
    /// Σ_x Δ_l(x)·|v_x|² / Σ_x |v_x|² with Δ_l(x) the bra-ket disagreement bit at site l; each entry
    /// lies in [0, 1] and the entries sum to ⟨n_XY⟩_v. On a single coherence the profile is its sharp
    /// bit pattern.</summary>
    public static double[] PerSiteLightProfile(ComplexVector vec)
    {
        if (vec is null) throw new ArgumentNullException(nameof(vec));
        int dim = vec.Count;
        int n = 0;
        while (n < 15 && (1 << (2 * n)) < dim) n++;
        if (n < 1 || (1 << (2 * n)) != dim)
            throw new ArgumentException($"vector length {dim} is not 4^N for some N >= 1.", nameof(vec));
        int d = 1 << n;

        var light = new double[n];
        double norm2 = 0.0;
        for (int x = 0; x < dim; x++)
        {
            Complex c = vec[x];
            double w = c.Real * c.Real + c.Imaginary * c.Imaginary;
            if (w == 0.0) continue;
            norm2 += w;
            int diff = (x / d) ^ (x % d);
            for (int l = 0; l < n; l++)
                if (((diff >> (n - 1 - l)) & 1) != 0) light[l] += w;
        }
        if (!double.IsFinite(norm2) || norm2 <= 0.0)
            throw new ArgumentException("the vector must be finite and nonzero.", nameof(vec));
        for (int l = 0; l < n; l++) light[l] /= norm2;
        return light;
    }

    /// <summary>Upper bound <c>2γ₀N</c> at uniform γ₀, following from 0 ≤ ⟨n_XY⟩_v ≤ N (per site,
    /// 2·Σ_l γ_l). This is a ceiling, not evidence that eigenvalues fill a quantized grid; the
    /// number-conserving XY/Heisenberg family reaches it, a generic Hermitian H need not.</summary>
    public double EigenmodeDecayRateCeiling(int n, double gammaZero)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be >= 1.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * n;
    }

    /// <summary>Recover the continuous eigenmode light expectation from its decay rate at uniform
    /// γ₀.</summary>
    public double AverageNXyFromEigenmodeDecayRate(double rate, double gammaZero)
    {
        if (!double.IsFinite(rate) || rate < 0)
            throw new ArgumentOutOfRangeException(nameof(rate), rate, "rate must be finite and >= 0.");
        if (!double.IsFinite(gammaZero) || gammaZero <= 0)
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "gamma must be finite and > 0.");
        return rate / (DissipatorCoefficient * gammaZero);
    }

    /// <summary>Basis-pair cell-cost sum under the bra Hamming complement at uniform γ₀:
    /// <c>2γ₀·n_diff + 2γ₀·(N − n_diff) = 2γ₀N</c>.</summary>
    public double HammingComplementCellCostSum(int blockSize, double gammaZero)
    {
        if (blockSize < 1)
            throw new ArgumentOutOfRangeException(nameof(blockSize), blockSize, "block size must be >= 1.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * blockSize;
    }

    /// <summary>True iff the Pi2 anchor's a_0 is the literal 2.0 of the proof, compared exactly:
    /// <c>Term(0)</c> = 2^(1−0) is exactly 2.0 in floating point.</summary>
    public bool DissipatorCoefficientMatchesLiteral() => DissipatorCoefficient == 2.0;

    public AbsorptionTheoremClaim(Pi2DyadicLadderClaim ladder)
        : base("Absorption Theorem: eigenmode decay −Re(λ) = 2·Σ_l γ_l·⟨Δ_l⟩_v; at uniform γ₀ the basis-pair cell cost is 2γ₀·n_diff and −Re(λ) = 2γ₀·⟨n_XY⟩_v; coefficient 2 = a_0",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md + " +
               "docs/ANALYTICAL_FORMULAS.md (AT and the entries in FNumberedDirectChildren) + " +
               "experiments/ABSORPTION_THEOREM_DISCOVERY.md + " +
               "experiments/IBM_ABSORPTION_THEOREM.md + " +
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "Absorption Theorem: basis-pair cell cost and eigenmode light expectation, site by site";

    public override string Summary =>
        $"Under uniform γ₀ the basis-pair dissipator diagonal / isolated-cell initial slope is 2γ₀·n_diff and " +
        $"the eigenmode decay is −Re(λ) = 2γ₀·⟨n_XY⟩_v with generally non-integer expectation; per site, " +
        $"−Re(λ) = 2·Σ_l γ_l·⟨Δ_l⟩_v for any Hermitian H; coefficient 2 = a_0; " +
        $"2γ₀N (2·Σ_l γ_l for a profile) is a ceiling, not an eigenvalue step ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            var inv = CultureInfo.InvariantCulture;
            yield return InspectableNode.RealScalar("dissipator coefficient (= a_0 = 2)", DissipatorCoefficient);
            yield return new InspectableNode("basis-pair cell reading",
                summary: "Under uniform γ₀, D[|A⟩⟨B|] contributes −2γ₀·n_diff(A,B)·|A⟩⟨B| (per site, −2·Σ_l γ_l·Δ_l). It is the diagonal cell cost and isolated-cell initial slope, a decay rate only where |A⟩⟨B| is an eigenmode.");
            yield return new InspectableNode("eigenmode reading",
                summary: "Under uniform γ₀, a normalized Liouvillian eigenmode v has −Re λ = 2γ₀·⟨n_XY⟩_v from the Hermitian-part Rayleigh quotient; Hamiltonian mixing permits non-integer expectations, and an integer one is no weight certificate (the {0,2}-coherence sits on −2γ₀ with histogram {0: ½, 2: ½} on the uniform XY chain at N = 2, 3 and on the Heisenberg chain at N = 2, above its Q*).");
            yield return PerSiteLawNode();
            yield return new InspectableNode("H-independence (any Hermitian H)",
                summary: "L_H = −i(H⊗I − I⊗Hᵀ) is anti-Hermitian for every Hermitian H (Hᵀ = H*, no reality assumption), so Herm(L) is the dephasing dissipator alone and H drops out of Re λ: complex H, Dzyaloshinskii-Moriya, transverse or Y terms, any graph (proof Step 1).");
            yield return new InspectableNode("Hamming-complement cell-cost sum",
                summary: "At uniform γ₀, n_diff maps to N − n_diff, so the two basis-cell costs sum to 2γ₀N; " +
                         $"at N = 3, γ₀ = 1 the sum is {HammingComplementCellCostSum(3, 1.0).ToString("G6", inv)}.");
            yield return new InspectableNode("Pi2 anchor drift check",
                summary: $"DissipatorCoefficientMatchesLiteral = {DissipatorCoefficientMatchesLiteral()} " +
                         $"(a_0 = {DissipatorCoefficient.ToString("G6", inv)}); the anchor records that the coefficients equal 2, " +
                         "the proof (§6) leaves open whether for one reason.");
            yield return new InspectableNode("F-numbered direct typed children",
                summary: string.Join(", ", FNumberedDirectChildren) +
                         ": the direct children that name themselves by the number (class name or Name opening " +
                         "with it, or the Name's title closing with it in parentheses); the others are not listed, " +
                         "and knowledge descendants AbsorptionTheoremClaim walks the whole set.");
            yield return new InspectableNode("numerical verification",
                summary: "1,342 active modes of the N = 2..5 Heisenberg chain (J = 1, γ₀ = 0.05), ratio mean 1.000000, CV 0.0000 to the printed digits (simulations/absorption_theorem_discovery.py); N = 5, all 1,024 modes, largest deviation 2.4·10⁻¹⁴ (simulations/path_d_bell_pair_absorption.py, experiments/COCKPIT_SCALING.md §11); per-site law on right and left eigenvectors with a complex H and non-uniform γ at N = 2, 3 gated in AbsorptionTheoremClaimTests.");
            yield return new InspectableNode("hardware reading and its fence",
                summary: "IBM Torino Q52 single qubit: excess/(2γ*) = 1.03 on the T2* baseline (experiments/IBM_ABSORPTION_THEOREM.md). " +
                         "At N = 1 γ* comes from the same envelope, so the ratio is 1 by construction up to the T1 each side " +
                         "subtracts (241.4 μs fitted, 221.2 μs calibrated): two fits to one decay, not the ladder. Beyond one " +
                         "site the prediction is additivity, a coherence differing at sites a and b decaying at 2(γ_a + γ_b), " +
                         "twice a single-site rate only at uniform γ₀; it needs N ≥ 2 and no registered ratio with an error bar " +
                         "reads it. The price-pair campaign (ibm_marrakesh 2026-07-04) registered its premise, local dephasing " +
                         "covariances in the clean run-3 session, while its rate test P1 was masked by coherent nearest-neighbour ZZ. " +
                         "The staircase null test (2026-07-26, two flights, B-BLOCK-INVALID both times) tested the other half, " +
                         "agreeing sites pay nothing, and the additivity of two spectators' relaxation, not a two-site coherence. " +
                         "Draft pre-registration, unfrozen: experiments/ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md.");
        }
    }

    /// <summary>The live Theorem 2 node: the per-site law stated, and evaluated on two sharp profiles
    /// at γ = (0.05, 0.10): |00⟩⟨11|, the proof's rung-2 example (an exact N = 2 eigenmode of any
    /// number-conserving H, the XY and Heisenberg chains included, so its cell cost is its decay rate), and |01⟩⟨00| (a cell cost, and not an eigenmode of the XY or Heisenberg
    /// chain at J ≠ 0, where |01⟩ is not an eigenstate).</summary>
    private InspectableNode PerSiteLawNode()
    {
        var gamma = new[] { 0.05, 0.10 };
        // |00⟩⟨11| is flat index 0·4 + 3; |01⟩⟨00| is flat index 1·4 + 0.
        var pair = PerSiteLightProfile(BasisVector(16, 3));
        var single = PerSiteLightProfile(BasisVector(16, 4));
        double pairRate = PerSiteEigenmodeDecayRate(gamma, pair);
        double singleRate = PerSiteEigenmodeDecayRate(gamma, single);
        var inv = CultureInfo.InvariantCulture;
        return new InspectableNode("the carrier is a vector (Theorem 2, the per-site law)",
            summary: "−Re λ_k = 2·Σ_l γ_l·⟨Δ_l⟩_k with ⟨Δ_l⟩_k = ⟨v_k|N_l|v_k⟩/‖v_k‖², N_l = (I − Z_l⊗Z_l)/2, " +
                     "for every right and left eigenvector and any Hermitian H; per-channel when a site is " +
                     "read as a dephasing channel. At γ = (0.05, 0.10) on sites (0, 1): |00⟩⟨11| has profile " +
                     $"({pair[0].ToString("G3", inv)}, {pair[1].ToString("G3", inv)}) and costs " +
                     $"{pairRate.ToString("G6", inv)}, both sites' prices, and as an exact N = 2 eigenmode of " +
                     "the XY and Heisenberg chains it decays at that rate; |01⟩⟨00| has profile " +
                     $"({single[0].ToString("G3", inv)}, {single[1].ToString("G3", inv)}) and cell cost " +
                     $"{singleRate.ToString("G6", inv)}, site 1's price alone, a decay rate only where H leaves |01⟩ " +
                     "and |00⟩ eigenstates, not on the XY or Heisenberg chain at J ≠ 0.",
            provenance: NodeProvenance.Live);
    }

    private static ComplexVector BasisVector(int dim, int index)
    {
        var v = ComplexVector.Build.Dense(dim);
        v[index] = Complex.One;
        return v;
    }

    private static void ValidateGamma(double gammaZero)
    {
        if (!double.IsFinite(gammaZero) || gammaZero < 0)
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero,
                "gamma must be finite and >= 0.");
    }
}
