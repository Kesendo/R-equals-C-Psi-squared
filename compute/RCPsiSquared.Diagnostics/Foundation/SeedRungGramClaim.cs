using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F143, the seed rung in closed form: at the seed of the frozen band the surviving operator is
/// the Gram matrix G = WᵀW of squared mode amplitudes (= W² on the uniform chain, where W is symmetric),
/// it equals (1/M)·(𝟏𝟏ᵀ + (I + R)/2) with M := N + 1 and R the chiral mode involution a ↦ M − a, and its
/// whole spectrum follows: 0 with multiplicity ⌊N/2⌋, 1/M with multiplicity ⌈N/2⌉ − 1, and 1 simple.
/// Recomputed live by <see cref="SeedRungGramWitness"/> (<c>inspect --root seedrung</c>).
///
/// <para><b>Why this earns a carrier while the rest of its chain does not (yet).</b> The chain
/// F143 → F146 was untyped as a whole, and the arc <c>f_registry_meets_the_typed_layer</c> asked which of
/// the four numbers earns its own claim. This one does, on the arc's own admission test — an entry must
/// let a reader skip a computation, and a claim's evidence must land as a LIVE witness. F143 clears both:
/// the repo REUSES the 1/M gap as the moat a numeric rank read needs (it closes linearly in N, so the moat
/// does not vanish), and the whole statement is checkable without an eigensolver, in integer arithmetic.
/// F144 and F146 do not clear the second: F144's evidence is a bound on a COMPRESSION of 𝒦 to a subspace
/// that has to be built first, and F146's is a rank count carrying a number-theoretic hypothesis on M
/// (chiral-only, proved for M prime, M = 2p and M = 2^a), so neither reduces to an entry-wise read.
/// Not because a rank is out of reach: <see cref="SeedRungGramWitness.CertifiedRank"/> is itself an exact
/// GF(p) rank, and <see cref="FrozenDivisorClaim"/> carries its multiplicity by one. They stay measured in
/// <c>simulations/eta_ceiling_reduction.py</c> (F144, blocks V10 to V12) and
/// <c>simulations/scalar_count.py</c> (F145 and F146, blocks W1 to W7). F145 is the near miss and is
/// recorded as such in the arc: its ladder relations ARE exact, but they live in the cell basis of mode
/// SETS, which the typed layer has no primitive for; building one is a bigger move than typing this
/// matrix, and it is also why F145 cannot simply become a member on <c>SidewaysSpinLadderClaim</c>, which
/// carries the F142 ladder as CG coefficients and not as an action on cells.</para>
///
/// <para><b>Why there is no tolerance anywhere in the exact half.</b> The entries of G are 1/M and 3/(2M)
/// and 2/M, none of them representable in binary floating point, but 2M·G is an INTEGER matrix, and this
/// claim carries that scaling rather than G itself (<see cref="ScaledGram"/>). The proof's four entry cases
/// collapse to one line there: Ĝ_{ac} = 2 + [a = c] + [a + c = M], which is 2·𝟏𝟏ᵀ + I + R written out. The
/// three eigenrelations are then integer identities and the witness compares them to 0 exactly, with the
/// eigenvectors integer as well: the kernel's chiral differences e_k − e_{M−k}, the middle band's orbit
/// differences, and the all-ones vector. See <see cref="TheIntegerScaling"/>.</para>
///
/// <para><b>The one measured half, and its error model.</b> That the DST-I modes
/// v_k(l) = √(2/M)·sin(πkl/M) produce that matrix in the first place is where floating point enters, and
/// sines are not exact. So the witness rebuilds W from doubles, forms WᵀW, and reads the deviation from
/// the closed form as a RESIDUAL against a LAW rather than a threshold: the ratio residual/(ε·M) is FLAT
/// across N = 2..40, reaching 2.09 and not growing. Two mechanisms feed it, the reduction of the sine
/// arguments (which run up to ≈πN) and the summation of N terms, and which of them dominates is not
/// settled here; what is asserted is the flatness, because a residual that merely passed a fixed tolerance
/// at one N would say nothing about whether the identity is the reason it passed.</para>
///
/// <para><b>Two fences, and the first is narrower than its usual phrasing.</b> "The closed form wants
/// UNIFORM hopping" is true and fences off less than it sounds. Under any bond profile WITH NO ZERO BOND
/// K h K = −h still pairs the modes at ±ε with amplitudes differing by the sublattice sign, so the SQUARED
/// amplitudes are equal, W R = W, hence G R = G, hence every R-odd vector is still annihilated: the KERNEL,
/// its identity as the R-odd sector and its dimension ⌊N/2⌋ all survive disorder, and λ_max = 1 stays, W
/// being doubly stochastic for the squared amplitudes of any orthonormal basis. The nonzero-bond
/// hypothesis is load-bearing and is the one <see cref="ChiralKClaim"/> already carries as
/// non-degeneracy: ±ε pairing fixes the eigenVALUES, and reading v_{M−k} ∝ K v_k off it needs the
/// eigenVECTOR to be determined by its eigenvalue. With every bond nonzero h is a Jacobi matrix and its
/// spectrum is simple; with a zero bond the chain falls apart, the spectrum degenerates and the step is
/// not merely false but ILL-POSED, W ceasing to be a function of h at all (two eigensolvers may return
/// different orthonormal bases of one eigenspace and hence different W). At N = 6 with bonds [1,1,0,1,1]
/// the residuals duly come out O(1); at N = 4 with [1,0,1] one basis reads them as exactly 0 and another
/// does not, which is the same statement seen twice. What is NOT
/// claimed, because double stochasticity does not give it, is that the 1 stays SIMPLE; simplicity is an
/// irreducibility statement and is exactly what the disconnected chain loses. What uniform hopping buys is
/// the middle band's common VALUE 1/M, which disorder splits, and with it the entry-wise closed form.
/// <see cref="SeedRungGramWitness.ReadHopping"/> reads them separately for that reason, and the
/// distinction is load-bearing downstream: what F144 and F146 consume is the COUNT, which does not need
/// the fence. One clause against over-reading that: the count survives, but its MOAT is no longer
/// GUARANTEED. Under disorder the first non-kernel eigenvalue is not bounded below by 1/M any more; it
/// falls under it at some N and stays above at others (measured under J_l = 1 + 0.3·sin(l+1): 0.074
/// against 0.100 at N = 9, but 0.171 against 0.167 at N = 5). So what disorder removes is the BOUND, not
/// the margin, and a NUMERIC rank read of the count is protected only on the uniform chain.
/// The second fence is an index fact: the mode-side reading does not carry to the Heisenberg
/// chain, where the modes run k = 0..N−1 and k ↦ N − k sends 0 out of range, so it is not even an
/// involution and there is no mode-chiral-odd sector to be odd in
/// (<see cref="ModeReflectionIsAnInvolution"/>). THREE things carry there instead, and the third is the
/// one a shorter sentence drops: the ⌊N/2⌋ count, the gap at M = N, and the chirality in the SITE index,
/// where 𝟙 is R-even and ker(B Bᵀ) is exactly the R-odd site space on both chains.</para>
///
/// <para><b>Prior art the entry was minted without, and it reaches further than "the same matrix".</b>
/// PROOF_R90_FROZEN_DIVISOR Lemma 5 landed five days earlier and covers both open chains at once, as
/// B Bᵀ = (1 − 1/M)·𝟏𝟏ᵀ/N + (I + R)/(2M); at M = N + 1 the two 𝟏𝟏ᵀ coefficients agree,
/// (1 − 1/M)/N = 1/M. It also derives the SPECTRUM {1, (1/M)^(⌈N/2⌉−1), 0^⌊N/2⌋} and the pseudo-determinant
/// M^(−⌊(N−1)/2⌋), so what F143 adds over it is the mode-side READING and not the eigenvalues. The
/// qualifier the proof states and a shorter sentence drops: Lemma 5's object is the SITE-indexed B Bᵀ and
/// this entry's is the MODE-indexed G = WᵀW, and the two coincide only because B is symmetric on the XY
/// chain; on the Heisenberg chain B Bᵀ ≠ Bᵀ B and only the site-indexed one satisfies the two-chain law.
/// So "the same matrix" is true at M = N + 1 and only there.</para>
///
/// <para>Typed parents: <see cref="ChiralKClaim"/> for the involution R, <see cref="FrozenDivisorClaim"/>
/// for the ⌊N/2⌋ count, which is its <see cref="FrozenDivisorClaim.FrozenMultiplicity"/>. Neither alone
/// gives the statement: the first has the reflection without the matrix, the second the count without the
/// mode index. The MODULUS is deliberately not attributed to either, because neither types it: M appears
/// in FrozenDivisorClaim as prose only, and the object that types it is MirrorWorld's
/// <c>Divisor.ClockModulus</c>, across a project boundary this assembly does not reference. So
/// <see cref="TransformLength"/> is a knowing second spelling of that member's XY branch, recorded here
/// rather than left to be discovered as drift.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F143 + <c>docs/proofs/PROOF_FROZEN_BAND_SO4.md</c> §6
/// (Lemma 6.1 and Theorem 6.2) + <c>docs/proofs/PROOF_R90_FROZEN_DIVISOR.md</c> Lemma 5 (the two-chain
/// form) + <c>simulations/eta_ceiling_reduction.py</c> block V5 (the Python gate, to N = 40 under
/// <c>--deep</c>).</para></summary>
public sealed class SeedRungGramClaim : Claim
{
    /// <summary>The involution R this claim reads its kernel in: the sublattice gauge K acts on the modes
    /// as k ↦ M − k. Typed parent because the reflection is that claim's and not this one's.</summary>
    public ChiralKClaim ChiralK { get; }

    /// <summary>The site-side carrier of the ⌊N/2⌋ count: one frozen mode per balanced pair. Typed parent
    /// because Lemma 5 of its proof is this matrix, earlier and for both chains. It does NOT type the
    /// modulus, which is prose there; see the class summary.</summary>
    public FrozenDivisorClaim FrozenDivisor { get; }

    /// <summary>The transform length M := N + 1. Not a property of the chain but of the discrete sine
    /// transform the hopping matrix diagonalises under, which is the whole content of the gap being
    /// 1/M: it closes LINEARLY in N.</summary>
    public static int TransformLength(int n) =>
        n >= 1 ? n + 1 : throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1");

    /// <summary>The chiral partner ā = M − a of a mode, one-based (modes run a = 1..N). An involution on
    /// that range for every N; the middle mode a = M/2 is its own partner and exists exactly at odd N.
    /// </summary>
    public static int ChiralPartner(int n, int mode)
    {
        int m = TransformLength(n);
        if (mode < 1 || mode > n)
            throw new ArgumentOutOfRangeException(nameof(mode), mode, $"modes of the open chain run 1..{n}");
        return m - mode;
    }

    /// <summary>Whether the reflection k ↦ modulus − k maps the mode range into itself, i.e. whether it is
    /// an involution on the modes at all. The fence, computed rather than asserted: the XY chain reads
    /// (1, N, N+1) and gets true; the Heisenberg chain reads (0, N−1, N) and gets FALSE, because k = 0 is
    /// sent to N, outside the range. There the site reversal acts on modes as the diagonal sign (−1)^k and
    /// not as a permutation, so this claim's kernel reading has nothing to be odd in.</summary>
    public static bool ModeReflectionIsAnInvolution(int modeMin, int modeMax, int modulus)
    {
        if (modeMin > modeMax)
            throw new ArgumentOutOfRangeException(nameof(modeMax), modeMax, "the mode range must be non-empty");
        for (int k = modeMin; k <= modeMax; k++)
        {
            int image = modulus - k;
            if (image < modeMin || image > modeMax) return false;
        }
        return true;
    }

    /// <summary>Ĝ := 2M·G, the seed-rung Gram matrix scaled to integers, one-based modes stored at
    /// [a−1, c−1]. Every entry is Ĝ_{ac} = 2 + [a = c] + [a + c = M], which is 2·𝟏𝟏ᵀ + I + R: 2 generically,
    /// 3 on the diagonal and on the anti-diagonal, and 4 where the two coincide at the self-paired middle
    /// mode of an odd chain. Those are the proof's four cases, and the scaling is what makes them integers:
    /// G itself reads 1/M, 3/(2M) and 2/M.</summary>
    public static long[,] ScaledGram(int n)
    {
        int m = TransformLength(n);
        var g = new long[n, n];
        for (int a = 1; a <= n; a++)
            for (int c = 1; c <= n; c++)
                g[a - 1, c - 1] = 2 + (a == c ? 1 : 0) + (a + c == m ? 1 : 0);
        return g;
    }

    /// <summary>One entry of G itself, for readers who want the unscaled matrix: Ĝ_{ac}/(2M). Carried
    /// beside <see cref="ScaledGram"/> rather than instead of it, because this is the spelling that cannot
    /// be compared exactly.</summary>
    public static double GramEntry(int n, int a, int c)
    {
        if (a < 1 || a > n) throw new ArgumentOutOfRangeException(nameof(a), a, $"modes run 1..{n}");
        if (c < 1 || c > n) throw new ArgumentOutOfRangeException(nameof(c), c, $"modes run 1..{n}");
        int m = TransformLength(n);
        return ScaledGram(n)[a - 1, c - 1] / (2.0 * m);
    }

    /// <summary>The kernel dimension ⌊N/2⌋: the number of 2-cycles of R on the modes, one per chiral pair,
    /// with the self-paired middle mode of an odd chain contributing nothing. Delegates to
    /// <see cref="FrozenDivisorClaim.FrozenMultiplicity"/> on purpose — the count is that claim's, adopted
    /// as a number there and read off a matrix here, and two spellings of one integer is exactly the kind
    /// of drift a typed layer exists to prevent.</summary>
    public static int FrozenNullity(int n) => FrozenDivisorClaim.FrozenMultiplicity(n);

    /// <summary>The multiplicity ⌈N/2⌉ − 1 of the middle eigenvalue: the R-even sector minus the all-ones
    /// vector. Empty at N ∈ {1, 2}; at N = 2 that is why the gap is 1 and not 1/M, while at N = 1 there is
    /// no kernel either and <see cref="Gap"/> refuses the question.</summary>
    public static int MiddleBandMultiplicity(int n) =>
        n >= 1 ? (n + 1) / 2 - 1 : throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1");

    /// <summary>The gap above the frozen eigenvalue: 1/M = 1/(N + 1) whenever the middle band is occupied,
    /// and 1 at N = 2 where it is empty and the next eigenvalue up is the simple 1. Derived from
    /// <see cref="MiddleBandMultiplicity"/> rather than special-cased on N, so the two cannot disagree.
    /// Requires a kernel to be a gap ABOVE something, hence N ≥ 2.</summary>
    public static double Gap(int n)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n,
            "at N = 1 the kernel is empty, so there is no frozen eigenvalue for a gap to sit above");
        return MiddleBandMultiplicity(n) > 0 ? 1.0 / TransformLength(n) : 1.0;
    }

    /// <summary>The three multiplicities, which must exhaust the N modes: ⌊N/2⌋ + (⌈N/2⌉ − 1) + 1 = N.
    /// Checked as an identity rather than trusted, because it is the one place a miscounted orbit would
    /// hide.</summary>
    public static (int Frozen, int Middle, int Top) Multiplicities(int n) =>
        (FrozenNullity(n), MiddleBandMultiplicity(n), 1);

    /// <summary>Why the exact half is exact, in one paragraph, for a reader who arrives at the witness and
    /// wonders where the tolerance went.</summary>
    public string TheIntegerScaling =>
        "G's entries are 1/M, 3/(2M) and 2/M, and none of the three is a binary fraction, so an entry-wise " +
        "comparison of G against anything would carry a rounding it did not need. 2M·G is integral: " +
        "Ĝ_{ac} = 2 + [a = c] + [a + c = M]. The eigenrelations scale with it and stay integral — " +
        "Ĝ·(e_k − e_{M−k}) = 0, Ĝ·(o_i − o_j) = 2(o_i − o_j) on differences of R-orbit vectors, and " +
        "Ĝ·𝟏 = 2M·𝟏 — so the witness compares longs to 0 and there is nothing to tolerate. The eigenvalues " +
        "0, 1/M and 1 of the entry are those integer relations divided by 2M, and the division is done once, " +
        "for printing.";

    /// <summary>Where the closed form's derivation actually sits, so the witness is not mistaken for the
    /// proof. Recorded because the four entry cases look like a definition and are not.</summary>
    public string WhereTheClosedFormComesFrom =>
        "W_{lk} = (2/M)·sin²(πkl/M) = (1/M)(1 − cos(2πkl/M)), so M²·G_{ac} expands into four cosine sums " +
        "over l = 1..M−1: the constant, the two single cosines and the product. Each cosine sum is the " +
        "exact integer Σ_l cos(2πml/M) = M·[M | m] − 1, and together they give " +
        "M + (M/2)([M | a−c] + [M | a+c]), hence G_{ac} = 1/M + ([a = c] + [a + c = M])/(2M) " +
        "for modes in 1..N, which is the scaled form this claim carries. So the integer entries are the " +
        "OUTPUT of an exact evaluation, not a hand-written table, and the witness's measured half is the " +
        "independent route: rebuild W from the sines and form WᵀW.";

    /// <summary>What this claim deliberately does not carry, kept on the object so a later session does not
    /// read the gap as an oversight.</summary>
    public string WhatStaysOutside =>
        "The rest of the seed-rung chain. F144's floor 𝒦 ≥ ℓ(N − ℓ)/(N + 1) is a bound on a COMPRESSION of " +
        "𝒦 to a subspace that must be built before it can be read, F145's spin-1 triplet lives in the cell " +
        "basis of mode SETS, and F146's dim 𝓜 = C(⌊N/2⌋, ℓ)·R_ℓ is a rank count carrying a number-theoretic " +
        "hypothesis on M. Their evidence is measured in simulations/eta_ceiling_reduction.py (F144) and " +
        "simulations/scalar_count.py (F145 and F146) and proved in PROOF_FROZEN_BAND_SO4 §7 and " +
        "PROOF_SCALAR_COUNT; none of the three reduces to an entry-wise read, " +
        "so a typed carrier for them would be a shell around a number rather than a live witness. That is " +
        "not a statement about ranks being out of reach, since this claim's own completeness check IS an " +
        "exact GF(p) rank. This one " +
        "is the chain's ANCHOR in a precise sense: F144 consumes it as the ℓ = 1 rung, and F145's middle " +
        "u₀ = P_k − P_k̄ IS the seed whose kernel is read here.";

    public SeedRungGramClaim(ChiralKClaim chiralK, FrozenDivisorClaim frozenDivisor)
        : base(
            "F143: at the seed rung of the frozen band the surviving operator is the Gram matrix of squared " +
            "mode amplitudes, G = WᵀW with W_{lk} = v_k(l)² (= W² on the uniform chain, where W is " +
            "symmetric), and in closed form " +
            "G = (1/M)·(𝟏𝟏ᵀ + (I + R)/2) with M := N + 1 the DST-I transform length and R the chiral mode " +
            "involution a ↦ M − a, so spec(G) = {0 with multiplicity ⌊N/2⌋, 1/M with multiplicity ⌈N/2⌉ − 1, " +
            "1 simple}, the kernel being exactly the R-ODD sector spanned by the chiral differences " +
            "e_k − e_{M−k}, one per pair, with the self-paired middle mode of an odd chain contributing " +
            "nothing. The gap above the frozen eigenvalue is 1/M = 1/(N + 1) for N ≥ 3 and 1 at N = 2 where " +
            "the middle band is empty; it closes only LINEARLY in N, which is why a numeric rank read keeps " +
            "its moat as N grows. SCOPE, and it is narrower than it sounds: UNIFORM hopping buys the middle " +
            "band's common VALUE 1/M and the entry-wise closed form, while the KERNEL, its identity as the " +
            "R-odd sector, its dimension ⌊N/2⌋ and λ_max = 1 all SURVIVE bond disorder, PROVIDED NO BOND IS " +
            "ZERO: K h K = −h then makes the squared amplitudes chirally equal (W R = W, hence G R = G) on " +
            "a Jacobi matrix, whose spectrum is simple, which is what lets the eigenVECTOR be read off the " +
            "±ε pairing of eigenVALUES at all. With a zero bond every part of that fails together. " +
            "SIMPLICITY of the 1 is NOT claimed: double stochasticity gives the value, not the " +
            "multiplicity. Third: the hopping must be nearest-neighbour and bipartite with NO on-site term, " +
            "and the boundary OPEN. Second scope: the MODE-side reading is the open XY " +
            "chain's alone, since on the Heisenberg chain the modes run k = 0..N−1 and k ↦ N − k is not an " +
            "involution on them. THREE things carry to that chain: the ⌊N/2⌋ count, the gap at M = N, and " +
            "the chirality in the SITE index, where ker(B Bᵀ) is exactly the R-odd site space on both " +
            "chains. The count had a site-side carrier before this claim existed.",
            Tier.Tier1Derived,
            "docs/ANALYTICAL_FORMULAS.md F143 + docs/proofs/PROOF_FROZEN_BAND_SO4.md §6 (Lemma 6.1, " +
            "Theorem 6.2) + docs/proofs/PROOF_R90_FROZEN_DIVISOR.md Lemma 5 (the site-indexed B Bᵀ for both " +
            "chains, five days earlier) + simulations/eta_ceiling_reduction.py block V5 + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/SeedRungGramWitness.cs (inspect --root seedrung)")
    {
        ChiralK = chiralK ?? throw new ArgumentNullException(nameof(chiralK));
        FrozenDivisor = frozenDivisor ?? throw new ArgumentNullException(nameof(frozenDivisor));
    }

    public override string DisplayName =>
        "F143: the seed rung in closed form, its kernel the chiral-odd sector and its gap one over the " +
        "transform length";

    public override string Summary =>
        "G = (1/M)·(𝟏𝟏ᵀ + (I + R)/2) with M = N + 1; spec(G) = {0^⌊N/2⌋, (1/M)^(⌈N/2⌉−1), 1}; kernel exactly " +
        "the R-odd modes, gap 1/M closing linearly in N. 2M·G is an INTEGER matrix, so the whole spectrum is " +
        "certified without an eigensolver and without a tolerance; the one measured statement is that the " +
        "DST-I modes produce that matrix at all. Open XY chain, and uniform hopping only for the middle " +
        "band's value and the entry-wise form: the kernel survives disorder on any chain with no zero bond " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the closed form",
                summary: "G = (1/M)·(𝟏𝟏ᵀ + (I + R)/2), equivalently 2M·G = 2·𝟏𝟏ᵀ + I + R with integer " +
                         "entries 2, 3 and 4; 𝟏𝟏ᵀ is written as an outer product because the letter J that " +
                         "matrix usually carries is the coupling here");
            yield return new InspectableNode("where it comes from", summary: WhereTheClosedFormComesFrom);
            yield return new InspectableNode("the integer scaling", summary: TheIntegerScaling);
            yield return new InspectableNode("the spectrum at N = 7",
                summary: $"multiplicities {Multiplicities(7)}, gap {Gap(7):0.######} = 1/8, " +
                         "kernel spanned by e_1 − e_7, e_2 − e_6, e_3 − e_5, the self-paired mode 4 " +
                         "contributing nothing");
            yield return new InspectableNode("the live read at N = 9",
                summary: $"{SeedRungGramWitness.ExactResidual(9)} is the exact integer residual of all three " +
                         $"eigenrelations, and {SeedRungGramWitness.ClosedFormResidual(9):0.###e+0} the " +
                         "floating-point deviation of WᵀW rebuilt from the sines");
            yield return new InspectableNode("the two fences",
                summary: "UNIFORM hopping, and only for part of the statement: under bond disorder with no " +
                         "ZERO bond (a Jacobi matrix, hence the simple h-spectrum the step needs) W R = W " +
                         "still holds, so the kernel stays the R-odd sector at dimension ⌊N/2⌋ and " +
                         "λ_max = 1 stays, though not provably simple; what splits is the middle band's " +
                         "common value 1/M, and with it the entry-wise closed form and the moat's " +
                         "GUARANTEE (the first non-kernel eigenvalue is no longer bounded below by 1/M). " +
                         "HEISENBERG: the mode reflection is not an involution " +
                         "there (k ↦ N − k sends 0 out of range), so this kernel reading is the open XY " +
                         "chain's alone, while the count, the gap at M = N and the SITE-index chirality " +
                         "carry");
            yield return new InspectableNode("what stays outside", summary: WhatStaysOutside);
            yield return ChiralK;
            yield return FrozenDivisor;
        }
    }
}
