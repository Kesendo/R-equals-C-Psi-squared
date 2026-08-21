using System.Collections.Generic;
using System.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The general (w_ket, w_bra) computational-basis coherence sub-block of the Z-dephasing Liouvillian
/// L = −i[H, ρ] + D[ρ] on an N-site chain (the XY hopping H = J·Σ(XX+YY) at J = q). Basis = |a⟩⟨b| with
/// popcount(a) = wKet, popcount(b) = wBra (ket configs outer, bra configs inner, both in ascending-mask order);
/// the diagonal is the Absorption-Theorem DIAGONAL −2·<see cref="DisagreementSum"/> = −2·Σ_s γ_s·[a_s ≠ b_s], ket
/// excitations hop −2iq (the −iHρ term) and bra excitations +2iq (the +iρH term), nearest-neighbour with Pauli
/// exclusion. q-linear: L(q) = A + q·C with A the real AT diagonal and C the pure-imaginary hopping.
///
/// <para><b>γ is a per-site PROFILE; uniform γ = 1 is its default, not its definition.</b> The diagonal entry of
/// the coherence |a⟩⟨b| is −2·Σ_s γ_s·[a_s ≠ b_s], the diagonal instance of the Absorption Theorem's vector form
/// (<see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>, 2026-05-29), so no law is added here; this
/// builder is catching up with one the repo already owned. Without a <c>gammaPerSite</c> argument the build is
/// γ_s ≡ 1 and the entry collapses to −2·n_diff: at uniform γ which sites disagree is invisible, only how many.
/// Under a profile it is a SUBSET SUM instead of a count.</para>
///
/// <para>What that costs SPECTRALLY is deliberately not stated here: the diagonal is not the spectrum, and the
/// answer depends on the block and on arg q. For the (0,1) edge block of a path,
/// docs/proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md owns it; F152 and F153 carry the wider statements.
/// Uniform-only and saying so where they live: <see cref="BuildReflectionSectorColumnMajor"/> and both emitters
/// of <see cref="WeightCoherenceSectorCsr"/>.</para>
///
/// <para><b>Site-index convention, and it is the OPPOSITE of the block-spectrum builder's.</b> Site s is BIT s
/// (site 0 = least significant bit): what the hopping loops mean by nearest-neighbour and what
/// <see cref="FieldEnergy"/> already means by <c>field[k]</c>, so the two per-site knobs agree.
/// <see cref="RCPsiSquared.Core.BlockSpectrum.PerBlockLiouvillianBuilder"/> reads site l as bit N−1−l, so a
/// profile handed to both must be REVERSED. Gated in <c>WeightCoherenceBlockGammaProfileTests</c>, which also
/// records why the mislabel hides: a uniform profile is its own reverse.</para>
///
/// <para>Promoted (verbatim physics) from the CLI's FoldCrossCommand.BuildBlock so the cross-fold partner block
/// (SE, w_{N−2}) and the (SE,DE) block are built by ONE shared builder. The partner-pairing carrier is
/// <see cref="BraComplementPermutation"/>: the branch-locus palindrome's bra bit-flip ρ[a,b] → ρ[a,b̄] maps the
/// (wKet, wBra) block to the (wKet, N−wBra) block, with S(a,b̄) = Σ_s γ_s − S(a,b) reflecting the AT diagonal
/// about −σ = −Σ_s γ_s, which is −N at uniform γ = 1. The similarity is ANTIUNITARY (entry-wise conjugation and
/// q → q̄) and holds at an arbitrary profile with the affine constant 2σ, gated in
/// <c>WeightCoherenceBlockGammaProfileTests</c>; it is not F1's linear palindrome, which shares the 2σ and
/// nothing else. See
/// experiments/F89_BRANCH_LOCUS_PALINDROME.md and the diabolic cross-fold (Move 4).</para></summary>
public static class WeightCoherenceBlock
{
    /// <summary>All n-bit masks with exactly w set bits, in ascending order.</summary>
    public static List<int> Configs(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }

    /// <summary>The (wKet, wBra) chain coherence block at complex coupling q (γ = 1), dim C(n,wKet)·C(n,wBra).
    /// Diagonal −2·n_diff; ket excitations hop −2iq, bra excitations +2iq (nearest-neighbour, Pauli-excluded).
    /// The pure-XY (Δ=0) case; delegates to the (q,Δ) overload.</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q) => Build(n, wKet, wBra, q, 0.0);

    /// <summary>The (wKet, wBra) XXZ-chain coherence block at complex coupling q (γ = 1) and real ZZ-anisotropy Δ,
    /// dim C(n,wKet)·C(n,wBra). The bond Hamiltonian is H = J·Σ(X_bX_{b+1}+Y_bY_{b+1}) + J·Δ·Σ Z_bZ_{b+1}, q = J.
    /// On top of the XY block (diagonal −2·n_diff; ket excitations hop −2iq, bra excitations +2iq, NN, Pauli-
    /// excluded), the Δ·ZZ term is a DIAGONAL Hermitian contribution, so it leaves the Absorption-Theorem
    /// diagonal untouched (−2·n_diff at uniform γ, −2·Σ_s γ_s·[a_s ≠ b_s] under a profile) and adds only the
    /// frequency −i·q·Δ·(zz(ket) − zz(bra)), with
    /// zz(c) = Σ_bond ⟨c|Z_bZ_{b+1}|c⟩ (<see cref="Zz"/>). Matches XxzCoherenceBlock.BuildFull's convention; at
    /// Δ=0 it reproduces the pure-XY block exactly. The Δ·ZZ term is EVEN under the global bit-flip
    /// (Z_bZ_{b+1} ↦ (−Z_b)(−Z_{b+1}) = Z_bZ_{b+1}, so zz(b̄) = zz(b)), which is exactly why the cross-fold
    /// antiunitary similarity (<see cref="BraComplementPermutation"/>, F89d) survives at every Δ.</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q, double delta) =>
        Build(n, wKet, wBra, q, delta, null, null);

    /// <summary>Σ_s γ_s·[a_s ≠ b_s], the γ-weighted count of bra-ket disagreements of ONE computational
    /// coherence; the block's diagonal ENTRY is −2 times it. It is not an eigenvalue and not a decay rate: the
    /// Absorption Theorem reads −Re λ = 2·Σ_l γ_l·⟨Δ_l⟩ over eigenmodes, with the expectation ⟨Δ_l⟩ ∈ [0,1],
    /// and the hopping is what turns the sharp bit into that expectation. Site s is bit s (see the class doc).
    /// A null profile means uniform γ ≡ 1 and takes the popcount route; a profile is summed in ascending site
    /// order.</summary>
    public static double DisagreementSum(int n, IReadOnlyList<double>? gammaPerSite, int a, int b)
    {
        int x = a ^ b;
        if (gammaPerSite is null) return BitOperations.PopCount((uint)x);
        if (gammaPerSite.Count != n)
            throw new ArgumentException($"gamma profile length {gammaPerSite.Count} != n={n}", nameof(gammaPerSite));
        double r = 0.0;
        for (int s = 0; s < n; s++)
            if (((x >> s) & 1) != 0) r += gammaPerSite[s];
        return r;
    }

    /// <summary>The (wKet, wBra) XXZ-chain coherence block at (q, Δ), an optional longitudinal field, and an
    /// optional per-site dephasing PROFILE γ_s (<c>gammaPerSite[s]</c> = the rate OF SITE s = bit s; null means
    /// uniform γ ≡ 1). The profile enters the Absorption-Theorem diagonal only, as −2·Σ_s γ_s·[a_s ≠ b_s]
    /// (<see cref="DisagreementSum"/>): dephasing is diagonal in the computational basis, so it cannot touch
    /// the hopping, and the q-linear split L(q) = A + q·C survives with the profile living entirely in A.
    /// Field and profile compose, gated by <c>TheFieldAndTheProfileCompose</c>.
    /// What a profile changes is that the diagonal entry stops being a function of n_diff alone; for what that
    /// costs spectrally see the class doc's pointer. Note that A is the Hermitian part of L only at REAL q,
    /// where C is anti-Hermitian, and this builder takes q complex.</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q, double delta, double[]? field,
                                   IReadOnlyList<double>? gammaPerSite)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var index = new Dictionary<(int, int), int>();
        var basis = new List<(int Ket, int Bra)>();
        foreach (var k in kets)
            foreach (var b in bras) { index[(k, b)] = basis.Count; basis.Add((k, b)); }
        int d = basis.Count;
        var l = new Complex[d, d];
        for (int col = 0; col < d; col++)
        {
            var (kc, bc) = basis[col];
            l[col, col] += new Complex(-2.0 * DisagreementSum(n, gammaPerSite, kc, bc), 0)
                         + (-Complex.ImaginaryOne) * q * (delta * (Zz(n, kc) - Zz(n, bc)));   // Δ·ZZ frequency
            if (field != null)
                l[col, col] += (-Complex.ImaginaryOne) * q * (FieldEnergy(n, field, kc) - FieldEnergy(n, field, bc));
            for (int s = 0; s < n; s++)
                if ((kc & (1 << s)) != 0)                                   // ket excitation hops (−2iq)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] += Complex.ImaginaryOne * -2.0 * q;
            for (int s = 0; s < n; s++)
                if ((bc & (1 << s)) != 0)                                   // bra excitation hops (+2iq)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] += Complex.ImaginaryOne * 2.0 * q;
        }
        return l;
    }

    /// <summary>The (wKet, wBra) XXZ-chain coherence block at (q, Δ) plus a per-site longitudinal Z-field
    /// Σ_k w_k Z_k (field[k] = w_k, the integrability-/symmetry-breaking disorder knob). The field is DIAGONAL and
    /// Hermitian, so like the Δ·ZZ term it leaves the Absorption-Theorem diagonal untouched (−2·n_diff at uniform
    /// γ) and adds
    /// only the frequency −i·q·(fe(ket) − fe(bra)), with fe(c) = Σ_k w_k·z_k (z_k = −1 if site k excited, +1 else,
    /// <see cref="FieldEnergy"/>). UNLIKE the Δ·ZZ term, the field is bit-flip-ODD (fe(c̄) = −fe(c)), so it BREAKS
    /// the cross-fold antiunitary similarity (the negative control in <see cref="WeightCoherenceBlockTests"/>) and
    /// breaks the S₂ reflection + conjugation symmetry — exactly the disorder needed to drive a dense coherence
    /// sector toward GinUE (the F89 Door-C filling-threshold test). field=null reproduces Build(n,wKet,wBra,q,Δ).</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q, double delta, double[]? field) =>
        Build(n, wKet, wBra, q, delta, field, null);

    /// <summary>fe(c) = Σ_k w_k·z_k, z_k = −1 if site k is excited (bit set), +1 otherwise — the longitudinal-field
    /// energy of a computational-basis config. Bit-flip-ODD: fe(c̄) = −fe(c) (each z_k flips sign).</summary>
    public static double FieldEnergy(int n, double[] w, int c)
    {
        double e = 0;
        for (int k = 0; k < n; k++) e += w[k] * (((c >> k) & 1) == 1 ? -1.0 : 1.0);
        return e;
    }

    /// <summary>zz(c) = Σ_{bond (b,b+1)} ⟨c|Z_bZ_{b+1}|c⟩ = Σ_b (+1 if bits b, b+1 are equal, −1 if they differ),
    /// the open-chain ZZ-bond sum of the computational-basis config c. Even under the global bit-flip, so
    /// zz(c̄) = zz(c) (each Z flips sign, the product is unchanged).</summary>
    public static int Zz(int n, int c)
    {
        int s = 0;
        for (int b = 0; b < n - 1; b++)
            s += (((c >> b) & 1) == ((c >> (b + 1)) & 1)) ? 1 : -1;
        return s;
    }

    /// <summary>The bra-complement permutation P: the basis index of |a⟩⟨b| in the (wKet, wBra) block ↦ the
    /// basis index of |a⟩⟨b̄| in the (wKet, n−wBra) block (b̄ = the n-site bitwise complement of b). The carrier
    /// of the cross-fold: since the disagreement sum obeys S(a,b̄) = Σ_s γ_s − S(a,b) (n − n_diff at uniform γ),
    /// conjugating L(wKet,wBra) by P, conjugating its ENTRIES and sending q → q̄ (the map is antiunitary, not
    /// a linear similarity) maps
    /// it onto L(wKet, n−wBra). A bijection because C(n,wBra) = C(n,n−wBra) and the ket weight is unchanged.
    /// Returns perm where perm[t] = the (wKet, n−wBra)-basis index that the (wKet, wBra)-basis index t maps to.</summary>
    public static int[] BraComplementPermutation(int n, int wKet, int wBra)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var brasC = Configs(n, n - wBra);
        var indexC = new Dictionary<int, int>();
        for (int j = 0; j < brasC.Count; j++) indexC[brasC[j]] = j;
        int full = (1 << n) - 1;
        var perm = new int[kets.Count * bras.Count];
        int t = 0;
        for (int ki = 0; ki < kets.Count; ki++)
            for (int bi = 0; bi < bras.Count; bi++)
            {
                int bBar = full ^ bras[bi];
                perm[t++] = ki * brasC.Count + indexC[bBar];     // same ket block, complemented bra index
            }
        return perm;
    }

    /// <summary>The ket-complement permutation Q: the basis index of |a⟩⟨b| in the (wKet, wBra) block ↦ the basis
    /// index of |ā⟩⟨b| in the (n−wKet, wBra) block (ā = the n-site bitwise complement of a). The KET-leg mirror of
    /// <see cref="BraComplementPermutation"/>: flipping the ket index is left-multiplication F·ρ (F = X^⊗N), so it
    /// maps (wKet, wBra) → (n−wKet, wBra) and, like the bra leg, flips the AT rate n_diff(ā,b) = n − n_diff(a,b),
    /// giving the SAME −2n affine reflection. Conjugating L(wKet,wBra) by Q and reflecting maps it onto
    /// L(n−wKet, wBra) (the ket-leg cross-fold). A bijection because C(n,wKet) = C(n,n−wKet) and the bra weight is
    /// unchanged. Returns perm where perm[t] = the (n−wKet, wBra)-basis index that the (wKet, wBra)-basis index t
    /// maps to.
    ///
    /// <para>Convention bridge: the D₄ proof docs (PROOF_PI_FACTORS_AS_R_TIMES_D, F118 / MirrorGroupD4Claim) name
    /// the spine V₄ = {I, F⊗F, I⊗F, F⊗I} ⊂ D₄ by multiplication side, calling left-mult F·ρ the "bra reflection";
    /// F89 names by the flipped INDEX, so this ket-complement Q (flips the ket) = left-mult F·ρ = the spine element
    /// 𝓕R = Π²·R. Its bra-leg partner P (<see cref="BraComplementPermutation"/>, flips the bra) = right-mult ρ·F =
    /// the spine element R, a factor of the F1 palindrome Π = R·D. Stating both by the flipped index keeps F89d's
    /// existing name; the spine docs' word for ρ·F is the opposite ("ket reflection").</para></summary>
    public static int[] KetComplementPermutation(int n, int wKet, int wBra)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var ketsC = Configs(n, n - wKet);
        var indexC = new Dictionary<int, int>();
        for (int j = 0; j < ketsC.Count; j++) indexC[ketsC[j]] = j;
        int full = (1 << n) - 1;
        var perm = new int[kets.Count * bras.Count];
        int t = 0;
        for (int ki = 0; ki < kets.Count; ki++)
            for (int bi = 0; bi < bras.Count; bi++)
            {
                int kBar = full ^ kets[ki];
                perm[t++] = indexC[kBar] * bras.Count + bi;      // complemented ket index, same bra block
            }
        return perm;
    }

    /// <summary>The site-reflection permutation R: the basis index of |a⟩⟨b| ↦ the index of |rev(a)⟩⟨rev(b)|
    /// (rev = the n-site bit-order reversal). An involution; commutes ENTRY-WISE with the uniform chain block
    /// at every q and Δ (the bond set is reflection-symmetric, each entry a single hop contribution). The
    /// general-(wKet,wBra) sibling of <see cref="F89PathKSeDeBlock.ReflectionPermutation"/>; carrier of the
    /// R-parity split the step-3 shell census uses as its LU-cost lever.</summary>
    public static int[] ReflectionPermutation(int n, int wKet, int wBra)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var ketIndex = new Dictionary<int, int>();
        var braIndex = new Dictionary<int, int>();
        for (int i = 0; i < kets.Count; i++) ketIndex[kets[i]] = i;
        for (int i = 0; i < bras.Count; i++) braIndex[bras[i]] = i;
        var perm = new int[kets.Count * bras.Count];
        int t = 0;
        foreach (var k in kets)
            foreach (var b in bras)
                perm[t++] = ketIndex[ReverseBits(n, k)] * bras.Count + braIndex[ReverseBits(n, b)];
        return perm;
    }

    static int ReverseBits(int n, int c)
    {
        int r = 0;
        for (int s = 0; s < n; s++)
            if (((c >> s) & 1) == 1) r |= 1 << (n - 1 - s);
        return r;
    }

    /// <summary>σ(c) = the sum of the OCCUPIED SITE INDICES of a computational-basis config (not the popcount).
    /// The staggered bipartite grading: a nearest-neighbour hop moves one excitation by one site, so it changes
    /// σ by exactly ±1 and therefore flips (−1)^σ. Under the site reversal σ(rev(c)) = w·(n−1) − σ(c), which is
    /// what makes <see cref="BipartiteGaugeCommutesWithReflection"/> a parity of (p+q)(n−1) and nothing else.</summary>
    public static int SiteIndexSum(int n, int c)
    {
        int s = 0;
        for (int k = 0; k < n; k++)
            if (((c >> k) & 1) == 1) s += k;
        return s;
    }

    /// <summary>The BIPARTITE GAUGE 𝒟 on the (wKet, wBra) block: the diagonal of signs
    /// 𝒟[t] = (−1)^(σ(a)+σ(b)) for the basis element |a⟩⟨b| at index t, σ = <see cref="SiteIndexSum"/>.
    /// Written 𝒟 in PROOF_CODIM1_BY_ADDITIVITY §7 ingredient (iv) and in experiments/F89_PATH_K_DIABOLIC.md;
    /// the same object is called T in <see cref="FoldResultantCertificate"/> and in the OpenArcs registry, and
    /// its p+q = 3 specialisation is BetaExoticPerNExclusionClaim's metric. It is NOT the fold antiunitary
    /// T = P·K, NOT the transpose leg, and NOT the bra/ket complement permutations above.
    ///
    /// <para>What it does, entry-wise and with no eigensolver: L(q) = A + q·C with A the real AT diagonal and
    /// C the pure-imaginary hopping, every entry of C a single nearest-neighbour hop. Conjugating by 𝒟 leaves
    /// A alone (𝒟 is diagonal) and negates C (a hop flips the bipartite parity on exactly one side), so
    /// 𝒟·L(q)·𝒟 = L(−q) at every complex q, and at REAL q that is also conj(L). Both are exact rearrangements
    /// of ±1 sign flips, so the residual is 0.0 bit-for-bit rather than machine zero.</para>
    ///
    /// <para>Scope, and both fences are gated: Δ = 0 and no longitudinal field. The Δ·ZZ term and the field
    /// contribute an IMAGINARY diagonal, which 𝒟 cannot touch, so the identity breaks there; that is the same
    /// Δ-scope PROOF_CODIM1 §7 puts on ingredient (iv). Nearest-neighbour hops on an open chain are load-bearing
    /// too: a next-nearest bond does not flip the bipartite parity.</para></summary>
    public static int[] BipartiteGaugeSigns(int n, int wKet, int wBra)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var signs = new int[kets.Count * bras.Count];
        int t = 0;
        foreach (var k in kets)
            foreach (var b in bras)
                signs[t++] = ((SiteIndexSum(n, k) + SiteIndexSum(n, b)) & 1) == 0 ? 1 : -1;
        return signs;
    }

    /// <summary>THE GAUGE CRITERION: 𝒟 commutes with the site reflection R exactly when (wKet + wBra)·(n − 1)
    /// is even, and anticommutes otherwise. One line: σ(rev(c)) = w·(n−1) − σ(c), so reversing both sides
    /// multiplies the gauge sign by (−1)^((wKet+wBra)(n−1)), the SAME factor at every basis index, which is why
    /// the alternative is commute-or-anticommute with nothing in between.
    ///
    /// <para>This is the clause PROOF_CODIM1_BY_ADDITIVITY §7's grading paragraph leaves open when it says
    /// P, Q, 𝒟, T "all commute with R up to a block-global sign" without saying when the sign is −1. Its
    /// consequence is the one that gets used: where it commutes, 𝒟 restricts to each R-parity sector, so each
    /// sector is a real matrix family and its spectrum is conjugation-closed; where it anticommutes, the gauge
    /// SWAPS the two sectors and supplies nothing about either. The criterion is SUFFICIENT, not necessary:
    /// nothing here excludes some other antiunitary closing a sector on the odd parity.</para></summary>
    public static bool BipartiteGaugeCommutesWithReflection(int n, int wKet, int wBra)
        => (((wKet + wBra) * (n - 1)) & 1) == 0;

    /// <summary>One R-parity sector of the (wKet,wBra) chain block at complex q (γ = 1, Δ = 0), assembled
    /// DIRECTLY in sector coordinates — the full block is never materialized, which is what fits the N=11
    /// census blocks under the managed-array LP64 wall (flat dim ≤ 46340). Returns (column-major flat
    /// matrix, sector dim). Sector basis: reflection fixed points e_f (even sector only, weight 1), then
    /// 2-cycle combinations (e_t ± e_{Rt})/√2 for orbit reps t &lt; Rt, in increasing t — the same convention
    /// as <see cref="F89PathKSeDeBlock.ROddBasis"/>. The basis is REAL orthonormal, so
    /// spec(full) = spec(even) ⊎ spec(odd) exactly and σ_min(full − s) = min over the two sectors.
    ///
    /// <para><b>Uniform γ ≡ 1, and here that is structural rather than an omission.</b> This route takes no
    /// <c>gammaPerSite</c> because the sector split is built out of <see cref="ReflectionPermutation"/>: the
    /// decomposition exists only when the chain reflection is a symmetry of L, and a γ profile that is not
    /// palindromic generically breaks it, so spec(full) = spec(even) ⊎ spec(odd) stops holding. The exceptions
    /// are the blocks whose achievable XOR masks are all palindromic, so that R commutes with any profile: dim 1,
    /// and n = 2 at (1,1). A PALINDROMIC
    /// profile would preserve it; that case is not implemented, and adding it means proving the commutation
    /// gated first rather than passing an array. Use the dense seven-argument <c>Build</c> for any profile. The same
    /// applies to <see cref="WeightCoherenceSectorCsr"/>, which TRANSCRIBES this rule rather than calling it,
    /// and so has to be threaded separately. Its full-basis emitter carries no such obstruction and is simply not
    /// threaded; see that file.</para></summary>
    public static (Complex[] A, int Dim) BuildReflectionSectorColumnMajor(int n, int wKet, int wBra, Complex q, bool odd)
    {
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        int nb = bras.Count;
        var ketIndex = new Dictionary<int, int>();
        var braIndex = new Dictionary<int, int>();
        for (int i = 0; i < kets.Count; i++) ketIndex[kets[i]] = i;
        for (int i = 0; i < bras.Count; i++) braIndex[bras[i]] = i;
        var perm = ReflectionPermutation(n, wKet, wBra);

        // orbit reps: fixed points (even only), then 2-cycle reps t < perm[t]; orbitOf = sector row or −1
        var reps = new List<int>();
        var orbitOf = new int[perm.Length];
        Array.Fill(orbitOf, -1);
        for (int i = 0; i < perm.Length; i++)
        {
            if (perm[i] == i) { if (!odd) { orbitOf[i] = reps.Count; reps.Add(i); } }
            else if (i < perm[i]) { orbitOf[i] = orbitOf[perm[i]] = reps.Count; reps.Add(i); }
        }
        int d = reps.Count;
        var a = new Complex[(long)d * d];
        double sSign = odd ? -1.0 : 1.0;
        double inv2 = 1.0 / Math.Sqrt(2.0);

        // the exact column action of L on full-basis index fullCol, scaled by weight (the same hop rule
        // as Build: diagonal −2·n_diff, ket hops −2iq, bra hops +2iq, nearest-neighbour, Pauli-excluded)
        void ApplyColumn(int fullCol, Complex weight, Dictionary<int, Complex> acc)
        {
            int kc = kets[fullCol / nb], bc = bras[fullCol % nb];
            AddTo(acc, fullCol, weight * new Complex(-2.0 * BitOperations.PopCount((uint)(kc ^ bc)), 0));
            for (int site = 0; site < n; site++)
                if ((kc & (1 << site)) != 0)
                    foreach (int s2 in new[] { site - 1, site + 1 })
                        if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                            AddTo(acc, ketIndex[(kc & ~(1 << site)) | (1 << s2)] * nb + braIndex[bc],
                                weight * Complex.ImaginaryOne * -2.0 * q);
            for (int site = 0; site < n; site++)
                if ((bc & (1 << site)) != 0)
                    foreach (int s2 in new[] { site - 1, site + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            AddTo(acc, ketIndex[kc] * nb + braIndex[(bc & ~(1 << site)) | (1 << s2)],
                                weight * Complex.ImaginaryOne * 2.0 * q);
        }

        var acc = new Dictionary<int, Complex>();
        for (int col = 0; col < d; col++)
        {
            acc.Clear();
            int t = reps[col];
            bool fixedPt = perm[t] == t;
            double wNorm = fixedPt ? 1.0 : inv2;
            ApplyColumn(t, wNorm, acc);
            if (!fixedPt) ApplyColumn(perm[t], sSign * wNorm, acc);

            foreach (var (fullRow, val) in acc)
            {
                int row = orbitOf[fullRow];
                if (row < 0) continue;                            // row absent from this sector (odd ∌ fixed)
                int rep = reps[row];
                // coefficient of e_fullRow inside the row's REAL sector basis vector (no conjugation)
                double coeff = perm[rep] == rep ? 1.0
                             : (fullRow == rep ? inv2 : sSign * inv2);
                a[(long)col * d + row] += coeff * val;
            }
        }
        return (a, d);
    }

    static void AddTo(Dictionary<int, Complex> acc, int key, Complex v) =>
        acc[key] = acc.TryGetValue(key, out var cur) ? cur + v : v;
}
