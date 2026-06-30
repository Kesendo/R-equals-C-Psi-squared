using System.Collections.Generic;
using System.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The general (w_ket, w_bra) computational-basis coherence sub-block of the Z-dephasing Liouvillian
/// L = −i[H, ρ] + D[ρ] on an N-site chain (γ = 1, the XY hopping H = J·Σ(XX+YY) at J = q). Basis = |a⟩⟨b| with
/// popcount(a) = wKet, popcount(b) = wBra (ket configs outer, bra configs inner, both in ascending-mask order);
/// the diagonal is −2·n_diff(a,b) (the Absorption-Theorem rate 2γ·n_diff, n_diff = popcount(a⊕b)), ket
/// excitations hop −2iq (the −iHρ term) and bra excitations +2iq (the +iρH term), nearest-neighbour with Pauli
/// exclusion. q-linear: L(q) = A + q·C with A the real AT diagonal and C the pure-imaginary hopping.
///
/// <para>Promoted (verbatim physics) from the CLI's FoldCrossCommand.BuildBlock so the cross-fold partner block
/// (SE, w_{N−2}) and the (SE,DE) block are built by ONE shared builder. The partner-pairing carrier is
/// <see cref="BraComplementPermutation"/>: the branch-locus palindrome's bra bit-flip ρ[a,b] → ρ[a,b̄] maps the
/// (wKet, wBra) block to the (wKet, N−wBra) block, with n_diff(a,b̄) = N − n_diff(a,b) reflecting the AT rate
/// about −N. See experiments/F89_BRANCH_LOCUS_PALINDROME.md and the diabolic cross-fold (Move 4).</para></summary>
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
    /// excluded), the Δ·ZZ term is a DIAGONAL Hermitian contribution, so it leaves the Absorption-Theorem rate
    /// Re λ = −2·n_diff untouched and adds only the frequency −i·q·Δ·(zz(ket) − zz(bra)), with
    /// zz(c) = Σ_bond ⟨c|Z_bZ_{b+1}|c⟩ (<see cref="Zz"/>). Matches XxzCoherenceBlock.BuildFull's convention; at
    /// Δ=0 it reproduces the pure-XY block exactly. The Δ·ZZ term is EVEN under the global bit-flip
    /// (Z_bZ_{b+1} ↦ (−Z_b)(−Z_{b+1}) = Z_bZ_{b+1}, so zz(b̄) = zz(b)), which is exactly why the cross-fold
    /// antiunitary similarity (<see cref="BraComplementPermutation"/>, F89d) survives at every Δ.</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q, double delta)
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
            l[col, col] += new Complex(-2.0 * BitOperations.PopCount((uint)(kc ^ bc)), 0)
                         + (-Complex.ImaginaryOne) * q * (delta * (Zz(n, kc) - Zz(n, bc)));   // Δ·ZZ frequency
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
    /// Hermitian, so like the Δ·ZZ term it leaves the Absorption-Theorem rate Re λ = −2·n_diff untouched and adds
    /// only the frequency −i·q·(fe(ket) − fe(bra)), with fe(c) = Σ_k w_k·z_k (z_k = −1 if site k excited, +1 else,
    /// <see cref="FieldEnergy"/>). UNLIKE the Δ·ZZ term, the field is bit-flip-ODD (fe(c̄) = −fe(c)), so it BREAKS
    /// the cross-fold antiunitary similarity (the negative control in <see cref="WeightCoherenceBlockTests"/>) and
    /// breaks the S₂ reflection + conjugation symmetry — exactly the disorder needed to drive a dense coherence
    /// sector toward GinUE (the F89 Door-C filling-threshold test). field=null reproduces Build(n,wKet,wBra,q,Δ).</summary>
    public static Complex[,] Build(int n, int wKet, int wBra, Complex q, double delta, double[]? field)
    {
        var l = Build(n, wKet, wBra, q, delta);
        if (field == null) return l;
        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        int col = 0;
        foreach (var kc in kets)
            foreach (var bc in bras)
            {
                l[col, col] += (-Complex.ImaginaryOne) * q * (FieldEnergy(n, field, kc) - FieldEnergy(n, field, bc));
                col++;
            }
        return l;
    }

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
    /// of the cross-fold: since n_diff(a,b̄) = n − n_diff(a,b), conjugating L(wKet,wBra) by P and reflecting maps
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
}
