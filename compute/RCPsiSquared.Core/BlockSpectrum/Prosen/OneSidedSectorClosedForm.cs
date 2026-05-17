using System.Numerics;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum.Prosen;

/// <summary>Closed-form Liouvillian eigenvalues for the one-sided joint-popcount sectors
/// (p_c = 0, p_r = m) of the chain XY + uniform Z-dephasing Liouvillian. These are the
/// sectors of operators <c>|ψ⟩⟨0|</c> where the ket has m fermions and the bra is the
/// vacuum — the "creation-only" coherences.
///
/// <para><b>Closed form.</b> On the (0, m) sector the dissipator reduces to a constant
/// diagonal: <c>D|ψ⟩⟨0| = −2γ·m·|ψ⟩⟨0|</c> for any m-fermion ψ (Hamming distance from
/// the empty bra is exactly m). The unitary commutator is <c>−i[H, |ψ⟩⟨0|] = −iH|ψ⟩⟨0|</c>
/// since <c>H|0⟩ = 0</c> for the chain-XY vacuum. Therefore</para>
/// <list type="bullet">
///   <item><c>L_(0,m) = −iH_m − 2γm·I</c> on the m-fermion ket-side Hilbert space.</item>
///   <item>Eigenvalues <c>λ_S = −2γm − i·ε_S = Σ_{k∈S} β_k</c> for <c>S ⊆ {1..N}, |S|=m</c>,
///         where <c>ε_S = Σ_{k∈S} ε_k</c> is the m-fermion energy and
///         <c>β_k = −2γ − i·ε_k</c> are the N <b>Prosen rapidities</b> of the model.</item>
///   <item>Sine-mode dispersion <c>ε_k = 2J·cos(πk/(N+1))</c> for OBC (k = 1..N), supplied
///         by <see cref="XyJordanWignerModes"/>.</item>
/// </list>
///
/// <para>This is the simplest sector family of the Prosen / Medvedyeva-Essler-Prosen
/// (2016, [arXiv:1606.09122](https://arxiv.org/abs/1606.09122)) third-quantization /
/// imaginary-Hubbard programme: subset sums of N rapidities. The (m, m̃) sectors with
/// m, m̃ ≥ 1 require the imaginary-U Bethe ansatz (not yet implemented in the repo); the
/// (m, 0) and (0, m̃) families are the closed-form leaves of that tree.</para>
///
/// <para><b>Combinatorial coverage at N=10.</b> Summing <c>C(N, m)</c> for m = 0..N gives
/// 2^N = 1024 eigenvalues per ket-side popcount family. With the complex-conjugate sister
/// family (m, 0) — eigenvalues <c>+iε_S − 2γm</c> from <c>−i[H, |0⟩⟨ψ|] = +i|0⟩⟨ψ|H</c> —
/// and the F1 mirrors at (N, N−m) and (N−m, N) which inherit the same closed form by F1
/// closure, four 2^N-sized eigenvalue families are accessible analytically: ~4096 of the
/// ~1M total eigenvalues at N=10, concentrated in the slow- and fast-decay extremes.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Closed-form algebraic
/// derivation: dissipator becomes constant on (0, m), unitary part is free-fermion, sums
/// over <c>C(N, m)</c> subsets enumerate the sector. Validated against
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> dense Evd at small N (multiset
/// match within 1e-10) in the test suite.</para>
///
/// <para>Anchor: <see cref="XyJordanWignerModes"/> (dispersion / sine modes) +
/// <see cref="F1.F1PalindromeIdentity"/> (the F1 closure that pairs (0, m) with (N, N−m));
/// MEP 2016 arXiv:1606.09122 (the imaginary-Hubbard programme this is the closed leaf of).</para>
/// </summary>
public sealed class OneSidedSectorClosedForm : Claim
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }

    /// <summary>Ket-side popcount; the sector is <c>(p_c = 0, p_r = M)</c>.</summary>
    public int M { get; }

    /// <summary>N Prosen rapidities <c>β_k = −2γ − i·ε_k</c> for k = 1..N (stored
    /// 0-indexed). The complete (0, m) sector eigenvalue family is the subset sums
    /// <c>{Σ_{k∈S} β_k : S ⊆ {1..N}, |S| = m}</c>.</summary>
    public Complex[] Rapidities { get; }

    /// <summary>Closed-form sector eigenvalues <c>λ_S = Σ_{k∈S} β_k</c> for all
    /// <c>S ⊆ {1..N}</c> with <c>|S| = M</c>, ordered by the canonical
    /// lexicographic subset enumeration. Length equals
    /// <c>C(N, M) = </c><see cref="SectorDim"/>.</summary>
    public Complex[] Eigenvalues { get; }

    public int SectorDim => Eigenvalues.Length;

    public static OneSidedSectorClosedForm Build(int N, int m, double J = 1.0, double gamma = 0.05)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (m < 0 || m > N) throw new ArgumentOutOfRangeException(nameof(m), m, $"m must be in [0, {N}].");
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");

        var modes = XyJordanWignerModes.Build(N, J);
        var beta = new Complex[N];
        for (int k = 1; k <= N; k++)
            beta[k - 1] = new Complex(-2.0 * gamma, -modes.Dispersion[k - 1]);

        var eigs = new List<Complex>(capacity: (int)Binomial(N, m));
        foreach (var subset in EnumerateSubsets(N, m))
        {
            Complex sum = Complex.Zero;
            foreach (var kIdx in subset) sum += beta[kIdx];
            eigs.Add(sum);
        }

        return new OneSidedSectorClosedForm(N, J, gamma, m, beta, eigs.ToArray());
    }

    /// <summary>Enumerate all m-element subsets of {0, 1, ..., N−1} in lexicographic order.
    /// Yields each subset as an int[] of length m (allocated fresh per yield — caller may
    /// retain). Empty subset yielded once if m = 0.</summary>
    private static IEnumerable<int[]> EnumerateSubsets(int N, int m)
    {
        if (m == 0) { yield return Array.Empty<int>(); yield break; }
        var idx = new int[m];
        for (int i = 0; i < m; i++) idx[i] = i;
        while (true)
        {
            yield return (int[])idx.Clone();
            int j = m - 1;
            while (j >= 0 && idx[j] == N - m + j) j--;
            if (j < 0) yield break;
            idx[j]++;
            for (int k = j + 1; k < m; k++) idx[k] = idx[k - 1] + 1;
        }
    }

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        k = Math.Min(k, n - k);
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    private OneSidedSectorClosedForm(int n, double j, double gamma, int m,
        Complex[] rapidities, Complex[] eigenvalues)
        : base($"One-sided sector closed form L_(0,{m}) on chain XY + uniform Z-dephasing " +
               $"(N={n}, J={j:G3}, γ={gamma:G3}): {eigenvalues.Length} eigenvalues as subset sums " +
               $"of N={n} Prosen rapidities β_k = −2γ − i·ε_k.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs (sine-mode dispersion) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/PerBlockLiouvillianBuilder.cs (per-block reference) + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (F1 closure that pairs (0,m) ↔ (N,N−m)); " +
               "Medvedyeva, Essler, Prosen 2016 arXiv:1606.09122 (imaginary-Hubbard programme, closed leaf).")
    {
        N = n; J = j; Gamma = gamma; M = m;
        Rapidities = rapidities;
        Eigenvalues = eigenvalues;
    }

    public override string DisplayName =>
        $"L_(0,{M}) closed form (N={N}, J={J:G3}, γ={Gamma:G3})";

    public override string Summary =>
        $"{SectorDim} eigenvalues, |Re λ|={Math.Abs(Eigenvalues[0].Real):F3} (constant), " +
        $"Im range [{Eigenvalues.Min(e => e.Imaginary):F3}, {Eigenvalues.Max(e => e.Imaginary):F3}] " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("J", J, "G4");
            yield return InspectableNode.RealScalar("γ", Gamma, "G4");
            yield return InspectableNode.RealScalar("m (ket popcount)", M);
            yield return InspectableNode.RealScalar("sector dim C(N,m)", SectorDim);
            for (int k = 0; k < Rapidities.Length; k++)
            {
                var b = Rapidities[k];
                yield return new InspectableNode($"β_{k + 1}",
                    summary: $"({b.Real:F5}, {b.Imaginary:F5}) = −2γ − i·ε_{k + 1}");
            }
        }
    }
}
