using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>An N-qubit Hamiltonian as a sum of weighted Pauli-string terms.
///
/// The dense 2^N × 2^N matrix is built lazily via <see cref="ToMatrix"/>. For typical
/// project N (≤ 8), this is up to 256 × 256 — tractable.
/// </summary>
public sealed record PauliHamiltonian(int N, IReadOnlyList<PauliTerm> Terms)
{
    /// <summary>Set of distinct Klein indices (bit_a, bit_b) across all non-identity terms.
    /// Identity terms are excluded (their Klein index (0,0) is trivial).</summary>
    public IReadOnlySet<(int BitA, int BitB)> KleinSet =>
        Terms.Where(t => t.KBody > 0).Select(t => t.KleinIndex).ToHashSet();

    /// <summary>True if all non-identity terms share the same Klein index. Empirical
    /// structural fact (k=2 full enumeration; k=3 sample): Klein-homogeneous Hamiltonians
    /// under matched-Klein dephasing are always F87 soft or truly, never F87 hard.</summary>
    public bool IsKleinHomogeneous => KleinSet.Count <= 1;

    /// <summary>Set of distinct full Z₂³ signatures (bit_a, bit_b, Y-par) across non-identity
    /// terms. At k=2 this has the same cardinality as <see cref="KleinSet"/>; at k≥3 it
    /// can be strictly finer (Y-parity becomes independent).</summary>
    public IReadOnlySet<(int BitA, int BitB, int YParity)> FullZ2SignatureSet =>
        Terms.Where(t => t.KBody > 0).Select(t => t.FullZ2Signature).ToHashSet();

    /// <summary>True if all non-identity terms share the same full Z₂³ signature. Strictly
    /// finer than <see cref="IsKleinHomogeneous"/> at k≥3; equivalent at k=2.</summary>
    public bool IsZ2Homogeneous => FullZ2SignatureSet.Count <= 1;

    /// <summary>Per-term Klein index list (in term order, non-identity terms only).</summary>
    public IReadOnlyList<(int BitA, int BitB)> PerTermKleinIndices =>
        Terms.Where(t => t.KBody > 0).Select(t => t.KleinIndex).ToList();

    /// <summary>Per-term full Z₂³ signature list (in term order, non-identity terms only).</summary>
    public IReadOnlyList<(int BitA, int BitB, int YParity)> PerTermFullZ2Signatures =>
        Terms.Where(t => t.KBody > 0).Select(t => t.FullZ2Signature).ToList();

    public ComplexMatrix ToMatrix()
    {
        int d = 1 << N;
        var H = Matrix<Complex>.Build.Dense(d, d);
        foreach (var term in Terms)
        {
            if (term.N != N)
                throw new ArgumentException($"term has N={term.N} letters; expected N={N}");
            H = H + term.Coefficient * PauliString.Build(term.Letters);
        }
        return H;
    }

    /// <summary>Uniform XY chain: H = (J/2) Σ_b (X_b X_{b+1} + Y_b Y_{b+1}). Throws
    /// <see cref="ArgumentOutOfRangeException"/> for N &lt; 1; the per-bond overload it
    /// forwards to also guards N &lt; 1 as defence-in-depth (the bare
    /// <c>new double[N - 1]</c> below would otherwise throw <see cref="OverflowException"/>
    /// for N = 0, which is a confusing surface for an invalid-argument condition).</summary>
    public static PauliHamiltonian XYChain(int N, double J)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        var bondJ = new double[N - 1];
        for (int b = 0; b < N - 1; b++) bondJ[b] = J;
        return XYChain(N, bondJ);
    }

    /// <summary>Non-uniform XY chain with per-bond coupling: H = Σ_b (J_b/2) (X_b X_{b+1} + Y_b Y_{b+1}).
    /// <paramref name="bondJ"/> must have length N − 1 (one coupling per nearest-neighbour bond).
    /// Used by per-bond J Builder paths in <c>BlockSpectrum/</c> for F100-territory experiments
    /// (palindromic J profiles, etc.). Scalar overload <see cref="XYChain(int, double)"/> calls
    /// this with a uniform list. Throws <see cref="ArgumentOutOfRangeException"/> for N &lt; 1
    /// before the length check so that N = 0 / empty bondJ does not slip through as a silent
    /// no-op Hamiltonian.</summary>
    public static PauliHamiltonian XYChain(int N, IReadOnlyList<double> bondJ)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        if (bondJ.Count != N - 1)
            throw new ArgumentException(
                $"bondJ length {bondJ.Count} != N - 1 = {N - 1}", nameof(bondJ));
        var terms = new List<PauliTerm>(2 * (N - 1));
        for (int b = 0; b < N - 1; b++)
        {
            Complex c = bondJ[b] / 2.0;
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
        }
        return new PauliHamiltonian(N, terms);
    }

    /// <summary>Uniform Heisenberg chain: H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}).
    /// Throws <see cref="ArgumentOutOfRangeException"/> for N &lt; 1 (the bare
    /// <c>new double[N - 1]</c> would otherwise raise <see cref="OverflowException"/> for
    /// N = 0, a confusing surface for an invalid-argument condition).</summary>
    public static PauliHamiltonian HeisenbergChain(int N, double J)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        var bondJ = new double[N - 1];
        for (int b = 0; b < N - 1; b++) bondJ[b] = J;
        return HeisenbergChain(N, bondJ);
    }

    /// <summary>Non-uniform Heisenberg chain with per-bond coupling:
    /// H = Σ_b (J_b/4) (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}). <paramref name="bondJ"/> must
    /// have length N − 1. Scalar overload <see cref="HeisenbergChain(int, double)"/> calls this
    /// with a uniform list. Throws <see cref="ArgumentOutOfRangeException"/> for N &lt; 1 before
    /// the length check so that N = 0 / empty bondJ does not slip through as a silent no-op
    /// Hamiltonian.</summary>
    public static PauliHamiltonian HeisenbergChain(int N, IReadOnlyList<double> bondJ)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        if (bondJ.Count != N - 1)
            throw new ArgumentException(
                $"bondJ length {bondJ.Count} != N - 1 = {N - 1}", nameof(bondJ));
        var terms = new List<PauliTerm>(3 * (N - 1));
        for (int b = 0; b < N - 1; b++)
        {
            Complex c = bondJ[b] / 4.0;
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Z, c));
        }
        return new PauliHamiltonian(N, terms);
    }

    /// <summary>Bilinear bond Hamiltonian H = Σ_bond Σ_term coeff · σ_la^i σ_lb^j on the given bonds.</summary>
    public static PauliHamiltonian Bilinear(int N, IReadOnlyList<Bond> bonds,
        IReadOnlyList<(PauliLetter La, PauliLetter Lb, Complex Coeff)> terms)
    {
        var allTerms = new List<PauliTerm>(bonds.Count * terms.Count);
        foreach (var bond in bonds)
            foreach (var (la, lb, coeff) in terms)
                allTerms.Add(PauliTerm.TwoSite(N, bond.Site1, la, bond.Site2, lb, coeff));
        return new PauliHamiltonian(N, allTerms);
    }
}
