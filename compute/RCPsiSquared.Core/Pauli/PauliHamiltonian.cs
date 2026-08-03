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

    /// <summary>True if all non-identity terms share the same Klein index. This does NOT
    /// settle the F87 class. Among the 21 content-bearing Klein-homogeneous two-letter pairs at N=3 the 6
    /// identity-free ones are 0/6 hard, but the 15 carrying an identity are 5/15 hard
    /// (IZ+ZI is Klein (0,1) and hard at N=3 and N=4), and at k=3 the 294-pair
    /// Z₂³-homogeneous sweep has 50 hard in the Klein-(0,1) cell under Z-dephasing.
    /// Among Klein-HOMOGENEOUS pairs, sitting in the dephase letter's own cell is
    /// necessary for hardness but not sufficient: XY+YX sits in Z's cell (0,1) and is
    /// soft. Without homogeneity not even necessary: 10 hard two-letter pairs at N=3
    /// have neither term in that cell. Use
    /// <c>Diagnostics.F87.PauliPairTrichotomy.Classify</c> for the verdict.</summary>
    public bool IsKleinHomogeneous => KleinSet.Count <= 1;

    /// <summary>Set of distinct full Z₂³ signatures (bit_a, bit_b, Y-par) across non-identity
    /// terms. Within ONE parity of <see cref="PauliTerm.KBody"/>, either one, this has the same
    /// cardinality as <see cref="KleinSet"/>; as soon as the terms differ in KBody parity
    /// it can be strictly finer (Y-parity becomes independent).</summary>
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
    /// Hamiltonian.
    ///
    /// <para>NO GAUGE TERM HERE, and the repository has a second builder that carries one:
    /// <c>Diagnostics.BlockSpectrumWitness.HeisenbergGraph</c> writes Σ_b (J_b/4)(XX+YY+ZZ−I).
    /// The −I is the occupation-versus-spin convention — per bond,
    /// (J/4)(ZZ−I) = J·n_a n_b − (J/2)(n_a + n_b) with n = (I−Z)/2 — so the two builders differ by
    /// the constant Σ_b J_b/4 and by nothing else. This one has tr H = 0, so the MEAN eigenvalue is 0
    /// (not a palindrome — the XXX spectrum is not symmetric about it), and the ferromagnetic vacuum
    /// sits at +Σ_b J_b/4; that one puts the vacuum at exactly 0. What that buys is not
    /// block-diagonality, which is a symmetry of any number-conserving H and holds either way, but the
    /// BIT-exactness of the witness's (0,1) generator identity: 6.0e-08 without the gauge at
    /// J = π·10⁸, exactly 0.0 with it. The two agree as mathematical objects ([c·I, ρ] = 0) but not as
    /// computed ones, which is the whole reason the gauge is there; anything reading H's own diagonal
    /// must know which side it is on.</para></summary>
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
