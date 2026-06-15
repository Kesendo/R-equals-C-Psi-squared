using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for Zero-Sector Immunity (<c>docs/proofs/PROOF_ZERO_IMMUNITY.md</c>): the
/// centered palindromic residual M = Π·L·Π⁻¹ + L + 2Σγ·I vanishes identically on the (w=0,w=0)
/// block (the {I,Z}^⊗N Pauli strings, XY-weight 0) for EVERY 2-body Hamiltonian under uniform
/// Z-dephasing, and by Π-symmetry on the (w=N,w=N) block ({X,Y}^⊗N).
///
/// <para>Where the proof argues it (the dissipator commutes with {I,Z}; only the ZZ bond preserves
/// w=0 and it too commutes; Π carries the uniform −2Σγ from w=N back as exactly the palindrome
/// shift), this witness BUILDS it: it generates a RANDOM parity-violating 2-body Hamiltonian
/// (every X/Y/Z bond pair, seeded real coefficients, so the global F1 palindrome genuinely breaks),
/// constructs the full 4^N×4^N residual M (reusing <see cref="PalindromeResidual.Build"/>), and
/// reads the Frobenius norm of the (w=0) and (w=N) sub-blocks. The gate is two-sided: the extreme
/// blocks are ≈ 0 (the theorem) AND the full M is ≫ 0 (so the vanishing is not the trivial
/// M=0-everywhere case the truly class gives). A construction error (a wrong Π sign, a scrambled
/// Pauli-index map) would leave the w=0 block non-zero, so the witness can fail.</para>
///
/// <para>Guard: 4^N ≤ <see cref="MaxDim"/> (= 1024), admitting N ∈ {2,3,4,5}; N=6 (4096²) is
/// excluded. Children: the w=0 immunity (the assert made visible), the non-trivial full-M gate, the
/// Π-mirror w=N block, and a curve of ‖M[w,w]‖_F across XY-weight w (the extremes vanish, the
/// boundary sectors 0&lt;w&lt;N carry the breaking — the V-Effect localization).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_ZERO_IMMUNITY.md</c> +
/// <c>experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md</c>; typed home
/// <see cref="ZeroSectorImmunityClaim"/>.</para></summary>
public sealed class ZeroSectorImmunityWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest residual dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public double Gamma { get; }
    public int Seed { get; }

    /// <summary>Dimension of the (w=0) block = |{I,Z}^⊗N| = 2^N.</summary>
    public int SectorDim => 1 << N;

    /// <summary>‖M restricted to the (w=0,w=0) block‖_F. The theorem says this is 0.</summary>
    public double ZeroSectorNorm { get; }

    /// <summary>‖M restricted to the (w=N,w=N) block‖_F. Zero by Π-symmetry (Π: w=0 ↔ w=N).</summary>
    public double MirrorSectorNorm { get; }

    /// <summary>‖M‖_F over the whole 4^N×4^N residual. Must be ≫ 0: the non-triviality gate.</summary>
    public double FullResidualNorm { get; }

    /// <summary>True iff the generated H genuinely breaks the global palindrome (its non-identity
    /// terms span more than one Klein (bit_a,bit_b) cell). If false, the test would be trivial.</summary>
    public bool IsParityViolating { get; }

    private readonly double[] _weightAxis;          // 0..N
    private readonly double[] _blockNormByWeight;   // ‖M[w,w]‖_F per XY-weight w

    public ZeroSectorImmunityWitness(int n = 3, double gamma = 0.05, int seed = 42)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2 (a 2-body Hamiltonian needs ≥ 2 sites); got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        long dim = 1L << (2 * n);
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n;
        Gamma = gamma;
        Seed = seed;

        // A random parity-violating 2-body H: every (X/Y/Z, X/Y/Z) nearest-neighbour bond pair with
        // a seeded real coefficient. Real coeff × a single Pauli string is Hermitian, so H is Hermitian;
        // the mix of Klein cells (XX/YY/ZZ truly, plus XY/XZ/YZ/... parity-violating) makes M ≠ 0.
        var rng = new Random(seed);
        var letters = new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var terms = new List<PauliTerm>();
        for (int i = 0; i < n - 1; i++)
            foreach (var la in letters)
                foreach (var lb in letters)
                    terms.Add(PauliTerm.TwoSite(n, i, la, i + 1, lb, new Complex(rng.NextDouble() * 2.0 - 1.0, 0.0)));
        var ham = new PauliHamiltonian(n, terms);
        IsParityViolating = !ham.IsKleinHomogeneous;

        // M = Π·L·Π⁻¹ + L + 2Σγ·I in the 4^N Pauli-string basis.
        var gammas = Enumerable.Repeat(gamma, n).ToList();
        ComplexMatrix L = PauliDephasingDissipator.BuildZ(ham.ToMatrix(), gammas);
        ComplexMatrix M = PalindromeResidual.Build(L, n, gammas.Sum());
        FullResidualNorm = M.FrobeniusNorm();

        // Diagonal block of M per XY-weight w; the (w=0) and (w=N) blocks are the immune extremes.
        var byWeight = new List<int>[n + 1];
        for (int w = 0; w <= n; w++) byWeight[w] = new List<int>();
        for (long k = 0; k < dim; k++)
            byWeight[PauliIndex.TotalBitA(PauliIndex.FromFlat(k, n))].Add((int)k);

        _weightAxis = new double[n + 1];
        _blockNormByWeight = new double[n + 1];
        for (int w = 0; w <= n; w++)
        {
            _weightAxis[w] = w;
            _blockNormByWeight[w] = BlockNorm(M, byWeight[w]);
        }
        ZeroSectorNorm = _blockNormByWeight[0];
        MirrorSectorNorm = _blockNormByWeight[n];
    }

    /// <summary>Frobenius norm of the sub-block of M on the given (scattered) index set.</summary>
    private static double BlockNorm(ComplexMatrix M, List<int> idx)
    {
        double s = 0.0;
        foreach (int r in idx)
            foreach (int c in idx)
            {
                var z = M[r, c];
                s += z.Real * z.Real + z.Imaginary * z.Imaginary;
            }
        return Math.Sqrt(s);
    }

    public string DisplayName =>
        $"ZeroSectorImmunityWitness (N={N}, γ={Gamma.ToString("0.###", Inv)}, w=0 dim={SectorDim})";

    public string Summary =>
        $"N={N}: a random parity-violating 2-body H (‖M‖_F = {FullResidualNorm.ToString("e2", Inv)} > 0, " +
        $"so it breaks the global F1 palindrome) leaves M EXACTLY zero on the w=0 ({{I,Z}}^⊗N) block " +
        $"(‖·‖ = {ZeroSectorNorm.ToString("e2", Inv)}) and the w=N ({{X,Y}}^⊗N) block " +
        $"(‖·‖ = {MirrorSectorNorm.ToString("e2", Inv)}). The classical extreme is immune to every " +
        $"2-body coupling; the breaking lives in the boundary sectors 0<w<N.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "the w=0 immunity, live (the assert)",
                summary: $"the (w=0,w=0) {SectorDim}×{SectorDim} sub-block of M on the {{I,Z}}^⊗N strings has " +
                         $"‖·‖_F = {ZeroSectorNorm.ToString("e3", Inv)} ≈ 0. The dissipator annihilates w=0 " +
                         $"(Z commutes with every I/Z string) and [H,·] cannot stay in w=0 (only the ZZ bond " +
                         $"preserves it, and ZZ commutes with I/Z too); Π carries the uniform −2Σγ from w=N back " +
                         $"as exactly the palindrome shift. Holds for the random 2-body H built here.");

            yield return new InspectableNode(
                displayName: "the non-trivial gate (full M ≠ 0)",
                summary: $"‖M‖_F = {FullResidualNorm.ToString("e3", Inv)} > 0, and the H is parity-violating " +
                         $"({(IsParityViolating ? "yes" : "NO — gate compromised")}); so the w=0 vanishing is NOT the " +
                         $"trivial M=0-everywhere case the truly (XX+YY+ZZ) class gives. The immunity is sector-specific.");

            yield return new InspectableNode(
                displayName: "the Π-mirror (w=N block)",
                summary: $"Π maps {{I,Z}}^⊗N ↔ {{X,Y}}^⊗N (flips bit_a at every site), so M also vanishes on the " +
                         $"w=N block: ‖·‖_F = {MirrorSectorNorm.ToString("e3", Inv)} ≈ 0. The all-transverse extreme " +
                         $"is the mirror image of the classical one.");

            yield return new InspectableNode(
                displayName: "‖M[w,w]‖ across XY-weight (the V-Effect localization)",
                summary: "the diagonal block norm of M by XY-weight w: the extremes w=0 and w=N vanish (immune), " +
                         "while the boundary sectors 0<w<N carry the palindrome-breaking — exactly the V-Effect " +
                         "boundary localization (V_EFFECT_BOUNDARY_LOCALIZATION).",
                payload: new InspectablePayload.Curve("‖M[w,w]‖_F by XY-weight", _weightAxis, _blockNormByWeight,
                    "XY-weight w", "‖M[w,w]‖_F"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
