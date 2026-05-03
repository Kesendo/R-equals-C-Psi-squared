using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F80;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Polynomial → Half → Mirror trio walked as a single chain at a given
/// chain configuration. Captures the trio's computable manifestations as a callable
/// diagnostic so the chain can be inspected on real chain systems and Hamiltonians.
///
/// <para><b>The trio</b> (named in <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>):
/// <list type="bullet">
///   <item><c>PolynomialFoundationClaim</c>: d²−2d=0 ↔ R=CΨ² selects the minimum-memory
///         dimension; d=1 algebraically excluded; solutions {0 = axis, 2 = count of
///         ±0.5 positions}.</item>
///   <item><c>HalfAsStructuralFixedPointClaim</c>: 1/2 number-anchor (three faces close at
///         d=2: V-Effect bridge, observation horizon ρ_mm = I/2, qubit dimensional anchor).</item>
///   <item><c>NinetyDegreeMirrorMemoryClaim</c>: 90° angle-anchor (the i in F80's
///         Spec(M) = ±2i · Spec(H_non-truly), the rotation back onto the mirror).</item>
/// </list></para>
///
/// <para><b>The four steps walked</b>:
/// <list type="number">
///   <item>Π eigenvalue partition: 4^N Pauli strings split into Π²-even (eigenvalue +1)
///         and Π²-odd (eigenvalue −1), each subspace of size 4^N / 2 (F88 bilinear apex 1/2).</item>
///   <item>Π scaled by 1/2: spectrum {+1, −1} shifted toward 0 by 0.5 = the structural
///         pair {+0.5, −0.5} on the d = 2 memory axis. Static anchor; see
///         <see cref="MemoryAxisChainResult.StructuralPair"/>.</item>
///   <item>Bloch realisation: at the qubit level (d = 2), pure-state populations differ
///         from the maximally mixed by ±0.5 (the (1 ± r_z)/2 form). Same static pair as
///         step 2; the 1/2 in ρ = (I + r·σ)/2 IS Π's eigenvalue transport onto the
///         memory axis.</item>
///   <item>F80 lift: chain Π²-odd Hamiltonian eigenvalues λ on the real axis lift to ±2iλ
///         on M's imaginary axis (the i factor is the 90° rotation, the 2 is chiral
///         pair-doubling). Bit-exact verified at N = 3..7 + k-body in
///         <c>docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md</c>.</item>
/// </list></para>
///
/// <para>Steps 1 and 4 produce data that depends on (N, H); steps 2 and 3 are static
/// anchors. <see cref="MemoryAxisChainResult.BilinearApexHolds"/> and
/// <see cref="MemoryAxisChainResult.F80SpectrumIsMirrorSymmetric"/> are the trio's
/// falsifiable structural assertions made directly checkable.</para>
/// </summary>
public sealed record MemoryAxisChainResult(
    int N,
    long Pi2EvenStringCount,
    long Pi2OddStringCount,
    IReadOnlyDictionary<double, int>? F80ImaginarySpectrum)
{
    /// <summary>F88 bilinear apex 1/2: the Π²-even and Π²-odd subspaces have equal
    /// cardinality (4^N / 2 each). Combinatorial identity, holds for any N ≥ 1.</summary>
    public bool BilinearApexHolds => Pi2EvenStringCount == Pi2OddStringCount;

    /// <summary>F80 chiral pair-doubling: every imaginary value ν in the spectrum has its
    /// mirror −ν with equal multiplicity. The 2 in 2i is exactly this pair-doubling
    /// (H's particle-hole symmetry lifted to M's imaginary axis). Returns true if no
    /// Hamiltonian was provided (vacuously satisfied).</summary>
    public bool F80SpectrumIsMirrorSymmetric =>
        F80ImaginarySpectrum == null
        || F80ImaginarySpectrum.All(kv =>
            F80ImaginarySpectrum.TryGetValue(Math.Round(-kv.Key, 10), out var mirrorMult)
            && mirrorMult == kv.Value);

    /// <summary>The static memory-axis pair {+0.5, −0.5}: Π's eigenvalues scaled by 1/2,
    /// Bloch deviation of pure-state populations from mixed I/2. The d = 2 manifestation
    /// of the polynomial's count-of-positions solution. Same pair appears at every layer.</summary>
    public static IReadOnlyList<double> StructuralPair { get; } = new[] { +0.5, -0.5 };
}

/// <summary>Compute the trio chain at a given chain configuration. The result captures
/// the data of step 1 (Π eigenspace counts) and step 4 (F80 imaginary spectrum, if a
/// chain Π²-odd Hamiltonian is provided), with the static ±0.5 anchor reachable via
/// <see cref="MemoryAxisChainResult.StructuralPair"/>.</summary>
public static class MemoryAxisChain
{
    /// <summary>Walk the trio chain. If <paramref name="hamiltonianTerms"/> is provided
    /// and contains chain Π²-odd 2-body bilinears, the F80 lift is computed via
    /// <see cref="BlochSignWalk.PredictMSpectrumImaginaryParts"/>; otherwise the F80
    /// spectrum field is null. <paramref name="dephaseLetter"/> selects which Π² axis
    /// defines the even/odd split (default Z, matching standard Z-dephasing).</summary>
    public static MemoryAxisChainResult Compute(
        ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm>? hamiltonianTerms = null,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        var (evenCount, oddCount) = Pi2Projection.Counts(chain.N, dephaseLetter);

        IReadOnlyDictionary<double, int>? f80Spec = null;
        if (hamiltonianTerms != null && hamiltonianTerms.Count > 0)
            f80Spec = BlochSignWalk.PredictMSpectrumImaginaryParts(chain, hamiltonianTerms);

        return new MemoryAxisChainResult(
            N: chain.N,
            Pi2EvenStringCount: evenCount,
            Pi2OddStringCount: oddCount,
            F80ImaginarySpectrum: f80Spec);
    }
}
