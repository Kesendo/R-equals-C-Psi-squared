using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.DZero;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The trio's reading of an arbitrary density matrix ρ. Decomposes ρ first
/// along the polynomial's two solutions (d=0 axis = static, d=2 axis = dynamic) via
/// <see cref="DZeroDecomposition"/>, then within the dynamic part along the Π²-axis
/// (Π²-even = residual non-static, Π²-odd = M_anti-related, where dynamics generators
/// live) via <see cref="Pi2Projection.Split"/>. The result answers the natural question
/// the trio raises on a state: "how much of this state is static (sitting on the d=0
/// axis), how much is memory (on the d=2 axis), and where on the memory axis does it
/// sit?"
///
/// <para>Frobenius-norm partition (orthogonal under Hilbert–Schmidt):
/// ‖ρ‖²_F = ‖ρ_d0‖²_F + ‖ρ_d2‖²_F = ‖ρ_d0‖²_F + ‖ρ_d2_even‖²_F + ‖ρ_d2_odd‖²_F.</para>
///
/// <para>Limit cases:
/// <list type="bullet">
///   <item>ρ = I/d (maximally mixed): entirely static; <c>StaticFraction = 1</c>.</item>
///   <item>ρ in a single sector P_n: entirely static (kernel of L); <c>StaticFraction = 1</c>.</item>
///   <item>Pure |0…0⟩ (computational basis vacuum): mixture of static (sector P_0)
///         + memory (off-diagonal content); both fractions positive.</item>
/// </list></para>
///
/// <para>Companion to <see cref="MemoryAxisChain"/>: the Chain version walks the trio
/// at the operator level (Π eigenspaces, F80 spectrum); this Rho version walks the trio
/// at the state level. Same trio, two views.</para>
/// </summary>
public sealed record MemoryAxisRhoResult(
    int N,
    double StaticNormSquared,
    double MemoryPi2EvenNormSquared,
    double MemoryPi2OddNormSquared,
    int KernelDimension)
{
    /// <summary>Total norm² = static + memory(even) + memory(odd). Equal to ‖ρ‖²_F.</summary>
    public double TotalNormSquared =>
        StaticNormSquared + MemoryPi2EvenNormSquared + MemoryPi2OddNormSquared;

    /// <summary>Memory norm² = ‖ρ_d2‖²_F = ‖ρ‖²_F − ‖ρ_d0‖²_F.</summary>
    public double MemoryNormSquared => MemoryPi2EvenNormSquared + MemoryPi2OddNormSquared;

    /// <summary>Static fraction: ‖ρ_d0‖² / ‖ρ‖². 1.0 at any kernel-of-L state, including
    /// the maximally mixed state I/d and any sector-projector mixture.</summary>
    public double StaticFraction =>
        TotalNormSquared > 0 ? StaticNormSquared / TotalNormSquared : 0;

    /// <summary>Memory fraction: ‖ρ_d2‖² / ‖ρ‖². The dynamics-carrying part. 0 at any
    /// kernel-of-L state.</summary>
    public double MemoryFraction =>
        TotalNormSquared > 0 ? MemoryNormSquared / TotalNormSquared : 0;

    /// <summary>Within the memory part (ρ_d2), the Π²-odd fraction. This is the
    /// M_anti-related content where the unitary commutator of any Π²-odd Hamiltonian
    /// component lives (F81). 1.0 means all dynamics is from Π²-odd sources; 0 means
    /// the dynamic content is purely Π²-even (richer F80-out-of-scope cluster).</summary>
    public double Pi2OddFractionWithinMemory =>
        MemoryNormSquared > 0 ? MemoryPi2OddNormSquared / MemoryNormSquared : 0;
}

/// <summary>Decompose ρ into the trio's three structural components: static (d=0 axis),
/// memory Π²-even (d=2 axis, Π²-fixed), memory Π²-odd (d=2 axis, Π²-anti-fixed,
/// M_anti-related). All decompositions use existing primitives — DZeroDecomposition for
/// the static/memory split, Pi2Projection.Split for the Π²-axis split — so the result
/// is the trio's framing applied to existing infrastructure, not new physics.</summary>
public static class MemoryAxisRho
{
    /// <summary>Decompose ρ. Builds a fresh <see cref="StationaryModesResult"/> for the
    /// chain; if you have one cached, use the overload taking it.</summary>
    public static MemoryAxisRhoResult Decompose(
        ComplexMatrix rho,
        ChainSystem chain,
        PauliLetter dephaseLetter = PauliLetter.Z,
        double tolerance = 1e-9)
    {
        var sm = StationaryModes.Compute(chain, tolerance: tolerance);
        return Decompose(rho, chain, sm, dephaseLetter);
    }

    /// <summary>Decompose ρ using a precomputed <see cref="StationaryModesResult"/>.</summary>
    public static MemoryAxisRhoResult Decompose(
        ComplexMatrix rho,
        ChainSystem chain,
        StationaryModesResult sm,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        var dzero = DZeroDecomposition.Decompose(rho, sm);
        var (rhoEven, rhoOdd) = Pi2Projection.Split(dzero.RhoD2, chain.N, dephaseLetter);

        double staticFro = dzero.RhoD0.FrobeniusNorm();
        double evenFro = rhoEven.FrobeniusNorm();
        double oddFro = rhoOdd.FrobeniusNorm();

        return new MemoryAxisRhoResult(
            N: chain.N,
            StaticNormSquared: staticFro * staticFro,
            MemoryPi2EvenNormSquared: evenFro * evenFro,
            MemoryPi2OddNormSquared: oddFro * oddFro,
            KernelDimension: dzero.KernelDimension);
    }
}
