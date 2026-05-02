using System.Numerics;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.ChainSystems;

/// <summary>An N-qubit chain with uniform Z-dephasing γ₀ at each site and a chosen H + topology.
///
/// <see cref="Bonds"/> is canonicalised at construction: derived from <see cref="Topology"/>
/// when no explicit list is given, otherwise the explicit list is used. Record equality
/// therefore compares the actual bond layout, not a sentinel.
///
/// <para><b>Convenience builders</b>: <see cref="BuildHamiltonian"/> and <see cref="BuildLiouvillian"/>
/// produce the dense matrices used by Lindblad / Symmetry / Observables. They build fresh on each
/// call (no caching) — the caller can hold the returned matrix.</para>
/// </summary>
public sealed record ChainSystem
{
    public int N { get; }
    public double J { get; }
    public double GammaZero { get; }
    public HamiltonianType HType { get; }
    public TopologyKind Topology { get; }
    public IReadOnlyList<Bond> Bonds { get; }

    public ChainSystem(int N, double J, double GammaZero,
        HamiltonianType HType = HamiltonianType.XY,
        TopologyKind Topology = TopologyKind.Chain)
    {
        this.N = N;
        this.J = J;
        this.GammaZero = GammaZero;
        this.HType = HType;
        this.Topology = Topology;
        Bonds = DefaultBonds(N, J, Topology);
    }

    private ChainSystem(ChainSystem other, IReadOnlyList<Bond> bonds)
    {
        N = other.N;
        J = other.J;
        GammaZero = other.GammaZero;
        HType = other.HType;
        Topology = other.Topology;
        Bonds = bonds;
    }

    public ChainSystem WithBonds(IReadOnlyList<Bond> bonds) => new(this, bonds);

    /// <summary>Build the 2^N × 2^N Hamiltonian matrix for this chain. Routes through
    /// <see cref="PauliHamiltonian.Bilinear"/> using <see cref="Bonds"/>, so any topology
    /// (Chain, Star, Ring, or custom <see cref="WithBonds"/> override) flows through correctly.</summary>
    public ComplexMatrix BuildHamiltonian()
    {
        var letterTerms = HType switch
        {
            HamiltonianType.XY => new (PauliLetter, PauliLetter, System.Numerics.Complex)[]
            {
                (PauliLetter.X, PauliLetter.X, J / 2.0),
                (PauliLetter.Y, PauliLetter.Y, J / 2.0),
            },
            HamiltonianType.Heisenberg => new (PauliLetter, PauliLetter, System.Numerics.Complex)[]
            {
                (PauliLetter.X, PauliLetter.X, J / 4.0),
                (PauliLetter.Y, PauliLetter.Y, J / 4.0),
                (PauliLetter.Z, PauliLetter.Z, J / 4.0),
            },
            _ => throw new ArgumentOutOfRangeException(nameof(HType)),
        };
        return PauliHamiltonian.Bilinear(N, Bonds, letterTerms).ToMatrix();
    }

    /// <summary>Build the d² × d² full Liouvillian L = −i[H, ·] + Σ_l γ₀ (Z_l ρ Z_l − ρ) for this
    /// chain under uniform Z-dephasing γ₀. For F1 palindrome verification, see
    /// <see cref="Symmetry.PalindromeResidual"/>.</summary>
    public ComplexMatrix BuildLiouvillian()
    {
        var H = BuildHamiltonian();
        var gammaList = Enumerable.Repeat(GammaZero, N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    /// <summary>Total dephasing rate Σ_l γ₀ (uniform γ₀ across N sites). Used by
    /// <see cref="Symmetry.PalindromeResidual"/> as the σ in Π·L·Π⁻¹ + L + 2σ·I = 0.</summary>
    public double SigmaGamma => N * GammaZero;

    private static IReadOnlyList<Bond> DefaultBonds(int N, double J, TopologyKind topology) => topology switch
    {
        TopologyKind.Chain => Enumerable.Range(0, N - 1)
            .Select(i => new Bond(i, i + 1, J))
            .ToArray(),
        TopologyKind.Star => Enumerable.Range(1, N - 1)
            .Select(i => new Bond(0, i, J))
            .ToArray(),
        TopologyKind.Ring => Enumerable.Range(0, N)
            .Select(i => new Bond(i, (i + 1) % N, J))
            .ToArray(),
        _ => throw new InvalidOperationException($"Unknown topology: {topology}"),
    };
}
