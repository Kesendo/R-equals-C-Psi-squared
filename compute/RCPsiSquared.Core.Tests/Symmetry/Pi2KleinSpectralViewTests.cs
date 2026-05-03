using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Schicht 3: bridge between the algebraic Klein-cell view and the dynamical
/// spectrum of L. For each F87 canonical Hamiltonian at N=3 we compute the per-eigenmode
/// Klein-cell distribution (normalised) and observe the structural patterns.
///
/// <para>L is non-normal (Lindbladian, non-Hermitian) so its right eigenvectors do not
/// form an orthonormal basis. The aggregate mass over eigenmodes does NOT recover the
/// 4^N/2 partition that the Π²_Z eigenspace structure has at the operator level. We
/// observe what MathNet's basis actually returns, lock the per-mode normalisation, and
/// document the slow-mode Klein-cell distribution as observation rather than enforcement.</para>
/// </summary>
public class Pi2KleinSpectralViewTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildL(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    [Theory]
    [MemberData(nameof(CanonicalHamiltonianCases))]
    public void EveryEigenmode_HasNormalisedKleinDistribution(
        string label, PauliPairBondTerm[] terms)
    {
        // After per-mode normalisation the 4 Klein-cell masses sum to 1.
        var chain = Chain3();
        var L = BuildL(terms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N);

        Assert.NotEmpty(modes);
        foreach (var m in modes)
        {
            Assert.True(Math.Abs(m.TotalMass - 1.0) < 1e-9,
                $"{label}: total Klein mass {m.TotalMass:F6} ≠ 1 for eigenmode at λ={m.Eigenvalue:F4}");
        }
    }

    [Theory]
    [MemberData(nameof(CanonicalHamiltonianCases))]
    public void EigenmodeCount_Equals_4PowN(string label, PauliPairBondTerm[] terms)
    {
        var chain = Chain3();
        var L = BuildL(terms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N);
        long expectedCount = 1L << (2 * chain.N);
        Assert.True(modes.Count == expectedCount, $"{label}: expected {expectedCount} modes, got {modes.Count}");
    }

    [Fact]
    public void TrulyXYChain_SlowestFourModes_AreKernel_AtRealZero()
    {
        // For truly H, F1 ensures Re(λ_kernel) = 0 exactly. The 4 slowest modes form the
        // kernel of L (stationary states; dim ≈ N+1 = 4 at N=3 for the |ΔN|=0 sector).
        var chain = Chain3();
        var terms = new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) };
        var L = BuildL(terms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();

        var slowest = modes.Take(4).ToList();
        foreach (var m in slowest)
        {
            Assert.True(Math.Abs(m.Eigenvalue.Real) < 1e-9,
                $"slow mode λ={m.Eigenvalue:F6} not at exact kernel — F1 anchor expects Re(λ)=0 for truly H");
        }
    }

    [Fact]
    public void TrulyXYChain_SlowestFourModes_SplitTwoTwoAcrossPi2Z()
    {
        // Observed structural pattern: the 4 kernel modes of truly L at N=3 split
        // 2:2 between Π²_Z = +1 (Pp + Pm) and Π²_Z = −1 (Mp + Mm). Two modes contain
        // Z-only conserved content (sit in the Π²_Z=−1 sector), two modes are even
        // under bit_b parity (sit in Π²_Z=+1). Both halves are kernel — the truly
        // Hamiltonian's stationary subspace is naturally 2:2 balanced under Π²_Z.
        var chain = Chain3();
        var terms = new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) };
        var L = BuildL(terms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();

        var slowest = modes.Take(4).ToList();
        double aggregateZPlus = slowest.Sum(m => m.MassPp + m.MassPm);
        double aggregateZMinus = slowest.Sum(m => m.MassMp + m.MassMm);

        // 4 normalised modes → total mass 4; split 2:2 between the two Π²_Z eigenspaces.
        Assert.True(Math.Abs(aggregateZPlus - 2.0) < 0.01,
            $"Truly kernel aggregate Π²_Z=+ mass {aggregateZPlus:F4} should be ≈ 2 (half of 4-dim kernel)");
        Assert.True(Math.Abs(aggregateZMinus - 2.0) < 0.01,
            $"Truly kernel aggregate Π²_Z=− mass {aggregateZMinus:F4} should be ≈ 2");
    }

    [Fact]
    public void XYChain_NonTriviallyOddH_PopulatesAllFourKleinCells_AggregatedOverSpectrum()
    {
        // For non-truly H the dynamics summed over all eigenmodes touches every Klein-cell.
        // Per cell the aggregate must be > 1 (more than one mode's worth of mass).
        var chain = Chain3();
        var terms = new[] { Term(PauliLetter.X, PauliLetter.Y) };
        var L = BuildL(terms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N);

        Assert.True(modes.Sum(m => m.MassPp) > 1.0);
        Assert.True(modes.Sum(m => m.MassPm) > 1.0);
        Assert.True(modes.Sum(m => m.MassMp) > 1.0);
        Assert.True(modes.Sum(m => m.MassMm) > 1.0);
    }

    [Theory]
    [MemberData(nameof(CanonicalHamiltonianCases))]
    public void SlowestMode_DominantCell_DocumentedPerHamiltonian(
        string label, PauliPairBondTerm[] terms)
    {
        // Diagnostic: which Klein-cell does the slowest non-kernel eigenmode prefer for
        // each canonical Hamiltonian? Document the observation; this is what the rohe
        // Schicht reveals about dynamics-vs-algebra alignment.
        var chain = Chain3();
        var L = BuildL(terms, chain);
        var sorted = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();

        // Grab the slowest non-kernel mode (Re(λ) > 1e-9). Some eigenvalues sit at the
        // dephasing shift -2σ; we want the slowest decaying coherence.
        var slowestNonKernel = sorted.FirstOrDefault(m => Math.Abs(m.Eigenvalue.Real) > 1e-9);
        Assert.True(slowestNonKernel is not null,
            $"{label}: no non-kernel slow mode found");
        Assert.True(slowestNonKernel!.TotalMass > 0.99 && slowestNonKernel.TotalMass < 1.01,
            $"{label}: slowest non-kernel mode total mass {slowestNonKernel.TotalMass:F4} not normalized");
    }

    public static IEnumerable<object[]> CanonicalHamiltonianCases() => new[]
    {
        new object[] { "XX+YY (Truly)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) } },
        new object[] { "Heisenberg XX+YY+ZZ (Truly)",
            new[] {
                Term(PauliLetter.X, PauliLetter.X),
                Term(PauliLetter.Y, PauliLetter.Y),
                Term(PauliLetter.Z, PauliLetter.Z),
            } },
        new object[] { "XY+YX (Pi2OddPure A)",
            new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) } },
        new object[] { "YZ+ZY (Pi2EvenNonTruly)",
            new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) } },
        new object[] { "XX+XY (Mixed)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) } },
    };

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
