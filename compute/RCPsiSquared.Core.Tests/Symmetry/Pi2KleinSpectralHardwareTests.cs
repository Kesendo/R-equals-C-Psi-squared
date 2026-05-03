using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Schicht 3 against the 2026-04-30 Marrakesh f83 hardware: for each
/// (H-class, fingerprint observable) pair the framework chose, compute the slow-mode
/// Klein-cell distribution and verify the observable's Klein-cell carries non-trivial
/// dynamic mass. This is the structural pre-condition for the Marrakesh fingerprint
/// observable to give a non-zero measurement under L-evolution; the actual measurement
/// depends on initial state and propagation, but the Klein-mass overlap is a necessary
/// condition that Schicht 3 must satisfy.
///
/// <para>Setup: N=3, J=1, γ_Z=0.05 (matching the Marrakesh γ_Z_eff). For each canonical
/// Hamiltonian we sort the L-eigenmodes by |Re(λ)| (slow → fast), take the slowest 16
/// (a quarter of all 64 modes — the long-lived dynamic content), and aggregate their
/// Klein-cell mass in the observable's cell.</para>
/// </summary>
public class Pi2KleinSpectralHardwareTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildL(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    private static double MassInCell(KleinSpectralMode m, int eigZ, int eigX) =>
        (eigZ, eigX) switch
        {
            (+1, +1) => m.MassPp,
            (+1, -1) => m.MassPm,
            (-1, +1) => m.MassMp,
            (-1, -1) => m.MassMm,
            _ => 0.0,
        };

    [Theory]
    [MemberData(nameof(F83FingerprintCases))]
    public void F83MarrakeshFingerprint_ObservableCell_HasSlowModeOverlap(
        string label, PauliPairBondTerm[] hTerms,
        PauliLetter[] observableLetters,
        double minSlowAggregateMass)
    {
        var chain = Chain3();
        var L = BuildL(hTerms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();

        int eigZ = PiOperator.SquaredEigenvalue(observableLetters, PauliLetter.Z);
        int eigX = PiOperator.SquaredEigenvalue(observableLetters, PauliLetter.X);

        // Aggregate slowest 16 modes' mass in the observable's Klein cell.
        double slowAggregate = modes.Take(16).Sum(m => MassInCell(m, eigZ, eigX));
        Assert.True(slowAggregate >= minSlowAggregateMass,
            $"{label}: slow-mode aggregate mass in observable cell ({eigZ:+#;-#}, {eigX:+#;-#}) " +
            $"= {slowAggregate:F4}, expected ≥ {minSlowAggregateMass}");
    }

    [Theory]
    [MemberData(nameof(F83FingerprintCases))]
    public void F83MarrakeshFingerprint_ObservableCell_AggregateOverFullSpectrum(
        string label, PauliPairBondTerm[] hTerms,
        PauliLetter[] observableLetters,
        double minSlowAggregateMass)
    {
        // Across the full 4^N spectrum, the observable's Klein cell aggregate must exceed
        // the slow-mode threshold (since slow modes are a subset of the full spectrum, the
        // total must dominate the slow-mode contribution).
        var chain = Chain3();
        var L = BuildL(hTerms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N);

        int eigZ = PiOperator.SquaredEigenvalue(observableLetters, PauliLetter.Z);
        int eigX = PiOperator.SquaredEigenvalue(observableLetters, PauliLetter.X);

        double totalAggregate = modes.Sum(m => MassInCell(m, eigZ, eigX));
        Assert.True(totalAggregate > minSlowAggregateMass,
            $"{label}: total aggregate mass {totalAggregate:F4} should exceed slow-mode threshold {minSlowAggregateMass}");
    }

    public static IEnumerable<object[]> F83FingerprintCases() => new[]
    {
        // Per ConfirmationsRegistry "f83_pi2_class_signature_marrakesh": the 4 H-classes
        // and their unique-fingerprint observables tested at path [4,5,6] on Marrakesh
        // 2026-04-30. We use the F87CanonicalWitness Hamiltonians as the closest match.
        // The minSlowAggregateMass thresholds are observed empirically from N=3 chain
        // computations and locked here as structural baselines.
        new object[]
        {
            "Truly XX+YY → ⟨Y₀ I Z₂⟩ in Pm",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) },
            new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.Z },
            0.5,
        },
        new object[]
        {
            "Pi2EvenNonTruly YZ+ZY → ⟨X₀ I X₂⟩ in Pp",
            new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) },
            new[] { PauliLetter.X, PauliLetter.I, PauliLetter.X },
            0.5,
        },
        new object[]
        {
            "Pi2OddPure XY+YX → ⟨X₀ I Z₂⟩ in Mm",
            new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) },
            new[] { PauliLetter.X, PauliLetter.I, PauliLetter.Z },
            0.5,
        },
        new object[]
        {
            "Mixed XX+XY → ⟨Z₀ I X₂⟩ in Mm",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) },
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.X },
            0.5,
        },
    };

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
