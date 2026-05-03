using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The raw view through Schicht 1 (Π²-Projection alone, no JW, no Bogoliubov,
/// no spectral test). What does Sein Blickwinkel reveal about the F87 canonical Hamiltonians
/// and the Marrakesh signature observables, before any further machinery is bolted on?
///
/// <para>These tests don't introduce new structure; they apply <see cref="Pi2Projection"/>
/// to known objects and lock what's seen. The findings are themselves the artifact.</para>
/// </summary>
public class Pi2ProjectionRawViewTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildH(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain) =>
        PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();

    private static double FrobNormSquared(ComplexMatrix m) =>
        (m.ConjugateTranspose() * m).Trace().Real;

    [Theory]
    [MemberData(nameof(CanonicalHamiltonianCases))]
    public void Pi2Decomposition_Of_F87_CanonicalHamiltonians(
        string label, PauliPairBondTerm[] terms,
        double expectedEvenNormSq, double expectedOddNormSq)
    {
        // Schicht 1 view: take H, split into (H_even, H_odd) under Π². No L, no Spec, no JW.
        var chain = Chain3();
        var H = BuildH(terms, chain);
        var (hEven, hOdd) = Pi2Projection.Split(H, chain.N);

        double evenSq = FrobNormSquared(hEven);
        double oddSq = FrobNormSquared(hOdd);
        Assert.True(Math.Abs(evenSq - expectedEvenNormSq) < 1e-9,
            $"{label}: ‖H_even‖² expected {expectedEvenNormSq}, got {evenSq}");
        Assert.True(Math.Abs(oddSq - expectedOddNormSq) < 1e-9,
            $"{label}: ‖H_odd‖² expected {expectedOddNormSq}, got {oddSq}");
        Assert.True((hEven + hOdd - H).FrobeniusNorm() < 1e-12,
            $"{label}: H_even + H_odd ≠ H");
    }

    public static IEnumerable<object[]> CanonicalHamiltonianCases() => new[]
    {
        // XX+YY at N=3 chain: 4 truly Π²-even bilinears (X₀X₁, X₁X₂, Y₀Y₁, Y₁Y₂),
        // each ‖σ‖²=2³=8. Total ‖H‖²=32. All in Π²-even subspace.
        new object[] { "XX+YY (Truly)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) },
            32.0, 0.0 },

        // Heisenberg at N=3: 6 truly Π²-even bilinears. ‖H‖²=48, all even.
        new object[] { "XX+YY+ZZ (Heisenberg, Truly)",
            new[]
            {
                Term(PauliLetter.X, PauliLetter.X),
                Term(PauliLetter.Y, PauliLetter.Y),
                Term(PauliLetter.Z, PauliLetter.Z),
            },
            48.0, 0.0 },

        // XY+YX at N=3: 4 Π²-odd bilinears. ‖H‖²=32, all odd.
        new object[] { "XY+YX (Pi2OddPure)",
            new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) },
            0.0, 32.0 },

        // YZ+ZY at N=3: 4 Π²-even non-truly bilinears. ‖H‖²=32, all even.
        // Indistinguishable from Truly under Π² alone — the truly criterion is a SECOND cut.
        new object[] { "YZ+ZY (Pi2EvenNonTruly)",
            new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) },
            32.0, 0.0 },

        // XX+XY at N=3: 2 Π²-even truly (XX×2 bonds, ‖²=16) + 2 Π²-odd (XY×2 bonds, ‖²=16).
        // Mass split exactly half/half between the two Π²-subspaces.
        new object[] { "XX+XY (Mixed)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) },
            16.0, 16.0 },
    };

    [Fact]
    public void RawView_Cannot_Distinguish_Truly_From_Pi2EvenNonTruly()
    {
        // The structural finding: Schicht 1 (Π²-Projection) alone groups Truly and
        // Pi2EvenNonTruly into the same shape — both have H_odd = 0 and H_even = H.
        // The truly-vs-non-truly distinction lives in a SECOND cut (the #Y / #Z parity
        // criterion in PauliPairBondTerm.IsTruly), invisible to Π² conjugation alone.
        var chain = Chain3();

        var truly = BuildH(new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.Y, PauliLetter.Y),
        }, chain);
        var (trulyEven, trulyOdd) = Pi2Projection.Split(truly, chain.N);

        var pi2EvenNonTruly = BuildH(new[]
        {
            Term(PauliLetter.Y, PauliLetter.Z),
            Term(PauliLetter.Z, PauliLetter.Y),
        }, chain);
        var (nonTrulyEven, nonTrulyOdd) = Pi2Projection.Split(pi2EvenNonTruly, chain.N);

        // Both have empty Π²-odd component. Schicht 1 sees no difference here.
        Assert.True(trulyOdd.FrobeniusNorm() < 1e-12);
        Assert.True(nonTrulyOdd.FrobeniusNorm() < 1e-12);
        // Both have the full mass on the Π²-even side.
        Assert.Equal(FrobNormSquared(truly), FrobNormSquared(trulyEven), precision: 9);
        Assert.Equal(FrobNormSquared(pi2EvenNonTruly), FrobNormSquared(nonTrulyEven), precision: 9);
    }

    [Theory]
    [MemberData(nameof(MarrakeshObservableCases))]
    public void Pi2Class_Of_Marrakesh_SignatureObservables(
        string label, PauliLetter[] letters, int expectedPi2Eigenvalue)
    {
        // The 2026-04-30 Marrakesh f83 4-class fingerprint test chose one observable per
        // F87 class. This test reads the Π²-eigenvalue of each chosen observable to lock
        // the framework's matching pattern: every signature observable lives in a definite
        // Π²-eigenspace.
        int eig = PiOperator.SquaredEigenvalue(letters);
        Assert.True(eig == expectedPi2Eigenvalue,
            $"{label}: expected Π²-eigenvalue {expectedPi2Eigenvalue}, got {eig}");
    }

    public static IEnumerable<object[]> MarrakeshObservableCases() => new[]
    {
        // Confirmation: f83_pi2_class_signature_marrakesh fingerprints (2026-04-30 N=3):
        new object[] { "⟨X₀ I Z₂⟩ for pi2_odd_pure",
            new[] { PauliLetter.X, PauliLetter.I, PauliLetter.Z }, -1 },
        new object[] { "⟨X₀ I X₂⟩ for pi2_even_nontruly",
            new[] { PauliLetter.X, PauliLetter.I, PauliLetter.X }, +1 },
        new object[] { "⟨Z₀ I X₂⟩ for mixed",
            new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.X }, -1 },
        new object[] { "⟨Y₀ I Z₂⟩ for truly",
            new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.Z }, +1 },
    };

    [Fact]
    public void Marrakesh_Observable_Pattern_OneClassEach_ButNotInBijectionWithHClasses()
    {
        // The 4 fingerprint observables split as 2 Π²-even + 2 Π²-odd (not 4 distinct
        // Π²-classes — there are only 2 Π²-classes total). The framework's f83 test
        // discriminates the 4 H classes via observables that are NOT Π²-class-distinct;
        // the discrimination lives in the SECOND cut (specific Pauli string within each
        // Π²-class), the same invariant that distinguishes Truly from Pi2EvenNonTruly.
        var observables = new (PauliLetter[] letters, string forClass)[]
        {
            (new[] { PauliLetter.X, PauliLetter.I, PauliLetter.Z }, "pi2_odd_pure"),
            (new[] { PauliLetter.X, PauliLetter.I, PauliLetter.X }, "pi2_even_nontruly"),
            (new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.X }, "mixed"),
            (new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.Z }, "truly"),
        };

        int evenCount = observables.Count(o => PiOperator.SquaredEigenvalue(o.letters) == +1);
        int oddCount = observables.Count(o => PiOperator.SquaredEigenvalue(o.letters) == -1);
        Assert.Equal(2, evenCount);
        Assert.Equal(2, oddCount);
    }

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
