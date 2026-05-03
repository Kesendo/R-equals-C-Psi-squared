using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The 2-axis Klein view of Schicht 1: each Pauli string carries a pair of
/// Π²-eigenvalues (Π²_Z, Π²_X). The 4 cells distinguish what Π²_Z alone could not — truly
/// bilinears (XX, YY, ZZ) sit in (+, +) while non-truly Π²-even bilinears (YZ, ZY) sit in
/// (−, +). F80's "universality across 4 Π²-odd cases" is itself a universality across 2
/// Klein-cells, (+, −) and (−, −).
///
/// <para>This is the framework's own 2-axis structural view, not an import from any standard
/// classification scheme. We're painting our own picture.</para>
/// </summary>
public class Pi2KleinViewTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildH(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain) =>
        PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();

    private static double FrobNormSquared(ComplexMatrix m) =>
        (m.ConjugateTranspose() * m).Trace().Real;

    [Theory]
    [MemberData(nameof(KleinViewCases))]
    public void KleinDecomposition_Of_F87_CanonicalHamiltonians(
        string label, PauliPairBondTerm[] terms,
        double expectedPpSq, double expectedMpSq, double expectedPmSq, double expectedMmSq)
    {
        // Schicht 1 with both axes: decompose H into 4 Klein-cells under (Π²_Z, Π²_X).
        var chain = Chain3();
        var H = BuildH(terms, chain);
        var klein = Pi2Projection.KleinSplit(H, chain.N);

        double ppSq = FrobNormSquared(klein.Pp);
        double mpSq = FrobNormSquared(klein.Mp);
        double pmSq = FrobNormSquared(klein.Pm);
        double mmSq = FrobNormSquared(klein.Mm);

        Assert.True(Math.Abs(ppSq - expectedPpSq) < 1e-9,
            $"{label}: ‖H_(+,+)‖² expected {expectedPpSq}, got {ppSq}");
        Assert.True(Math.Abs(mpSq - expectedMpSq) < 1e-9,
            $"{label}: ‖H_(−,+)‖² expected {expectedMpSq}, got {mpSq}");
        Assert.True(Math.Abs(pmSq - expectedPmSq) < 1e-9,
            $"{label}: ‖H_(+,−)‖² expected {expectedPmSq}, got {pmSq}");
        Assert.True(Math.Abs(mmSq - expectedMmSq) < 1e-9,
            $"{label}: ‖H_(−,−)‖² expected {expectedMmSq}, got {mmSq}");

        // Round-trip: 4 cells sum to H exactly.
        var sum = klein.Pp + klein.Mp + klein.Pm + klein.Mm;
        Assert.True((sum - H).FrobeniusNorm() < 1e-12,
            $"{label}: 4-cell sum ≠ H");
    }

    public static IEnumerable<object[]> KleinViewCases() => new[]
    {
        // Truly Heisenberg XX+YY+ZZ on N=3: 6 truly bilinears (XX×2, YY×2, ZZ×2).
        // All three letter pairs are in cell (+, +). ‖H‖² = 6·8 = 48 → all in Pp.
        new object[] { "XX+YY+ZZ (Heisenberg) → all (+, +)",
            new[]
            {
                Term(PauliLetter.X, PauliLetter.X),
                Term(PauliLetter.Y, PauliLetter.Y),
                Term(PauliLetter.Z, PauliLetter.Z),
            },
            48.0, 0.0, 0.0, 0.0 },

        // XX+YY truly: 4 bilinears all in (+, +). ‖H‖² = 32.
        new object[] { "XX+YY (Truly) → all (+, +)",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) },
            32.0, 0.0, 0.0, 0.0 },

        // YZ+ZY non-truly Π²-even: bit_a sum = 1+0 or 0+1 = 1 → Π²_X=−. bit_b sum = 2 → Π²_Z=+.
        // 4 bilinears all in cell (+, −) = Pm. Distinct from Truly (Pp).
        new object[] { "YZ+ZY (Pi2EvenNonTruly) → all (+, −)",
            new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) },
            0.0, 0.0, 32.0, 0.0 },

        // XY+YX bond-flip-Paar A: bit_a=1+1=0 → Π²_X=+. bit_b=0+1=1 → Π²_Z=−.
        // 4 bilinears all in cell (−, +) = Mp.
        new object[] { "XY+YX (Pi2OddPure subgroup A) → all (−, +)",
            new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) },
            0.0, 32.0, 0.0, 0.0 },

        // XZ+ZX bond-flip-Paar B: bit_a=1 → Π²_X=−. bit_b=1 → Π²_Z=−.
        // 4 bilinears all in cell (−, −) = Mm. F80's "universality" treats XY/YX and
        // XZ/ZX as one class; the Klein view sees them as two distinct cells.
        new object[] { "XZ+ZX (Pi2OddPure subgroup B) → all (−, −)",
            new[] { Term(PauliLetter.X, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.X) },
            0.0, 0.0, 0.0, 32.0 },

        // XX+XY mixed: XX in (+, +) = Pp, XY in (−, +) = Mp. Half/half between Pp and Mp.
        new object[] { "XX+XY (Mixed across (+,+) and (−,+))",
            new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) },
            16.0, 16.0, 0.0, 0.0 },
    };

    [Fact]
    public void KleinView_Distinguishes_Truly_From_Pi2EvenNonTruly()
    {
        // What Π²_Z alone could not see: the Klein view places truly bilinears in cell (+, +)=Pp
        // and non-truly Π²-even bilinears in cell (+, −)=Pm. The orthogonal axis Π²_X is the
        // discriminator (X-axis flips between truly and non-truly Π²-even).
        var chain = Chain3();

        var truly = BuildH(new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.Y, PauliLetter.Y),
        }, chain);
        var trulyKlein = Pi2Projection.KleinSplit(truly, chain.N);

        var pi2EvenNonTruly = BuildH(new[]
        {
            Term(PauliLetter.Y, PauliLetter.Z),
            Term(PauliLetter.Z, PauliLetter.Y),
        }, chain);
        var nonTrulyKlein = Pi2Projection.KleinSplit(pi2EvenNonTruly, chain.N);

        // Truly lives in Pp; Pi2EvenNonTruly lives in Pm. Mutually exclusive cells.
        Assert.True(FrobNormSquared(trulyKlein.Pp) > 1e-9);
        Assert.True(FrobNormSquared(trulyKlein.Pm) < 1e-12);

        Assert.True(FrobNormSquared(nonTrulyKlein.Pm) > 1e-9);
        Assert.True(FrobNormSquared(nonTrulyKlein.Pp) < 1e-12);
    }

    [Fact]
    public void F80_Universality_Spans_Two_Klein_Cells_Not_One()
    {
        // F80 universality: all 4 Π²-odd 2-body bilinears (XY, YX, XZ, ZX) give the same
        // M spectrum. The Klein view shows these split into TWO cells: XY/YX in (+, −) and
        // XZ/ZX in (−, −). F80's universality is therefore a universality across two
        // Klein-cells, not one.
        var chain = Chain3();
        var subgroupA = BuildH(new[] { Term(PauliLetter.X, PauliLetter.Y) }, chain);
        var subgroupB = BuildH(new[] { Term(PauliLetter.X, PauliLetter.Z) }, chain);

        var kleinA = Pi2Projection.KleinSplit(subgroupA, chain.N);
        var kleinB = Pi2Projection.KleinSplit(subgroupB, chain.N);

        // Subgroup A entirely in Pm = (Π²_Z=+1, Π²_X=−1). Wait, let me re-derive: XY has
        // bit_b=0+1=1 → Π²_Z=−1. bit_a=1+1=0 → Π²_X=+1. So XY is in Mp = (−, +).
        // Hmm — that contradicts my hand-calculation earlier. Let me check.
        // Per-letter bit_a for X is 1, for Y is 1. Sum = 2 mod 2 = 0. So Π²_X = +1.
        // Per-letter bit_b for X is 0, for Y is 1. Sum = 1. So Π²_Z = −1.
        // So XY is in cell (Π²_Z=−1, Π²_X=+1) = Mp. Subgroup A → Mp.
        Assert.True(FrobNormSquared(kleinA.Mp) > 1e-9, "subgroup A (XY) should be in Mp = (−, +)");
        Assert.True(FrobNormSquared(kleinA.Mm) < 1e-12);

        // XZ: bit_a sum = 1+0=1 → Π²_X=−1. bit_b sum = 0+1=1 → Π²_Z=−1. Cell Mm = (−, −).
        Assert.True(FrobNormSquared(kleinB.Mm) > 1e-9, "subgroup B (XZ) should be in Mm = (−, −)");
        Assert.True(FrobNormSquared(kleinB.Mp) < 1e-12);
    }

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
