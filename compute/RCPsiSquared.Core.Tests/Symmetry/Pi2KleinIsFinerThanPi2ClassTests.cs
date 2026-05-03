using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The 2-axis Klein view is genuinely finer than the 1-axis Pi2Class.
/// Pi2Class.Mixed (which means "Π²_Z mixed across bilinears") splits into TWO
/// distinct Klein sub-cases depending on whether Π²_X is uniform or also mixed.
/// This test exhibits two Mixed Hamiltonians that are indistinguishable under
/// Pi2Class but live in different Klein-cell distributions.
///
/// <para>The structural finding: Pi2Class encodes one parity axis (Z); the Klein
/// view encodes both (Z and X). Going from Pi2Class to Klein-cells is a true 2-fold
/// refinement on the Mixed class, and reveals analogous sub-structure in other
/// classes when crossed with longer term lists.</para>
/// </summary>
public class Pi2KleinIsFinerThanPi2ClassTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildH(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain) =>
        PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();

    private static double FrobNormSquared(ComplexMatrix m) =>
        (m.ConjugateTranspose() * m).Trace().Real;

    [Fact]
    public void TwoMixedHamiltonians_HaveSamePi2Class_ButDifferentKleinDistributions()
    {
        // H_A = XX + XY: bit_a uniform = 0 (both X-or-Y count = 2, mod 2 = 0).
        //   XX in Pp = (+, +), XY in Mp = (−, +). Klein: 50/50 across Π²_Z, uniform Π²_X = +.
        //
        // H_B = XX + YZ: bit_a per term varies.
        //   XX in Pp = (+, +). YZ: bit_a sum = 1, bit_b sum = 0. Cell = (+, −) = Pm.
        //   Klein: 50/50 across Π²_X, uniform Π²_Z = +.
        //
        // Both Hamiltonians are Pi2Class.Mixed because they have BOTH Π²-odd and Π²-even
        // bilinears... wait, let me re-derive. XY has Π²_Z = bit_b mod 2 = 1 → Π²-odd.
        // YZ has Π²_Z = 1+1 mod 2 = 0 → Π²-even. So H_A (XX truly + XY Π²-odd) IS
        // Pi2Class.Mixed. H_B (XX truly + YZ Π²-even-non-truly) is also a non-trivial
        // class, but per the Pi2ClassOf algorithm both Π²-even, no Π²-odd → Pi2EvenNonTruly.
        //
        // So H_A is Mixed, H_B is Pi2EvenNonTruly — they're already distinguished by Pi2Class.
        // Let me find a genuinely tricky case.
        var chain = Chain3();

        // H_A: XX (Pp truly) + XY (Mp Π²-odd). Pi2Class = Mixed. Klein: Pp + Mp.
        var hA = BuildH(new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.X, PauliLetter.Y),
        }, chain);
        var kleinA = Pi2Projection.KleinSplit(hA, chain.N);
        var pi2A = new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.X, PauliLetter.Y),
        }.Pi2ClassOf();

        // H_C: XX (Pp truly) + XZ (Mm Π²-odd). Pi2Class = Mixed. Klein: Pp + Mm — DIFFERENT!
        var hC = BuildH(new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.X, PauliLetter.Z),
        }, chain);
        var kleinC = Pi2Projection.KleinSplit(hC, chain.N);
        var pi2C = new[]
        {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.X, PauliLetter.Z),
        }.Pi2ClassOf();

        // Both are Pi2Class.Mixed:
        Assert.Equal(Pi2Class.Mixed, pi2A);
        Assert.Equal(Pi2Class.Mixed, pi2C);

        // But Klein distributions differ:
        // H_A occupies (Pp, Mp); H_C occupies (Pp, Mm). The Π²_X axis distinguishes them.
        Assert.True(FrobNormSquared(kleinA.Pp) > 1e-9 && FrobNormSquared(kleinA.Mp) > 1e-9,
            "H_A should occupy Pp and Mp");
        Assert.True(FrobNormSquared(kleinA.Pm) < 1e-12 && FrobNormSquared(kleinA.Mm) < 1e-12,
            "H_A should NOT touch Pm or Mm");

        Assert.True(FrobNormSquared(kleinC.Pp) > 1e-9 && FrobNormSquared(kleinC.Mm) > 1e-9,
            "H_C should occupy Pp and Mm");
        Assert.True(FrobNormSquared(kleinC.Pm) < 1e-12 && FrobNormSquared(kleinC.Mp) < 1e-12,
            "H_C should NOT touch Pm or Mp");
    }

    [Fact]
    public void Pi2Class_Pi2OddPure_AlsoSplitsIntoTwoKleinSubgroups()
    {
        // F80 universality: 4 Π²-odd bilinears (XY, YX, XZ, ZX) all give Pi2Class.Pi2OddPure.
        // But they split across two Klein-cells: (XY, YX) in Mp, (XZ, ZX) in Mm.
        // F80 universality is therefore over 2 Klein-cells, not one.
        var chain = Chain3();

        var hXY = BuildH(new[] { Term(PauliLetter.X, PauliLetter.Y) }, chain);
        var hXZ = BuildH(new[] { Term(PauliLetter.X, PauliLetter.Z) }, chain);

        var kleinXY = Pi2Projection.KleinSplit(hXY, chain.N);
        var kleinXZ = Pi2Projection.KleinSplit(hXZ, chain.N);

        // Both are Pi2Class.Pi2OddPure:
        Assert.Equal(Pi2Class.Pi2OddPure,
            new[] { Term(PauliLetter.X, PauliLetter.Y) }.Pi2ClassOf());
        Assert.Equal(Pi2Class.Pi2OddPure,
            new[] { Term(PauliLetter.X, PauliLetter.Z) }.Pi2ClassOf());

        // But XY lives in Mp and XZ in Mm — different Klein subgroups.
        Assert.True(FrobNormSquared(kleinXY.Mp) > 1e-9);
        Assert.True(FrobNormSquared(kleinXY.Mm) < 1e-12);

        Assert.True(FrobNormSquared(kleinXZ.Mm) > 1e-9);
        Assert.True(FrobNormSquared(kleinXZ.Mp) < 1e-12);
    }

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
