using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>X- and Y-dephasing on the block engines through the exact letter turn
/// (<see cref="LetterTurn"/>): the turn is a proper rotation carrying the jump to Z, written exactly on
/// Pauli letters, and the engines' spectrum of the turned Hamiltonian equals the dense eigensolver of
/// the J-dephased Lindbladian of the original one. Test Hamiltonians are built in the Z frame and
/// turned back by the inverse turn (the Hadamard is its own inverse, the quarter turn about X has
/// order four), so they conserve the jump letter's magnetisation by construction. The sign of each turned letter
/// is pinned only by the exact conjugation test: an improper sign assignment (H mapped to its transpose) leaves the
/// spectrum unchanged and would pass the spectral tests.</summary>
public class BlockSpectrumLetterTurnTests
{
    private readonly ITestOutputHelper _out;
    public BlockSpectrumLetterTurnTests(ITestOutputHelper output) => _out = output;

    private static readonly Complex I = Complex.ImaginaryOne;
    private static ComplexMatrix M(PauliLetter l) => PauliString.Build(new[] { l });

    /// <summary>Entry-for-entry equality (MathNet's == compares references).</summary>
    private static bool Same(ComplexMatrix a, ComplexMatrix b)
    {
        if (a.RowCount != b.RowCount || a.ColumnCount != b.ColumnCount) return false;
        for (int i = 0; i < a.RowCount; i++)
            for (int j = 0; j < a.ColumnCount; j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    // V = √2·U with entries in {0, ±1, ±i}, so V·P·V† = 2·(sign)·P' is an exact comparison.
    private static ComplexMatrix V(PauliLetter jump) => jump switch
    {
        PauliLetter.X => Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 1 }, { 1, -1 } }),           // Hadamard
        PauliLetter.Y => Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, -I }, { -I, 1 } }),          // 1 − iX
        _ => throw new ArgumentException(),
    };

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void TheTurn_IsTheConjugationByAProperClifford(PauliLetter jump)
    {
        var v = V(jump);
        foreach (var p in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
        {
            var (image, sign) = LetterTurn.TurnToZ(p, jump);
            Assert.True(Same(v * M(p) * v.ConjugateTranspose(), 2.0 * sign * M(image)), $"{p} under the turn for {jump}");
        }
        Assert.Equal((PauliLetter.Z, 1), LetterTurn.TurnToZ(jump, jump));
        Assert.Equal((PauliLetter.I, 1), LetterTurn.TurnToZ(PauliLetter.I, jump));
    }

    [Fact]
    public void TheZTurn_IsTheIdentity()
    {
        foreach (var p in new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
            Assert.Equal((p, 1), LetterTurn.TurnToZ(p, PauliLetter.Z));
    }

    // The Heisenberg bond XX + YY + ZZ is invariant under every proper rotation, so its turned
    // matrix equals the original entry for entry.
    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void TheHeisenbergChain_IsFixedByTheTurn(PauliLetter jump)
    {
        foreach (int N in new[] { 3, 4 })
        {
            var H = PauliHamiltonian.HeisenbergChain(N, new[] { 0.7, 1.3, 0.9 }.Take(N - 1).ToArray());
            Assert.True(Same(LetterTurn.Turn(H, jump).ToMatrix(), H.ToMatrix()));
        }
    }

    private static PauliHamiltonian InverseTurn(PauliHamiltonian HZ, PauliLetter jump) =>
        jump == PauliLetter.X ? LetterTurn.Turn(HZ, jump)
            : LetterTurn.Turn(LetterTurn.Turn(LetterTurn.Turn(HZ, jump), jump), jump);

    // Z-frame Hamiltonians: an XY chain with a Z field on one site and an XXZ chain with a uniform Z
    // field (both fail the pairing guard, so every sector is solved), an XXZ chain without a field
    // (passes it and is not fixed by the turn, so the pairing runs on a turned Hamiltonian), and a
    // Heisenberg chain with weighted bonds.
    private static PauliHamiltonian ZFrame(string name, int N)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            double J = 0.8 + 0.3 * b;
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, J));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, J));
            if (name != "XY + Z on site 0") terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Z, name.StartsWith("XXZ") ? 0.6 * J : J));
        }
        if (name == "XY + Z on site 0") terms.Add(PauliTerm.SingleSite(N, 0, PauliLetter.Z, 0.7));
        if (name == "XXZ + Z field") for (int l = 0; l < N; l++) terms.Add(PauliTerm.SingleSite(N, l, PauliLetter.Z, 0.4));
        return new PauliHamiltonian(N, terms);
    }

    // No exact route connects two eigensolves; the deviation (largest gap of a greedy one-to-one
    // matching) is read against eps·‖L‖_F at three scales a hundredfold apart, the error model of
    // BlockSpectrumPairingGuardTests: bounded and flat across the scales where the spectra agree.
    private const double Eps = 2.220446049250313e-16;
    private static readonly double[] Scales = { 1e-2, 1.0, 1e2 };
    // Measured on this machine's MKL: 2.9 to 8.9 where the spectra agree (every case, size and scale below),
    // 4.4e14 for the wrong channel, flat across the scales.
    private const double RoundingRatioCeiling = 100.0;
    private const double SpreadCeiling = 10.0;
    private const double StructuralRatioFloor = 1e10;

    private static double GreedyMatchDistance(IReadOnlyList<Complex> a, IReadOnlyList<Complex> b)
    {
        Assert.Equal(a.Count, b.Count);
        var taken = new bool[b.Count];
        double worst = 0.0;
        foreach (var z in a)
        {
            int best = -1; double bestDistance = double.MaxValue;
            for (int j = 0; j < b.Count; j++)
            {
                if (taken[j]) continue;
                double distance = (z - b[j]).Magnitude;
                if (distance < bestDistance) { bestDistance = distance; best = j; }
            }
            taken[best] = true; worst = Math.Max(worst, bestDistance);
        }
        return worst;
    }

    private void AssertAgree(IReadOnlyList<double> ratios, string context)
    {
        _out.WriteLine($"{context}: deviation / (eps·‖L‖_F) = " + string.Join(" | ", ratios.Select(r => r.ToString("G4"))));
        foreach (double r in ratios) Assert.True(r <= RoundingRatioCeiling, $"{context}: ratio {r:G4}");
        double spread = ratios.Max(r => Math.Max(r, 1.0)) / ratios.Min(r => Math.Max(r, 1.0));
        Assert.True(spread <= SpreadCeiling, $"{context}: spread {spread:G4}");
    }

    private static PauliHamiltonian Scale(PauliHamiltonian H, double s) =>
        new(H.N, H.Terms.Select(t => t with { Coefficient = s * t.Coefficient }).ToList());

    [Theory]
    [InlineData("XY + Z on site 0", PauliLetter.X)]
    [InlineData("XY + Z on site 0", PauliLetter.Y)]
    [InlineData("XXZ + Z field", PauliLetter.X)]
    [InlineData("XXZ + Z field", PauliLetter.Y)]
    [InlineData("XXZ", PauliLetter.X)]
    [InlineData("XXZ", PauliLetter.Y)]
    [InlineData("Heisenberg", PauliLetter.Y)]
    public void PerBlock_MatchesTheFullEigensolver_UnderXAndYDephasing(string name, PauliLetter jump)
    {
        foreach (int N in new[] { 3, 4 })
        {
            var HX = InverseTurn(ZFrame(name, N), jump);
            Assert.True(Same(LetterTurn.Turn(HX, jump).ToMatrix(), ZFrame(name, N).ToMatrix()), "the inverse turn returns the Z frame");
            var ratios = new List<double>();
            foreach (double s in Scales)
            {
                var H = Scale(HX, s);
                var gamma = Enumerable.Range(0, N).Select(l => s * (0.3 + 0.1 * l)).ToArray();
                var L = PauliDephasingDissipator.Build(H.ToMatrix(), gamma, jump);
                var full = L.Evd().EigenValues.ToArray();
                var perBlock = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gamma, jump);
                ratios.Add(GreedyMatchDistance(perBlock, full) / (Eps * L.FrobeniusNorm()));
            }
            AssertAgree(ratios, $"{name}, {jump}-dephasing, N={N}");
        }
    }

    // Control: the same per-block spectrum against the dense Z-dephased Lindbladian of the unturned
    // H, the wrong channel, disagrees structurally (the ratio sits near 1/eps, not near rounding).
    [Fact]
    public void PerBlock_DisagreesWithTheWrongChannel()
    {
        const int N = 3;
        var HX0 = InverseTurn(ZFrame("XXZ + Z field", N), PauliLetter.X);
        var ratios = new List<double>();
        foreach (double s in Scales)
        {
            var HX = Scale(HX0, s);
            var gamma = new[] { 0.3 * s, 0.4 * s, 0.5 * s };
            var LZ = PauliDephasingDissipator.Build(HX.ToMatrix(), gamma, PauliLetter.Z);
            var wrong = LZ.Evd().EigenValues.ToArray();
            var perBlock = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(HX, gamma, PauliLetter.X);
            ratios.Add(GreedyMatchDistance(perBlock, wrong) / (Eps * LZ.FrobeniusNorm()));
        }
        _out.WriteLine("wrong channel: " + string.Join(" | ", ratios.Select(r => r.ToString("G4"))));
        foreach (double r in ratios) Assert.True(r >= StructuralRatioFloor, $"ratio {r:G4}");
        double spread = ratios.Max() / ratios.Min();
        Assert.True(spread <= SpreadCeiling, $"spread {spread:G4}");
    }

    // The F71 refinement returns L's spectrum for reflection-symmetric H and rates.
    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void F71Refined_MatchesTheFullEigensolver_UnderXAndYDephasing(PauliLetter jump)
    {
        foreach (int N in new[] { 3, 4 })
        {
            var terms = new List<PauliTerm>();
            for (int b = 0; b < N - 1; b++)
                foreach (var p in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
                    terms.Add(PauliTerm.TwoSite(N, b, p, b + 1, p, 0.9));
            for (int l = 0; l < N; l++) terms.Add(PauliTerm.SingleSite(N, l, PauliLetter.Z, 0.4));
            var HX = InverseTurn(new PauliHamiltonian(N, terms), jump);
            var gamma = Enumerable.Repeat(0.35, N).ToArray();
            var L = PauliDephasingDissipator.Build(HX.ToMatrix(), gamma, jump);
            var full = L.Evd().EigenValues.ToArray();
            var refined = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(HX, gamma, jump);
            double ratio = GreedyMatchDistance(refined, full) / (Eps * L.FrobeniusNorm());
            _out.WriteLine($"F71, {jump}, N={N}: ratio {ratio:G4}");
            Assert.True(ratio <= RoundingRatioCeiling, $"F71, {jump}, N={N}: ratio {ratio:G4}");
        }
    }

    // In the frame of H the guard reads [H, Z^⊗N] = 0 under X-dephasing and [H, X^⊗N] = 0 under
    // Y-dephasing: computed directly on the X- or Y-frame matrix, the commutator agrees with the guard
    // the engine applies to the turned matrix, on a Hamiltonian that passes and one that fails.
    [Theory]
    [InlineData(PauliLetter.X, PauliLetter.Z)]
    [InlineData(PauliLetter.Y, PauliLetter.X)]
    public void ThePairingGuard_ReadsTheFrameOfH(PauliLetter jump, PauliLetter guardLetter)
    {
        const int N = 4;
        var P = PauliString.Build(Enumerable.Repeat(guardLetter, N).ToArray());
        foreach (var (name, passes) in new[] { ("XXZ", true), ("XY + Z on site 0", false) })
        {
            var Hj = InverseTurn(ZFrame(name, N), jump).ToMatrix();
            bool commutes = Same(Hj * P, P * Hj);
            bool guard = LiouvillianBlockSpectrum.CommutesWithXN(LetterTurn.Turn(InverseTurn(ZFrame(name, N), jump), jump).ToMatrix(), N);
            Assert.Equal(passes, commutes);
            Assert.Equal(commutes, guard);
        }
    }

    // The XY chain under X-dephasing: XX + YY turns to ZZ + YY, and YY breaks popcount in the
    // turned frame, so the engine refuses it exactly as it refuses such an H under Z.
    [Fact]
    public void AHamiltonianThatBreaksTheJumpPopcount_IsRefused()
    {
        var H = PauliHamiltonian.XYChain(3, 1.0);
        var ex = Assert.Throws<ArgumentException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, new[] { 0.3, 0.3, 0.3 }, PauliLetter.X));
        Assert.Contains("popcount", ex.Message);
        var exF71 = Assert.Throws<ArgumentException>(() =>
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, new[] { 0.3, 0.3, 0.3 }, PauliLetter.X));
        Assert.Contains("popcount", exF71.Message);
    }

    [Fact]
    public void TheIdentity_IsNotADephasingLetter()
    {
        var H = PauliHamiltonian.HeisenbergChain(3, new[] { 1.0, 1.0 });
        Assert.Throws<ArgumentException>(() => LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, new[] { 0.3, 0.3, 0.3 }, PauliLetter.I));
        Assert.Throws<ArgumentException>(() => F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, new[] { 0.3, 0.3, 0.3 }, PauliLetter.I));
    }
}
