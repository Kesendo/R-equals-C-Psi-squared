using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Stage 3b: 8-cell extension of <see cref="Pi2KleinSpectralView"/>.
/// Verifies that the 8-cell mass distribution is well-formed, normalised,
/// reduces to the 4-cell view when summed across y_par, and reveals the
/// y_par-axis content per eigenmode.
/// </summary>
public class Pi2KleinEightSpectralViewTests
{
    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildL(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);

    [Fact]
    public void EveryEigenmode_Has8CellDistributionSummingTo1()
    {
        // After per-mode normalisation the 8 cell masses sum to 1.
        var chain = Chain3();
        var L = BuildL(new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) }, chain);
        var modes = Pi2KleinEightSpectralView.ComputeFor(L, chain.N);

        Assert.NotEmpty(modes);
        foreach (var m in modes)
        {
            Assert.True(Math.Abs(m.TotalMass - 1.0) < 1e-9,
                $"total 8-cell mass {m.TotalMass:F6} ≠ 1 for eigenmode at λ={m.Eigenvalue:F4}");
        }
    }

    [Fact]
    public void EigenmodeCount_Equals_4PowN()
    {
        var chain = Chain3();
        var L = BuildL(new[] { Term(PauliLetter.X, PauliLetter.X) }, chain);
        var modes = Pi2KleinEightSpectralView.ComputeFor(L, chain.N);
        long expected = 1L << (2 * chain.N);
        Assert.Equal(expected, modes.Count);
    }

    [Fact]
    public void KleinMassAggregation_RecoversPi2KleinSpectralViewMasses()
    {
        // The 8-cell view summed across y_par must equal the 4-cell view's MassPp/Pm/Mp/Mm.
        // This is the structural consistency check: 8-cell is a refinement, not a replacement.
        var chain = Chain3();
        var L = BuildL(new[] {
            Term(PauliLetter.X, PauliLetter.X),
            Term(PauliLetter.Y, PauliLetter.Y),
        }, chain);

        var eightCellModes = Pi2KleinEightSpectralView.ComputeFor(L, chain.N);
        var fourCellModes = Pi2KleinSpectralView.ComputeFor(L, chain.N);

        Assert.Equal(eightCellModes.Count, fourCellModes.Count);
        for (int i = 0; i < eightCellModes.Count; i++)
        {
            var e = eightCellModes[i];
            var f = fourCellModes[i];
            // KleinMass(bit_a, bit_b) maps to fourCell's MassPp/Pm/Mp/Mm convention.
            // Per Pi2KleinSpectralView lines 66-69:
            //   Pp = (eigZ +1, eigX +1) = (bit_b=0, bit_a=0)
            //   Pm = (eigZ +1, eigX −1) = (bit_b=0, bit_a=1)
            //   Mp = (eigZ −1, eigX +1) = (bit_b=1, bit_a=0)
            //   Mm = (eigZ −1, eigX −1) = (bit_b=1, bit_a=1)
            Assert.True(Math.Abs(e.KleinMass(bitA: 0, bitB: 0) - f.MassPp) < 1e-9,
                $"mode {i}: KleinMass(0,0)={e.KleinMass(0, 0):F6} ≠ Pp={f.MassPp:F6}");
            Assert.True(Math.Abs(e.KleinMass(bitA: 1, bitB: 0) - f.MassPm) < 1e-9,
                $"mode {i}: KleinMass(1,0)={e.KleinMass(1, 0):F6} ≠ Pm={f.MassPm:F6}");
            Assert.True(Math.Abs(e.KleinMass(bitA: 0, bitB: 1) - f.MassMp) < 1e-9,
                $"mode {i}: KleinMass(0,1)={e.KleinMass(0, 1):F6} ≠ Mp={f.MassMp:F6}");
            Assert.True(Math.Abs(e.KleinMass(bitA: 1, bitB: 1) - f.MassMm) < 1e-9,
                $"mode {i}: KleinMass(1,1)={e.KleinMass(1, 1):F6} ≠ Mm={f.MassMm:F6}");
        }
    }

    [Fact]
    public void YParityZeroPlusOneMass_SumsTo1PerMode()
    {
        // y_par axis is a clean partition: YParityZeroMass + YParityOneMass = 1 per
        // mode (after normalisation). Y-par non-conservation by L doesn't affect this
        // partition identity; it only means an eigenmode can have non-zero mass in
        // BOTH y_par=0 and y_par=1.
        var chain = Chain3();
        var L = BuildL(new[] { Term(PauliLetter.X, PauliLetter.Y) }, chain);
        var modes = Pi2KleinEightSpectralView.ComputeFor(L, chain.N);

        foreach (var m in modes)
        {
            double sum = m.YParityZeroMass + m.YParityOneMass;
            Assert.True(Math.Abs(sum - 1.0) < 1e-9,
                $"y_par partition: y0={m.YParityZeroMass:F6} + y1={m.YParityOneMass:F6} = {sum:F6} ≠ 1 at λ={m.Eigenvalue:F4}");
        }
    }

    [Fact]
    public void XYChain_TrulyHamiltonian_KernelModes_AreYParityZeroPure()
    {
        // F107 says: truly H ⟹ M=0 ⟹ kernel modes are linear combinations of truly
        // Pauli strings, which per F85 have #Y even (= y_par=0). So for XX+YY truly H,
        // the 4 kernel modes should be y_par=0-pure (YParityZeroMass ≈ 1).
        //
        // Note: this is an empirical test of F107's IMPLICATION for L's eigenmode
        // structure, not a direct restatement of F107. The 4 kernel modes are
        // stationary states living in L's null space; their decomposition over the
        // Pauli basis lands on truly bilinears (XX, YY, ZZ, II), each with #Y even,
        // hence y_par=0.
        var chain = Chain3();
        var L = BuildL(new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) }, chain);
        var sorted = Pi2KleinEightSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();

        var kernel = sorted.Take(4).ToList();
        foreach (var m in kernel)
        {
            Assert.True(m.YParityZeroMass > 0.99,
                $"kernel mode λ={m.Eigenvalue:F6}: y_par=0 mass {m.YParityZeroMass:F4} < 0.99 (F107 expects truly kernel y_par-pure)");
            Assert.True(m.YParityOneMass < 0.01,
                $"kernel mode λ={m.Eigenvalue:F6}: y_par=1 mass {m.YParityOneMass:F4} > 0.01 (F107 expects ≈ 0)");
        }
    }

    [Fact]
    public void DominantCellLabel_RenderableHumanLabel()
    {
        // Smoke check: DominantCellLabel returns a non-empty string for every mode and
        // matches KleinEightCellClaim's naming convention.
        var chain = Chain3();
        var L = BuildL(new[] { Term(PauliLetter.X, PauliLetter.Z) }, chain);
        var modes = Pi2KleinEightSpectralView.ComputeFor(L, chain.N);

        foreach (var m in modes)
        {
            var label = m.DominantCellLabel;
            Assert.False(string.IsNullOrEmpty(label));
            Assert.Contains(label, new[]
            {
                "Mother/no-Y", "Mother/odd-Y",
                "Z-Klein/no-Y", "Z-Klein/odd-Y",
                "X-Klein/no-Y", "X-Klein/odd-Y",
                "Y-Klein/no-Y (paradox)", "Y-Klein/odd-Y (canonical)",
            });
        }
    }

    [Fact]
    public void WrongDimension_ThrowsArgumentException()
    {
        var chain = Chain3();
        long d2 = 1L << (2 * chain.N);
        // Build a deliberately wrong-sized matrix (e.g., d² + 1).
        var wrongL = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.Dense(
            (int)(d2 + 1), (int)(d2 + 1));
        Assert.Throws<ArgumentException>(() => Pi2KleinEightSpectralView.ComputeFor(wrongL, chain.N));
    }
}
