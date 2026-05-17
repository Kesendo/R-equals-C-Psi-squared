using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.BlockSpectrum.Prosen;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.Prosen;

/// <summary>Witnesses for <see cref="OneSidedSectorClosedForm"/>: closed-form eigenvalues
/// of the chain XY + uniform Z-dephasing Liouvillian on the (p_c = 0, p_r = m) sector,
/// expressed as subset sums of N Prosen rapidities β_k = −2γ − i·ε_k. Cross-validated
/// against <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> dense Evd at small N;
/// F1-mirror identity verified across the sector pair (0, m) ↔ (N, N−m).</summary>
public class OneSidedSectorClosedFormTests
{
    private readonly ITestOutputHelper _out;
    public OneSidedSectorClosedFormTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(3, 1)]
    [InlineData(3, 2)]
    [InlineData(3, 3)]
    [InlineData(4, 1)]
    [InlineData(4, 2)]
    [InlineData(4, 3)]
    [InlineData(5, 2)]
    [InlineData(5, 3)]
    [InlineData(6, 3)]
    public void Build_ClosedFormSpectrum_MatchesPerBlockBuilderEvd(int N, int m)
    {
        // For (p_c = 0, p_r = m), the closed-form eigenvalues {Σ_{k∈S} β_k : |S|=m} must
        // equal the eigenvalues of the per-block dense Liouvillian as a multiset to FP
        // precision. XY-chain Liouvillian spectra are highly degenerate so we use the
        // greedy nearest-neighbour multiset matcher.
        const double gamma = 0.05;
        const double J = 1.0;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);
        var sector = decomp.SectorRanges.First(s => s.PCol == 0 && s.PRow == m);
        var flat = new int[sector.Size];
        for (int k = 0; k < sector.Size; k++) flat[k] = decomp.Permutation[sector.Offset + k];

        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaArr, flat);
        var blockEigs = block.Evd().EigenValues.ToArray();

        var closed = OneSidedSectorClosedForm.Build(N, m, J, gamma);
        Assert.Equal(blockEigs.Length, closed.SectorDim);

        MultisetAssert.NearestNeighbourEqual(
            closed.Eigenvalues, blockEigs, tolerance: 1e-10, context: $"N={N}, m={m}");
        _out.WriteLine($"N={N}, (0,{m}): {closed.SectorDim} eigenvalues, closed-form vs dense Evd match within 1e-10");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void RapiditiesMatchSineModeDispersion(int N)
    {
        // β_k = −2γ − i·ε_k where ε_k comes from XyJordanWignerModes.
        const double gamma = 0.07;
        const double J = 1.2;
        var modes = XyJordanWignerModes.Build(N, J);
        var closed = OneSidedSectorClosedForm.Build(N, m: 1, J, gamma);

        Assert.Equal(N, closed.Rapidities.Length);
        for (int k = 0; k < N; k++)
        {
            Assert.Equal(-2.0 * gamma, closed.Rapidities[k].Real, 12);
            Assert.Equal(-modes.Dispersion[k], closed.Rapidities[k].Imaginary, 12);
        }
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    public void XGlobalConjugationPair_LinksZeroMSectorToNMinusMNSectorAsEqualSpectra(int N, int m)
    {
        // X⊗N global charge conjugation maps |ψ⟩⟨φ| ↦ |¬ψ⟩⟨¬φ| (bitwise complement),
        // sending (p_c, p_r) → (N − p_c, N − p_r). For chain XY: X⊗N commutes with H
        // (X_l X_{l+1} unchanged; Y_l Y_{l+1} acquires (−1)² = +1). For Z-dephasing: each
        // X⊗N Z_l X⊗N = −Z_l, so each L_l → −L_l, leaving the dissipator Z_l·Z_l unchanged.
        // Therefore L commutes with X⊗N as a superoperator, and the paired sectors have
        // EQUAL spectra (not F1-mirrored — F1 is a separate, Π-mediated identity that
        // does not restrict to a simple sector permutation).
        const double gamma = 0.05;
        const double J = 1.0;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);

        var sectorPaired = decomp.SectorRanges.First(s => s.PCol == N && s.PRow == N - m);
        var flatPaired = new int[sectorPaired.Size];
        for (int k = 0; k < sectorPaired.Size; k++) flatPaired[k] = decomp.Permutation[sectorPaired.Offset + k];
        var blockPaired = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaArr, flatPaired);
        var pairedEigs = blockPaired.Evd().EigenValues.ToArray();

        var closed = OneSidedSectorClosedForm.Build(N, m, J, gamma);

        Assert.Equal(closed.Eigenvalues.Length, pairedEigs.Length);
        MultisetAssert.NearestNeighbourEqual(
            closed.Eigenvalues, pairedEigs, tolerance: 1e-10, context: $"N={N}, m={m} X⊗N pair");
        _out.WriteLine($"N={N}, m={m}: (0,{m}) closed form equals ({N},{N - m}) dense Evd (X⊗N pair) within 1e-10");
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    public void F1Mirror_OfZeroMEigenvalues_AppearInFullSpectrum(int N, int m)
    {
        // F1 (Π-conjugation theorem): each λ in the full spectrum pairs with −λ − 2·Σγ.
        // The F1 mirrors of (0, m) eigenvalues need NOT live in any single sector — Π acts
        // on Pauli strings, not on the (p_c, p_r) computational-basis labels, so F1 mirrors
        // distribute across sectors. This test verifies the F1 prediction at the full-
        // spectrum level: every −λ − 2·Σγ for λ ∈ (0, m) closed form appears somewhere in
        // the union of all dense per-block Evds.
        const double gamma = 0.05;
        const double J = 1.0;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();
        double sumGamma = N * gamma;

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);

        var allEigs = new List<Complex>();
        foreach (var sec in decomp.SectorRanges)
        {
            if (sec.Size == 0) continue;
            var flat = new int[sec.Size];
            for (int k = 0; k < sec.Size; k++) flat[k] = decomp.Permutation[sec.Offset + k];
            var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaArr, flat);
            allEigs.AddRange(block.Evd().EigenValues.ToArray());
        }

        var closed = OneSidedSectorClosedForm.Build(N, m, J, gamma);
        foreach (var lam in closed.Eigenvalues)
        {
            var f1Mirror = new Complex(-lam.Real - 2.0 * sumGamma, -lam.Imaginary);
            double minDist = allEigs.Min(e => (e - f1Mirror).Magnitude);
            Assert.True(minDist < 1e-10,
                $"F1 mirror of (0,{m}) eigenvalue λ={lam} predicted at {f1Mirror} but " +
                $"nearest in full spectrum is at distance {minDist:G3}");
        }
        _out.WriteLine($"N={N}, m={m}: all {closed.SectorDim} F1 mirrors of (0,{m}) closed form appear in the full spectrum within 1e-10");
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 3)]
    public void ConjugateSister_MZeroSectorIsComplexConjugateOfZeroMSector(int N, int m)
    {
        // The sister sector (p_c = m, p_r = 0) corresponds to operators |0⟩⟨ψ| with ψ
        // having m fermions. Its Liouvillian is the complex conjugate of the (0, m) sector
        // (since the Hamiltonian acts on the opposite side and the dissipator is the same).
        // Hence its eigenvalues are the complex conjugates of (0, m)'s.
        const double gamma = 0.05;
        const double J = 1.0;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);
        var sectorSister = decomp.SectorRanges.First(s => s.PCol == m && s.PRow == 0);
        var flatSister = new int[sectorSister.Size];
        for (int k = 0; k < sectorSister.Size; k++) flatSister[k] = decomp.Permutation[sectorSister.Offset + k];
        var blockSister = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaArr, flatSister);
        var sisterEigs = blockSister.Evd().EigenValues.ToArray();

        var closed = OneSidedSectorClosedForm.Build(N, m, J, gamma);
        var predictedSister = closed.Eigenvalues.Select(Complex.Conjugate).ToArray();

        MultisetAssert.NearestNeighbourEqual(
            predictedSister, sisterEigs, tolerance: 1e-10, context: $"N={N}, m={m} sister sector");
        _out.WriteLine($"N={N}, m={m}: (0,{m}) closed form conjugated equals ({m},0) dense Evd within 1e-10");
    }

    [Fact]
    public void Build_AtN10_ProducesAllOneSidedSectorsClosedForm()
    {
        // Reconnaissance: at N=10 the full set of (0, m) sectors for m=0..10 collectively
        // has Σ C(10, m) = 2^10 = 1024 eigenvalues, available in closed form without any
        // diagonalization. Combined with sister (m, 0) sectors (conjugates), F1 mirrors at
        // (10, 10−m), and conjugates of mirrors, four 1024-eigenvalue families are accessible
        // analytically — bypassing the dense-Evd-infeasible (5, 5) bottleneck for these
        // specific spectral regions.
        const int N = 10;
        const double J = 1.0;
        const double gamma = 0.05;

        long totalClosedForm = 0;
        for (int m = 0; m <= N; m++)
        {
            var closed = OneSidedSectorClosedForm.Build(N, m, J, gamma);
            totalClosedForm += closed.SectorDim;
            _out.WriteLine($"  m={m}: C(10,{m})={closed.SectorDim} eigenvalues, " +
                           $"Re λ = {-2 * gamma * m:F3} (constant)");
        }
        _out.WriteLine($"Total (0, m) closed-form eigenvalues at N=10: {totalClosedForm} (= 2^N = {1L << N})");
        Assert.Equal(1L << N, totalClosedForm);
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var closed = OneSidedSectorClosedForm.Build(N: 4, m: 2);
        Assert.Equal(Tier.Tier1Derived, closed.Tier);
    }

    [Fact]
    public void Build_RejectsInvalidArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OneSidedSectorClosedForm.Build(N: 0, m: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => OneSidedSectorClosedForm.Build(N: 4, m: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => OneSidedSectorClosedForm.Build(N: 4, m: 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => OneSidedSectorClosedForm.Build(N: 4, m: 2, gamma: -0.01));
    }
}
