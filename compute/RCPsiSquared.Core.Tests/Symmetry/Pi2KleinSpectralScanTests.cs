using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Schicht 3 scan across N = 1..N_max for the truly XY chain (XX+YY) under uniform
/// Z-dephasing. Locks the cross-N pattern of the Klein-cell decomposition of L's spectrum.
///
/// <para>Observed N-invariant structural locks:</para>
/// <list type="bullet">
///   <item><b>Kernel dim = N + 1</b>: span{I, Z_total, Z_total², ..., Z_total^N} commutes
///     with the truly XY Hamiltonian and the Z-dephasing dissipator. Z_total has N+1
///     distinct eigenvalues, giving N+1 independent kernel operators.</item>
///   <item><b>Kernel mass lives in Π²_X = +1 only</b>: kernel modes are I,Z-only operators
///     (bit_a = 0 always → Π²_X = +1). Pp and Mp carry full kernel mass; Pm and Mm carry zero.</item>
/// </list>
///
/// <para>Observed N-dependent locks:</para>
/// <list type="bullet">
///   <item><b>Slow non-kernel modes at small N (2, 3) preserve Π²_X-bilinear apex</b>:
///     the Pp ≈ Mp and Pm ≈ Mm balance (p ≈ 1/2) holds.</item>
///   <item><b>Slow non-kernel modes at N ≥ 4 concentrate in Π²_X = −1 only</b>: Pp + Mp ≈ 0,
///     all the mass goes to Pm + Mm. The bilinear apex still holds within Π²_X = −1
///     (Pm ≈ Mm), but the Π²_X = +1 axis is empty in the slowest non-kernel window.</item>
/// </list>
///
/// <para>Computational limit: dense eigendecomp of L is 4^N × 4^N. N=6 takes ~7 minutes
/// per eigendecomp; N=7 is borderline (~30+ minutes); N=8 (4 billion entries, 32 GB) is
/// out of reach without sparse Krylov methods. Scan stops at N=5 for a practical CI test.
/// N=1 is the trivial degenerate case (kernel is full operator space).</para>
/// </summary>
public class Pi2KleinSpectralScanTests
{
    private static ComplexMatrix BuildL(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    private static PauliPairBondTerm[] TrulyXYTerms => new[]
    {
        new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
        new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
    };

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void TrulyXYChain_KernelDim_EqualsNPlus1(int N)
    {
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var L = BuildL(TrulyXYTerms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N);

        int kernelDim = modes.Count(m => Math.Abs(m.Eigenvalue.Real) < 1e-9);
        Assert.True(kernelDim == N + 1, $"N={N}: expected kernel dim {N + 1}, got {kernelDim}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void TrulyXYChain_KernelMass_LivesInPi2XPlus_Only(int N)
    {
        // Kernel is span{I, Z_total^k} → I,Z-only operators → bit_a = 0 → Π²_X = +1.
        // All kernel mass in Pp + Mp; zero in Pm + Mm.
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var L = BuildL(TrulyXYTerms, chain);
        var kernel = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .Where(m => Math.Abs(m.Eigenvalue.Real) < 1e-9)
            .ToList();

        double aggregatePm = kernel.Sum(m => m.MassPm);
        double aggregateMm = kernel.Sum(m => m.MassMm);
        double total = kernel.Sum(m => m.TotalMass);

        Assert.True(aggregatePm < 1e-6,
            $"N={N}: kernel Pm mass {aggregatePm:E3} should be ≈ 0 (kernel is I,Z-only)");
        Assert.True(aggregateMm < 1e-6,
            $"N={N}: kernel Mm mass {aggregateMm:E3} should be ≈ 0");
        Assert.True(Math.Abs(total - (N + 1)) < 1e-6,
            $"N={N}: kernel total mass {total:F4} should equal N+1 = {N + 1}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void TrulyXYChain_SmallN_SlowNonKernel_PreservesPi2XApex(int N)
    {
        // At small N (2, 3), slow non-kernel modes preserve the Π²_X-axis bilinear apex:
        // Pp ≈ Mp and Pm ≈ Mm. The framework's recurring 1/2 manifests both axes.
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var L = BuildL(TrulyXYTerms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();
        var slowNonKernel = modes.Skip(N + 1).Take(4 * (N + 1)).ToList();

        double aggPp = slowNonKernel.Sum(m => m.MassPp);
        double aggPm = slowNonKernel.Sum(m => m.MassPm);
        double aggMp = slowNonKernel.Sum(m => m.MassMp);
        double aggMm = slowNonKernel.Sum(m => m.MassMm);

        double rPp = aggPp / Math.Max(1e-12, aggPp + aggMp);
        double rPm = aggPm / Math.Max(1e-12, aggPm + aggMm);

        Assert.True(Math.Abs(rPp - 0.5) < 0.15,
            $"N={N}: Π²_X=+ p = Pp/(Pp+Mp) = {rPp:F4} should sit near apex 1/2");
        Assert.True(Math.Abs(rPm - 0.5) < 0.15,
            $"N={N}: Π²_X=− p = Pm/(Pm+Mm) = {rPm:F4} should sit near apex 1/2");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void TrulyXYChain_LargerN_SlowNonKernel_ConcentratesInOnePi2XAxis(int N)
    {
        // At N ≥ 4 the slowest 4·(N+1) non-kernel modes concentrate in Π²_X = −1 only:
        // Pp + Mp ≈ 0, with all mass in Pm + Mm. The Π²_X = +1 axis is empty in this slow
        // window. Within Π²_X = −1, the Π²_Z bilinear apex still holds (Pm ≈ Mm).
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var L = BuildL(TrulyXYTerms, chain);
        var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
            .OrderBy(m => Math.Abs(m.Eigenvalue.Real))
            .ToList();
        var slowNonKernel = modes.Skip(N + 1).Take(4 * (N + 1)).ToList();

        double aggPp = slowNonKernel.Sum(m => m.MassPp);
        double aggPm = slowNonKernel.Sum(m => m.MassPm);
        double aggMp = slowNonKernel.Sum(m => m.MassMp);
        double aggMm = slowNonKernel.Sum(m => m.MassMm);

        // Π²_X = +1 sub-axis is essentially empty.
        Assert.True(aggPp + aggMp < 0.5,
            $"N={N}: Π²_X=+1 aggregate {aggPp + aggMp:F4} should be ≈ 0 in slow non-kernel window");

        // Within Π²_X = −1, Π²_Z bilinear apex still holds.
        double rPm = aggPm / Math.Max(1e-12, aggPm + aggMm);
        Assert.True(Math.Abs(rPm - 0.5) < 0.15,
            $"N={N}: within Π²_X=−1, p = Pm/(Pm+Mm) = {rPm:F4} should sit near apex 1/2");
    }

    [Theory]
    [InlineData(2, 1.0, false)]
    [InlineData(3, 1.5, true)]
    [InlineData(4, 2.0, false)]
    [InlineData(5, 2.5, true)]
    [InlineData(6, 3.0, false)]
    [InlineData(7, 3.5, true)]
    [InlineData(8, 4.0, false)]
    public void HalfIntegerMirrorRegime_ClassifiedByNParity(int N, double expectedWXY, bool expectedHalfIntegerRegime)
    {
        // Tom's half-integer-mirror family: w_XY = N/2. Odd N ⇒ half-integer (1/2, 3/2, ...);
        // even N ⇒ integer (1, 2, 3, ...). The framework's recurring 0.5 is the universal
        // fractional part across odd N. At N=8 (even, w_XY=4): integer-mirror regime — the
        // operator-level analysis would proceed without the "no modes on axis" benefit that
        // odd N gives.
        double wXY = N / 2.0;
        Assert.Equal(expectedWXY, wXY);
        bool actualHalfInteger = N % 2 == 1;
        Assert.Equal(expectedHalfIntegerRegime, actualHalfInteger);
    }
}
