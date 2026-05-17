using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Native C# verification of F96's bit-exact Dyson matrix elements
/// (M_3 = −4 for |01⟩, M_5 = −20 for |11⟩, U_2 = 3/4 for |01⟩, U_4 = 3/2 for
/// |11⟩) on the same setup as F94. Reuses Core/Pauli + Core/ChainSystems +
/// Core/States/PartialTrace primitives.
///
/// <para>Mirrors <c>simulations/_born_rule_subdominant_dyson.py</c>. Where F94's
/// native test verifies the dominant outcome's sym3 element (= 8 on |00⟩), F96's
/// verifies (a) the singly-subdominant sym3 element (= −4 on |01⟩) and U_h²
/// element (= 3/4 on |01⟩), and (b) the doubly-subdominant sym5 element
/// (= −20 on |11⟩) and L_h⁴ element (= 3/2 on |11⟩). From these the F96 slopes
/// −16/9 and −8/3 follow algebraically.</para>
/// </summary>
public class F96NativeDysonVerificationTests
{
    private const int N = 4;

    private static ComplexVector Build_0P0P()
    {
        var psi = ComplexVector.Build.Dense(1 << N);
        const double half = 0.5;
        psi[0b0000] = half;
        psi[0b0001] = half;
        psi[0b0100] = half;
        psi[0b0101] = half;
        return psi;
    }

    private static ComplexMatrix ApplyLH(ComplexMatrix H, ComplexMatrix rho) =>
        -Complex.ImaginaryOne * (H * rho - rho * H);

    private static ComplexMatrix ApplyLDis(ComplexMatrix rho, IReadOnlyList<ComplexMatrix> ZSites)
    {
        var result = Matrix<Complex>.Build.Dense(rho.RowCount, rho.ColumnCount);
        foreach (var Zl in ZSites)
            result += Zl * rho * Zl - rho;
        return result;
    }

    private static ComplexMatrix BuildHeisenbergRing() =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.0,
                        HType: HamiltonianType.Heisenberg,
                        Topology: TopologyKind.Ring).BuildHamiltonian();

    private static IReadOnlyList<ComplexMatrix> BuildZSites() =>
        Enumerable.Range(0, N)
            .Select(l => PauliString.SiteOp(N, l, PauliLetter.Z))
            .ToArray();

    /// <summary>sym3 = L_H² L_dis + L_H L_dis L_H + L_dis L_H² applied to ρ_0.</summary>
    private static ComplexMatrix BuildSym3RhoZero(ComplexMatrix H, IReadOnlyList<ComplexMatrix> ZSites)
    {
        var rho0 = DensityMatrix.FromStateVector(Build_0P0P());

        var a = ApplyLDis(rho0, ZSites);
        a = ApplyLH(H, a);
        a = ApplyLH(H, a);

        var b = ApplyLH(H, rho0);
        b = ApplyLDis(b, ZSites);
        b = ApplyLH(H, b);

        var c = ApplyLH(H, rho0);
        c = ApplyLH(H, c);
        c = ApplyLDis(c, ZSites);

        return a + b + c;
    }

    /// <summary>L_H^k applied to ρ_0 (no dissipator). For F96 we need k = 2 and k = 4
    /// to extract U_2 and U_4 matrix elements on the unitary Taylor expansion.</summary>
    private static ComplexMatrix BuildLHPowerRhoZero(ComplexMatrix H, int k)
    {
        var rho = DensityMatrix.FromStateVector(Build_0P0P());
        for (int i = 0; i < k; i++) rho = ApplyLH(H, rho);
        return rho;
    }

    /// <summary>sym5 = sum over 5 orderings of (L_H L_H L_H L_H L_dis), applied to ρ_0.
    /// γ¹·J⁴·t⁵ Dyson term — the leading non-vanishing γ¹ contribution for the
    /// doubly-subdominant |11⟩ outcome (where sym3 vanishes by parity).</summary>
    private static ComplexMatrix BuildSym5RhoZero(ComplexMatrix H, IReadOnlyList<ComplexMatrix> ZSites)
    {
        var rho0 = DensityMatrix.FromStateVector(Build_0P0P());

        // Position of L_dis in the ordering (0 = first applied, 4 = last applied).
        // For each position, apply L_dis at that step and L_H at the others.
        var total = Matrix<Complex>.Build.Dense(rho0.RowCount, rho0.ColumnCount);
        for (int pos = 0; pos < 5; pos++)
        {
            var x = rho0.Clone();
            for (int step = 0; step < 5; step++)
            {
                if (step == pos) x = ApplyLDis(x, ZSites);
                else x = ApplyLH(H, x);
            }
            total += x;
        }
        return total;
    }

    private static double PairElement_01(ComplexMatrix x) =>
        PairElement(x, row: 1, col: 1);

    private static double PairElement_11(ComplexMatrix x) =>
        PairElement(x, row: 3, col: 3);

    private static double PairElement(ComplexMatrix x, int row, int col)
    {
        var reduced = PartialTrace.Of(x, N, new[] { 0, 2 });
        var val = reduced[row, col];
        if (Math.Abs(val.Imaginary) > 1e-10)
            throw new InvalidOperationException(
                $"Non-real pair element at ({row},{col}): {val}");
        return val.Real;
    }

    [Fact]
    public void NativeDerivation_M3_SingleFlipped_EqualsMinusFour()
    {
        // ⟨01|_pair Tr_{1,3}[sym3 · ρ_0] |01⟩_pair = -4 bit-exact.
        // This is F96.M3_SingleFlipped, derived natively in C#.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3 = BuildSym3RhoZero(H, ZSites);
        double m3 = PairElement_01(sym3);
        Assert.Equal((double)F96BornSubdominantSlopesPi2Inheritance.M3_SingleFlipped,
            m3, precision: 10);
    }

    [Fact]
    public void NativeDerivation_U2_SingleFlipped_EqualsThreeQuarters()
    {
        // ⟨01|_pair Tr_{1,3}[L_h² · ρ_0] |01⟩_pair = 3/4 bit-exact.
        // F96 stores this as U2_SingleFlipped_TimesFour = 3 (implicit denom 4).
        var H = BuildHeisenbergRing();
        var lh2 = BuildLHPowerRhoZero(H, 2);
        double u2 = PairElement_01(lh2);
        double u2Expected = (double)F96BornSubdominantSlopesPi2Inheritance.U2_SingleFlipped_TimesFour / 4.0;
        Assert.Equal(u2Expected, u2, precision: 10);
        Assert.Equal(0.75, u2, precision: 10);
    }

    [Fact]
    public void NativeDerivation_M5_DoubleFlipped_EqualsMinusTwenty()
    {
        // ⟨11|_pair Tr_{1,3}[sym5 · ρ_0] |11⟩_pair = -20 bit-exact.
        // F96.M5_DoubleFlipped.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym5 = BuildSym5RhoZero(H, ZSites);
        double m5 = PairElement_11(sym5);
        Assert.Equal((double)F96BornSubdominantSlopesPi2Inheritance.M5_DoubleFlipped,
            m5, precision: 9);
    }

    [Fact]
    public void NativeDerivation_U4_DoubleFlipped_EqualsThreeHalves()
    {
        // ⟨11|_pair Tr_{1,3}[L_h⁴ · ρ_0] |11⟩_pair = 3/2 bit-exact.
        // F96 stores this as U4_DoubleFlipped_TimesTwo = 3 (implicit denom 2).
        var H = BuildHeisenbergRing();
        var lh4 = BuildLHPowerRhoZero(H, 4);
        double u4 = PairElement_11(lh4);
        double u4Expected = (double)F96BornSubdominantSlopesPi2Inheritance.U4_DoubleFlipped_TimesTwo / 2.0;
        Assert.Equal(u4Expected, u4, precision: 10);
        Assert.Equal(1.5, u4, precision: 10);
    }

    [Fact]
    public void NativeDerivation_M3_DoubleFlipped_IsZero()
    {
        // ⟨11|_pair Tr_{1,3}[sym3 · ρ_0] |11⟩_pair = 0 (sym3 vanishes for |11⟩;
        // this is why M_5 + U_4 take over for the doubly-subdominant slope).
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3 = BuildSym3RhoZero(H, ZSites);
        double m3_11 = PairElement_11(sym3);
        Assert.True(Math.Abs(m3_11) < 1e-10,
            $"|sym3 on |11⟩| must vanish; got {m3_11}");
    }

    [Fact]
    public void NativeDerivation_U2_DoubleFlipped_IsZero()
    {
        // ⟨11|_pair Tr_{1,3}[L_h² · ρ_0] |11⟩_pair = 0 (P_u(|11⟩) starts at t⁴,
        // not t²; this forces F96 to use the next-higher unitary order U_4).
        var H = BuildHeisenbergRing();
        var lh2 = BuildLHPowerRhoZero(H, 2);
        double u2_11 = PairElement_11(lh2);
        Assert.True(Math.Abs(u2_11) < 1e-10,
            $"|L_h² on |11⟩| must vanish; got {u2_11}");
    }

    [Fact]
    public void NativeDerivation_SlopeSingleFlipped_EqualsMinus16Over9()
    {
        // Combining native M_3 and U_2 via F96's universal slope formula
        // slope_i = M_{2k+1} / ((2k+1) · U_{2k}) at k = 1:
        // slope = -4 / (3 · 3/4) = -16/9. Self-verifies F96.SlopeSingleFlipped.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3 = BuildSym3RhoZero(H, ZSites);
        var lh2 = BuildLHPowerRhoZero(H, 2);
        double m3 = PairElement_01(sym3);
        double u2 = PairElement_01(lh2);
        double slope = m3 / (3.0 * u2);
        Assert.Equal(-16.0 / 9.0, slope, precision: 10);
    }

    [Fact]
    public void NativeDerivation_SlopeDoubleFlipped_EqualsMinus8Over3()
    {
        // slope_|11⟩ = M_5 / (5 · U_4) at k = 2:
        // slope = -20 / (5 · 3/2) = -8/3. Self-verifies F96.SlopeDoubleFlipped.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym5 = BuildSym5RhoZero(H, ZSites);
        var lh4 = BuildLHPowerRhoZero(H, 4);
        double m5 = PairElement_11(sym5);
        double u4 = PairElement_11(lh4);
        double slope = m5 / (5.0 * u4);
        Assert.Equal(-8.0 / 3.0, slope, precision: 10);
    }

    [Fact]
    public void NativeDerivation_CrossOutcomeRatio_EqualsMinus16Over3()
    {
        // F96's cross-outcome universality drift check: M_3 / U_2 = -16/3 identically
        // for both dominant (|00⟩: 8/(-3/2)) and singly-subdominant (|01⟩: -4/(3/4)).
        // Both ratios verified natively in C# from the same H and ρ_0.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3 = BuildSym3RhoZero(H, ZSites);
        var lh2 = BuildLHPowerRhoZero(H, 2);

        double m3_dominant = PartialTrace.Of(sym3, N, new[] { 0, 2 })[0, 0].Real;       // 8
        double u2_dominant = PartialTrace.Of(lh2, N, new[] { 0, 2 })[0, 0].Real;        // -3/2
        double m3_subdom = PairElement_01(sym3);                                          // -4
        double u2_subdom = PairElement_01(lh2);                                           // 3/4

        double ratio_dominant = m3_dominant / u2_dominant;
        double ratio_subdom = m3_subdom / u2_subdom;

        Assert.Equal(-16.0 / 3.0, ratio_dominant, precision: 10);
        Assert.Equal(-16.0 / 3.0, ratio_subdom, precision: 10);
        Assert.Equal(ratio_dominant, ratio_subdom, precision: 10);
    }
}
