using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Native C# verification of F94's bit-exact Dyson constants. Builds
/// the Heisenberg ring + Z-dephasing apparatus from Core/Pauli + Core/ChainSystems
/// primitives, applies the F94 sym3 expansion to ρ_0 = |0+0+⟩⟨0+0+|, partial-traces
/// to pair (0, 2), and verifies the matrix elements against F94's asserted constants.
///
/// <para>Mirrors the Python script <c>simulations/_born_rule_tier1_derivation.py</c>
/// and its enumeration sibling <c>_born_rule_sym3_decomposition.py</c>, but as a
/// self-contained C# test — no external Python dependency. If these tests pass,
/// the F94 typed-claim constants are bit-exact from-scratch verified.</para>
/// </summary>
public class F94NativeDysonVerificationTests
{
    private const int N = 4;
    private const int Dim = 1 << N;  // 16

    /// <summary>|0+0+⟩ on N=4 qubits, site 0 = msb (per Core's Pauli convention).
    /// Components: |0⟩_0 ⊗ |+⟩_1 ⊗ |0⟩_2 ⊗ |+⟩_3. Non-zero amplitudes 1/2 at
    /// indices {0000, 0001, 0100, 0101} = {0, 1, 4, 5}.</summary>
    private static ComplexVector Build_0P0P()
    {
        var psi = ComplexVector.Build.Dense(Dim);
        const double half = 0.5;
        psi[0b0000] = half;
        psi[0b0001] = half;
        psi[0b0100] = half;
        psi[0b0101] = half;
        return psi;
    }

    /// <summary>Apply L_H = −i [H, ρ] to a density matrix.</summary>
    private static ComplexMatrix ApplyLH(ComplexMatrix H, ComplexMatrix rho) =>
        -Complex.ImaginaryOne * (H * rho - rho * H);

    /// <summary>Apply single-site Z-dephasing L_l[ρ] = Z_l ρ Z_l − ρ.</summary>
    private static ComplexMatrix ApplyLDisSite(ComplexMatrix Zl, ComplexMatrix rho) =>
        Zl * rho * Zl - rho;

    /// <summary>Apply full uniform-γ Z-dephasing L'_dis = Σ_l (Z_l ρ Z_l − ρ) at γ = 1.</summary>
    private static ComplexMatrix ApplyLDis(ComplexMatrix rho, IReadOnlyList<ComplexMatrix> ZSites)
    {
        var result = Matrix<Complex>.Build.Dense(rho.RowCount, rho.ColumnCount);
        foreach (var Zl in ZSites)
            result += ApplyLDisSite(Zl, rho);
        return result;
    }

    /// <summary>Heisenberg ring N=4 at J=1: H = (1/4) Σ_b (XX + YY + ZZ).
    /// Returned via Core's ChainSystem builder for the exact same normalization
    /// as Python's heisenberg_ring(N, J=1.0).</summary>
    private static ComplexMatrix BuildHeisenbergRing() =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.0,
                        HType: HamiltonianType.Heisenberg,
                        Topology: TopologyKind.Ring).BuildHamiltonian();

    private static IReadOnlyList<ComplexMatrix> BuildZSites() =>
        Enumerable.Range(0, N)
            .Select(l => PauliString.SiteOp(N, l, PauliLetter.Z))
            .ToArray();

    /// <summary>Sym3 = L_H² L'_dis + L_H L'_dis L_H + L'_dis L_H² applied to ρ_0.
    /// The γ¹·J² coefficient of L³ in the time-Taylor expansion of e^{Lt} ρ_0.</summary>
    private static ComplexMatrix BuildSym3RhoZero(ComplexMatrix H, IReadOnlyList<ComplexMatrix> ZSites)
    {
        var rho0 = DensityMatrix.FromStateVector(Build_0P0P());

        // Ordering 1: L_H L_H L_dis ρ_0
        var a = ApplyLDis(rho0, ZSites);
        a = ApplyLH(H, a);
        a = ApplyLH(H, a);

        // Ordering 2: L_H L_dis L_H ρ_0
        var b = ApplyLH(H, rho0);
        b = ApplyLDis(b, ZSites);
        b = ApplyLH(H, b);

        // Ordering 3: L_dis L_H L_H ρ_0
        var c = ApplyLH(H, rho0);
        c = ApplyLH(H, c);
        c = ApplyLDis(c, ZSites);

        return a + b + c;
    }

    [Fact]
    public void NativeDerivation_Sym3PairElement_Equals_F94_Constant()
    {
        // ⟨00|_pair Tr_{1,3}[sym3 · ρ_0] |00⟩_pair = 8 bit-exact.
        // This is F94.Sym3PartialTraceInteger, derived natively in C#.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3Rho0 = BuildSym3RhoZero(H, ZSites);
        var reduced = PartialTrace.Of(sym3Rho0, N, new[] { 0, 2 });

        double matrixElement = reduced[0, 0].Real;
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.Sym3PartialTraceInteger,
            matrixElement, precision: 10);
        Assert.True(Math.Abs(reduced[0, 0].Imaginary) < 1e-10,
            $"Matrix element must be real; got imag = {reduced[0, 0].Imaginary}");
    }

    [Fact]
    public void NativeDerivation_Coefficient_Is4Over3()
    {
        // Combining Sym3 = 8 with Taylor 3! = 6 gives the F94 coefficient 4/3.
        var H = BuildHeisenbergRing();
        var ZSites = BuildZSites();
        var sym3Rho0 = BuildSym3RhoZero(H, ZSites);
        var reduced = PartialTrace.Of(sym3Rho0, N, new[] { 0, 2 });
        double coefficient = reduced[0, 0].Real
            / F94BornDeviationFourThirdsPi2Inheritance.TaylorThreeFactorial;
        Assert.Equal(4.0 / 3.0, coefficient, precision: 10);
    }

    [Fact]
    public void NativeDerivation_PunitaryAt00IsOne()
    {
        // ⟨00|_pair Tr_{1,3}[ρ_0] |00⟩_pair = 1 (Bell+-like sites 0 and 2 are |0⟩ deterministically).
        var rho0 = DensityMatrix.FromStateVector(Build_0P0P());
        var reduced = PartialTrace.Of(rho0, N, new[] { 0, 2 });
        Assert.Equal(1.0, reduced[0, 0].Real, precision: 12);
    }

    // ────────────────────────────────────────────────────────────────────
    // F94 enumeration: 32 surviving (b1, b2, s, ord, c1, c2) sextuples
    // ────────────────────────────────────────────────────────────────────

    private static readonly (int A, int B)[] _bonds = { (0, 1), (1, 2), (2, 3), (3, 0) };
    private static readonly PauliLetter[] _components = { PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    /// <summary>(J/4) · σ_c_a σ_c_b for one bond, one Heisenberg component, at J = 1.</summary>
    private static ComplexMatrix BondComponentHamiltonian(int siteA, int siteB, PauliLetter c)
    {
        var sA = PauliString.SiteOp(N, siteA, c);
        var sB = PauliString.SiteOp(N, siteB, c);
        return 0.25 * (sA * sB);
    }

    /// <summary>Evaluate the |00⟩_pair element of one (b1, b2, s, ord, c1, c2) sextuple in sym3.
    /// Returns the real part (the imaginary part is asserted ≤ 1e-12).</summary>
    private static double EvaluateSextuple(
        ComplexMatrix rho0, IReadOnlyList<ComplexMatrix> ZSites,
        int b1, int b2, int s, int ordering, PauliLetter c1, PauliLetter c2)
    {
        var Hb1 = BondComponentHamiltonian(_bonds[b1].A, _bonds[b1].B, c1);
        var Hb2 = BondComponentHamiltonian(_bonds[b2].A, _bonds[b2].B, c2);
        var Zs = ZSites[s];

        ComplexMatrix x = ordering switch
        {
            1 => ApplyLH(Hb1, ApplyLH(Hb2, ApplyLDisSite(Zs, rho0))),
            2 => ApplyLH(Hb1, ApplyLDisSite(Zs, ApplyLH(Hb2, rho0))),
            3 => ApplyLDisSite(Zs, ApplyLH(Hb1, ApplyLH(Hb2, rho0))),
            _ => throw new ArgumentOutOfRangeException(nameof(ordering)),
        };
        var reduced = PartialTrace.Of(x, N, new[] { 0, 2 });
        var val = reduced[0, 0];
        if (Math.Abs(val.Imaginary) > 1e-10)
            throw new InvalidOperationException(
                $"Non-real pair element at (b1={b1}, b2={b2}, s={s}, ord={ordering}, c1={c1}, c2={c2}): {val}");
        return val.Real;
    }

    private record SextupleResult(int B1, int B2, int Site, int Ordering, PauliLetter C1, PauliLetter C2, double Value);

    private static List<SextupleResult> EnumerateSurvivors()
    {
        var rho0 = DensityMatrix.FromStateVector(Build_0P0P());
        var ZSites = BuildZSites();
        var survivors = new List<SextupleResult>();
        const double tol = 1e-10;
        for (int b1 = 0; b1 < _bonds.Length; b1++)
            for (int b2 = 0; b2 < _bonds.Length; b2++)
                for (int s = 0; s < N; s++)
                    for (int ordering = 1; ordering <= 3; ordering++)
                        foreach (var c1 in _components)
                            foreach (var c2 in _components)
                            {
                                double v = EvaluateSextuple(rho0, ZSites, b1, b2, s, ordering, c1, c2);
                                if (Math.Abs(v) > tol)
                                    survivors.Add(new SextupleResult(b1, b2, s, ordering, c1, c2, v));
                            }
        return survivors;
    }

    [Fact]
    public void NativeEnumeration_SurvivorCount_Equals32()
    {
        // F94.SurvivingDysonDiagrams = 32 verified by direct enumeration of all
        // 4·4·4·3·3·3 = 1728 (b1, b2, s, ord, c1, c2) sextuples in sym3.
        var survivors = EnumerateSurvivors();
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.SurvivingDysonDiagrams,
            survivors.Count);
    }

    [Fact]
    public void NativeEnumeration_AllSurvivorsContribute_OneQuarter()
    {
        // Every surviving diagram contributes exactly 1/4 (= (J/4)² · 4 raw) in
        // the J = γ = 1 normalization. This is the pure-counting form of F94:
        // no signs, no cancellation; 32 × (1/4) = 8.
        var survivors = EnumerateSurvivors();
        foreach (var s in survivors)
            Assert.Equal(0.25, s.Value, precision: 10);
        double total = survivors.Sum(s => s.Value);
        Assert.Equal((double)F94BornDeviationFourThirdsPi2Inheritance.Sym3PartialTraceInteger,
            total, precision: 10);
    }

    [Fact]
    public void NativeEnumeration_CellCountsMatchF94Constants()
    {
        // Cell A: ord = 1, (X, X) → 8 diagrams
        // Cell B: ord = 2, (X, X) → 16 diagrams
        // Cell C: ord = 2, (Y, Y) → 8 diagrams
        // 8 + 16 + 8 = 32 (and 8/16/8 matches the F94 CellA/B/C constants).
        var survivors = EnumerateSurvivors();
        int cellA = survivors.Count(s => s.Ordering == 1 && s.C1 == PauliLetter.X && s.C2 == PauliLetter.X);
        int cellB = survivors.Count(s => s.Ordering == 2 && s.C1 == PauliLetter.X && s.C2 == PauliLetter.X);
        int cellC = survivors.Count(s => s.Ordering == 2 && s.C1 == PauliLetter.Y && s.C2 == PauliLetter.Y);
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.CellA_Ord1XX_AdjKeptSide, cellA);
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.CellB_Ord2XX_SelfOrAdjKeptSide, cellB);
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.CellC_Ord2YY_Self, cellC);
        Assert.Equal(F94BornDeviationFourThirdsPi2Inheritance.SurvivingDysonDiagrams,
            cellA + cellB + cellC);
    }

    [Fact]
    public void NativeEnumeration_NoSurvivorsOutsideThreeCells()
    {
        // Beyond the 3 cells above, no (ordering, c1, c2) combination has any
        // surviving diagram. Component-pair rule: only (X,X) and (Y,Y) survive
        // (no cross, no Z,Z); ordering rule: ord = 3 (L_dis last) is empty.
        var survivors = EnumerateSurvivors();
        foreach (var s in survivors)
        {
            bool inCellA = s.Ordering == 1 && s.C1 == PauliLetter.X && s.C2 == PauliLetter.X;
            bool inCellB = s.Ordering == 2 && s.C1 == PauliLetter.X && s.C2 == PauliLetter.X;
            bool inCellC = s.Ordering == 2 && s.C1 == PauliLetter.Y && s.C2 == PauliLetter.Y;
            Assert.True(inCellA || inCellB || inCellC,
                $"Survivor outside the 3 F94 cells: ord={s.Ordering}, c1={s.C1}, c2={s.C2}");
        }
    }
}
