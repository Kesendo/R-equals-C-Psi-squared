using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using GM = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;   // not M: the class has a method M(k)

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>The levels of one spin I = TwoI/2, level index k ↔ m = I − k (index 0 is m = I), and its
/// spin operators in two forms. Exact, in the gauge Φ(A) = G·A·G⁻¹ with G = diag(σ_k), σ₀ = 1,
/// σ_k = σ_{k−1}·a(m_k), a(m) = √(I(I + 1) − m(m + 1)) = ⟨m + 1|I₊|m⟩: there I₊ is ones on the
/// superdiagonal and I₋ carries a(m)² on the subdiagonal, so every entry is rational at every I
/// (the kind of similarity <c>SodiumChainWitness</c> uses for its √3 entries, and
/// <c>SodiumTwoMirrorsWitness</c> for its control spins). Float, in the physical basis, where I₊
/// carries a(m) itself.</summary>
public sealed class SpinSpace
{
    private SpinSpace(int twoI) => TwoI = twoI;

    public static SpinSpace Of(int twoI)
    {
        if (twoI < 1) throw new ArgumentOutOfRangeException(nameof(twoI), twoI, "twice the spin, at least 1.");
        return new SpinSpace(twoI);
    }

    public int TwoI { get; }
    public int Levels => TwoI + 1;
    public BigRational Spin => new(TwoI, 2);

    /// <summary>I(I + 1).</summary>
    public BigRational SpinSquared => Spin * (Spin + 1);

    /// <summary>m at level index k.</summary>
    public BigRational M(int k) => new(TwoI - 2 * k, 2);

    /// <summary>a(m)² = I(I + 1) − m(m + 1) for m = M(k), the squared entry ⟨m + 1|I₊|m⟩ that takes
    /// level k to level k − 1; defined for k = 1 … Levels − 1.</summary>
    public BigRational RaiseSquared(int k)
    {
        if (k < 1 || k >= Levels) throw new ArgumentOutOfRangeException(nameof(k));
        var m = M(k);
        return SpinSquared - m * (m + 1);
    }

    /// <summary>σ_k² = a(m₁)²·…·a(m_k)², the squared gauge weight of level k (σ₀ = 1).</summary>
    public BigRational GaugeWeightSquared(int k)
    {
        var w = BigRational.One;
        for (int j = 1; j <= k; j++) w *= RaiseSquared(j);
        return w;
    }

    public GaussianRational[,] GaugedIz()
    {
        var m = GM.Zero(Levels);
        for (int k = 0; k < Levels; k++) m[k, k] = M(k);
        return m;
    }

    public GaussianRational[,] GaugedIplus()
    {
        var m = GM.Zero(Levels);
        for (int k = 1; k < Levels; k++) m[k - 1, k] = GaussianRational.One;
        return m;
    }

    public GaussianRational[,] GaugedIminus()
    {
        var m = GM.Zero(Levels);
        for (int k = 1; k < Levels; k++) m[k, k - 1] = RaiseSquared(k);
        return m;
    }

    /// <summary>C = Σ_a ad(I_a)² = ad(I_z)² + ½(ad(I₊)ad(I₋) + ad(I₋)ad(I₊)), ad(A)ρ = Aρ − ρA, the Casimir
    /// of the rotations ρ ↦ UρU† on the row-major vec of the d × d matrices, in the gauge: ad of a gauged
    /// matrix is Φ·ad(A)·Φ⁻¹, so every polynomial identity between C and a gauged superoperator holds in
    /// the physical basis as well. The matrices of one spin hold the spherical-tensor ranks 0 … 2I, each
    /// once, and C is k(k + 1) on rank k.</summary>
    public GaussianRational[,] GaugedRotationCasimir()
    {
        var one = GM.Identity(Levels);
        GaussianRational[,] Ad(GaussianRational[,] a) => GM.Subtract(GM.Sandwich(a, one), GM.Sandwich(one, a));
        var z = Ad(GaugedIz());
        var plus = Ad(GaugedIplus());
        var minus = Ad(GaugedIminus());
        return GM.Add(GM.Multiply(z, z), GM.Scale(GM.Add(GM.Multiply(plus, minus), GM.Multiply(minus, plus)), new BigRational(1, 2)));
    }

    /// <summary>The nullity of <see cref="GaugedRotationCasimir"/> − r(r + 1)·1 for every rank r = 0 … 2I, exactly (the
    /// Casimir is real in the gauge; r is a rank here, not a level index). Nullities 2r + 1, which add up to d², make C diagonalizable with exactly these
    /// eigenspaces, one per rank.</summary>
    public int[] RotationCasimirNullities()
    {
        var c = GaugedRotationCasimir();
        int n = c.GetLength(0);
        var nullities = new int[TwoI + 1];
        for (int rank = 0; rank <= TwoI; rank++)
        {
            var a = new BigRational[n, n];
            for (int r = 0; r < n; r++)
                for (int col = 0; col < n; col++)
                {
                    if (!c[r, col].Im.IsZero) throw new InvalidOperationException("the Casimir has an imaginary entry.");
                    a[r, col] = c[r, col].Re;
                }
            for (int i = 0; i < n; i++) a[i, i] -= rank * (rank + 1);
            nullities[rank] = BigRationalLinearAlgebra.Nullspace(a).Count;
        }
        return nullities;
    }

    public Matrix<Complex> Iz()
    {
        var m = Matrix<Complex>.Build.Dense(Levels, Levels);
        for (int k = 0; k < Levels; k++) m[k, k] = ToDouble(M(k));
        return m;
    }

    public Matrix<Complex> Iplus()
    {
        var m = Matrix<Complex>.Build.Dense(Levels, Levels);
        for (int k = 1; k < Levels; k++) m[k - 1, k] = Math.Sqrt(ToDouble(RaiseSquared(k)));
        return m;
    }

    public Matrix<Complex> Iminus() => Iplus().Transpose();

    public static double ToDouble(BigRational r) => (double)r.Numerator / (double)r.Denominator;
}
