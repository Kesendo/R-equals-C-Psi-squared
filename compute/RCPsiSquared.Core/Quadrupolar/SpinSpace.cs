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
