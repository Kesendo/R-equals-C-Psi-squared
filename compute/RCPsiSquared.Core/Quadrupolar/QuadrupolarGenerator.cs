using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>L = −i[H, ·] + Σ_m J_|m|·D[T₂ₘ] for one spin I, on the row-major vec of its d × d
/// matrices (index d·i + j for |i⟩⟨j|, level index 0 = m = I). Views: <see cref="Gauged"/>, exact
/// for every I in the gauge Φ(ρ) = GρG⁻¹ of <see cref="SpinSpace"/> (it keeps the spectrum, ranks,
/// kernels and every cell, not transpose or adjoint relations); <see cref="Physical"/>, exact in the
/// physical basis where that is rational; <see cref="Float"/>, the physical basis in floating
/// point, built from the float spin matrices (a separate route from the exact views).</summary>
public sealed class QuadrupolarGenerator
{
    public QuadrupolarGenerator(SpinSpace space, QuadrupolarRates rates, SecularHamiltonian hamiltonian)
    {
        Space = space ?? throw new ArgumentNullException(nameof(space));
        Rates = rates ?? throw new ArgumentNullException(nameof(rates));
        Hamiltonian = hamiltonian ?? throw new ArgumentNullException(nameof(hamiltonian));
    }

    public SpinSpace Space { get; }
    public QuadrupolarRates Rates { get; }
    public SecularHamiltonian Hamiltonian { get; }
    public int Dim => Space.Levels * Space.Levels;

    private static readonly BigRational Half = new(1, 2);
    private static readonly int[] OffOrders = { 1, -1, 2, -2 };
    private GaussianRational[,]? _dissipator;

    /// <summary>The gauged dissipator, built once per generator (the census asks for every order block).</summary>
    private GaussianRational[,] DissipatorCache => _dissipator ??= BuildGaugedDissipator();

    /// <summary>J₀·D[T₂₀], a diagonal superoperator, the same in both bases.</summary>
    public GaussianRational[,] J0Part()
    {
        var a = QuadrupolarTensors.Root6T20(Space);
        return M.Scale(Dissipate(a, a), Rates.J0 / 6);
    }

    /// <summary>J₀·D[T₂₀] in the physical basis: identical to <see cref="J0Part"/> (diagonal).</summary>
    public GaussianRational[,] PhysicalJ0Part() => J0Part();

    /// <summary>The dissipator in the gauge: real, exact, every I (a copy; the generator keeps its own).</summary>
    public GaussianRational[,] GaugedDissipator() => (GaussianRational[,])DissipatorCache.Clone();

    private GaussianRational[,] BuildGaugedDissipator()
    {
        var acc = J0Part();
        foreach (int q in OffOrders)
        {
            var j = Rates.ByOrder(q);
            if (j.IsZero) continue;
            var t = QuadrupolarTensors.Gauged(Space, q);
            var tDagger = M.Scale(QuadrupolarTensors.Gauged(Space, -q), QuadrupolarTensors.Parity(q));
            acc = M.Add(acc, M.Scale(Dissipate(t, tDagger), j));
        }
        return acc;
    }

    /// <summary>H in the gauge (the lock's I_x becomes (I₊ + I₋)/2 with the gauged ladders). An extra
    /// physical matrix has no gauge image here: the view refuses it.</summary>
    public GaussianRational[,] GaugedHamiltonian()
    {
        if (Hamiltonian.ExtraPhysical is not null)
            throw new InvalidOperationException("an extra physical matrix is honoured by the physical views only.");
        return GaugedSecularHamiltonian();
    }

    public GaussianRational[,] Gauged() => M.Add(CommutatorPart(GaugedHamiltonian()), DissipatorCache);

    /// <summary>L in the physical basis, exactly: the gauged view taken back through Φ⁻¹, every
    /// nonzero entry multiplied by the root of <see cref="GaugeFactorSquared"/>; where that square is not
    /// a rational square the physical basis is not rational and the view throws. The view accepts 2I ≤ 3
    /// (a sweep to 2I = 24 found no other case), and there the factor is 1 on every nonzero entry: at spin
    /// 1/2 and 1 σ is constant or geometric, at 3/2 the entries whose factor is √3/2 vanish because T₂,±₁ is
    /// zero on the central transition. So there the physical view is the gauged one plus the extra physical
    /// matrix (gated); the factor's direction is gated through the float view from 2I = 4 on.</summary>
    public GaussianRational[,] Physical()
    {
        var lg = M.Add(CommutatorPart(GaugedSecularHamiltonian()), DissipatorCache);
        var l = M.Zero(Dim);
        for (int r = 0; r < Dim; r++)
            for (int c = 0; c < Dim; c++)
            {
                if (lg[r, c].IsZero) continue;
                if (!TryRationalSqrt(GaugeFactorSquared(Space, r, c), out var f))
                    throw new NotSupportedException($"the physical basis is not rational at 2I = {Space.TwoI} with these terms.");
                l[r, c] = lg[r, c] * f;
            }
        if (Hamiltonian.ExtraPhysical is not null) l = M.Add(l, CommutatorPart(Hamiltonian.ExtraPhysical));
        return l;
    }

    /// <summary>The square of σ_u·σ_y/(σ_v·σ_x), the factor that takes the gauged entry [d·x + y, d·u + v] to
    /// the physical one (L_phys = Φ⁻¹·L_gauged·Φ, Φ scaling the cell |i⟩⟨j| by σ_i/σ_j).</summary>
    public static BigRational GaugeFactorSquared(SpinSpace s, int row, int col)
    {
        int d = s.Levels;
        int x = row / d, y = row % d, u = col / d, v = col % d;
        return s.GaugeWeightSquared(u) * s.GaugeWeightSquared(y) / (s.GaugeWeightSquared(v) * s.GaugeWeightSquared(x));
    }

    /// <summary>L in the physical basis in floating point, built from the float spin matrices:
    /// vec(AρB) = (A ⊗ Bᵀ)·vec(ρ) on the row-major vec.</summary>
    public Matrix<Complex> Float()
    {
        int d = Space.Levels;
        var one = Matrix<Complex>.Build.DenseIdentity(d);
        var iz = Space.Iz();
        double q = SpinSpace.ToDouble(Hamiltonian.Q), w = SpinSpace.ToDouble(Hamiltonian.Omega);
        var h = (q / 3) * (3 * iz * iz - SpinSpace.ToDouble(Space.SpinSquared) * one) + w * iz;
        Matrix<Complex> TdT(int m) { var t = QuadrupolarTensors.Physical(Space, m); return t.ConjugateTranspose() * t; }
        h += SpinSpace.ToDouble(Hamiltonian.K1) * (TdT(1) - TdT(-1)) + SpinSpace.ToDouble(Hamiltonian.K2) * (TdT(2) - TdT(-2));
        h += SpinSpace.ToDouble(Hamiltonian.Nu) * (Space.Iplus() + Space.Iminus()) / 2;
        if (Hamiltonian.ExtraPhysical is not null)
        {
            var e = Matrix<Complex>.Build.Dense(d, d);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    e[i, j] = new Complex(SpinSpace.ToDouble(Hamiltonian.ExtraPhysical[i, j].Re), SpinSpace.ToDouble(Hamiltonian.ExtraPhysical[i, j].Im));
            h += e;
        }
        var minusI = new Complex(0, -1);
        var l = minusI * (h.KroneckerProduct(one) - one.KroneckerProduct(h.Transpose()));
        foreach (int m in new[] { 0, 1, -1, 2, -2 })
        {
            double j = SpinSpace.ToDouble(Rates.ByOrder(m));
            if (j == 0) continue;
            var t = QuadrupolarTensors.Physical(Space, m);
            var td = t.ConjugateTranspose();
            var tdt = td * t;
            l += j * (t.KroneckerProduct(td.Transpose()) - 0.5 * (tdt.KroneckerProduct(one) + one.KroneckerProduct(tdt.Transpose())));
        }
        return l;
    }

    /// <summary>s = −Tr D/d², the mean rate; the candidate centre of a palindrome is −s.</summary>
    public BigRational MeanRate()
    {
        var dg = DissipatorCache;
        var tr = BigRational.Zero;
        for (int i = 0; i < Dim; i++) tr += dg[i, i].Re;
        return -tr / (Space.Levels * Space.Levels);
    }

    /// <summary>The coherence order m_i − m_j = j − i of the cell at vec index d·i + j.</summary>
    public static int OrderOf(int levels, int vecIndex) => vecIndex % levels - vecIndex / levels;

    /// <summary>The vec indices of the cells of coherence order p, in increasing order.</summary>
    public static int[] OrderCells(int levels, int p)
    {
        var cells = new List<int>();
        for (int i = 0; i < levels; i++)
        {
            int j = i + p;
            if (j >= 0 && j < levels) cells.Add(levels * i + j);
        }
        return cells.ToArray();
    }

    /// <summary>The gauged dissipator restricted to the cells of order p, as a real matrix (the
    /// dissipator conserves the order; an imaginary entry would be a construction error).</summary>
    public BigRational[,] DissipatorBlock(int p)
    {
        var dg = DissipatorCache;
        var cells = OrderCells(Space.Levels, p);
        var b = new BigRational[cells.Length, cells.Length];
        for (int r = 0; r < cells.Length; r++)
            for (int c = 0; c < cells.Length; c++)
            {
                var z = dg[cells[r], cells[c]];
                if (!z.Im.IsZero) throw new InvalidOperationException("the dissipator has an imaginary entry.");
                b[r, c] = z.Re;
            }
        return b;
    }

    /// <summary>−i[H, ·] as a superoperator on the row-major vec.</summary>
    public static GaussianRational[,] CommutatorPart(GaussianRational[,] h)
    {
        var one = M.Identity(h.GetLength(0));
        return M.Scale(M.Subtract(M.Sandwich(h, one), M.Sandwich(one, h)), -GaussianRational.I);
    }

    private GaussianRational[,] GaugedSecularHamiltonian()
    {
        var h = M.Add(M.Scale(QuadrupolarTensors.Root6T20(Space), Hamiltonian.Q / 3),
                      M.Scale(Space.GaugedIz(), Hamiltonian.Omega));
        if (!Hamiltonian.K1.IsZero) h = M.Add(h, M.Scale(LambPiece(1), Hamiltonian.K1));
        if (!Hamiltonian.K2.IsZero) h = M.Add(h, M.Scale(LambPiece(2), Hamiltonian.K2));
        if (!Hamiltonian.Nu.IsZero)
            h = M.Add(h, M.Scale(M.Add(Space.GaugedIplus(), Space.GaugedIminus()), Hamiltonian.Nu * Half));
        return h;
    }

    /// <summary>T₂q†T₂q − T₂,−q†T₂,−q in the gauge (diagonal, so the same in both bases).</summary>
    private GaussianRational[,] LambPiece(int q)
    {
        GaussianRational[,] TdT(int m) => M.Scale(M.Multiply(QuadrupolarTensors.Gauged(Space, -m), QuadrupolarTensors.Gauged(Space, m)),
                                                  QuadrupolarTensors.Parity(m));
        return M.Subtract(TdT(q), TdT(-q));
    }

    /// <summary>ρ ↦ TρT† − ½{T†T, ρ}, with the image of T† given (in the gauge it is not the adjoint).</summary>
    private static GaussianRational[,] Dissipate(GaussianRational[,] t, GaussianRational[,] tDagger)
    {
        var one = M.Identity(t.GetLength(0));
        var tdt = M.Multiply(tDagger, t);
        return M.Subtract(M.Sandwich(t, tDagger), M.Scale(M.Add(M.Sandwich(tdt, one), M.Sandwich(one, tdt)), Half));
    }

    private static bool TryRationalSqrt(BigRational x, out BigRational root)
    {
        root = BigRational.Zero;
        if (x.Sign < 0) return false;
        var n = IntegerSqrt(x.Numerator);
        var dd = IntegerSqrt(x.Denominator);
        if (n * n != x.Numerator || dd * dd != x.Denominator) return false;
        root = new BigRational(n, dd);
        return true;
    }

    private static BigInteger IntegerSqrt(BigInteger v)
    {
        if (v.Sign <= 0) return BigInteger.Zero;
        var x = (BigInteger)Math.Sqrt((double)v);
        while (x * x > v) x--;
        while ((x + 1) * (x + 1) <= v) x++;
        return x;
    }
}
