using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>The computing gate for <see cref="ShiftedGeneratorSectorwisePClaim"/>: on an actual
/// Liouvillian, per Π² sector, the phase-normalized P = √p_x·Π is an involution that anticommutes
/// with the shifted generator L' = L + σI.
///
/// <para><b>An exact route, so exact comparisons.</b> L is built in the 4^N Pauli-string basis
/// from the Pauli algebra itself: −i[H, σ_b] for one Hamiltonian term σ_t is (ph(σ_tσ_b) −
/// ph(σ_bσ_t))·σ_{t·b}, a unit phase difference in {0, ±2, ±2i} times the coupling, and two
/// distinct terms never map σ_b to the same string, so every off-diagonal entry is ONE product of
/// a coupling with a Gaussian-integer unit and carries no rounding for any real coupling. The
/// dephasing diagonal −2Σ_{l∈XY(b)} γ_l is a sum of rates, exact for the dyadic rates used here.
/// Π (<see cref="PiOperator.BuildFull"/>) is a signed permutation with phases in {1, i}, so every
/// conjugation below moves entries and multiplies them by units. Every assertion is ==.</para>
///
/// <para>The couplings are generic (J = 0.7, 1.1, Δ = 1.3, −0.4, non-uniform rates), so nothing
/// holds by a coincidence of the inputs. The rates are dyadic for a stated reason: the diagonal
/// identity σ − 2A = −(σ − 2(σ − A)) is formed from float sums of rates, so with non-dyadic rates
/// its residual would be a function of those sums' rounding alone; dyadic rates license == there,
/// they do not make the identity easier to satisfy. Each clause has a control that can fail it:
/// the involution clause, because the unnormalized Π squares to −I on the p_x = −1 sector (the
/// phase √−1 = i is what makes P an involution there); the anticommutation and sector clauses,
/// because a longitudinal Z field, which breaks F1 and [L, Π²] = 0 on this chain, breaks both. The builder is checked once against the repo's own vec-form route
/// (<see cref="LindbladianBuilder"/> + <see cref="PauliBasis.VecToPauliBasisTransform"/>) at
/// dyadic inputs, where that route is exact too, so the two are the same object and not two
/// conventions.</para></summary>
public class ShiftedGeneratorSectorwisePGateTests
{
    public sealed record BondSpec(int A, int B, double J, double Delta);

    private static readonly Complex MinusI = new(0.0, -1.0);

    /// <summary>σ_p·σ_q = phase·σ_r for single letters.</summary>
    private static (PauliLetter Letter, Complex Phase) Mul(PauliLetter p, PauliLetter q)
    {
        if (p == PauliLetter.I) return (q, Complex.One);
        if (q == PauliLetter.I) return (p, Complex.One);
        if (p == q) return (PauliLetter.I, Complex.One);
        var i = Complex.ImaginaryOne;
        return (p, q) switch
        {
            (PauliLetter.X, PauliLetter.Y) => (PauliLetter.Z, i),
            (PauliLetter.Y, PauliLetter.X) => (PauliLetter.Z, -i),
            (PauliLetter.Y, PauliLetter.Z) => (PauliLetter.X, i),
            (PauliLetter.Z, PauliLetter.Y) => (PauliLetter.X, -i),
            (PauliLetter.Z, PauliLetter.X) => (PauliLetter.Y, i),
            (PauliLetter.X, PauliLetter.Z) => (PauliLetter.Y, -i),
            _ => throw new InvalidOperationException(),
        };
    }

    private static List<(PauliLetter[] Letters, double Coefficient)> Terms(
        int n, IEnumerable<BondSpec> bonds, (int Site, PauliLetter Letter, double H)? field)
    {
        var terms = new List<(PauliLetter[], double)>();
        foreach (var bond in bonds)
            foreach (var (letter, c) in new[] { (PauliLetter.X, bond.J), (PauliLetter.Y, bond.J), (PauliLetter.Z, bond.J * bond.Delta) })
            {
                var letters = Enumerable.Repeat(PauliLetter.I, n).ToArray();
                letters[bond.A] = letter;
                letters[bond.B] = letter;
                terms.Add((letters, c));
            }
        if (field is { } f)
        {
            var letters = Enumerable.Repeat(PauliLetter.I, n).ToArray();
            letters[f.Site] = f.Letter;
            terms.Add((letters, f.H));
        }
        return terms;
    }

    /// <summary>L in the 4^N Pauli-string coefficient basis, straight from the Pauli algebra.</summary>
    private static Complex[,] BuildPauliL(int n, IEnumerable<BondSpec> bonds, double[] gamma,
        (int Site, PauliLetter Letter, double H)? field = null)
    {
        int d2 = 1 << (2 * n);
        var l = new Complex[d2, d2];
        var terms = Terms(n, bonds, field);
        for (int b = 0; b < d2; b++)
        {
            var lb = PauliIndex.FromFlat(b, n);
            foreach (var (lt, c) in terms)
            {
                Complex tb = Complex.One, bt = Complex.One;
                var target = new PauliLetter[n];
                for (int s = 0; s < n; s++)
                {
                    var (r, p1) = Mul(lt[s], lb[s]);
                    var (_, p2) = Mul(lb[s], lt[s]);
                    tb *= p1;
                    bt *= p2;
                    target[s] = r;
                }
                Complex comm = tb - bt;
                if (comm == Complex.Zero) continue;
                int row = (int)PauliIndex.ToFlat(target);
                Assert.True(l[row, b] == Complex.Zero, "two terms hit one entry: the exact route would not hold");
                l[row, b] = MinusI * c * comm;
            }
            double rate = 0.0;
            for (int s = 0; s < n; s++)
                if (lb[s].BitA() == 1) rate += gamma[s];
            l[b, b] += -2.0 * rate;
        }
        return l;
    }

    private static Complex[,] Pi(int n)
    {
        var sparse = PiOperator.BuildFull(n, PauliLetter.Z);
        int d2 = sparse.RowCount;
        var pi = new Complex[d2, d2];
        for (int i = 0; i < d2; i++)
            for (int j = 0; j < d2; j++)
                pi[i, j] = sparse[i, j];
        return pi;
    }

    private static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        int n = a.GetLength(0), m = b.GetLength(1), k = a.GetLength(1);
        var c = new Complex[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                Complex s = Complex.Zero;
                for (int t = 0; t < k; t++)
                    if (a[i, t] != Complex.Zero && b[t, j] != Complex.Zero) s += a[i, t] * b[t, j];
                c[i, j] = s;
            }
        return c;
    }

    private static Complex[,] Dagger(Complex[,] a)
    {
        int n = a.GetLength(0);
        var c = new Complex[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                c[j, i] = Complex.Conjugate(a[i, j]);
        return c;
    }

    private static Complex[,] Shift(Complex[,] l, double sigma)
    {
        int n = l.GetLength(0);
        var c = (Complex[,])l.Clone();
        for (int i = 0; i < n; i++) c[i, i] += sigma;
        return c;
    }

    private static Complex[,] Sub(Complex[,] a, int[] rows, int[] cols)
    {
        var c = new Complex[rows.Length, cols.Length];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 0; j < cols.Length; j++)
                c[i, j] = a[rows[i], cols[j]];
        return c;
    }

    private static bool AllZero(Complex[,] a)
    {
        foreach (var x in a) if (x != Complex.Zero) return false;
        return true;
    }

    private static bool IsIdentity(Complex[,] a, Complex scale)
    {
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (a[i, j] != (i == j ? scale : Complex.Zero)) return false;
        return true;
    }

    private static Complex[,] AntiCommutator(Complex[,] p, Complex[,] l)
    {
        var a = Mul(p, l);
        var b = Mul(l, p);
        int n = a.GetLength(0);
        var c = new Complex[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                c[i, j] = a[i, j] + b[i, j];
        return c;
    }

    private static int[] Sector(int n, int px)
    {
        int d2 = 1 << (2 * n);
        return Enumerable.Range(0, d2)
            .Where(k => PiOperator.SquaredEigenvalue(PauliIndex.FromFlat(k, n), PauliLetter.Z) == px)
            .ToArray();
    }

    public static IEnumerable<object[]> Chains() => new[]
    {
        new object[] { 2, new[] { new BondSpec(0, 1, 0.7, 1.3) }, new[] { 0.25, 0.75 } },
        new object[] { 3, new[] { new BondSpec(0, 1, 0.7, 1.3), new BondSpec(1, 2, 1.1, -0.4) }, new[] { 0.25, 0.75, 0.5 } },
    };

    [Theory]
    [MemberData(nameof(Chains))]
    public void PerPi2Sector_PhaseNormalizedPi_IsAnInvolutionAnticommutingWithTheShiftedGenerator(
        int n, BondSpec[] bonds, double[] gamma)
    {
        var l = BuildPauliL(n, bonds, gamma);
        var shifted = Shift(l, gamma.Sum());
        var pi = Pi(n);

        // The F1 premise the claim rests on: Π L' Π⁻¹ = −L', entry by entry.
        var conj = Mul(Mul(pi, shifted), Dagger(pi));
        for (int i = 0; i < conj.GetLength(0); i++)
            for (int j = 0; j < conj.GetLength(1); j++)
                Assert.True(conj[i, j] == -shifted[i, j], $"Π L' Π⁻¹ ≠ −L' at ({i},{j})");

        foreach (int px in new[] { +1, -1 })
        {
            int[] inSector = Sector(n, px);
            int[] outside = Enumerable.Range(0, shifted.GetLength(0)).Except(inSector).ToArray();
            Assert.Equal(shifted.GetLength(0) / 2, inSector.Length);

            // [L, Π²] = 0: L' does not leak out of the sector, so the sector is a resolved block.
            Assert.True(AllZero(Sub(shifted, outside, inSector)), $"L' leaks out of the p_x = {px} sector");

            var piSector = Sub(pi, inSector, inSector);
            // Unnormalized, Π squares to p_x·I on the sector: −I on the odd one.
            Assert.True(IsIdentity(Mul(piSector, piSector), px), $"Π² ≠ {px}·I on the p_x = {px} sector");

            Complex phase = ShiftedGeneratorSectorwisePClaim.PhaseForPi2Character(px);
            var p = new Complex[inSector.Length, inSector.Length];
            for (int i = 0; i < inSector.Length; i++)
                for (int j = 0; j < inSector.Length; j++)
                    p[i, j] = phase * piSector[i, j];
            var lSector = Sub(shifted, inSector, inSector);

            Assert.True(IsIdentity(Mul(p, p), Complex.One), $"P² ≠ I on the p_x = {px} sector");
            Assert.True(AllZero(AntiCommutator(p, lSector)), $"{{P, L'}} ≠ 0 on the p_x = {px} sector");
        }
    }

    [Fact]
    public void ALongitudinalZField_BreaksTheAnticommutation()
    {
        // Control: the same N = 2 chain with h·Z on site 0. F1 fails for this H, and so does the
        // sectorwise statement; a gate that could not see this would pass anything.
        var bonds = new[] { new BondSpec(0, 1, 0.7, 1.3) };
        double[] gamma = { 0.25, 0.75 };
        var shifted = Shift(BuildPauliL(2, bonds, gamma, field: (0, PauliLetter.Z, 0.3)), gamma.Sum());
        var pi = Pi(2);
        var anti = AntiCommutator(pi, shifted);
        Assert.False(AllZero(anti), "a Z field should break Π L' Π⁻¹ = −L'");
        int[] even = Sector(2, +1);
        int[] odd = Sector(2, -1);
        Assert.False(AllZero(Sub(shifted, odd, even)), "a Z field should make L' leak between Π² sectors");
        // The X field on the same site keeps it: Π maps the X term into itself with the sign F1 needs.
        var shiftedX = Shift(BuildPauliL(2, bonds, gamma, field: (0, PauliLetter.X, 0.3)), gamma.Sum());
        Assert.True(AllZero(AntiCommutator(pi, shiftedX)));
    }

    [Fact]
    public void ThePauliAlgebraBuilder_IsTheRepoLiouvillian()
    {
        // Dyadic inputs keep the vec-form route exact as well, so the comparison is ==. The vec
        // transform stores σᵀ, so its Pauli-basis L carries the sign (−1)^(n_Y(a)+n_Y(b))
        // relative to the natural one (see PauliBasis.VecToPauliBasisTransform).
        const int n = 3;
        var bonds = new[] { new BondSpec(0, 1, 0.75, 1.25), new BondSpec(1, 2, 1.5, -0.5) };
        double[] gamma = { 0.25, 1.0, 0.5625 };
        var direct = BuildPauliL(n, bonds, gamma);

        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        foreach (var (letters, c) in Terms(n, bonds, null))
            h += (Complex)c * PauliString.Build(letters);
        var cOps = Enumerable.Range(0, n)
            .Select(s => (Complex)Math.Sqrt(gamma[s]) * PauliString.SiteOp(n, s, PauliLetter.Z))
            .ToList();
        var lVec = LindbladianBuilder.Build(h, cOps);
        var t = PauliBasis.VecToPauliBasisTransform(n);
        var viaVec = t.ConjugateTranspose() * lVec * t / Math.Pow(2, n);

        int d2 = 1 << (2 * n);
        for (int a = 0; a < d2; a++)
        {
            int ya = PauliIndex.FromFlat(a, n).Count(x => x == PauliLetter.Y);
            for (int b = 0; b < d2; b++)
            {
                int yb = PauliIndex.FromFlat(b, n).Count(x => x == PauliLetter.Y);
                Complex expected = ((ya + yb) & 1) == 0 ? direct[a, b] : -direct[a, b];
                Assert.True(viaVec[a, b] == expected, $"entry ({a},{b}): vec route {viaVec[a, b]}, Pauli route {expected}");
            }
        }
    }
}
