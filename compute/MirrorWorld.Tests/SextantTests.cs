using System.Numerics;
using MirrorWorld;
using static MirrorWorldTests.LetterCubeTests;

namespace MirrorWorldTests;

// From-below guard for the sextant (built 2026-10-14 beside docs/THE_ONE_SQUARE.md §9): a setting read against the
// letter cube, every move asked the relation form of §9's two questions on the full generator, exactly. Judged by the
// verdicts §9 and its gates already hold (Z dephasing, loss, a thermal bath, F1's Pi_Z), a jump with a trace, the
// sign of L_H against a dense commutator, and every verdict against the generator built densely in the matrix-unit
// basis at N = 2.
public class SextantTests
{
    static readonly World W = new();

    [Fact]
    public void The_Sextant_Owns_The_Verdict_And_The_Search_And_Inherits_The_Cube()
    {
        var sextant = new Sextant(new LetterCube(W));
        Assert.Equal(new[] { "verdict", "search" }, sextant.Own);
        Assert.Equal(new[] { "coordinates", "height", "step", "arrow", "moves", "x", "y", "z" }, sextant.Inherited);
    }

    // ---- the two questions (§9, gates M1, M2, M5, M8, M9) ----

    static Sextant.Setting Heisenberg(int n, IEnumerable<(int, LetterCube.Gauss[,], long)> jumps)
        => Sextant.Setting.Chain(n, new[] { ("XX", 1L), ("YY", 1L), ("ZZ", 1L) }, jumps);

    static IEnumerable<(int, LetterCube.Gauss[,], long)> Every(int n, LetterCube.Gauss[,] jump, long rate)
        => Sextant.Setting.OnEverySite(n, jump, rate);

    static Dictionary<string, Sextant.Reading> Table(Sextant.Setting s)
        => LetterCube.TableMoves(s.N).Select(m => Sextant.Ask(s, m)).ToDictionary(r => r.Move.Split(' ')[0] == "R," ? "R" : r.Move, r => r);

    // §9 under Z dephasing on the Heisenberg chain (gates M1, M2, M5): R and the shift by Y^N reflect with c = 2 sigma,
    // the shift by Z^N, the three turns, the quarter-turn about Z and the swap X <-> Y are copies, the transpose flips
    // L_H and keeps the dissipator, the quarter-turns about X and Y and the swaps that move Z are neither, and S_XY (the
    // turn by X^N then the shift by Z^N) is a copy. The transpose followed by R is F1's Pi_Z: S L S^-1 = -L - 2 sigma.
    [Fact]
    public void Under_Z_Dephasing_The_Table_Reads_As_Section_Nine()
    {
        const int n = 3;
        var t = Table(Heisenberg(n, Every(n, LetterCube.LetterMatrix('Z'), 1)));
        Assert.Equal(Sextant.Verdict.Reflection, t["R"].Verdict);
        Assert.Equal($"{2 * n}", t["R"].C);
        Assert.Equal(Sextant.Verdict.Reflection, t["the shift by Y^N (rho.Y^N)"].Verdict);
        foreach (var copy in new[] { "the shift by Z^N (rho.Z^N)", "the turn by X^N", "the turn by Y^N", "the turn by Z^N",
                                     "the quarter-turn about Z", "the letter swap X <-> Y", "S_XY: rho -> X^N rho Y^N" })
            Assert.Equal(Sextant.Verdict.Copy, t[copy].Verdict);
        Assert.Equal(Sextant.OnHamiltonian.Flips, t["the transpose"].OnH);
        Assert.Equal(Sextant.OnDissipator.Keeps, t["the transpose"].OnD);
        foreach (var broken in new[] { "the quarter-turn about X", "the quarter-turn about Y", "the letter swap X <-> Z", "the letter swap Y <-> Z" })
            Assert.Equal(Sextant.Verdict.Neither, t[broken].Verdict);
        var moves1 = LetterCube.TableMoves(1);
        var piZ = moves1.Single(m => m.Name == "the transpose").Sites[0].Then(moves1.Single(m => m.Name.StartsWith("R,")).Sites[0], "Pi_Z");
        var pi = Sextant.Ask(Heisenberg(n, Every(n, LetterCube.LetterMatrix('Z'), 1)), LetterCube.Move.Uniform(piZ, n));
        Assert.Equal(Sextant.Verdict.PiReflection, pi.Verdict);
        Assert.Equal((Sextant.OnHamiltonian.Flips, Sextant.OnDissipator.Reflects), (pi.OnH, pi.OnD));
        Assert.Equal($"{2 * n}", pi.C);
    }

    // §9 under amplitude damping (gates M8, M9): a copy keeps the arrow (the turn by Z^N, the quarter-turn about Z),
    // the transpose keeps the damping and flips L_H, the turns by X^N and Y^N and every shift break, and S_XY turns
    // the arrow round: a reflection with c = sum of the rates, at loss alone and with gain beside it.
    [Fact]
    public void Under_Damping_A_Copy_Keeps_The_Arrow_And_S_XY_Turns_It_Round()
    {
        const int n = 3;
        var loss = Table(Heisenberg(n, Every(n, LetterCube.Lowering, 2)));
        Assert.Equal(Sextant.Verdict.Copy, loss["the turn by Z^N"].Verdict);
        Assert.Equal(Sextant.Verdict.Copy, loss["the quarter-turn about Z"].Verdict);
        Assert.Equal(Sextant.OnDissipator.Keeps, loss["the transpose"].OnD);
        Assert.Equal(Sextant.OnHamiltonian.Flips, loss["the transpose"].OnH);
        foreach (var broken in new[] { "the turn by X^N", "the turn by Y^N", "R", "the shift by Y^N (rho.Y^N)", "the shift by Z^N (rho.Z^N)" })
            Assert.Equal(Sextant.OnDissipator.Breaks, loss[broken].OnD);
        Assert.Equal(Sextant.Verdict.Reflection, loss["S_XY: rho -> X^N rho Y^N"].Verdict);
        Assert.Equal($"{n * 2}", loss["S_XY: rho -> X^N rho Y^N"].C);
        var jumps = Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1));
        var thermal = Table(Heisenberg(n, jumps));
        Assert.Equal(Sextant.Verdict.Reflection, thermal["S_XY: rho -> X^N rho Y^N"].Verdict);
        Assert.Equal($"{n * 4}", thermal["S_XY: rho -> X^N rho Y^N"].C);
        // at equal loss and gain the arrows cancel: the shift by Z^N reflects and the turns by X^N and Y^N are copies
        var hot = Table(Heisenberg(n, Every(n, LetterCube.Lowering, 2).Concat(Every(n, LetterCube.Raising, 2))));
        Assert.Equal(Sextant.Verdict.Reflection, hot["the shift by Z^N (rho.Z^N)"].Verdict);
        Assert.Equal(Sextant.Verdict.Copy, hot["the turn by X^N"].Verdict);
        Assert.Equal(Sextant.Verdict.Neither, hot["R"].Verdict);
    }

    // The search finds S_XY by itself: searching the whole closure on the thermal Heisenberg chain gives 4 copies and
    // 4 reflections, the reflections exactly S_XY composed with the copies, and nothing of Pi's shape; a longitudinal Z
    // field (odd under X^N) leaves no reflection at all.
    [Fact]
    public void The_Search_Finds_The_Thermal_Reflection_By_Itself()
    {
        const int n = 3;
        var jumps = Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)).ToList();
        var found = Sextant.Search(Heisenberg(n, jumps));
        var sxy = LetterCube.TableMoves(1).Single(m => m.Name.StartsWith("S_XY")).Sites[0].Key;
        Assert.Contains(found, f => f.Reading.Verdict == Sextant.Verdict.Reflection && f.Map.Key == sxy);
        Assert.All(found.Where(f => f.Reading.Verdict == Sextant.Verdict.Reflection), f => Assert.Equal($"{n * 4}", f.Reading.C));
        var copies = found.Where(f => f.Reading.Verdict == Sextant.Verdict.Copy).Select(f => f.Map).ToList();
        var reflections = found.Where(f => f.Reading.Verdict == Sextant.Verdict.Reflection).Select(f => f.Map.Key).ToHashSet();
        Assert.Equal(4, copies.Count);
        Assert.Equal(4, reflections.Count);
        Assert.DoesNotContain(found, f => f.Reading.Verdict == Sextant.Verdict.PiReflection);
        var sxyMap = LetterCube.TableMoves(1).Single(m => m.Name.StartsWith("S_XY")).Sites[0];
        Assert.True(reflections.SetEquals(copies.Select(c => c.Then(sxyMap).Key)));
        var loss = Sextant.Search(Heisenberg(n, Every(n, LetterCube.Lowering, 2)));
        Assert.Equal(4, loss.Count(f => f.Reading.Verdict == Sextant.Verdict.Copy));
        Assert.Equal(4, loss.Count(f => f.Reading.Verdict == Sextant.Verdict.Reflection));
        var withField = Sextant.Setting.Chain(n, new[] { ("XX", 1L), ("YY", 1L), ("ZZ", 1L), ("Z", 1L) }, jumps);
        Assert.DoesNotContain(Sextant.Search(withField), f => f.Reading.Verdict == Sextant.Verdict.Reflection);
    }

    // The verdicts judged densely on two sites, in the matrix-unit basis, compared exactly (every entry a small dyadic):
    // S L S^-1 against L, against -L^dagger - c and against -L - c, with L and its adjoint built from H and the jumps here (the adjoint Lindbladian
    // i[H, X] + sum r (J^dagger X J - {J^dagger J, X}/2)), for the 14 named moves and all 192 closure maps, on Heisenberg with loss and gain and on Heisenberg with Z
    // dephasing, where F1's Pi_Z is among the Pi-reflections.
    [Fact]
    public void The_Verdicts_Agree_With_The_Dense_Generator()
    {
        const int n = 2;
        int piSeen = 0;
        foreach (var jumps in new[] { Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)).ToList(),
                                      Every(n, LetterCube.LetterMatrix('Z'), 1).ToList() })
        {
        var setting = Heisenberg(n, jumps);
        var H = new Complex[4, 4];
        foreach (var (term, c) in setting.H) Acc(H, DenseString(term.ToString(n)), c);
        var Js = jumps.Select(j => (Site(n, j.Item1, ToComplex(j.Item2)), (double)j.Item3)).ToList();
        Complex[,] L(Complex[,] rho, bool adjoint)
        {
            var outm = new Complex[4, 4];
            var hr = Mul(H, rho); var rh = Mul(rho, H);
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
                    outm[i, j] = (adjoint ? Complex.ImaginaryOne : -Complex.ImaginaryOne) * (hr[i, j] - rh[i, j]);
            foreach (var (J, r) in Js)
            {
                var Jd = Dag4(J); var JdJ = Mul(Jd, J);
                var sand = adjoint ? Mul(Jd, Mul(rho, J)) : Mul(J, Mul(rho, Jd));
                var l = Mul(JdJ, rho); var rr = Mul(rho, JdJ);
                for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) outm[i, j] += r * (sand[i, j] - 0.5 * (l[i, j] + rr[i, j]));
            }
            return outm;
        }
        foreach (var move in LetterCube.TableMoves(n).Concat(LetterCube.Closure().Select(m => LetterCube.Move.Uniform(m, n, m.Describe()))))
        {
            // the move on a matrix: through the Pauli expansion, rho = sum_s c_s s / 4, each s sent to phase * image
            Complex[,] Apply(Complex[,] rho, bool inverse)
            {
                var outm = new Complex[4, 4];
                foreach (var s in LetterCube.Strings(n))
                {
                    var S = DenseString(s.ToString(n));
                    var coeff = Trace(Mul(S, rho)) / 4;
                    if (coeff == Complex.Zero) continue;
                    var (img, ph) = move.Apply(s);
                    var phase = new Complex(ph.Re, ph.Im);
                    if (inverse)
                    {
                        // the preimage of s: the string t with image s
                        var t = LetterCube.Strings(n).Single(u => move.Apply(u).Image == s);
                        var (_, pt) = move.Apply(t);
                        Acc(outm, DenseString(t.ToString(n)), coeff / new Complex(pt.Re, pt.Im));
                    }
                    else Acc(outm, DenseString(img.ToString(n)), coeff * phase);
                }
                return outm;
            }
            var reading = Sextant.Ask(setting, move);
            bool copy = true, refl = true, piRefl = true;
            double c = (double)reading.MinusTwoTraceTimesFour / (4 * reading.Dim);
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    var e = new Complex[4, 4]; e[i, j] = 1;
                    var sls = Apply(L(Apply(e, inverse: true), false), inverse: false);
                    var l = L(e, false);
                    var la = L(e, true);
                    for (int a = 0; a < 4; a++)
                        for (int b = 0; b < 4; b++)
                        {
                            if (sls[a, b] - l[a, b] != Complex.Zero) copy = false;
                            if (sls[a, b] + la[a, b] + c * e[a, b] != Complex.Zero) refl = false;
                            if (sls[a, b] + l[a, b] + c * e[a, b] != Complex.Zero) piRefl = false;
                        }
                }
            Assert.Equal(copy, reading.Verdict == Sextant.Verdict.Copy);
            Assert.Equal(refl && !copy, reading.Verdict == Sextant.Verdict.Reflection);
            Assert.Equal(piRefl && !copy && !refl, reading.Verdict == Sextant.Verdict.PiReflection);
            if (piRefl && !copy && !refl) piSeen++;
        }
        }
        Assert.True(piSeen > 0);                                         // the Pi branch is checked in the positive too
    }

    // The sign of L_H, which no verdict can see (every verdict is even under L_H -> -L_H): the stored Hamiltonian part
    // applied to every string on two sites equals 4(-i[H, sigma]) built densely here, compared exactly.
    [Fact]
    public void The_Hamiltonian_Part_Has_The_Right_Sign()
    {
        const int n = 2;
        var setting = Sextant.Setting.Chain(n, new[] { ("XX", 1L), ("YY", 2L), ("ZZ", 3L), ("X", 1L) }, Array.Empty<(int, LetterCube.Gauss[,], long)>());
        var H = new Complex[4, 4];
        foreach (var (term, c) in setting.H) Acc(H, DenseString(term.ToString(n)), c);
        foreach (var s in LetterCube.Strings(n))
        {
            var S = DenseString(s.ToString(n));
            var hs = Mul(H, S); var sh = Mul(S, H);
            var expected = new Complex[4, 4];
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) expected[i, j] = 4 * -Complex.ImaginaryOne * (hs[i, j] - sh[i, j]);
            var got = new Complex[4, 4];
            foreach (var ((row, col), v) in setting.HamiltonianTimesFour)
                if (col == s) Acc(got, DenseString(row.ToString(n)), new Complex(v.Re, v.Im));
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) Assert.True(expected[i, j] == got[i, j]);
        }
    }

    // A jump with a trace hides part of its dissipator in a commutator: 1 + iZ at rate 1 with a field Z on every site is
    // pure Z dephasing, so every map of the closure gets the same verdict as on pure Z dephasing, while the readings
    // part by part do not (the turn by X^N flips L_H, breaks the dissipator and is a copy of the whole).
    [Fact]
    public void A_Jump_With_A_Trace_Is_Read_On_The_Whole_Generator()
    {
        const int n = 2;
        var onePlusIZ = new[,] { { new LetterCube.Gauss(1, 1), LetterCube.Gauss.Zero }, { LetterCube.Gauss.Zero, new LetterCube.Gauss(1, -1) } };
        var hidden = Sextant.Setting.Chain(n, new[] { ("Z", 1L) }, Every(n, onePlusIZ, 1));
        var plain = Sextant.Setting.Chain(n, Array.Empty<(string, long)>(), Every(n, LetterCube.LetterMatrix('Z'), 1));
        foreach (var m in LetterCube.Closure())
        {
            var move = LetterCube.Move.Uniform(m, n);
            Assert.Equal(Sextant.Ask(plain, move).Verdict, Sextant.Ask(hidden, move).Verdict);
        }
        var turnX = Sextant.Ask(hidden, LetterCube.TableMoves(n).Single(m => m.Name == "the turn by X^N"));
        Assert.Equal(Sextant.Verdict.Copy, turnX.Verdict);
        Assert.Equal((Sextant.OnHamiltonian.Flips, Sextant.OnDissipator.Breaks), (turnX.OnH, turnX.OnD));
    }

    static void Acc(Complex[,] m, Complex[,] add, Complex k)
    {
        for (int i = 0; i < m.GetLength(0); i++) for (int j = 0; j < m.GetLength(1); j++) m[i, j] += k * add[i, j];
    }

    static Complex Trace(Complex[,] m) { Complex t = 0; for (int i = 0; i < m.GetLength(0); i++) t += m[i, i]; return t; }

    static Complex[,] Dag4(Complex[,] a)
    {
        int d = a.GetLength(0); var r = new Complex[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) r[i, j] = Complex.Conjugate(a[j, i]);
        return r;
    }

    static Complex[,] Site(int n, int l, Complex[,] one)
    {
        Complex[,] m = l == 0 ? one : Letter('I');
        for (int k = 1; k < n; k++) m = Kron(m, k == l ? one : Letter('I'));
        return m;
    }

}
