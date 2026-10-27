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
        Assert.Equal(new[] { "verdict", "search", "decision" }, sextant.Own);
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

    // One bit watched along Z with a field along X, Y or Z (§9's corner-shift row; the one bit with a field of
    // reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md). Read on the cube by hand first: the watching puts I
    // and Z on the quiet edge and X and Y on the far edge, and a field along P leaves I and P at rest and moves the
    // other two into each other. The X field moves Y and Z, one on each edge, and R (I <-> X, Y <-> Z) flips the height
    // and carries each pair onto itself: a reflection. The Y field moves X and Z, and the shift by Y (I <-> Y, X <-> Z)
    // does the same. The Z field moves X and Y, both on the far edge, and leaves I and Z at rest, both on the quiet
    // edge: a map that flips the height would carry the moving pair onto the resting one, so nothing reflects, and the
    // shift by Z is a copy. Pi_Z (the transpose, then R) needs a real field: it reflects under X, not under Y.
    [Fact]
    public void One_Watched_Bit_With_A_Field_Reads_As_The_Cube_By_Hand()
    {
        var moves1 = LetterCube.TableMoves(1);
        var piZ = LetterCube.Move.Uniform(moves1.Single(m => m.Name == "the transpose").Sites[0].Then(moves1.Single(m => m.Name.StartsWith("R,")).Sites[0], "Pi_Z"), 1);
        Sextant.Setting Bit(char? field) => new(1, field is { } f ? new[] { (PauliString.Parse(f.ToString()), 1L) } : Array.Empty<(PauliString, long)>(),
                                               Every(1, LetterCube.LetterMatrix('Z'), 1));
        const string shiftY = "the shift by Y^N (rho.Y^N)", shiftZ = "the shift by Z^N (rho.Z^N)";

        var none = Table(Bit(null));
        Assert.Equal(Sextant.Verdict.Reflection, none["R"].Verdict);
        Assert.Equal(Sextant.Verdict.Reflection, none[shiftY].Verdict);

        var x = Table(Bit('X'));
        Assert.Equal((Sextant.OnHamiltonian.Keeps, Sextant.OnDissipator.Reflects, Sextant.Verdict.Reflection), (x["R"].OnH, x["R"].OnD, x["R"].Verdict));
        Assert.Equal("2", x["R"].C);
        Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.Verdict.Neither), (x[shiftY].OnH, x[shiftY].Verdict));
        Assert.Equal(Sextant.Verdict.PiReflection, Sextant.Ask(Bit('X'), piZ).Verdict);

        var y = Table(Bit('Y'));
        Assert.Equal((Sextant.OnHamiltonian.Keeps, Sextant.Verdict.Reflection), (y[shiftY].OnH, y[shiftY].Verdict));
        Assert.Equal("2", y[shiftY].C);
        Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.Verdict.Neither), (y["R"].OnH, y["R"].Verdict));
        Assert.Equal(Sextant.Verdict.Neither, Sextant.Ask(Bit('Y'), piZ).Verdict);

        var z = Table(Bit('Z'));
        Assert.Equal(Sextant.Verdict.Copy, z[shiftZ].Verdict);
        Assert.Equal(Sextant.Verdict.Neither, z["R"].Verdict);
        Assert.Equal(Sextant.Verdict.Neither, z[shiftY].Verdict);

        Assert.Equal(Sextant.Palindrome.Witnessed, Sextant.Decide(Bit('X')).Kind);
        Assert.Equal(Sextant.Palindrome.Witnessed, Sextant.Decide(Bit('Y')).Kind);
        Assert.Equal(Sextant.Palindrome.Obstructed, Sextant.Decide(Bit('Z')).Kind);
        Assert.DoesNotContain(Sextant.Search(Bit('Z')), f => f.Reading.Verdict is Sextant.Verdict.Reflection or Sextant.Verdict.PiReflection);
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

    // F158 on a tilted axis: -3XX + 3YX under dephasing along X + Z at (2, 1). H uses one axis per site, (-1, 1, 0) on
    // site 0 and X on site 1, and a letter perpendicular to it and to the watched axis (1, 0, 1) anticommutes with
    // both, so the colour (X + Y - Z)⊗Y commutes with H and anticommutes with every jump. By a second route, dense
    // matrices built letter by letter that never see the masks: G·H = H·G and G·J = -J·G exactly, and the left
    // multiplication by G carries the generator (J² = 2, so σ = 6) to -L^dagger - 12 exactly, while a candidate the
    // search rejects, (X - Z)⊗Y, perpendicular to the watching but not to H's axis on site 0, does not. No theodolite
    // reading sees it (its factors are monomial, X + Y - Z is not). A Z field on site 0 asks of that site's colour to
    // commute with Z as well, a YY bond to be perpendicular to Y too; neither leaves a colour, and the odd traces then
    // say no, exactly. A Heisenberg bond asks one colour of a whole component: the normal of two sites' axes, or a
    // field's axis carried to a site without a field.
    [Fact]
    public void A_Colour_On_The_Circle_Witnesses_A_Tilted_Axis()
    {
        const int n = 2;
        var tilt = new[,] { { LetterCube.Gauss.One, LetterCube.Gauss.One }, { LetterCube.Gauss.One, -LetterCube.Gauss.One } };
        var jumps = new[] { (0, tilt, 2L), (1, tilt, 1L) };
        var s = Sextant.Setting.Chain(n, new[] { ("XX", -3L), ("YX", 3L) }, jumps);
        Assert.Equal(new[] { (1L, 1L, -1L), (0L, 1L, 0L) }, Sextant.CircleColour(s));
        var d = Sextant.Decide(s);
        Assert.Equal(Sextant.Palindrome.Witnessed, d.Kind);
        Assert.Contains("G = (X + Y - Z)⊗Y", d.Detail);
        var t = new Theodolite.Setting(s.N, s.H, s.Jumps);
        foreach (var rel in new[] { Theodolite.Relation.Reflection, Theodolite.Relation.Pi })
        {
            Assert.Null(Theodolite.Read(t, rel));
            Assert.Null(Theodolite.ReadRelabelled(t, rel));
            Assert.Null(Theodolite.ReadRotated(t, rel));
            Assert.Null(Theodolite.ReadSectors(t, rel));
        }
        Assert.Contains("c = 12", d.Detail);

        int dd = 1 << n, dim = dd * dd;
        var H = new Complex[dd, dd];
        foreach (var (term, c) in s.H) Acc(H, DenseString(term.ToString(n)), c);
        var g0 = new Complex[2, 2]; Acc(g0, Letter('X'), 1); Acc(g0, Letter('Y'), 1); Acc(g0, Letter('Z'), -1);
        var G = Kron(g0, Letter('Y'));
        var j = new Complex[2, 2]; Acc(j, Letter('X'), 1); Acc(j, Letter('Z'), 1);
        Assert.True(Same(Mul(G, H), Mul(H, G)));
        foreach (var J in new[] { Site(n, 0, j), Site(n, 1, j) })
            Assert.True(Same(Mul(G, J), Neg(Mul(J, G))));

        // row-major vec: A·rho·B -> A ⊗ B^T
        var id = Site(n, 0, Letter('I'));
        var L = new Complex[dim, dim];
        Acc(L, Kron(H, id), -Complex.ImaginaryOne);
        Acc(L, Kron(id, Transpose(H)), Complex.ImaginaryOne);
        foreach (var (J, rate) in new[] { (Site(n, 0, j), 2), (Site(n, 1, j), 1) })
        {
            Acc(L, Kron(J, Conj(J)), rate);
            Acc(L, Kron(id, id), -2 * rate);                                // J² = 2: -½{J², ·} = -2·id
        }
        var S = Kron(G, id);
        var minus = Dag4(L);
        Acc(minus, Kron(id, id), 12);
        Assert.True(Same(Mul(S, L), Neg(Mul(minus, S))));
        var r0 = new Complex[2, 2]; Acc(r0, Letter('X'), 1); Acc(r0, Letter('Z'), -1);
        var rejected = Kron(Kron(r0, Letter('Y')), id);
        Assert.False(Same(Mul(rejected, L), Neg(Mul(minus, rejected))));            // the check can fail

        // a Heisenberg chain asks one colour of all its sites: under axes (2, 2, -1), (2, -1, -2), (2, 2, -1) it is
        // their normal (5, -2, 6), read from the axes of two sites; with one axis (1, 1, 0) everywhere and fields along
        // (1, -1, -2) on sites 0 and 1 only, it is that field's axis, carried to site 2
        LetterCube.Gauss[,] Axis(long x, long y, long z) => new[,] { { new LetterCube.Gauss(z, 0), new LetterCube.Gauss(x, -y) }, { new LetterCube.Gauss(x, y), new LetterCube.Gauss(-z, 0) } };
        var heis = new[] { ("XX", 1L), ("YY", 1L), ("ZZ", 1L) };
        var twoAxes = Sextant.Setting.Chain(3, heis, new[] { (0, Axis(2, 2, -1), 2L), (1, Axis(2, -1, -2), 2L), (2, Axis(2, 2, -1), 3L) });
        Assert.Equal(Enumerable.Repeat((5L, -2L, 6L), 3), Sextant.CircleColour(twoAxes));
        var ring = LetterCube.ChainTerms(new[] { "XX", "YY", "ZZ" }, 3).Select(t => (t, 1L))
            .Concat(new[] { "XIX", "YIY", "ZIZ" }.Select(t => (PauliString.Parse(t), 1L)))
            .Concat(new[] { ("XII", 1L), ("YII", -1L), ("ZII", -2L), ("IXI", 2L), ("IYI", -2L), ("IZI", -4L) }.Select(f => (PauliString.Parse(f.Item1), f.Item2)));
        var carried = new Sextant.Setting(3, ring, Enumerable.Range(0, 3).Select(l => (l, Axis(1, 1, 0), (long)(l + 1))));
        Assert.Equal(Enumerable.Repeat((1L, -1L, -2L), 3), Sextant.CircleColour(carried));

        var field = new Sextant.Setting(n, s.H.Append((PauliString.Parse("ZI"), 1L)), jumps);
        var yy = Sextant.Setting.Chain(n, new[] { ("XX", -3L), ("YX", 3L), ("YY", 1L) }, jumps);
        foreach (var broken in new[] { field, yy })
        {
            Assert.Null(Sextant.CircleColour(broken));
            Assert.Equal(Sextant.Palindrome.Obstructed, Sextant.Decide(broken).Kind);
        }
    }

    // F138 off the letters (the open arc f138_clause_two_sweep, its tilted half; TiltedAxes holds the draws and the
    // rule). On F138's bond classes under dephasing along integer axes off the letters the cube's rule reads the
    // palindrome: a bond with coupling matrix W carries a colour a to W·a, so one colour per component, perpendicular
    // to every watched axis and along every field, for XXZ (Δ = 2) and XX + YY along Z or perpendicular to it; for ZZ
    // one colour per site, all along Z or all perpendicular to it.
    // Draws with a fixed seed on the three-site chain and ring and the four-site star are read against the sextant's
    // decision: every one is decided exactly, the rule agrees with every one, both answers occur, and F138's clauses
    // read with directions disagree with some. Two rows pin the two ways the clauses miss: three coplanar axes on a
    // Heisenberg chain pair (their normal (4, -3, -2) is the colour), which clause 1 forbids; two axes whose normal
    // (1, 0, 2) is neither along Z nor perpendicular to it do not pair on an XX + YY ring, which the clauses allow.
    [Fact]
    public void Off_The_Letters_The_Cube_Rule_Reads_F138s_Classes()
    {
        bool Read(TiltedAxes.Row row, string label)
        {
            var kind = Sextant.Decide(row.Setting()).Kind;
            Assert.True(kind is Sextant.Palindrome.Witnessed or Sextant.Palindrome.Obstructed, $"{label}: {kind}");
            Assert.True(row.Rule() == (kind == Sextant.Palindrome.Witnessed), $"{label}: rule {row.Rule()}, decision {kind}, {row}");
            return kind == Sextant.Palindrome.Witnessed;
        }
        var rnd = new Random(20261027);
        var three = TiltedAxes.Graphs(3);
        var star = TiltedAxes.Graphs(4).Single(g => g.Name == "star").Edges;
        int yes = 0, no = 0, clausesMiss = 0;
        for (int trial = 0; trial < 86; trial++)
        {
            var row = trial < 80 ? TiltedAxes.Draw(rnd, 3, three[trial % 2].Edges) : TiltedAxes.Draw(rnd, 4, star);   // the star costs more
            bool pairs = Read(row, $"trial {trial}");
            if (pairs) yes++; else no++;
            if (row.Clauses() != pairs) clausesMiss++;
        }
        Assert.True(yes > 0 && no > 0, $"yes {yes}, no {no}");
        Assert.True(clausesMiss > 0);

        var none = new (long X, long Y, long Z)[3];
        var coplanar = new TiltedAxes.Row(3, three[0].Edges, TiltedAxes.Classes[0], new (long, long, long)[] { (2, 2, 1), (1, 2, -1), (1, 0, 2) }, new[] { 2L, 1, 2 }, none);
        Assert.True(Read(coplanar, "coplanar"));
        Assert.False(coplanar.Clauses());
        Assert.Equal(Enumerable.Repeat((4L, -3L, -2L), 3), Sextant.CircleColour(coplanar.Setting()));
        var tiltedNormal = new TiltedAxes.Row(3, three[1].Edges, TiltedAxes.Classes[2], new (long, long, long)[] { (2, 0, -1), (2, 2, -1), (2, 2, -1) }, new[] { 3L, 2, 3 }, none);
        Assert.False(Read(tiltedNormal, "tilted normal"));
        Assert.True(tiltedNormal.Clauses());
        // ZZ takes one colour per site: axes perpendicular to Z with Z fields pair (every colour Z); a field on the
        // (1, 0, 1)-watched site 0 pairs along Y (perpendicular to Z and to the axis) and breaks along X + Y
        var zz = TiltedAxes.Classes[3];
        Assert.True(Read(new TiltedAxes.Row(3, three[0].Edges, zz, new (long, long, long)[] { (1, 1, 0), (1, -1, 0), (1, 2, 0) },
            new[] { 1L, 2, 3 }, new (long, long, long)[] { (0, 0, 1), (0, 0, -2), (0, 0, 0) }), "zz along"));
        var watched = new (long, long, long)[] { (1, 0, 1), (1, 0, 1), (1, 0, 1) };
        Assert.True(Read(new TiltedAxes.Row(3, three[0].Edges, zz, watched, new[] { 1L, 2, 1 },
            new (long, long, long)[] { (0, 1, 0), (0, 0, 0), (0, 0, 0) }), "zz across"));
        Assert.False(Read(new TiltedAxes.Row(3, three[0].Edges, zz, watched, new[] { 1L, 2, 1 },
            new (long, long, long)[] { (1, 1, 0), (0, 0, 0), (0, 0, 0) }), "zz across, field off"));
    }

    static bool Same(Complex[,] a, Complex[,] b)
    {
        for (int i = 0; i < a.GetLength(0); i++) for (int k = 0; k < a.GetLength(1); k++) if (a[i, k] != b[i, k]) return false;
        return true;
    }

    static Complex[,] Neg(Complex[,] a) { var r = (Complex[,])a.Clone(); for (int i = 0; i < r.GetLength(0); i++) for (int k = 0; k < r.GetLength(1); k++) r[i, k] = -r[i, k]; return r; }

    static Complex[,] Conj(Complex[,] a) { var r = (Complex[,])a.Clone(); for (int i = 0; i < r.GetLength(0); i++) for (int k = 0; k < r.GetLength(1); k++) r[i, k] = Complex.Conjugate(r[i, k]); return r; }

    static Complex[,] Transpose(Complex[,] a) { int d = a.GetLength(0); var r = new Complex[d, d]; for (int i = 0; i < d; i++) for (int k = 0; k < d; k++) r[i, k] = a[k, i]; return r; }

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

    // ---- the no (the odd centred traces) and the decision ----

    // The sextant's three exact answers on settings §9 and F137 already hold: the thermal Heisenberg chain is WITNESSED
    // (S_XY's class); the warm DNA model (transverse X fields, ZZ, loss, gain and Z) is OBSTRUCTED at the third power,
    // so no element of any kind reflects it; DM bonds under loss, which F137 pairs by another route, are WITNESSED by the
    // theodolite, a site-local sandwich outside the closure. The cooled chain XX + 2XY + ZZ at rates (1, 2, 1) is
    // WITNESSED by the theodolite after the chain's reversal (its carrier, the reversal times a turn about Z); at
    // (1, 2, 3) the seventh trace obstructs, as experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md has it. YZ with Z fields
    // under X dephasing on one site and gain on the other, at (1, 3), stays SYMMETRIC-MOD-P (no instrument reads an
    // element there), and read only to the third power, OPEN. With no jumps the identity reflects (WITNESSED).
    [Fact]
    public void The_Sextant_Answers_Yes_No_And_Symmetric_Unseen()
    {
        const int n = 3;
        var thermal = Heisenberg(n, Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)));
        Assert.Equal(Sextant.Palindrome.Witnessed, Sextant.Decide(thermal).Kind);
        var warm = Sextant.Setting.Chain(n, new[] { ("X", -1L), ("ZZ", 2L) },
            Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)).Concat(Every(n, LetterCube.LetterMatrix('Z'), 1)));
        var dWarm = Sextant.Decide(warm);
        Assert.Equal(Sextant.Palindrome.Obstructed, dWarm.Kind);
        Assert.Equal(3, Sextant.OddTraceObstruction(warm, 9)!.Value.Power);
        var dm = Sextant.Setting.Chain(n, new[] { ("XY", 1L), ("YX", -1L) }, Every(n, LetterCube.Lowering, 1));
        var dDm = Sextant.Decide(dm);
        Assert.Equal(Sextant.Palindrome.Witnessed, dDm.Kind);
        Assert.StartsWith("the theodolite", dDm.Detail);
        Assert.DoesNotContain(Sextant.Search(dm), f => f.Reading.Verdict is Sextant.Verdict.Reflection or Sextant.Verdict.PiReflection);
        Assert.Null(Sextant.OddTraceObstruction(dm, 64));                // the theodolite's yes meets no exact no
        // -3YY - YZ - 2XX under loss at (1, 0) is outside F158 (loss is no letter): the theodolite reads it in the frame
        // carrying X onto Z on the unwatched site
        var offAxis = Sextant.Setting.Chain(2, new[] { ("YY", -3L), ("YZ", -1L), ("XX", -2L) },
            new[] { (0, LetterCube.Lowering, 1L), (1, LetterCube.Lowering, 0L) });
        var dOff = Sextant.Decide(offAxis);
        Assert.Equal(Sextant.Palindrome.Witnessed, dOff.Kind);
        Assert.Contains("in the frame turning the axes ZX", dOff.Detail);
        // an unwatched site drops its jump and the end count decides (F158 section (a)): -3XY with fields -2X under Z
        // dephasing at (0, 3), X conserved on site 0 and no element of the closure reflecting, is witnessed there; so
        // is the same setting with a loss jump at rate zero on site 0, the same generator
        var unwatched = Sextant.Setting.Chain(2, new[] { ("XY", -3L), ("X", -2L) },
            new[] { (0, LetterCube.LetterMatrix('Z'), 0L), (1, LetterCube.LetterMatrix('Z'), 3L) });
        var dUnwatched = Sextant.Decide(unwatched);
        Assert.Equal(Sextant.Palindrome.Witnessed, dUnwatched.Kind);
        Assert.StartsWith("F158", dUnwatched.Detail);
        Assert.DoesNotContain(Sextant.Search(unwatched), f => f.Reading.Verdict is Sextant.Verdict.Reflection or Sextant.Verdict.PiReflection);
        var idleLoss = Sextant.Setting.Chain(2, new[] { ("XY", -3L), ("X", -2L) },
            new[] { (0, LetterCube.Lowering, 0L), (1, LetterCube.LetterMatrix('Z'), 3L) });
        Assert.StartsWith("F158", Sextant.Decide(idleLoss).Detail);
        // dephasing along X + Z, an axis between the letters: no element of the closure; F158 holds on any axis and a
        // colour on the circle witnesses it (A_Colour_On_The_Circle_Witnesses_A_Tilted_Axis)
        var tilt = new[,] { { LetterCube.Gauss.One, LetterCube.Gauss.One }, { LetterCube.Gauss.One, -LetterCube.Gauss.One } };
        var tilted = Sextant.Setting.Chain(2, new[] { ("XX", -3L), ("YX", 3L) }, new[] { (0, tilt, 2L), (1, tilt, 1L) });
        Assert.StartsWith("F158, a colour on the circle", Sextant.Decide(tilted).Detail);
        Assert.DoesNotContain(Sextant.Search(tilted), f => f.Reading.Verdict is Sextant.Verdict.Reflection or Sextant.Verdict.PiReflection);
        var unseen = Sextant.Setting.Chain(2, new[] { ("YZ", 1L), ("Z", -3L) },
            new[] { (0, LetterCube.LetterMatrix('X'), 1L), (1, LetterCube.Raising, 3L) });
        Assert.Equal(Sextant.Palindrome.SymmetricModP, Sextant.Decide(unseen).Kind);
        Assert.Equal(Sextant.Palindrome.Open, Sextant.Decide(unseen, maxPower: 3).Kind);
        Sextant.Setting Cooled(params long[] rates) => Sextant.Setting.Chain(n, new[] { ("XX", 1L), ("XY", 2L), ("ZZ", 1L) },
            Enumerable.Range(0, n).Select(l => (l, LetterCube.Lowering, rates[l])));
        var mirrored = Sextant.Decide(Cooled(1, 2, 1));
        Assert.Equal(Sextant.Palindrome.Witnessed, mirrored.Kind);
        Assert.Contains("relabelled l -> [2, 1, 0]", mirrored.Detail);
        Assert.Equal(7, Sextant.OddTraceObstruction(Cooled(1, 2, 3), 9)!.Value.Power);
        var bare = Sextant.Setting.Chain(n, new[] { ("XY", 1L), ("X", 1L), ("ZX", 1L) }, Array.Empty<(int, LetterCube.Gauss[,], long)>());
        Assert.Equal(Sextant.Palindrome.Witnessed, Sextant.Decide(bare).Kind);
    }

    // The DM reading by an exact second route: on two sites the generator is built densely from H and the jumps, scaled
    // to 4 dim M with Gaussian-integer entries, and its odd powers up to the dimension are traced in BigInteger
    // arithmetic: every one is exactly zero for DM bonds under loss (a second route beside the theodolite's element,
    // which the decision reports),
    // and the warm DNA model's third trace is the integer whose residue the obstruction reads.
    [Fact]
    public void The_DM_Symmetry_Is_Exact_On_Two_Sites()
    {
        const int n = 2;
        (System.Numerics.BigInteger Re, System.Numerics.BigInteger Im)[,] Scaled(Sextant.Setting s)
        {
            int d = 1 << n, dim = d * d;
            var H = new Complex[d, d];
            foreach (var (term, c) in s.H) Acc(H, DenseString(term.ToString(n)), c);
            var L = new Complex[dim, dim];
            for (int a = 0; a < dim; a++)
            {
                var e = new Complex[d, d]; e[a / d, a % d] = 1;
                var hr = Mul(H, e); var rh = Mul(e, H);
                var o = new Complex[d, d];
                for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) o[i, j] = -Complex.ImaginaryOne * (hr[i, j] - rh[i, j]);
                foreach (var (site, jump, rate) in s.Jumps)
                {
                    var J = Site(n, site, ToComplex(jump)); var Jd = Dag4(J); var JdJ = Mul(Jd, J);
                    var sand = Mul(J, Mul(e, Jd)); var l = Mul(JdJ, e); var r = Mul(e, JdJ);
                    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) o[i, j] += rate * (sand[i, j] - 0.5 * (l[i, j] + r[i, j]));
                }
                for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) L[i * d + j, a] = o[i, j];
            }
            Complex tr = 0; for (int i = 0; i < dim; i++) tr += L[i, i];
            var M = new (System.Numerics.BigInteger, System.Numerics.BigInteger)[dim, dim];
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    var v = 4 * dim * L[i, j] - (i == j ? 4 * tr : Complex.Zero);   // 4 dim (L + c/2), integral
                    Assert.True(v.Real == Math.Round(v.Real) && v.Imaginary == Math.Round(v.Imaginary));
                    M[i, j] = (new System.Numerics.BigInteger(v.Real), new System.Numerics.BigInteger(v.Imaginary));
                }
            return M;
        }
        (System.Numerics.BigInteger, System.Numerics.BigInteger)[,] MulG((System.Numerics.BigInteger Re, System.Numerics.BigInteger Im)[,] a, (System.Numerics.BigInteger Re, System.Numerics.BigInteger Im)[,] b)
        {
            int d = a.GetLength(0);
            var c = new (System.Numerics.BigInteger, System.Numerics.BigInteger)[d, d];
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                {
                    System.Numerics.BigInteger re = 0, im = 0;
                    for (int k = 0; k < d; k++)
                    {
                        re += a[i, k].Re * b[k, j].Re - a[i, k].Im * b[k, j].Im;
                        im += a[i, k].Re * b[k, j].Im + a[i, k].Im * b[k, j].Re;
                    }
                    c[i, j] = (re, im);
                }
            return c;
        }
        (System.Numerics.BigInteger Re, System.Numerics.BigInteger Im) Trace((System.Numerics.BigInteger Re, System.Numerics.BigInteger Im)[,] m)
        {
            System.Numerics.BigInteger re = 0, im = 0;
            for (int i = 0; i < m.GetLength(0); i++) { re += m[i, i].Re; im += m[i, i].Im; }
            return (re, im);
        }
        var dm = Sextant.Setting.Chain(n, new[] { ("XY", 1L), ("YX", -1L) }, Every(n, LetterCube.Lowering, 1));
        var M = Scaled(dm);
        var P = M;
        for (int m = 1; m <= 16; m++)
        {
            if (m > 1) P = MulG(P, M);
            if (m % 2 == 1) Assert.Equal((System.Numerics.BigInteger.Zero, System.Numerics.BigInteger.Zero), Trace(P));
        }
        Assert.Equal(Sextant.Palindrome.Witnessed, Sextant.Decide(dm).Kind);
        var warm = Sextant.Setting.Chain(n, new[] { ("X", -1L), ("ZZ", 2L) },
            Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)).Concat(Every(n, LetterCube.LetterMatrix('Z'), 1)));
        var W = Scaled(warm);
        var t3 = Trace(MulG(MulG(W, W), W));
        Assert.NotEqual(System.Numerics.BigInteger.Zero, t3.Re);
        var (power, prime) = Sextant.OddTraceObstruction(warm, 15)!.Value;
        Assert.Equal(3, power);
        Assert.NotEqual(System.Numerics.BigInteger.Zero, System.Numerics.BigInteger.Remainder(t3.Re, prime));
    }

    // The traces by a second route: Tr(M^m), M = L + c/2, built densely in complex numbers on two sites from H and the
    // jumps, is not zero at m = 3 for the warm DNA model and zero at m = 1, 3, 5 for the thermal Heisenberg chain, as
    // the modular reading says (all entries small dyadics, compared exactly).
    [Fact]
    public void The_Odd_Traces_Agree_With_The_Dense_Generator()
    {
        const int n = 2;
        Complex OddTrace(Sextant.Setting s, int m)
        {
            int d = 1 << n, dim = d * d;
            var H = new Complex[d, d];
            foreach (var (term, c) in s.H) Acc(H, DenseString(term.ToString(n)), c);
            var L = new Complex[dim, dim];
            for (int a = 0; a < dim; a++)
            {
                var e = new Complex[d, d]; e[a / d, a % d] = 1;
                var hr = Mul(H, e); var rh = Mul(e, H);
                var out1 = new Complex[d, d];
                for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) out1[i, j] = -Complex.ImaginaryOne * (hr[i, j] - rh[i, j]);
                foreach (var (site, jump, rate) in s.Jumps)
                {
                    var J = Site(n, site, ToComplex(jump)); var Jd = Dag4(J); var JdJ = Mul(Jd, J);
                    var sand = Mul(J, Mul(e, Jd)); var l = Mul(JdJ, e); var r = Mul(e, JdJ);
                    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) out1[i, j] += rate * (sand[i, j] - 0.5 * (l[i, j] + r[i, j]));
                }
                for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) L[i * d + j, a] = out1[i, j];
            }
            Complex tr = 0; for (int i = 0; i < dim; i++) tr += L[i, i];
            for (int i = 0; i < dim; i++) L[i, i] -= tr / dim;          // M = L + c/2, c = -2 Tr L / dim
            var P = L;
            for (int k = 1; k < m; k++) P = Mul(P, L);
            Complex t = 0; for (int i = 0; i < dim; i++) t += P[i, i];
            return t;
        }
        var warm = Sextant.Setting.Chain(n, new[] { ("X", -1L), ("ZZ", 2L) },
            Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)).Concat(Every(n, LetterCube.LetterMatrix('Z'), 1)));
        Assert.True(OddTrace(warm, 1) == Complex.Zero);
        Assert.True(OddTrace(warm, 3) != Complex.Zero);
        Assert.Equal(3, Sextant.OddTraceObstruction(warm, 9)!.Value.Power);
        var thermal = Heisenberg(n, Every(n, LetterCube.Lowering, 3).Concat(Every(n, LetterCube.Raising, 1)));
        foreach (int m in new[] { 1, 3, 5 }) Assert.True(OddTrace(thermal, m) == Complex.Zero);
        Assert.Null(Sextant.OddTraceObstruction(thermal, 5));
    }

    // Five palindromes at N = 3 from a random sweep that no element of the one-site closure reflects: Pauli-letter jumps,
    // so F158 decides them, and EndCount's verdict is exact for each (by an element or by the two lifted counts); the
    // sextant reports them WITNESSED through it. Where EndCount names no element, a product with a turned letter is one:
    // X⊗Z⊗(Y − 2X) in the first setting and Z⊗(2X − Y)⊗Z in the third, every factor anticommuting with its site's jumps,
    // the product commuting with H, its square a nonzero scalar, so it is invertible and lies in F158's far kernel.
    [Fact]
    public void Five_Palindromes_Outside_The_Closure_Are_Witnessed_Through_The_End_Count()
    {
        var cases = new (string[] Terms, long[] Coeffs, string[] Jumps, long[] Rates)[]
        {
            (new[] { "IXZ", "IZX", "IZY", "XZI", "ZXI", "ZYI" }, new[] { 1L, 2, -1, 1, 2, -1 }, new[] { "IYI" }, new[] { 1L }),
            (new[] { "IYX", "IZY", "IZZ", "YXI", "ZYI", "ZZI" }, new[] { 2L, 2, -2, 2, 2, -2 }, new[] { "XII" }, new[] { 1L }),
            (new[] { "IXZ", "IYZ", "XZI", "YZI" }, new[] { 2L, -1, 2, -1 }, new[] { "YII", "IIX" }, new[] { 3L, 1 }),
            (new[] { "IXY", "IXZ", "XYI", "XZI" }, new[] { 2L, 1, 2, 1 }, new[] { "ZII", "YII", "IXI", "IIX" }, new[] { 2L, 3, 2, 1 }),
            (new[] { "IXY", "IXZ", "IYX", "XYI", "XZI", "YXI" }, new[] { 1L, -1, -2, 1, -1, -2 }, new[] { "XII", "IYI", "IIX" }, new[] { 3L, 3, 1 }),
        };
        var counts = new List<EndCount>();
        foreach (var (t, c, j, r) in cases)
        {
            var h = t.Zip(c, (a, b) => (a, b)).ToList();
            var jumps = j.Zip(r, (s, rate) => { int l = s.IndexOfAny(new[] { 'X', 'Y', 'Z' }); return (l, LetterCube.LetterMatrix(s[l]), rate); });
            var setting = new Sextant.Setting(3, h.Select(x => (PauliString.Parse(x.a), x.b)), jumps);
            Assert.DoesNotContain(Sextant.Search(setting), f => f.Reading.Verdict is Sextant.Verdict.Reflection or Sextant.Verdict.PiReflection);
            var d = Sextant.Decide(setting);
            Assert.Equal(Sextant.Palindrome.Witnessed, d.Kind);
            Assert.StartsWith("F158", d.Detail);
            var count = new EndCount(new World(), 3, h, j);
            Assert.True(EndCount.IsExact(count.Verdict()) && EndCount.IsPalindrome(count.Verdict()));
            counts.Add(count);
        }
        foreach (var (index, element) in new[] { (0, new[] { ("XZY", 1L), ("XZX", -2L) }), (2, new[] { ("ZXZ", 2L), ("ZYZ", -1L) }) })
        {
            var check = counts[index].CheckElement(element);
            Assert.True(check.AllLit && check.CommutesWithH);
            Assert.Equal(new System.Numerics.BigInteger(5), check.SquareScalar);
        }
        // the control: a turned letter at the wrong angle is not in the far kernel
        Assert.False(counts[0].CheckElement(new[] { ("XZY", 1L), ("XZX", 2L) }).CommutesWithH);
    }

    // The end count's no, and the sites: H = ZZ + XI with X watched on site 0 and Z on site 1 has the ends 2 and 0, so
    // F158's necessity rules every element out; the same H with the two watchings swapped pairs (the string XY is a
    // colouring). A dock that read the sites reversed would answer each with the other's verdict.
    [Fact]
    public void The_End_Count_Says_No_And_Reads_The_Sites_In_Order()
    {
        var h = new[] { (PauliString.Parse("ZZ"), 1L), (PauliString.Parse("XI"), 1L) };
        var X = LetterCube.LetterMatrix('X'); var Z = LetterCube.LetterMatrix('Z');
        var broken = Sextant.Decide(new Sextant.Setting(2, h, new[] { (0, X, 1L), (1, Z, 1L) }));
        Assert.Equal(Sextant.Palindrome.Obstructed, broken.Kind);
        Assert.StartsWith("F158", broken.Detail);
        var swapped = Sextant.Decide(new Sextant.Setting(2, h, new[] { (1, X, 1L), (0, Z, 1L) }));
        Assert.Equal(Sextant.Palindrome.Witnessed, swapped.Kind);
        var e = new EndCount(new World(), 2, new[] { ("ZZ", 1L), ("XI", 1L) }, new[] { "XI", "IZ" });
        Assert.Equal((2, 0), e.UpperCounts());
        Assert.Throws<ArgumentException>(() => new Sextant.Setting(2, new[] { (PauliString.Parse("XXX"), 1L) }, new[] { (0, X, 1L) }));
    }
}
