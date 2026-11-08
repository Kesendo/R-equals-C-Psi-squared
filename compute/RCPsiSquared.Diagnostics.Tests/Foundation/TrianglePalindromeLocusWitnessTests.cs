using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="TrianglePalindromeLocusWitness"/> (sections "The triangle's palindromic
/// locus at N = 3" and "The triangle with planar fields", and the paragraph "The sheet and blindness, at N = 3 and N = 4",
/// of <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>): the rows the proof, its gates and
/// <c>EndCountTests</c> pin, read on the witness's own exact route, the grids with a planted mutation per class or
/// condition, and the two settings read alike where they meet.</summary>
public class TrianglePalindromeLocusWitnessTests
{
    static BigRational[] T(string s) => s.Split(',').Select(p =>
    {
        var q = p.Split('/');
        return new BigRational(long.Parse(q[0]), q.Length > 1 ? long.Parse(q[1]) : 1);
    }).ToArray();

    [Theory]
    [InlineData("ZXY", "3,6,-2", "0,1,1", 1, 1, true, 'C')]     // the reciprocal family's pinned row
    [InlineData("ZXY", "3,6,-2", "0,0,0", 10, 10, true, 'A')]   // at h = 0 it is coloured
    [InlineData("ZXY", "3,6,-3", "0,1,1", 1, 0, false, null)]   // off the locus it breaks
    [InlineData("ZXY", "2,2,7", "0,1,1", 1, 1, true, 'B')]      // the swap branch
    [InlineData("XXX", "1,2,3", "1,2,3", 1, 1, true, 'A')]      // a colouring
    public void Pinned_Rows(string word, string j, string h, int near, int far, bool pal, char? cls)
    {
        var r = TrianglePalindromeLocusWitness.ReadRow(word, T(j), T(h));
        Assert.Equal((near, far, pal, cls), (r.Near, r.Far, r.Palindrome, r.InClass));
        Assert.True(r.Agrees);
    }

    /// <summary>The near and far counts are n3_planar_field_locus_gate.py G7's, there from the block equations over
    /// Q(i): the E row and the corner pair with near = far = 1, each perturbation leaves near 1 and far 0.</summary>
    [Theory]
    [InlineData("1,2,6", "8,0", "-5,3", 1, 1, true, 'E')]       // the E row (w = 1, xi = 8)
    [InlineData("1,2,6", "8,0", "-5,4", 1, 0, false, null)]     // its field perturbed
    [InlineData("1,2,7", "8,0", "-5,3", 1, 0, false, null)]     // its coupling perturbed
    [InlineData("2,-2,0", "1,4", "4,1", 1, 1, true, 'E')]       // the corner J_uv = 0, f_u.f_v = 2 J_du^2
    [InlineData("2,-2,0", "1,4", "4,-1", 1, 0, false, null)]    // the corner perturbed (f_u.f_v = 0)
    [InlineData("2,-2,0", "1,4", "-1,4", 1, 0, false, null)]    // the corner perturbed (f_u.f_v = 15)
    [InlineData("-4,-1,-1", "6,0", "-3,3", 1, 0, false, null)]  // R2 = R3 = 0 with R1 = 11: E needs its first relation
    public void Pinned_Planar_Rows(string j, string fu, string fv, int near, int far, bool pal, char? cls)
    {
        var r = TrianglePalindromeLocusWitness.ReadPlanarRow(T(j), T(fu), T(fv));
        Assert.Equal((near, far, pal, cls), (r.Near, r.Far, r.Palindrome, r.InClass));
        Assert.True(r.Agrees);
    }

    [Theory]
    [InlineData("3,6,-2", "3,4", "5,0", true, 'C')]       // C at an angle neither parallel nor perpendicular
    [InlineData("3,6,-2", "3,4", "4,4", false, null)]     // unequal lengths on e2 = 0
    [InlineData("2,2,3", "3,4", "0,5", true, 'B')]        // the swap class at a general angle
    [InlineData("0,2,3", "3,4", "5,0", false, null)]      // d at an end: equal lengths do not suffice
    [InlineData("0,2,3", "1,2", "-2,-4", true, 'A')]      // d at an end: parallel fields pair
    public void Planar_Classes_At_General_Angles(string j, string fu, string fv, bool pal, char? cls)
    {
        var r = TrianglePalindromeLocusWitness.ReadPlanarRow(T(j), T(fu), T(fv));
        Assert.Equal((pal, cls), (r.Palindrome, r.InClass));
        Assert.True(r.Agrees);
    }

    [Fact]
    public void The_Grid_Agrees_And_Every_Class_Occurs()
    {
        var g = TrianglePalindromeLocusWitness.ReadGrid();
        Assert.Equal(27 * 8 * 6, g.Rows);
        Assert.Equal(0, g.Disagreements);
        Assert.True(g.A > 0 && g.B > 0 && g.C > 0, $"A {g.A}, B {g.B}, C {g.C}");
    }

    [Fact]
    public void The_Planar_Grid_Agrees_And_Every_Class_Occurs()
    {
        var g = TrianglePalindromeLocusWitness.ReadPlanarGrid();
        Assert.Equal(TrianglePalindromeLocusWitness.PlanarCouplings.Length * TrianglePalindromeLocusWitness.PlanarFieldPairs.Length, g.Rows);
        Assert.Equal(0, g.Disagreements);
        Assert.True(g.A > 0 && g.B > 0 && g.C > 0 && g.E > 0, $"A {g.A}, B {g.B}, C {g.C}, E {g.E}");
        Assert.True(g.Palindromic < g.Rows, "the grid holds broken rows too");
    }

    [Theory]
    [InlineData('A')]
    [InlineData('B')]
    [InlineData('C')]
    public void Leaving_A_Class_Out_Breaks_The_Grid(char without)
    {
        Assert.True(TrianglePalindromeLocusWitness.ReadGrid(without).Disagreements > 0);
    }

    [Theory]
    [InlineData('A')]
    [InlineData('B')]
    [InlineData('C')]
    [InlineData('E')]
    public void Leaving_A_Class_Out_Breaks_The_Planar_Grid(char without)
    {
        Assert.True(TrianglePalindromeLocusWitness.ReadPlanarGrid(without).Disagreements > 0);
    }

    /// <summary>Where the settings meet (no field on d, letters X or Y on u and v) the two rules give the same verdict.
    /// The two ends are one kernel routine read twice, so their equality checks only the planar wiring (f_u on site 1,
    /// f_v on site 2, nothing on d); the content is that the letter rule and the planar rule agree with it.</summary>
    [Fact]
    public void The_Two_Rules_Agree_Where_The_Settings_Meet()
    {
        int shared = 0;
        foreach (var (w, j, h) in TrianglePalindromeLocusWitness.Grid())
        {
            if (!h[0].IsZero || w[1] == 'Z' || w[2] == 'Z') continue;
            BigRational[] V(char letter, BigRational x) => letter == 'X' ? new[] { x, BigRational.Zero } : new[] { BigRational.Zero, x };
            var lr = TrianglePalindromeLocusWitness.ReadRow(w, j, h);
            var pr = TrianglePalindromeLocusWitness.ReadPlanarRow(j, V(w[1], h[1]), V(w[2], h[2]));
            Assert.Equal((lr.Near, lr.Far, lr.Palindrome), (pr.Near, pr.Far, pr.Palindrome));
            Assert.True(lr.Agrees && pr.Agrees);
            shared++;
        }
        Assert.True(shared > 0);
    }

    /// <summary>A field on d: the two ends exactly, no pairing rule. The counts are a second route's (sympy over Q(i), the
    /// commutant on the 64 strings): all fields parallel, B with f_d on the swap axis and D (d at an end, f_d ⊥ f_v,
    /// f_d + f_u + f_v = 0) pair, each perturbation breaks, and the u, v singlet (J_du = J_dv, f_u = f_v) holds a far
    /// element but a second near one.</summary>
    [Theory]
    [InlineData("1,2,3", "1,0", "2,0", "3,0", 1, 1)]      // A, all three parallel
    [InlineData("2,2,3", "1,1", "3,4", "4,3", 1, 1)]      // B, f_d on the swap axis
    [InlineData("2,2,3", "1,0", "3,4", "4,3", 1, 0)]      // B, f_d off the axis
    [InlineData("1,0,3", "1,0", "-1,-2", "0,2", 1, 1)]    // D: d-u-v, f_d = X, f_v = 2Y, the sum zero
    [InlineData("1,0,3", "1,0", "-1,-2", "0,3", 1, 0)]    // D perturbed
    [InlineData("2,2,3", "1,2", "3,1", "3,1", 2, 1)]      // the singlet: a far element, two near ones
    [InlineData("1,2,6", "1,0", "8,0", "-5,3", 1, 0)]     // the E row with a field on d
    public void A_Field_On_D_Is_Read_Without_A_Pairing_Rule(string j, string fd, string fu, string fv, int near, int far)
    {
        var r = TrianglePalindromeLocusWitness.ReadPlanarRow(T(j), T(fd), T(fu), T(fv));
        Assert.Equal((near, far, false), (r.Near, r.Far, r.RuleStated));
        Assert.Null(r.InClass);
    }

    /// <summary>The far element's operator Schmidt ranks across d|uv, u|dv, v|du, exact over Q; the expected ranks
    /// are a second route's (sympy, the kernel vector over Q(i) reshaped). Only A's is a product; B's is one across
    /// d|uv alone (letters times SWAP_uv), C's, D's and E's across no cut.</summary>
    [Theory]
    [InlineData("1,2,3", "1,0", "2,0", "3,0", 1, 1, 1)]   // A: X X X
    [InlineData("2,2,3", "1,1", "3,4", "4,3", 1, 4, 4)]   // B: (m·σ)_d ⊗ SWAP·(m·σ ⊗ m·σ)
    [InlineData("1,0,3", "1,0", "-1,-2", "0,2", 2, 3, 2)] // D: XXX + c·XIY + c'·YYY
    [InlineData("3,6,-2", "0,0", "3,4", "5,0", 2, 4, 4)]  // C
    [InlineData("2,-2,0", "0,0", "1,4", "4,1", 2, 4, 4)]  // E's corner
    [InlineData("1,2,6", "0,0", "8,0", "-5,3", 2, 4, 4)]  // E
    public void The_Far_Element_Is_A_Product_Only_In_A(string j, string fd, string fu, string fv, int dUv, int uDv, int vDu)
    {
        var loc = TrianglePalindromeLocusWitness.ReadFarLocality(T(j),
            TrianglePalindromeLocusWitness.PlanarFields(T(fd), T(fu), T(fv)));
        Assert.NotNull(loc);
        Assert.Equal((dUv, uDv, vDu), (loc!.DUv, loc.UDv, loc.VDu));
        Assert.Equal(dUv == 1 && uDv == 1 && vDu == 1, loc.Product);
    }

    /// <summary>The far element's kind patterns (d, u, v; 1 for X or Y), against simulations/framework/diagnostics/
    /// far_element_shape.py, a separate route over Q(i): A's product, B's letters times SWAP and C at three rows of
    /// e₂ = 0 keep the pair (u, v) on one parity, D keeps v lit, E uses every pattern with d lit.</summary>
    [Theory]
    [InlineData("1,2,3", "1,0", "2,0", "3,0", "111")]                // A
    [InlineData("2,2,3", "1,1", "3,4", "4,3", "100,111")]            // B
    [InlineData("3,6,-2", "0,0", "3,4", "5,0", "100,111")]           // C
    [InlineData("4,12,-3", "0,0", "7,1", "5,5", "100,111")]          // C at t = 3, another angle
    [InlineData("-2,6,3", "0,0", "2,1", "1,2", "100,111")]           // C at t = -3
    [InlineData("1,0,3", "1,0", "-1,-2", "0,2", "101,111")]          // D
    [InlineData("1,2,6", "0,0", "8,0", "-5,3", "100,101,110,111")]   // E
    public void The_Kind_Patterns_Of_The_Far_Element(string j, string fd, string fu, string fv, string patterns)
    {
        var p = TrianglePalindromeLocusWitness.ReadFarPatterns(T(j), TrianglePalindromeLocusWitness.PlanarFields(T(fd), T(fu), T(fv)));
        Assert.Equal(patterns.Split(','), p);
    }

    [Fact]
    public void The_Locality_Reading_Is_Made_Only_Where_The_Row_Pairs_Through_One_Far_Element()
    {
        // the singlet row: one far element but two near ones, so the row does not pair and no mirror exists at all
        Assert.Null(TrianglePalindromeLocusWitness.ReadFarLocality(T("2,2,3"),
            TrianglePalindromeLocusWitness.PlanarFields(T("1,2"), T("3,1"), T("3,1"))));
        // a broken row: no far element
        Assert.Null(TrianglePalindromeLocusWitness.ReadFarLocality(T("1,2,6"),
            TrianglePalindromeLocusWitness.PlanarFields(T("1,0"), T("8,0"), T("-5,3"))));
    }

    [Fact]
    public void The_Witness_Reads_A_Field_On_D()
    {
        var w = new TrianglePalindromeLocusWitness(null, "1,0,3", null, "-1,-2", "0,2", "1,0");
        var r = w.Read();
        Assert.True(r.Palindrome && !r.RuleStated);
        Assert.Equal(new TrianglePalindromeLocusWitness.FarLocality(2, 3, 2), w.ReadFarLocality());
        Assert.Contains("no site-local mirror", w.Summary);
    }

    [Fact]
    public void A_Pairing_Grid_Refuses_A_Row_Without_A_Pairing_Rule()
    {
        var unruled = TrianglePalindromeLocusWitness.ReadPlanarRow(T("1,0,3"), T("1,0"), T("-1,-2"), T("0,2"));
        var tally = typeof(TrianglePalindromeLocusWitness).GetMethod("Tally",
            System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)!;
        var ex = Assert.Throws<System.Reflection.TargetInvocationException>(() => tally.Invoke(null, new object[] { new[] { unruled } }));
        Assert.IsType<ArgumentException>(ex.InnerException);
    }

    [Fact]
    public void A_Disconnected_Row_Is_Refused()
    {
        Assert.Throws<ArgumentException>(() => new TrianglePalindromeLocusWitness("ZXY", "0,0,3", "0,1,1"));
        Assert.Throws<ArgumentException>(() => new TrianglePalindromeLocusWitness(null, "0,0,3", null, "8,0", "-5,3"));
    }

    /// <summary>The sheet and blindness (the paragraph of that name in the proof's section on blind dephased seats),
    /// on the witness's own route. The far counts by height are the second route of simulations/blind_seat_sheet_gate.py
    /// (its table of far elements by height, sympy over Q(i) on the Hermitian W = |1⟩⟨0|_d ⊗ M + h.c. with the pair's
    /// parity kept or flipped): C and B hold a far element of odd height, B also at a second coupling with f_d on its
    /// line; unequal lengths, a seat that is not blind, d at an end, a field on d off B's axis and one at C's couplings hold none; E and D
    /// hold an element of no definite height (all 1, odd and even 0) at a seat that is not blind. Outside the setting
    /// the rule is not read: parallel and antiparallel fields (A's colouring, at a seat that is not blind), the singlet
    /// at J_du = J_dv with parallel fields of one length and f_d not along them (off B's line), and a disconnected row (u isolated, unequal lengths)
    /// each hold a far element of odd height that the rule would not predict, which is why the fences are there.</summary>
    [Theory]
    [InlineData("3,6,-2", "0,0", "3,4", "5,0", 1, 1, 0, 1, true)]     // C
    [InlineData("3,6,-2", "0,0", "3,4", "4,4", 1, 0, 0, 0, false)]    // C, unequal lengths
    [InlineData("2,2,3", "0,0", "3,4", "0,5", 1, 1, 0, 1, true)]      // B
    [InlineData("1,2,3", "0,0", "3,4", "5,0", 0, 0, 0, 0, false)]     // not blind
    [InlineData("1,2,6", "0,0", "8,0", "-5,3", 0, 0, 0, 1, false)]    // E: an element of no definite height
    [InlineData("0,2,3", "0,0", "3,4", "5,0", 0, 0, 0, 0, false)]     // d at an end
    [InlineData("3,3,0", "0,0", "3,4", "5,0", 1, 1, 0, 1, true)]      // B at J_uv = 0
    [InlineData("2,2,3", "1,1", "3,4", "4,3", 1, 1, 0, 1, true)]      // B, f_d on the line of f_u + f_v
    [InlineData("2,2,3", "1,0", "3,4", "4,3", 1, 0, 0, 0, false)]     // B, f_d off it
    [InlineData("3,6,-2", "8,4", "3,4", "5,0", 1, 0, 0, 0, false)]    // C, f_d on the bisector: none
    [InlineData("1,0,3", "1,0", "-1,-2", "0,2", 0, 0, 0, 1, false)]   // D: no definite height
    [InlineData("1,1,-5", "1,1", "3,4", "4,3", 1, 1, 0, 1, true)]     // B at a second coupling, f_d on the line
    [InlineData("2,2,-1", "0,0", "3,4", "5,0", 1, 1, 0, 1, true)]     // B and C at once: the swap branch alone
    public void The_Sheet_And_Blindness_At_Pinned_Rows(string j, string fd, string fu, string fv,
        int blind, int odd, int even, int all, bool predicted)
    {
        var r = TrianglePalindromeLocusWitness.ReadSheet(T(j), T(fd), T(fu), T(fv));
        Assert.Equal((blind, odd, even, all, true, predicted), (r.Blind, r.FarOdd, r.FarEven, r.FarAll, r.InSetting, r.Predicted));
        Assert.True(r.Agrees);
    }

    [Theory]
    [InlineData("1,2,3", "0,0", "1,2", "2,4", 1, 0, 1)]      // parallel fields: A's colouring, at a seat that is not blind
    [InlineData("1,2,3", "0,0", "1,0", "-1,0", 1, 0, 1)]     // antiparallel fields: A again
    [InlineData("2,2,3", "1,2", "3,1", "3,1", 1, 0, 1)]      // the singlet at J_du = J_dv, f_d not along the fields
    [InlineData("0,1,0", "0,0", "3,4", "2,0", 1, 1, 2)]      // disconnected: u isolated, unequal lengths
    public void Outside_The_Setting_The_Rule_Is_Not_Read(string j, string fd, string fu, string fv, int odd, int even, int all)
    {
        var r = TrianglePalindromeLocusWitness.ReadSheet(T(j), T(fd), T(fu), T(fv));
        Assert.Equal((odd, even, all, false, false), (r.FarOdd, r.FarEven, r.FarAll, r.InSetting, r.Predicted));
    }

    /// <summary>The Krylov count of seat d and the factor (J_du − J_dv)·e₂ of the proof name the same rows, on every
    /// coupling of the planar grid (all connected): F157 at N = 3.</summary>
    [Fact]
    public void Blindness_By_Krylov_Rank_Is_The_Factor()
    {
        foreach (var j in TrianglePalindromeLocusWitness.PlanarCouplings)
        {
            var e2 = j[0] * j[1] + j[0] * j[2] + j[1] * j[2];
            Assert.Equal(((j[0] - j[1]) * e2).IsZero, TrianglePalindromeLocusWitness.BlindDimension(j) > 0);
        }
    }

    /// <summary>The sheet grid's tallies, counted by hand from its composition: the planar grid's 13 couplings and the
    /// 14 of its 16 field pairs that are neither parallel, antiparallel nor zero, each with no field on d and three on
    /// it, 13·14·4 = 728 rows, 546 with a field on d. A far element of odd height at 20: the three blind couplings
    /// ((2, 2, 3), (3, 6, −2), (2, 5, −10/7)) on the five pairs of equal lengths with no field on d, 15, where the closed
    /// form is read, and B's coupling on those pairs with f_d on the line of f_u + f_v, 5. Read, not counted: elements
    /// of other heights on as many rows as the planar grid has E rows, none of even height.</summary>
    [Fact]
    public void The_Sheet_Grid_Agrees_And_The_Closed_Form_Holds()
    {
        var g = TrianglePalindromeLocusWitness.ReadSheetGrid();
        Assert.Equal((728, 546, 20, 0), (g.Rows, g.WithFieldOnD, g.OddElements, g.Disagreements));
        Assert.Equal(8, g.OtherHeights);   // the planar grid's E rows: (1, 2, 6) with three pairs, its swap, double and corners
        Assert.Equal(TrianglePalindromeLocusWitness.ReadPlanarGrid().E, g.OtherHeights);
        Assert.Equal(0, g.EvenElements);
        Assert.Equal((15, 0), (g.ClosedForms, g.ClosedFormFailures));
    }

    /// <summary>The live node the proof and the claim point to, at the census row and at B with a field on d.</summary>
    [Fact]
    public void The_Live_Node_Reads_The_Sheet_And_Blindness()
    {
        string Node(TrianglePalindromeLocusWitness w) =>
            w.Children.Single(c => c.DisplayName == "the sheet and blindness").Summary;
        var c = Node(new TrianglePalindromeLocusWitness(null, "3,6,-2", null, "3,4", "5,0"));
        Assert.Contains("blind of dimension 1", c);
        Assert.Contains("odd height 1, of even height 0, all 1", c);
        Assert.Contains("The rule (odd height exactly when d is blind", c);
        Assert.DoesNotContain("The rule with a field on d", c);
        Assert.Contains("predicts at least one and holds", c);
        Assert.Contains("lies in the far end and spans it", c);
        Assert.Contains("728 rows", c);
        var b = Node(new TrianglePalindromeLocusWitness(null, "2,2,3", null, "3,4", "4,3", "1,1"));
        Assert.Contains("The rule with a field on d", b);
        Assert.Contains("predicts at least one and holds", b);
        Assert.DoesNotContain("closed form lies", b);
    }

    /// <summary>Each condition of the rule carries weight, in both branches, and the disagreements are counted by hand
    /// from the grid's composition (13 couplings, 14 field pairs of which 5 have equal lengths, 3 blind couplings of
    /// which 1 is B's, and f_d on the line of f_u + f_v at one of its three values). Dropping blindness, or B's couplings
    /// with a field on d: 13·5 − 15 = 50 rows with no field on d and 13·5 − 5 = 60 with one, 110. Dropping the equal
    /// lengths: 3·14 − 15 = 27 and 14 − 5 = 9, 36. Dropping the axis of f_d: 5·3 − 5 = 10.</summary>
    [Theory]
    [InlineData('b', 110)]
    [InlineData('l', 36)]
    [InlineData('a', 10)]
    public void Dropping_A_Condition_Breaks_The_Sheet_Grid(char without, int disagreements)
    {
        Assert.Equal(disagreements, TrianglePalindromeLocusWitness.ReadSheetGrid(without).Disagreements);
    }

    /// <summary>The closed form: on C a far element of rank two across d|uv spanning the one-dimensional far end, at the
    /// census row and at another angle; on B rank one.</summary>
    [Theory]
    [InlineData("3,6,-2", "3,4", "5,0", 2)]      // the census row, ζ_swap/ζ_prod = 16i/15
    [InlineData("-2,6,3", "1,2", "2,1", 2)]      // C at another chart point and angle
    [InlineData("2,2,3", "3,4", "0,5", 1)]       // B
    [InlineData("2,2,-1", "3,4", "5,0", 1)]      // B and C at once: B's branch, never a division by J_du - J_dv
    public void The_Closed_Form_Is_The_Far_Element(string j, string fu, string fv, int rank)
    {
        var c = TrianglePalindromeLocusWitness.ReadClosedForm(T(j), T(fu), T(fv));
        Assert.NotNull(c);
        Assert.True(c!.InFarEnd && c.SpansFarEnd && c.SquareIsScalar, $"{c}");
        Assert.Equal(rank, c.RankDUv);
    }

    /// <summary>No closed form at unequal lengths, at antiparallel fields, or at a seat that is not blind.</summary>
    [Fact]
    public void The_Closed_Form_Needs_Equal_Lengths_An_Angle_And_A_Blind_Seat()
    {
        Assert.Null(TrianglePalindromeLocusWitness.ClosedForm(T("3,6,-2"), T("3,4"), T("4,4")));
        Assert.Null(TrianglePalindromeLocusWitness.ClosedForm(T("3,6,-2"), T("3,4"), T("-3,-4")));
        Assert.Null(TrianglePalindromeLocusWitness.ClosedForm(T("1,2,3"), T("3,4"), T("5,0")));
    }

    /// <summary>The census row's numbers in the closed form: n_prod = (J_du + J_uv)·f_u + J_uv·f_v = (−7, 4) and,
    /// ζ_swap/ζ_prod being 16i/15, n_swap = (16/15)·(−4, −7). Only the swap part reaches the strings X_d ⊗ 1 and
    /// Y_d ⊗ 1 (SWAP′X′X′ holds 1/2 of the identity), so they carry n_swap / 2. That the element is the far element is
    /// the test above; this one pins the arithmetic of the construction at the row the proof quotes.</summary>
    [Fact]
    public void The_Census_Row_Has_Its_Ratio()
    {
        var w = TrianglePalindromeLocusWitness.ClosedForm(T("3,6,-2"), T("3,4"), T("5,0"))!;
        // n_prod = (J_du + J_uv)·f_u + J_uv·f_v = (3, 4) − 2·(5, 0) = (−7, 4); n_swap = (16/15)·(−4, −7)
        // X_d ⊗ 1 carries n_swap,x / 2 and Y_d ⊗ 1 carries n_swap,y / 2 (SWAP′X′X′ holds 1/2 of the identity)
        Assert.Equal(new BigRational(-32, 15), w[16 * 1 + 0]);
        Assert.Equal(new BigRational(-56, 15), w[16 * 2 + 0]);
        Assert.Equal(12, w.Count(x => !x.IsZero));   // twelve strings, as the gate reads at this row
    }

    [Fact]
    public void The_Two_Settings_Are_Not_Mixed()
    {
        Assert.Throws<ArgumentException>(() => new TrianglePalindromeLocusWitness("ZXY", "1,2,6", null, "8,0", "-5,3"));
        var w = new TrianglePalindromeLocusWitness(null, "1,2,6", null, "8,0", "-5,3");
        Assert.True(w.Read().Palindrome && w.Read().InClass == 'E');
    }
}
