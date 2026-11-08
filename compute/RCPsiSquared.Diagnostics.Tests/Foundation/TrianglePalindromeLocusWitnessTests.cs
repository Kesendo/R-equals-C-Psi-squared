using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="TrianglePalindromeLocusWitness"/> (sections "The triangle's palindromic
/// locus at N = 3" and "The triangle with planar fields" of <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>):
/// the rows the proof, its gates and <c>EndCountTests</c> pin, read on the witness's own exact route, both grids with
/// a planted mutation per class, and the two settings read alike where they meet.</summary>
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

    /// <summary>A field on d: the two ends exactly, no rule. The counts are a second route's (sympy over Q(i), the
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
    public void A_Field_On_D_Is_Read_Without_A_Rule(string j, string fd, string fu, string fv, int near, int far)
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
    public void A_Grid_Refuses_A_Row_Without_A_Rule()
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

    [Fact]
    public void The_Two_Settings_Are_Not_Mixed()
    {
        Assert.Throws<ArgumentException>(() => new TrianglePalindromeLocusWitness("ZXY", "1,2,6", null, "8,0", "-5,3"));
        var w = new TrianglePalindromeLocusWitness(null, "1,2,6", null, "8,0", "-5,3");
        Assert.True(w.Read().Palindrome && w.Read().InClass == 'E');
    }
}
