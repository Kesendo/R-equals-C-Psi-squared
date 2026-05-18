using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F1DepolResidualClosedFormPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F1DepolResidualClosedFormPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F1DepolResidualClosedFormPi2Inheritance Build() =>
        new F1DepolResidualClosedFormPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // F1 depol closed form is Tier1Derived (PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM
        // + bit-exact verification N=2..5). Both coefficients (16/9, 16) reduce to
        // trivial algebra over the Pi2-Foundation primitives, so the inheritance
        // claim sits at the same Tier as its parent (no algebra gap).
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void DSquared_IsExactlyFour_FromLadder()
    {
        // d² = 4 from Pi2DyadicLadderClaim.Term(-1) = a_{-1}. (Cross-check
        // available via Pi2OperatorSpaceMirrorClaim.PairAt(1).OperatorSpace
        // which feeds DSquaredMinusOne below.)
        var f = Build();
        Assert.Equal(4.0, f.DSquared, precision: 14);
    }

    [Fact]
    public void DSquaredMinusOne_IsExactlyThree()
    {
        // d² − 1 = 3 (non-identity Paulis per qubit). Same anchor F5 uses for
        // its scalar 2/3 denominator.
        var f = Build();
        Assert.Equal(3.0, f.DSquaredMinusOne, precision: 14);
    }

    [Fact]
    public void PerPauliDepolarizingRate_IsExactly4Over3()
    {
        // d²/(d²−1) = 4/3. One power of d higher than F5's d/(d²−1) = 2/3
        // because the residual squares the rate before per-site assembly.
        var f = Build();
        Assert.Equal(4.0 / 3.0, f.PerPauliDepolarizingRate, precision: 14);
    }

    [Fact]
    public void LocalCoefficient_IsExactly16Over9()
    {
        // (d²/(d²−1))² = (4/3)² = 16/9. Squaring the per-Pauli rate is what
        // separates F1's squared Frobenius residual coefficient from F5's
        // linear scalar.
        var f = Build();
        Assert.Equal(16.0 / 9.0, f.LocalCoefficient, precision: 14);
    }

    [Fact]
    public void CrossSiteCoefficient_IsExactly16()
    {
        // (d²)² = 4² = 16. Squared because the cooperative tensor-assembly
        // |tr(M_l)|² · 4^(N−2) for l ≠ l′ already squares the trace, and
        // |tr(M_l)| ∝ d² at γ=1.
        var f = Build();
        Assert.Equal(16.0, f.CrossSiteCoefficient, precision: 14);
    }

    [Fact]
    public void LiveLocalCoefficient_MatchesParentConstant()
    {
        // Drift guard: the Pi2-derived live value must agree bit-exact with the
        // parent closed-form constant on F1DepolResidualClosedForm. If either
        // side moves, this test surfaces it.
        var f = Build();
        Assert.Equal(F1DepolResidualClosedForm.LocalCoefficient, f.LiveLocalCoefficient, precision: 14);
    }

    [Fact]
    public void LiveCrossSiteCoefficient_MatchesParentConstant()
    {
        var f = Build();
        Assert.Equal(F1DepolResidualClosedForm.CrossSiteCoefficient, f.LiveCrossSiteCoefficient, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullLadder()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F1DepolResidualClosedFormPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
    }

    [Fact]
    public void Constructor_RejectsNullMirror()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F1DepolResidualClosedFormPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_AllRequiredSources()
    {
        var f = Build();
        Assert.Contains("PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim", f.Anchor);
        Assert.Contains("F5DepolarizingErrorPi2Inheritance", f.Anchor);
    }

    [Fact]
    public void ExtraChildren_ContainAllInspectionNodes()
    {
        var f = Build();
        // Walk Children (ExtraChildren are folded into the public Children
        // sequence by the Claim base class) and collect display names.
        var names = f.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(names, n => n.Contains("DSquared (= a_{-1} = d²)"));
        Assert.Contains(names, n => n.Contains("DSquaredMinusOne"));
        Assert.Contains(names, n => n.Contains("PerPauliDepolarizingRate"));
        Assert.Contains(names, n => n.Contains("LocalCoefficient"));
        Assert.Contains(names, n => n.Contains("CrossSiteCoefficient"));
    }

    [Fact]
    public void Reconnaissance_EmitsDecomposition()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F1 depol-residual closed form: ‖M(depol)‖² = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²]");
        _out.WriteLine("    (Tier 1 derived; H-independent; γ_Z-independent; per-site only)");
        _out.WriteLine("");
        _out.WriteLine($"    d²                           = a_(-1)         = {f.DSquared}");
        _out.WriteLine($"    d² − 1                       = a_(-1) − 1     = {f.DSquaredMinusOne}");
        _out.WriteLine($"    per-Pauli rate (= d²/(d²−1)) = 4/3            = {f.PerPauliDepolarizingRate:F6}");
        _out.WriteLine($"    local  (= (d²/(d²−1))² )     = 16/9           = {f.LocalCoefficient:F6}");
        _out.WriteLine($"    cross  (= (d²)²)             = 16             = {f.CrossSiteCoefficient:F6}");
        _out.WriteLine("");
        _out.WriteLine("    F5 sibling: linear scalar uses d/(d²−1) = 2/3");
        _out.WriteLine("    F1 depol:   squared residual uses (d²/(d²−1))² = 16/9 (squared, one d higher)");
        _out.WriteLine("    Shared anchor: DSquaredMinusOne = 3 from Pi2-Foundation primitive");
    }
}
