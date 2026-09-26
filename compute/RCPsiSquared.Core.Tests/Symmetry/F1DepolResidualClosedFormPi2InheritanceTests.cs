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
        // + machine-precision verification N=2..5). The centered coefficients (16/9, 0) reduce to
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
    public void CrossSiteCoefficient_IsExactlyZeroAfterCentering()
    {
        // Computed from the per-site kernel in integers: centered at its spectral mean the kernel
        // is (2, 2, −2, −2)/3 on (I, X, Y, Z), traceless, so the cross-site term vanishes exactly.
        var f = Build();
        var centered = f.ResidualCoefficients(centered: true);
        Assert.True(centered.Cross == 0.0, $"centered cross {centered.Cross:R}");
        Assert.True(centered.Local == 16.0 / 9.0, $"centered local {centered.Local:R}");
        Assert.True(f.CrossSiteCoefficient == 0.0);
    }

    [Fact]
    public void BareResidual_CarriesTheCrossSiteSixteen()
    {
        // The control that makes the zero above a finding about centering: the same computation at
        // σ = 0 (bare kernel (−4, −4, −8, −8)/3, trace −8) returns the bare formula's (16/9, 16).
        var bare = Build().ResidualCoefficients(centered: false);
        Assert.True(bare.Cross == 16.0, $"bare cross {bare.Cross:R}");
        Assert.True(bare.Local == 16.0 / 9.0, $"bare local {bare.Local:R}");
    }

    [Fact]
    public void LiveLocalCoefficient_MatchesParentConstant()
    {
        // Drift guard: the Pi2-derived live value must agree to floating-point tolerance with the
        // parent closed-form constant on F1DepolResidualClosedForm. If either
        // side moves, this test surfaces it.
        var f = Build();
        Assert.Equal(F1DepolResidualClosedForm.LocalCoefficient, f.LiveLocalCoefficient, precision: 14);
    }

    [Fact]
    public void LiveCrossSiteCoefficient_MatchesParentConstant()
    {
        // Exact: the kernel computation and the parent's constant are the same number by two routes.
        var f = Build();
        Assert.True(f.LiveCrossSiteCoefficient == F1DepolResidualClosedForm.CrossSiteCoefficient);
        Assert.True(f.ResidualCoefficients(centered: true).Local == F1DepolResidualClosedForm.LocalCoefficient);
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
        _out.WriteLine("    F1 centered depol-residual: ‖M_F1(depol)‖² = 4^(N−1)·(16/9)·Σγ²");
        _out.WriteLine("    (Tier 1 derived; H-independent; γ_Z-independent; per-site only)");
        _out.WriteLine("");
        _out.WriteLine($"    d²                           = a_(-1)         = {f.DSquared}");
        _out.WriteLine($"    d² − 1                       = a_(-1) − 1     = {f.DSquaredMinusOne}");
        _out.WriteLine($"    per-Pauli rate (= d²/(d²−1)) = 4/3            = {f.PerPauliDepolarizingRate:F6}");
        _out.WriteLine($"    local  (= (d²/(d²−1))² )     = 16/9           = {f.LocalCoefficient:F6}");
        _out.WriteLine($"    cross  (centered trace zero)  = 0              = {f.CrossSiteCoefficient:F6}");
        _out.WriteLine("");
        _out.WriteLine("    F5 sibling: linear scalar uses d/(d²−1) = 2/3");
        _out.WriteLine("    F1 depol:   squared residual uses (d²/(d²−1))² = 16/9 (squared, one d higher)");
        _out.WriteLine("    Shared anchor: DSquaredMinusOne = 3 from Pi2-Foundation primitive");
    }
}
