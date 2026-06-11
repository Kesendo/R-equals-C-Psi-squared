using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class AntilinearTriangleClaimTests
{
    private static AntilinearTriangleClaim MakeClaim()
    {
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        var f114 = new CommutatorDConjugationSign(kleinV4);
        var mirror = new MirrorGroupD4Claim(
            new KleinEightCellClaim(new KleinFourCellClaim()),
            f114,
            kleinV4);
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f112 = new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));
        return new AntilinearTriangleClaim(mirror, f114, f112);
    }

    // ------------------------------------------------------------------
    // Claim metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = MakeClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Anchor_ReferencesProofAndVerifier()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_ANTILINEAR_TRIANGLE.md", claim.Anchor);
        Assert.Contains("antilinear_triangle.py", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.MirrorGroup);
        Assert.NotNull(claim.F114);
        Assert.NotNull(claim.F112);
        // The θ-leg parent is the same F114 instance the mirror group carries.
        Assert.Same(claim.F114, claim.MirrorGroup.F114);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        var f114 = new CommutatorDConjugationSign(kleinV4);
        var mirror = new MirrorGroupD4Claim(
            new KleinEightCellClaim(new KleinFourCellClaim()), f114, kleinV4);
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f112 = new LindbladBitBPiBalance(
            new F108Part1Pi2EvenAlwaysPalindromic(part2),
            new LindbladBitAPiBalance(part2));

        Assert.Throws<ArgumentNullException>(() => new AntilinearTriangleClaim(null!, f114, f112));
        Assert.Throws<ArgumentNullException>(() => new AntilinearTriangleClaim(mirror, null!, f112));
        Assert.Throws<ArgumentNullException>(() => new AntilinearTriangleClaim(mirror, f114, null!));
    }

    [Fact]
    public void NotAnIZ2AxisClaim_CubeMapCountsUnchanged()
    {
        // Pinned invariant: the Cubic3Claims count in KleinEightCellClaimRegistrationTests
        // stays 1; this claim is cross-axis structural like MirrorGroupD4Claim.
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(AntilinearTriangleClaim)));
    }

    // ------------------------------------------------------------------
    // Self-check battery (N = 2 dense operators, built in the ctor)
    // ------------------------------------------------------------------

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(10, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Battery_TriangleClosure_IsExact()
    {
        var claim = MakeClaim();
        var triangle = claim.Cases.Single(c => c.Name.Contains("† = θ∘conj"));
        Assert.True(triangle.Passes, $"triangle closure dev out of tolerance: {triangle.Actual}");
    }

    [Fact]
    public void Battery_PauliAction_AllSixteenStrings()
    {
        var claim = MakeClaim();
        var pauli = claim.Cases.Single(c => c.Name.StartsWith("Pauli action", StringComparison.Ordinal));
        Assert.Equal("16/16", pauli.Actual);
    }

    [Fact]
    public void Battery_TransportSigns_AllFourMaps()
    {
        var claim = MakeClaim();
        foreach (var rhs in new[] { "+L_H", "−L_{Hᵀ}", "−L_{H̄}", "+L_{H†}" })
        {
            var c = claim.Cases.Single(x =>
                x.Name.StartsWith("transport ", StringComparison.Ordinal) && x.Name.EndsWith(rhs, StringComparison.Ordinal));
            Assert.True(c.Passes, $"transport case '{c.Name}' failed: {c.Actual}");
        }
    }

    [Fact]
    public void Battery_WordReversal_CoversBothSigns()
    {
        var claim = MakeClaim();
        var words = claim.Cases.Single(c => c.Name.StartsWith("word reversal", StringComparison.Ordinal));
        Assert.Equal("6/6", words.Actual);
    }

    [Fact]
    public void Battery_AntilinearDouble_IsOrderSixteen()
    {
        var claim = MakeClaim();
        var dbl = claim.Cases.Single(c => c.Name.Contains("antilinear double"));
        Assert.Equal("16 elements, 8 antilinear", dbl.Expected);
        Assert.Equal("16 elements, 8 antilinear", dbl.Actual);
    }

    [Fact]
    public void Summary_CarriesTheEngineAndTheDouble()
    {
        var claim = MakeClaim();
        Assert.Contains("ℓ(μ)·m(μ)", claim.Summary);
        Assert.Contains("antilinear double", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }

    // ------------------------------------------------------------------
    // Direct mathematical spot-check, independent of the claim's battery:
    // the transport-law sign table on a fixed 1-qubit H and ρ.
    // ------------------------------------------------------------------

    [Fact]
    public void TransportLaw_SignTable_DirectCheck()
    {
        var h = Matrix<Complex>.Build.Dense(2, 2);
        h[0, 0] = new Complex(0.3, -0.7); h[0, 1] = new Complex(1.1, 0.4);
        h[1, 0] = new Complex(-0.2, 0.9); h[1, 1] = new Complex(-0.5, 0.6);
        var rho = Matrix<Complex>.Build.Dense(2, 2);
        rho[0, 0] = new Complex(0.8, 0.1); rho[0, 1] = new Complex(-0.4, 0.5);
        rho[1, 0] = new Complex(0.2, -0.6); rho[1, 1] = new Complex(-0.9, 0.3);

        var minusI = new Complex(0.0, -1.0);
        ComplexMatrix L(ComplexMatrix g, ComplexMatrix r) => (g * r - r * g).Multiply(minusI);

        // θ: θ∘L_H∘θ = −L_{Hᵀ}
        var theta = L(h, rho.Transpose()).Transpose();
        Assert.True(MaxAbsDiff(theta, L(h.Transpose(), rho).Multiply(-Complex.One)) <= 1e-12);
        // conj: conj∘L_H∘conj = −L_{H̄}
        var conj = L(h, rho.Conjugate()).Conjugate();
        Assert.True(MaxAbsDiff(conj, L(h.Conjugate(), rho).Multiply(-Complex.One)) <= 1e-12);
        // †: †∘L_H∘† = +L_{H†} (the two signs cancel)
        var dag = L(h, rho.ConjugateTranspose()).ConjugateTranspose();
        Assert.True(MaxAbsDiff(dag, L(h.ConjugateTranspose(), rho)) <= 1e-12);

        // The signs are load-bearing: the opposite sign visibly fails on the same data.
        Assert.True(MaxAbsDiff(theta, L(h.Transpose(), rho)) > 1e-6);
        Assert.True(MaxAbsDiff(conj, L(h.Conjugate(), rho)) > 1e-6);
        Assert.True(MaxAbsDiff(dag, L(h.ConjugateTranspose(), rho).Multiply(-Complex.One)) > 1e-6);
    }

    private static double MaxAbsDiff(ComplexMatrix a, ComplexMatrix b)
    {
        double m = 0.0;
        for (int i = 0; i < a.RowCount; i++)
            for (int j = 0; j < a.ColumnCount; j++)
            {
                double v = (a[i, j] - b[i, j]).Magnitude;
                if (v > m) m = v;
            }
        return m;
    }
}
