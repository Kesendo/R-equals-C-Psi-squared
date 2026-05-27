using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CommutatorDConjugationSignTests
{
    private const double Tol = 1e-10;

    private static CommutatorDConjugationSign MakeClaim() =>
        new CommutatorDConjugationSign(new Pi2KleinV4DephaseSwapGroup());

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
    public void Anchor_ReferencesF114Sources()
    {
        var claim = MakeClaim();
        Assert.Contains("F114", claim.Anchor);
        Assert.Contains("_m_level_sign_functional_explore.py", claim.Anchor);
        Assert.Contains("Pi2KleinV4DephaseSwapGroup", claim.Anchor);
    }

    [Fact]
    public void KleinV4_ParentIsWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.KleinV4);
        Assert.Equal(Tier.Tier1Derived, claim.KleinV4.Tier);
    }

    [Fact]
    public void Constructor_RejectsNullKleinV4()
    {
        Assert.Throws<ArgumentNullException>(
            () => new CommutatorDConjugationSign(null!));
    }

    // ------------------------------------------------------------------
    // Single-term ε(σ) at N = 1, 2, 3: closed-form check
    // ------------------------------------------------------------------

    [Theory]
    [InlineData(PauliLetter.I, DConjugationSign.Zero)]
    [InlineData(PauliLetter.X, DConjugationSign.Minus)]   // n_Y = 0 even
    [InlineData(PauliLetter.Z, DConjugationSign.Minus)]   // n_Y = 0 even
    [InlineData(PauliLetter.Y, DConjugationSign.Plus)]    // n_Y = 1 odd
    public void Compute_SingleSite_MatchesClosedForm(PauliLetter letter, DConjugationSign expected)
    {
        var term = new PauliTerm(new[] { letter }, Complex.One);
        Assert.Equal(expected, CommutatorDConjugationSign.Compute(term));
    }

    [Theory]
    // N=2: bit-exact match against the empirical Welle 15 Task A polish observations.
    [InlineData(new[] { PauliLetter.X, PauliLetter.Z }, DConjugationSign.Minus)] // n_Y = 0
    [InlineData(new[] { PauliLetter.Z, PauliLetter.X }, DConjugationSign.Minus)] // n_Y = 0
    [InlineData(new[] { PauliLetter.Y, PauliLetter.Z }, DConjugationSign.Plus)]  // n_Y = 1
    [InlineData(new[] { PauliLetter.Z, PauliLetter.Y }, DConjugationSign.Plus)]  // n_Y = 1
    [InlineData(new[] { PauliLetter.X, PauliLetter.X }, DConjugationSign.Minus)] // n_Y = 0
    [InlineData(new[] { PauliLetter.Y, PauliLetter.Y }, DConjugationSign.Minus)] // n_Y = 2 even
    [InlineData(new[] { PauliLetter.Z, PauliLetter.Z }, DConjugationSign.Minus)] // n_Y = 0
    [InlineData(new[] { PauliLetter.X, PauliLetter.Y }, DConjugationSign.Plus)]  // n_Y = 1
    [InlineData(new[] { PauliLetter.I, PauliLetter.I }, DConjugationSign.Zero)]
    // N=3: multi-site cases.
    [InlineData(new[] { PauliLetter.Y, PauliLetter.Y, PauliLetter.Y }, DConjugationSign.Plus)]  // n_Y = 3
    [InlineData(new[] { PauliLetter.Z, PauliLetter.Z, PauliLetter.Z }, DConjugationSign.Minus)] // n_Y = 0
    [InlineData(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, DConjugationSign.Plus)]  // n_Y = 1
    public void Compute_MultiSite_MatchesClosedForm(PauliLetter[] letters, DConjugationSign expected)
    {
        var term = new PauliTerm(letters, Complex.One);
        Assert.Equal(expected, CommutatorDConjugationSign.Compute(term));
    }

    [Fact]
    public void Compute_Term_RejectsNull()
    {
        Assert.Throws<ArgumentNullException>(
            () => CommutatorDConjugationSign.Compute((PauliTerm)null!));
    }

    // ------------------------------------------------------------------
    // Linear-combination ε(H): well-defined + Mixed cases
    // ------------------------------------------------------------------

    [Fact]
    public void Compute_LinearCombination_XZ_Plus_ZX_IsMinus()
    {
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.X }, Complex.One),
        };
        Assert.Equal(DConjugationSign.Minus, CommutatorDConjugationSign.Compute(terms));
    }

    [Fact]
    public void Compute_LinearCombination_YZ_Plus_ZY_IsPlus()
    {
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.Y, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Y }, Complex.One),
        };
        Assert.Equal(DConjugationSign.Plus, CommutatorDConjugationSign.Compute(terms));
    }

    [Fact]
    public void Compute_Heisenberg_XX_YY_ZZ_IsMinus()
    {
        // XX (n_Y=0), YY (n_Y=2), ZZ (n_Y=0): all even-parity n_Y → all ε = -1 per term.
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.X }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Y, PauliLetter.Y }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.Z }, Complex.One),
        };
        Assert.Equal(DConjugationSign.Minus, CommutatorDConjugationSign.Compute(terms));
    }

    [Fact]
    public void Compute_MixedParity_IsMixed()
    {
        // X (n_Y=0, ε=-1) + Y (n_Y=1, ε=+1): mixed parity → no single ε(H).
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.I }, Complex.One),
            new PauliTerm(new[] { PauliLetter.I, PauliLetter.Y }, Complex.One),
        };
        Assert.Equal(DConjugationSign.Mixed, CommutatorDConjugationSign.Compute(terms));
    }

    [Fact]
    public void Compute_EmptyList_IsZero()
    {
        Assert.Equal(DConjugationSign.Zero,
            CommutatorDConjugationSign.Compute(Array.Empty<PauliTerm>()));
    }

    [Fact]
    public void Compute_AllIdentityTerms_IsZero()
    {
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.I, PauliLetter.I }, Complex.One),
            new PauliTerm(new[] { PauliLetter.I, PauliLetter.I }, new Complex(2.0, 0)),
        };
        Assert.Equal(DConjugationSign.Zero, CommutatorDConjugationSign.Compute(terms));
    }

    [Fact]
    public void Compute_LinearCombination_RejectsNull()
    {
        Assert.Throws<ArgumentNullException>(
            () => CommutatorDConjugationSign.Compute((IReadOnlyList<PauliTerm>)null!));
    }

    // ------------------------------------------------------------------
    // IsWellDefined
    // ------------------------------------------------------------------

    [Fact]
    public void IsWellDefined_SameParity_IsTrue()
    {
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.Z, PauliLetter.X }, Complex.One),
        };
        Assert.True(CommutatorDConjugationSign.IsWellDefined(terms));
    }

    [Fact]
    public void IsWellDefined_MixedParity_IsFalse()
    {
        var terms = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.I }, Complex.One),
            new PauliTerm(new[] { PauliLetter.I, PauliLetter.Y }, Complex.One),
        };
        Assert.False(CommutatorDConjugationSign.IsWellDefined(terms));
    }

    [Fact]
    public void IsWellDefined_ZeroCase_IsTrue()
    {
        Assert.True(CommutatorDConjugationSign.IsWellDefined(Array.Empty<PauliTerm>()));
    }

    // ------------------------------------------------------------------
    // Matrix-level verification: D · L_σ · D = ε · L_σ bit-exact at N = 2
    // ------------------------------------------------------------------

    [Theory]
    [InlineData(new[] { PauliLetter.X, PauliLetter.Z })] // n_Y=0 even → ε=-1
    [InlineData(new[] { PauliLetter.Y, PauliLetter.Z })] // n_Y=1 odd → ε=+1
    [InlineData(new[] { PauliLetter.Y, PauliLetter.Y })] // n_Y=2 even → ε=-1
    [InlineData(new[] { PauliLetter.X, PauliLetter.Y })] // n_Y=1 odd → ε=+1
    public void MatrixLevel_DLsigmaD_EqualsEpsilonLsigma_AtN2(PauliLetter[] letters)
    {
        const int N = 2;
        var term = new PauliTerm(letters, Complex.One);
        var expectedSign = CommutatorDConjugationSign.Compute(term);
        Assert.True(
            expectedSign == DConjugationSign.Plus || expectedSign == DConjugationSign.Minus,
            "test design: expected non-trivial ε on these Pauli strings");

        var L = BuildLsigmaInPauliBasis(letters);
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N);

        var DLD = D * L * D;
        var expected = expectedSign == DConjugationSign.Plus ? L : -L;
        var residual = (DLD - expected).FrobeniusNorm();

        Assert.True(L.FrobeniusNorm() > 1e-6, "test design: L_σ should be non-trivial");
        Assert.True(residual < Tol,
            $"F114 violation: ‖D·L·D − ε·L‖_F = {residual:E3} (ε = {expectedSign}, tol {Tol:E1})");
    }

    /// <summary>Build L_σ = −i[σ, ·] in the 4^N Pauli basis directly via
    /// [L_σ]_(α,β) = Tr(σ_α · [σ, σ_β]) / 2^N, matching the script
    /// simulations/_m_level_sign_functional_explore.py construction.</summary>
    private static ComplexMatrix BuildLsigmaInPauliBasis(PauliLetter[] sigmaLetters)
    {
        int N = sigmaLetters.Length;
        int d = 1 << N;
        long d2 = 1L << (2 * N);
        var sigma = PauliString.Build(sigmaLetters);
        var L = Matrix<Complex>.Build.Dense((int)d2, (int)d2);
        for (long j = 0; j < d2; j++)
        {
            var betaLetters = PauliIndex.FromFlat(j, N);
            var sigmaBeta = PauliString.Build(betaLetters);
            var commutator = -Complex.ImaginaryOne * (sigma * sigmaBeta - sigmaBeta * sigma);
            for (long i = 0; i < d2; i++)
            {
                var alphaLetters = PauliIndex.FromFlat(i, N);
                var sigmaAlpha = PauliString.Build(alphaLetters);
                L[(int)i, (int)j] = (sigmaAlpha * commutator).Trace() / d;
            }
        }
        return L;
    }
}
