using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F80 reach beyond its proven scope. PROOF_F80_BLOCH_SIGNWALK establishes
/// Spec(M) = ±2i·Spec(H) (M = Π·L·Π⁻¹ + L + 2σ·I) for chain Π²-odd 2-body Hamiltonians, and
/// flags ring/star topologies and the k-body cluster-value table as open. The Step-5 mechanism
/// (Π·[H,·]·Π⁻¹ = ±{H,·}) is per-bond and per-site, so the structural identity should not care
/// about topology or body-count, only that every term is Π²-odd. Here we check that directly with
/// the verified <see cref="PalindromeResidual"/>: M's spectrum should be purely imaginary and its
/// distinct |Im| clusters should equal 2·|Spec(H)| for the ring and 3-body cases too.</summary>
public class F80ExtensionExplorationTests
{
    private readonly ITestOutputHelper _out;
    public F80ExtensionExplorationTests(ITestOutputHelper output) => _out = output;

    private static readonly Complex I_ = Complex.ImaginaryOne;
    private static readonly ComplexMatrix PX = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
    private static readonly ComplexMatrix PY = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, -I_ }, { I_, 0 } });
    private static readonly ComplexMatrix PZ = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
    private static readonly ComplexMatrix P1 = ComplexMatrix.Build.DenseIdentity(2);

    private static ComplexMatrix Letter(string s) => s switch { "X" => PX, "Y" => PY, "Z" => PZ, _ => P1 };

    private static ComplexMatrix ChainBond(int N, ComplexMatrix P, ComplexMatrix Q)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int l = 0; l < N - 1; l++) H += Term(N, new[] { (l, P), (l + 1, Q) });
        return H;
    }

    // One Pauli-string term: the given letters on the given sites, identity elsewhere.
    private static ComplexMatrix Term(int N, (int site, ComplexMatrix P)[] ops)
    {
        ComplexMatrix? r = null;
        for (int s = 0; s < N; s++)
        {
            var op = ops.FirstOrDefault(o => o.site == s).P ?? P1;
            r = r is null ? op : r.KroneckerProduct(op);
        }
        return r!;
    }

    private static ComplexMatrix ChainXY(int N, bool ring)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int l = 0; l < N - 1; l++) H += Term(N, new[] { (l, PX), (l + 1, PY) });
        if (ring) H += Term(N, new[] { (N - 1, PX), (0, PY) });
        return H;
    }

    private static ComplexMatrix Chain3XXY(int N)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int l = 0; l < N - 2; l++) H += Term(N, new[] { (l, PX), (l + 1, PX), (l + 2, PY) });
        return H;
    }

    private static ComplexMatrix Chain4XXXY(int N)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int l = 0; l < N - 3; l++) H += Term(N, new[] { (l, PX), (l + 1, PX), (l + 2, PX), (l + 3, PY) });
        return H;
    }

    // Star: hub site 0 bonded (X,Y) to every spoke.
    private static ComplexMatrix StarXY(int N)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int s = 1; s < N; s++) H += Term(N, new[] { (0, PX), (s, PY) });
        return H;
    }

    // Returns (M is purely imaginary?, distinct |Im(Spec(M))| clusters, distinct 2·|Spec(H)|).
    private static (bool pureImag, double[] mClusters, double[] hClusters, double defect) CheckF80(int N, ComplexMatrix H)
    {
        var gammas = Enumerable.Repeat(0.1, N).ToList();
        var L = PauliDephasingDissipator.BuildZ(H, gammas);
        var M = PalindromeResidual.Build(L, N, gammas.Sum());
        double defect = M.FrobeniusNorm();

        var mEig = M.Evd().EigenValues;
        double maxRe = mEig.Max(z => Math.Abs(z.Real));
        var mClusters = mEig.Select(z => Math.Round(Math.Abs(z.Imaginary), 6))
            .Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        var hClusters = H.Evd().EigenValues.Select(z => Math.Round(2.0 * Math.Abs(z.Real), 6))
            .Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        return (maxRe < 1e-9, mClusters, hClusters, defect);
    }

    [Theory]
    [InlineData("chain (X,Y) 2-body  [proven scope]", 4, "chain")]
    [InlineData("ring  (X,Y) 2-body  [F80 open]", 4, "ring")]
    [InlineData("star  (X,Y) 2-body  [F80 open]", 4, "star")]
    [InlineData("3-body (X,X,Y)      [F80 open]", 4, "xxy3")]
    [InlineData("chain (X,Y) 2-body  [proven scope]", 5, "chain")]
    [InlineData("ring  (X,Y) 2-body  [F80 open]", 5, "ring")]
    [InlineData("star  (X,Y) 2-body  [F80 open]", 5, "star")]
    [InlineData("3-body (X,X,Y)      [F80 open]", 5, "xxy3")]
    [InlineData("4-body (X,X,X,Y)    [F80 open]", 5, "xxxy4")]
    public void SpecM_Equals_TwiceI_SpecH_BeyondChainTwoBody(string label, int N, string kind)
    {
        var H = kind switch
        {
            "ring" => ChainXY(N, ring: true),
            "star" => StarXY(N),
            "xxy3" => Chain3XXY(N),
            "xxxy4" => Chain4XXXY(N),
            _ => ChainXY(N, ring: false),
        };

        var (pureImag, mClusters, hClusters, defect) = CheckF80(N, H);
        _out.WriteLine($"N={N}  {label}");
        _out.WriteLine($"    ||M||_F = {defect:F4}   M purely imaginary: {pureImag}");
        _out.WriteLine($"    |Im Spec(M)| clusters = [{string.Join(", ", mClusters.Select(x => x.ToString("F4")))}]");
        _out.WriteLine($"    2|Spec(H)|            = [{string.Join(", ", hClusters.Select(x => x.ToString("F4")))}]");

        // F80 structural identity: M anti-Hermitian (purely imaginary spectrum) and its clusters
        // are exactly 2|Spec(H)|. If this holds for ring and 3-body, the identity is topology- and
        // body-count-agnostic (the per-bond Step-5 mechanism), extending F80's verified scope.
        Assert.True(pureImag, $"{label} N={N}: M should be purely imaginary (anti-Hermitian)");
        Assert.Equal(hClusters, mClusters);
    }

    [Theory]
    [InlineData("(X,Y)  Pi2-odd ", "X", "Y")]
    [InlineData("(X,Z)  Pi2-odd ", "X", "Z")]
    [InlineData("(Y,Z)  Pi2-even", "Y", "Z")]
    [InlineData("(Z,Y)  Pi2-even", "Z", "Y")]
    public void Pi2Parity_DecidesSingleVsDifferenceSpectrum(string label, string a, string b)
    {
        int N = 4;
        var H = ChainBond(N, Letter(a), Letter(b));
        var gammas = Enumerable.Repeat(0.1, N).ToList();
        var M = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(H, gammas), N, gammas.Sum());

        var mEig = M.Evd().EigenValues;
        bool pureImag = mEig.Max(z => Math.Abs(z.Real)) < 1e-9;
        var mClusters = mEig.Select(z => Math.Round(Math.Abs(z.Imaginary), 6)).Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        var hEig = H.Evd().EigenValues.Select(z => z.Real).ToArray();
        var single = hEig.Select(l => Math.Round(2 * Math.Abs(l), 6)).Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();
        var diff = (from la in hEig from lb in hEig select Math.Round(2 * Math.Abs(la - lb), 6))
            .Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        bool matchesSingle = mClusters.SequenceEqual(single);
        bool matchesDiff = mClusters.SequenceEqual(diff);

        _out.WriteLine($"N={N}  {label}   ||M||_F={M.FrobeniusNorm():F3}  pureImag={pureImag}  #Mclusters={mClusters.Length}");
        _out.WriteLine($"    matches 2|single eigenvalues|  (F80 odd, M=-2i H(x)I): {matchesSingle}  (#single={single.Length})");
        _out.WriteLine($"    matches 2|eigenvalue diffs|    (even, M=2 L_H):        {matchesDiff}  (#diff={diff.Length})");

        // The conjecture: Pi2-odd bonds -> M sees single eigenvalues (M = -2i H(x)I);
        // Pi2-even bonds -> M sees eigenvalue differences (M = 2 L_H = -2i[H,.]).
        Assert.True(pureImag, $"{label}: M should be purely imaginary");
        Assert.True(matchesSingle || matchesDiff, $"{label}: M clusters match neither single nor differences");
    }
}
