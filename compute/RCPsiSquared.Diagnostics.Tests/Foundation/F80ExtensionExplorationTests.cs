using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
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

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void StarDispersion_ClosedForm_EvenIntegers(int N)
    {
        // H_star = sum_{spokes s} X_hub Y_s = X_hub (x) (sum_s Y_s). The m = N-1 spoke Y's commute
        // and each has spectrum {+-1}, so sum_s Y_s has spectrum {m - 2j : j=0..m}, and the X_hub
        // factor (spectrum {+-1}) gives Spec(H_star) = {+-(m - 2j)}. By the F80 identity the cluster
        // values are 2|m - 2j| , the even integers 2m, 2m-4, 2m-8, ... (the star's "dispersion").
        var H = StarXY(N);
        var hAbs = H.Evd().EigenValues.Select(z => Math.Round(Math.Abs(z.Real), 6))
            .Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        int m = N - 1;
        var predicted = Enumerable.Range(0, m + 1).Select(j => (double)Math.Abs(m - 2 * j))
            .Where(x => x > 1e-6).Distinct().OrderBy(x => x).ToArray();

        _out.WriteLine($"N={N} star: |Spec(H)| = [{string.Join(", ", hAbs)}]   predicted |m-2j| = [{string.Join(", ", predicted)}]   clusters = 2x");
        Assert.Equal(predicted, hAbs);
    }

    [Fact]
    public void FourFoldKleinFingerprint_Neither_OneSide_OtherSide_Both()
    {
        // Make the 4-fold Klein V4 structure of the mirror-defect M visible as "which side of rho
        // does M touch". M acting as H(x)I changes only the left operator-index of rho; I(x)H^T only
        // the right; the commutator [H,.] both. The two Klein cells are exactly (touches-left?) x
        // (touches-right?). The D-twist of the vec basis is diagonal, so "which index changes" is
        // invariant (only the absolute left=ket / right=bra naming carries a vec_F transpose, which
        // is why we report "left/right index" rather than hard-claiming ket vs bra).
        int N = 3, d = 1 << N, d2 = d * d;
        var gammas = Enumerable.Repeat(0.1, N).ToList();
        double sigma = gammas.Sum();
        var T = PauliBasis.VecToPauliBasisTransform(N);
        double inv2N = 1.0 / (1 << N);

        (double mNorm, double leftCh, double rightCh) Side(ComplexMatrix H)
        {
            var mPauli = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(H, gammas), N, sigma);
            var mVec = (T * mPauli * T.ConjugateTranspose()).Multiply((Complex)inv2N);  // back to vec basis
            double total = 0, left = 0, right = 0;
            for (int row = 0; row < d2; row++)
                for (int col = 0; col < d2; col++)
                {
                    var z = mVec[row, col];
                    double w = z.Real * z.Real + z.Imaginary * z.Imaginary;
                    if (w < 1e-18) continue;
                    total += w;
                    if (row / d != col / d) left += w;   // left operator-index changed
                    if (row % d != col % d) right += w;  // right operator-index changed
                }
            return (Math.Sqrt(total), left, right);
        }

        var mother = ChainBond(N, PX, PX) + ChainBond(N, PY, PY);  // XX+YY  truly      -> Pp
        var fatherA = ChainBond(N, PX, PY);                         // XY     Pi2-odd a  -> Mp
        var fatherB = ChainBond(N, PX, PZ);                         // XZ     Pi2-odd b  -> Mm
        var child = ChainBond(N, PY, PZ);                           // YZ     Pi2-even   -> Pm

        var (mN, mL, mR) = Side(mother);
        var (aN, aL, aR) = Side(fatherA);
        var (bN, bL, bR) = Side(fatherB);
        var (cN, cL, cR) = Side(child);

        _out.WriteLine($"Mother  XX+YY (Pp): ||M||={mN:F3}  left={mL:F3}  right={mR:F3}  -> touches NEITHER side (silence)");
        _out.WriteLine($"FatherA XY    (Mp): ||M||={aN:F3}  left={aL:F3}  right={aR:F3}  -> ONE side only");
        _out.WriteLine($"FatherB XZ    (Mm): ||M||={bN:F3}  left={bL:F3}  right={bR:F3}  -> the OPPOSITE side only");
        _out.WriteLine($"Child   YZ    (Pm): ||M||={cN:F3}  left={cL:F3}  right={cR:F3}  -> BOTH sides ([H,.])");

        // The 4 Klein cells = (touches-left?) x (touches-right?): neither / one / other / both.
        Assert.True(mN < 1e-9, "truly -> M=0 (Mother, touches neither)");
        Assert.True(aL < 1e-9 && aR > 1e-6, "XY -> one side only (Father_a, Mp)");
        Assert.True(bR < 1e-9 && bL > 1e-6, "XZ -> the opposite side only (Father_b, Mm)");
        Assert.True(cL > 1e-6 && cR > 1e-6, "YZ -> both sides (Child, Pm)");
        // The two fathers sit on opposite sides , that opposition is the bit_a / Klein distinction
        // the 3-role Trinity collapsed.
        Assert.True(aR > 1e-6 && bL > 1e-6 && aL < 1e-9 && bR < 1e-9, "the two fathers are on opposite sides");
    }

    [Fact]
    public void DefectLadder_NormSquared_IsSideTimesBondsTimes4PowerN()
    {
        // The defect magnitude is a ladder: ||M||^2 = side * #bonds * 4^N, with side = 4 (Pi2-odd,
        // one-sided, M = -2i H(x)I) or 8 (Pi2-even, two-sided, M = -2i[H,.]). From the F80 relation
        // ||M||^2 = 4||H||^2 2^N with ||H||^2 = #bonds * 2^N for orthogonal Pauli bonds. Each bond
        // adds a fixed quantum (4*4^N odd / 8*4^N even); the child is sqrt(2)x the father per bond.
        (string label, ComplexMatrix H, int bonds, int side, int N)[] cases =
        {
            ("N=3 chain XY (odd, 1-side)",  ChainBond(3, PX, PY),    2, 4, 3),  // sqrt(512)
            ("N=3 chain YZ (even, 2-side)", ChainBond(3, PY, PZ),    2, 8, 3),  // sqrt(1024)
            ("N=4 ring  XY (odd, 1-side)",  ChainXY(4, ring: true),  4, 4, 4),  // sqrt(4096)
            ("N=5 chain XY (odd, 1-side)",  ChainBond(5, PX, PY),    4, 4, 5),  // sqrt(16384)
        };
        foreach (var (label, H, bonds, side, N) in cases)
        {
            var gammas = Enumerable.Repeat(0.1, N).ToList();
            var M = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(H, gammas), N, gammas.Sum());
            double normSq = M.FrobeniusNorm();
            normSq *= normSq;
            double predicted = (double)side * bonds * Math.Pow(4, N);
            _out.WriteLine($"{label,-28}: ||M||^2 = {normSq,7:F0} = {side}*{bonds}*4^{N}   ||M|| = sqrt({normSq:F0}) = {Math.Sqrt(normSq):F3}");
            Assert.Equal(predicted, normSq, precision: 2);
        }
    }

    [Fact]
    public void MixedLetter_M_IsLinear_SoItIsTheSumOfPerParityPieces()
    {
        // A mixed-letter chain: one Pi2-odd family (X,Y) plus one Pi2-even family (Y,Z).
        int N = 4;
        var Hxy = ChainBond(N, PX, PY);  // Pi2-odd  -> M sees single eigenvalues
        var Hyz = ChainBond(N, PY, PZ);  // Pi2-even -> M sees eigenvalue differences
        var gammas = Enumerable.Repeat(0.1, N).ToList();
        double sigma = gammas.Sum();

        var Mmixed = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(Hxy + Hyz, gammas), N, sigma);
        var Mxy = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(Hxy, gammas), N, sigma);
        var Myz = PalindromeResidual.Build(PauliDephasingDissipator.BuildZ(Hyz, gammas), N, sigma);

        // M is linear in H: the dissipator + 2sigma terms cancel (F1), so M(Hxy+Hyz) = M(Hxy)+M(Hyz)
        // bit-exact. A mixed-letter Hamiltonian therefore needs no new sign-walk , its mirror-defect
        // is the sum of the per-bond pieces, each fixed by that bond's Pi2-parity.
        double residual = (Mmixed - (Mxy + Myz)).FrobeniusNorm();
        _out.WriteLine($"||M(Hxy+Hyz) - (M(Hxy)+M(Hyz))||_F = {residual:E3}  (linearity of M in H)");
        Assert.True(residual < 1e-9, $"M should be linear in H; got residual {residual:E3}");

        // The mixed spectrum is still purely imaginary (odd part single energies + even part
        // differences, summed as one operator) but richer than either piece alone.
        var mEig = Mmixed.Evd().EigenValues;
        bool pureImag = mEig.Max(z => Math.Abs(z.Real)) < 1e-9;
        int mixedClusters = mEig.Select(z => Math.Round(Math.Abs(z.Imaginary), 6)).Where(x => x > 1e-6).Distinct().Count();
        int xyClusters = Mxy.Evd().EigenValues.Select(z => Math.Round(Math.Abs(z.Imaginary), 6)).Where(x => x > 1e-6).Distinct().Count();
        _out.WriteLine($"mixed M purely imaginary: {pureImag}   #clusters mixed={mixedClusters} vs odd-only={xyClusters}");
        Assert.True(pureImag, "mixed M should still be purely imaginary");
    }
}
