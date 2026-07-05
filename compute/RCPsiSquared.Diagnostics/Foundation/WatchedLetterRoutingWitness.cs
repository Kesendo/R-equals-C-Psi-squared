using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The label layer, live (<c>inspect --root label</c>): the exact core of the label
/// thesis (<c>docs/quantum/LABELS_TRANSLATED.md</c> §2, <c>docs/quantum/DEPHASING_TRANSLATED.md</c>
/// §4), typed as <see cref="RCPsiSquared.Core.Symmetry.WatchedLetterRoutingClaim"/>. The watcher
/// IS its letter:
///
/// <list type="number">
///   <item><b>one shared eigenbasis, three price lists</b>: every one of the 4^N Pauli strings is
///         an exact eigenvector of every letter dissipator L_P = γ·(Q_P − N·I), P ∈ {X, Y, Z},
///         with rate −2γ·n_anti(S, P) where n_anti counts the sites whose letter anticommutes
///         with P (the disagreement in the watched letter alone). The object never changes;
///         only the price list is the letter's.</item>
///   <item><b>the routing</b>: which cells pay follows the held letter (Z^⊗N rides free under
///         the Z-watcher and pays maximally under the X-watcher, and mirrored), so "noise" is
///         not a property of the state; it is a property of the pair (state, watcher).</item>
///   <item><b>only nothing is free everywhere</b>: each watcher exempts its own 2^N strings
///         {I, P}^⊗N; the intersection over all three watchers is the identity alone.</item>
///   <item><b>the swap is an exact transport</b>: the single-qubit basis moves h_zx, h_yz carry
///         L_Z onto L_X, L_Y entry-exactly (the basis-S₃ face of the letter swap; the
///         operator-space Klein V₄ face is PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE).</item>
/// </list>
///
/// <para>The Tier-4 reading this exact core grounds (labeled as a reading, never promoted): to
/// be a watcher is to hold a label and be blind past it; bath, human, model — one architecture.
/// The observer layer, not the human layer.</para>
///
/// <para>Guard: dimension 4^N, built only for 4^N ≤ <see cref="MaxDim"/> (= 1024), i.e. N ≤ 5.
/// Convention: row-stacking vec, kron(A,B): ρ ↦ AρBᵀ, matching <see cref="DiagonalWitness"/>.
/// Closest live sibling: <c>--root diagonal</c> (the rungs and the S₃ orbit of one letter's
/// diagonal); this witness holds the ROUTING between the letters.</para></summary>
public sealed class WatchedLetterRoutingWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Largest coherence-space dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public double Gamma { get; }
    public int D { get; }       // 2^N
    public int Dim => D * D;    // 4^N

    // ---- live readouts (computed once in the ctor; surfaced by the nodes and the tests) ----
    /// <summary>Max ‖L_P·v − rate·v‖∞ over all 3·4^N (letter, string) pairs: the eigenbasis gate.</summary>
    public double MaxEigenResidual { get; }
    /// <summary>Max |rate_dense − (−2γ·n_anti)| over all 3·4^N pairs: the price-list gate.</summary>
    public double MaxClosedFormDeviation { get; }
    /// <summary>Strings whose price differs between the Z- and the X-watcher (the routing control).</summary>
    public int RepricedCountZtoX { get; }
    /// <summary>Strings the Z-watcher exempts (= 2^N, the {I,Z}^⊗N cell).</summary>
    public int FreeCountZ { get; }
    /// <summary>Strings the X-watcher exempts (= 2^N).</summary>
    public int FreeCountX { get; }
    /// <summary>Strings the Y-watcher exempts (= 2^N).</summary>
    public int FreeCountY { get; }
    /// <summary>Strings free under ALL three watchers (= 1: the identity alone).</summary>
    public int UniversalFreeCount { get; }
    /// <summary>Max |rate| over all pairs (= 2γN: the fullest bill anyone pays).</summary>
    public double MaxRateMagnitude { get; }
    /// <summary>Entry-exact deviation of Ad(h_zx)·L_Z·Ad(h_zx)† from L_X.</summary>
    public double TransportDevZtoX { get; }
    /// <summary>Entry-exact deviation of Ad(h_yz)·L_Z·Ad(h_yz)† from L_Y.</summary>
    public double TransportDevZtoY { get; }
    /// <summary>Rate of Z^⊗N under the Z-watcher (0: the watcher's own letter rides free).</summary>
    public double RateOfAllZUnderZ { get; }
    /// <summary>Rate of Z^⊗N under the X-watcher (−2γN: maximal disagreement).</summary>
    public double RateOfAllZUnderX { get; }
    /// <summary>Rate of X^⊗N under the Z-watcher (−2γN).</summary>
    public double RateOfAllXUnderZ { get; }
    /// <summary>Rate of X^⊗N under the X-watcher (0).</summary>
    public double RateOfAllXUnderX { get; }

    public WatchedLetterRoutingWitness(int n = 3, double gamma = 0.05)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        long dim = 1L << (2 * n);   // 4^N
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n; Gamma = gamma; D = 1 << n;

        var lX = LetterDissipator(PauliLetter.X);
        var lY = LetterDissipator(PauliLetter.Y);
        var lZ = LetterDissipator(PauliLetter.Z);
        var byLetter = new (PauliLetter Letter, ComplexMatrix L)[]
            { (PauliLetter.X, lX), (PauliLetter.Y, lY), (PauliLetter.Z, lZ) };

        int total = 1 << (2 * n);
        double eigenResidual = 0.0, closedFormDev = 0.0, maxRate = 0.0;
        int repriced = 0, freeZ = 0, freeX = 0, freeY = 0, universalFree = 0;
        double allZunderZ = double.NaN, allZunderX = double.NaN,
               allXunderZ = double.NaN, allXunderX = double.NaN;
        int allZIndex = StringIndexOfUniform(PauliLetter.Z);
        int allXIndex = StringIndexOfUniform(PauliLetter.X);

        for (int s = 0; s < total; s++)
        {
            var letters = LettersOf(s);
            var v = VecOf(PauliString.Build(letters));
            var rates = new double[3];
            for (int p = 0; p < byLetter.Length; p++)
            {
                var (letter, l) = byLetter[p];
                var lv = l.Multiply(v);
                int m = MaxComponent(v);
                Complex rate = lv[m] / v[m];
                eigenResidual = Math.Max(eigenResidual, (lv - v.Multiply(rate)).AbsoluteMaximum().Magnitude);
                double closedForm = -2.0 * Gamma * AntiCount(letters, letter);
                closedFormDev = Math.Max(closedFormDev, (rate - new Complex(closedForm, 0)).Magnitude);
                rates[p] = rate.Real;
                maxRate = Math.Max(maxRate, rate.Magnitude);
            }
            // Free = rate 0; the rates are exact multiples of 2γ, so 1e-12 separates cleanly.
            bool fX = Math.Abs(rates[0]) < 1e-12, fY = Math.Abs(rates[1]) < 1e-12, fZ = Math.Abs(rates[2]) < 1e-12;
            if (fX) freeX++;
            if (fY) freeY++;
            if (fZ) freeZ++;
            if (fX && fY && fZ) universalFree++;
            if (Math.Abs(rates[0] - rates[2]) > 1e-12) repriced++;   // X-watcher vs Z-watcher
            if (s == allZIndex) { allZunderZ = rates[2]; allZunderX = rates[0]; }
            if (s == allXIndex) { allXunderZ = rates[2]; allXunderX = rates[0]; }
        }

        MaxEigenResidual = eigenResidual;
        MaxClosedFormDeviation = closedFormDev;
        RepricedCountZtoX = repriced;
        FreeCountZ = freeZ; FreeCountX = freeX; FreeCountY = freeY;
        UniversalFreeCount = universalFree;
        MaxRateMagnitude = maxRate;
        RateOfAllZUnderZ = allZunderZ; RateOfAllZUnderX = allZunderX;
        RateOfAllXUnderZ = allXunderZ; RateOfAllXUnderX = allXunderX;

        var adZX = AdUnitary(Hadamard1());
        var adYZ = AdUnitary(RxHalfPi1());
        TransportDevZtoX = MaxAbsDiff(adZX * lZ * adZX.ConjugateTranspose(), lX);
        TransportDevZtoY = MaxAbsDiff(adYZ * lZ * adYZ.ConjugateTranspose(), lY);
    }

    public string DisplayName =>
        $"the watched letter (label routing, N={N}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary =>
        $"N={N}: {1 << (2 * N)} Pauli strings, ONE shared eigenbasis, three price lists (rate −2γ·n_anti); "
      + $"eigen residual {MaxEigenResidual:0.0e+00}, closed-form dev {MaxClosedFormDeviation:0.0e+00}; "
      + $"repriced(Z→X) = {RepricedCountZtoX}; universally free = {UniversalFreeCount} (the identity alone); "
      + $"swap transport dev {Math.Max(TransportDevZtoX, TransportDevZtoY):0.0e+00}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the shared eigenbasis, three price lists (the gate)",
                summary: $"all {3 * (1 << (2 * N))} (letter, string) pairs: L_P·S = −2γ·n_anti(S,P)·S exactly "
                       + $"(eigen residual {MaxEigenResidual:0.0e+00}, closed-form dev {MaxClosedFormDeviation:0.0e+00}); "
                       + "the object never changes, only the price list is the letter's",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the routing (which cells pay follows the held letter)",
                summary: $"Z^⊗N: rate {RateOfAllZUnderZ.ToString("0.###", Inv)} under the Z-watcher, "
                       + $"{RateOfAllZUnderX.ToString("0.###", Inv)} under the X-watcher (X^⊗N mirrored: "
                       + $"{RateOfAllXUnderZ.ToString("0.###", Inv)} / {RateOfAllXUnderX.ToString("0.###", Inv)}); "
                       + $"{RepricedCountZtoX}/{1 << (2 * N)} strings repriced by the Z→X swap: 'noise' is a "
                       + "property of the pair (state, watcher), not of the state",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("only nothing rides free everywhere",
                summary: $"each watcher exempts its own 2^N = {1 << N} strings ({{I,P}}^⊗N: Z {FreeCountZ}, "
                       + $"X {FreeCountX}, Y {FreeCountY}); free under all three = {UniversalFreeCount} "
                       + "(the identity alone); everything that exists pays someone",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the letter swap is an exact transport",
                summary: $"Ad(h_zx)·L_Z·Ad† = L_X (dev {TransportDevZtoX:0.0e+00}), Ad(h_yz)·L_Z·Ad† = L_Y "
                       + $"(dev {TransportDevZtoY:0.0e+00}): the basis-S₃ face; the operator-space Klein V₄ "
                       + "face is PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE (typed parent)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("the reading (Tier 4, labeled)",
                summary: "to be a watcher is to hold a label and be blind past it; the label layer is the "
                       + "observer layer (bath, human, model — one architecture); labels are canvases, "
                       + "perspective-true, transported raw. Sources: docs/quantum/LABELS_TRANSLATED.md §2, "
                       + "docs/quantum/DEPHASING_TRANSLATED.md, docs/quantum/THE_LABEL_MAP.md");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    // ============================ builders (DiagonalWitness conventions) ============================

    /// <summary>L_P = γ·(Q_P − N·I), Q_P = Σ_l kron(P_l, P_lᵀ): the full letter dissipator.</summary>
    private ComplexMatrix LetterDissipator(PauliLetter letter)
    {
        var acc = Matrix<Complex>.Build.Dense(Dim, Dim);
        for (int l = 0; l < N; l++)
        {
            var p = PauliString.SiteOp(N, l, letter);
            acc += p.KroneckerProduct(p.Transpose());
        }
        return (acc - Matrix<Complex>.Build.DenseIdentity(Dim).Multiply(N)).Multiply(Gamma);
    }

    /// <summary>Row-stacking vec: v[i·d + j] = S[i, j].</summary>
    private ComplexVector VecOf(ComplexMatrix s)
    {
        var v = ComplexVector.Build.Dense(Dim);
        for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
                v[i * D + j] = s[i, j];
        return v;
    }

    private PauliLetter[] LettersOf(int s)
    {
        var letters = new PauliLetter[N];
        for (int site = 0; site < N; site++)
            letters[site] = (PauliLetter)((s >> (2 * site)) & 3);
        return letters;
    }

    private int StringIndexOfUniform(PauliLetter letter)
    {
        int s = 0;
        for (int site = 0; site < N; site++) s |= (int)letter << (2 * site);
        return s;
    }

    /// <summary>Sites whose letter anticommutes with the watched letter (∉ {I, P}).</summary>
    private static int AntiCount(PauliLetter[] letters, PauliLetter watched)
    {
        int k = 0;
        foreach (var l in letters)
            if (l != PauliLetter.I && l != watched) k++;
        return k;
    }

    private static int MaxComponent(ComplexVector v)
    {
        int m = 0; double best = -1.0;
        for (int i = 0; i < v.Count; i++)
            if (v[i].Magnitude > best) { best = v[i].Magnitude; m = i; }
        return m;
    }

    private ComplexMatrix AdUnitary(ComplexMatrix u1)
    {
        var u = u1;
        for (int s = 1; s < N; s++) u = u.KroneckerProduct(u1);   // U = u1^⊗N
        return u.KroneckerProduct(u.Conjugate());                  // Ad_U: vec(UρU†) = (U⊗U*)vec(ρ)
    }

    private static ComplexMatrix Hadamard1()
    {
        double s = 1.0 / Math.Sqrt(2.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { s, s }, { s, -s } });
    }

    private static ComplexMatrix RxHalfPi1()
    {
        double c = Math.Cos(Math.PI / 4.0), sn = Math.Sin(Math.PI / 4.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
            { { c, new Complex(0, -sn) }, { new Complex(0, -sn), c } });
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
