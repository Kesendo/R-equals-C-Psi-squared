using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The three-ladder hinge, recomputed live (<c>inspect --root ladders</c>): the disagreement
/// RUNG k = popcount(i⊕j), the F87 GIRTH ℓ, and the F120 MOMENT j are not three orthogonal axes — they
/// are the two factors of one F87-hardness coefficient on M = A + γQ, hinged by <b>Q</b>.
///
/// <para>On the 4^N coherence space: A = −i[H,·] carries the girth/moment side (A's closed walks on H's
/// hopping graph; the moments t_j(l) = Tr(Z_l H^j), girth ℓ = their onset). Q = Σ_l Z_l⊗Z_l is the rung
/// side: diagonal, Q_x = N − 2k(x). The F87 hardness coefficient (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §4)
/// is P_{m,1} = m·Tr(Q·A^{m−1}) = the rung Q weighting A's closed walks; the rung k(x) literally weights
/// each coherence's walk count (A^{m−1})_xx, and the supertrace factorizes the result into the girth
/// moments. So Q's SPECTRUM is the rung ladder, Q's ACTION is the girth-moment projector — the same
/// operator is both ladders, and that is why the three meet.</para>
///
/// <para>Five live nodes: (1) the hinge Q; (2) the bridge identity P_{m,1} = m·Tr(Q·A^{m−1}) = the girth
/// moments (m=3, the cell-free face); (3) the per-rung decomposition (how each rung k contributes); (4)
/// the rung is essential (remove Q → the moment projection is gone); (5) the general rung (m=5, not an
/// m=3 accident). The illustrative H is a uniform-Z-field XY chain (c_l = 1 per site, so the deg-1 face
/// fires). Verifier: <c>simulations/_three_ladders_bridge.py</c>. Convention: row-stacking vec,
/// kron(A,B): ρ ↦ AρBᵀ. Guard: 4^N ≤ <see cref="MaxDim"/> (N ≤ 5).</para></summary>
public sealed class LadderHingeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Largest coherence-space dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public int D { get; }       // 2^N
    public int Dim => D * D;     // 4^N coherence space

    private readonly ComplexMatrix _H, _A, _Q;   // H (2^N), A = −i[H,·], Q = Σ Z_l⊗Z_l (4^N)
    private readonly int[] _rung;                 // popcount(i⊕j) per coherence index m = i*d + j

    public LadderHingeWitness(int n = 3)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2; got {n}");
        long dim = 1L << (2 * n);
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n; D = 1 << n;
        _H = BuildH(N);
        _A = HamiltonianSuper(_H, N);        // A = −i[H,·]
        _Q = DephasingDiagonal(N);           // Q = Σ_l Z_l⊗Z_l (diagonal, Q_x = N − 2k)
        _rung = RungLabels(N);
    }

    public string DisplayName => $"the three-ladder hinge (N={N})";

    public string Summary
    {
        get
        {
            double dev3 = BridgeIdentityDev(3);
            double dev5 = BridgeIdentityDev(5);
            return $"N={N}: P_m1 = m·Tr(Q·A^(m-1)) = girth moments — rung×girth-walks, hinged by Q; "
                 + $"bridge dev m=3 {dev3:0.0e+00}, m=5 {dev5:0.0e+00}; the rung is essential";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheHinge();
            yield return TheBridgeIdentity();
            yield return ThePerRungDecomposition();
            yield return TheRungIsEssential();
            yield return TheGeneralRung();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    // ============================ public live readouts ============================
    /// <summary>The F87 deg-1 hardness coefficient as the rung weighting A's closed walks:
    /// P_{m,1} = m·Tr(Q·A^{m−1}) = m·Σ_x (N−2k(x))·(A^{m−1})_xx.</summary>
    public double RungWeightedCoefficient(int m) => m * (_Q * MatrixPow(_A, m - 1)).Trace().Real;

    /// <summary>The same coefficient as the girth moments: m·(−1)^k·Σ_l Σ_j (−1)^j C(2k,j) t_j(l) t_{2k−j}(l),
    /// k=(m−1)/2, t_j(l) = Tr(Z_l H^j) (the supertrace factorization, proof §4).</summary>
    public double GirthMomentCoefficient(int m)
    {
        int twoK = m - 1, k = twoK / 2;
        var towers = new IReadOnlyList<double>[twoK + 1];   // towers[j][l] = t_j(l) = Tr(Z_l H^j), j ≥ 1
        for (int j = 1; j <= twoK; j++) towers[j] = MomentTowerPumpChannelClaim.MomentTower(_H, j, N);
        double total = 0.0;
        for (int l = 0; l < N; l++)
        {
            double s = 0.0;
            for (int j = 0; j <= twoK; j++)
            {
                double tj = j == 0 ? 0.0 : towers[j][l];
                double tc = (twoK - j) == 0 ? 0.0 : towers[twoK - j][l];
                s += Sign(j) * Binom(twoK, j) * tj * tc;   // t_0 = Tr(Z_l) = 0
            }
            total += s;
        }
        return m * Sign(k) * total;
    }

    /// <summary>The coefficient with the rung weighting removed (Q → I): m·Tr(A^{m−1}). It does NOT equal
    /// the girth-moment form — the rung is essential to the projection.</summary>
    public double RungLessCoefficient(int m) => m * MatrixPow(_A, m - 1).Trace().Real;

    /// <summary>|rung-weighted − girth-moment| at moment m: the bridge identity dev (the gate).</summary>
    public double BridgeIdentityDev(int m) => Math.Abs(RungWeightedCoefficient(m) - GirthMomentCoefficient(m));

    /// <summary>How each rung k contributes to P_{m,1}: m·Σ_{x: k(x)=k} (N−2k)·(A^{m−1})_xx.</summary>
    public IReadOnlyDictionary<int, double> PerRungContribution(int m)
    {
        var apow = MatrixPow(_A, m - 1);
        var acc = new SortedDictionary<int, double>();
        for (int x = 0; x < Dim; x++)
        {
            int k = _rung[x];
            double contrib = m * (N - 2 * k) * apow[x, x].Real;
            acc[k] = acc.TryGetValue(k, out var v) ? v + contrib : contrib;
        }
        return acc;
    }

    // ============================ nodes ============================
    private IInspectable TheHinge()
    {
        return new InspectableNode("the hinge: Q is both ladders at once",
            summary: "M = A + γQ on the 4^N coherences. Q = Σ_l Z_l⊗Z_l is diagonal, Q_x = N − 2k(x): its "
                   + "SPECTRUM is the rung ladder k (the −2γk dissipative reading). Its ACTION (Σ Z_l⊗Z_l) "
                   + "projects A = −i[H,·]'s closed walks onto the Z-weighted girth moments t_j(l) = Tr(Z_l H^j) "
                   + "(girth ℓ = their onset). One operator, both ladders — the hinge.");
    }

    private IInspectable TheBridgeIdentity()
    {
        double rung = RungWeightedCoefficient(3);
        double girth = GirthMomentCoefficient(3);
        double dev = BridgeIdentityDev(3);
        return new InspectableNode("the bridge identity: P_{m,1} = m·Tr(Q·A^{m−1}) = the girth moments (m=3)",
            summary: $"the F87 hardness coefficient (PROOF_F87 §4) IS the rung weighting A's closed walks: "
                   + $"3·Tr(Q·A²) = {rung:0.###} (rung-weighted walks) = the girth-moment form {girth:0.###} "
                   + $"(= 6·4^N·Σc_l², the cell-free ℓ=1 face); dev = {dev:0.0e+00}. The rung k(x) literally "
                   + "weights each coherence's closed-walk count (A²)_xx.");
    }

    private IInspectable ThePerRungDecomposition()
    {
        var per = PerRungContribution(3);
        string body = string.Join(", ", per.Select(kv => $"k={kv.Key}:{kv.Value:+0.0;-0.0}"));
        return new InspectableNode("the per-rung decomposition (how the rung ladder enters)",
            summary: $"each rung k contributes 3·Σ_x(N−2k)(A²)_xx with weight (N−2k): {body}. The walk-diagonal "
                   + "(A²)_xx ≤ 0, so rungs below N/2 (weight +) contribute negative, above N/2 (weight −) "
                   + "positive — the rung k grades the girth-walk weight, made concrete.");
    }

    private IInspectable TheRungIsEssential()
    {
        double rungLess = RungLessCoefficient(3);
        double girth = GirthMomentCoefficient(3);
        return new InspectableNode("the rung is essential (remove Q → the moment projection is gone)",
            summary: $"drop the rung weighting (Q → I): 3·Tr(A²) = {rungLess:0.###}, nowhere near the girth-moment "
                   + $"form {girth:0.###}. Q's specific structure (Σ Z_l⊗Z_l) is what projects A's walks onto the "
                   + "single-site-Z moments t_ℓ; without it there is no bridge. The rung is not incidental.");
    }

    private IInspectable TheGeneralRung()
    {
        double rung = RungWeightedCoefficient(5);
        double girth = GirthMomentCoefficient(5);
        double dev = BridgeIdentityDev(5);
        return new InspectableNode("the general rung (m=5): not an m=3 accident",
            summary: $"5·Tr(Q·A⁴) = {rung:0.#} (rung-weighted walks) = the supertrace moment bilinear {girth:0.#} "
                   + $"(dev {dev:0.0e+00}): the rung-weights-the-girth-walks bridge holds at every rung, the "
                   + "coefficient factorizing into the d-leg girth moments at each power of A.");
    }

    // ============================ builders (reused conventions) ============================
    /// <summary>Illustrative H: a uniform Z field (c_l = 1 per site, so the deg-1 girth face fires)
    /// plus 0.5·XY hopping on the chain. H = Σ_l Z_l + 0.5·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}).</summary>
    private static ComplexMatrix BuildH(int n)
    {
        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        for (int l = 0; l < n; l++)
            h += PauliString.SiteOp(n, l, PauliLetter.Z).Multiply(Complex.One);              // Σ_l Z_l
        for (int b = 0; b < n - 1; b++)
        {
            var xx = PauliString.SiteOp(n, b, PauliLetter.X) * PauliString.SiteOp(n, b + 1, PauliLetter.X);
            var yy = PauliString.SiteOp(n, b, PauliLetter.Y) * PauliString.SiteOp(n, b + 1, PauliLetter.Y);
            h += (xx + yy).Multiply(new Complex(0.5, 0));
        }
        return h;
    }

    private static ComplexMatrix HamiltonianSuper(ComplexMatrix h, int n)
    {
        int d = 1 << n;
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var lh = h.KroneckerProduct(idH) - idH.KroneckerProduct(h.Transpose());   // [H,·]
        return lh.Multiply(new Complex(0, -1));                                     // A = −i[H,·]
    }

    private static ComplexMatrix DephasingDiagonal(int n)
    {
        int d = 1 << n, d2 = d * d;
        var acc = Matrix<Complex>.Build.Dense(d2, d2);
        for (int l = 0; l < n; l++)
        {
            var z = PauliString.SiteOp(n, l, PauliLetter.Z);
            acc += z.KroneckerProduct(z.Transpose());   // Z_lᵀ = Z_l ⟹ Σ_l Z_l⊗Z_l
        }
        return acc;
    }

    private static int[] RungLabels(int n)
    {
        int d = 1 << n;
        var r = new int[d * d];
        for (int i = 0; i < d; i++)
            for (int jj = 0; jj < d; jj++)
                r[i * d + jj] = System.Numerics.BitOperations.PopCount((uint)(i ^ jj));
        return r;
    }

    // ============================ helpers ============================
    private static ComplexMatrix MatrixPow(ComplexMatrix a, int p)
    {
        var r = Matrix<Complex>.Build.DenseIdentity(a.RowCount);
        for (int i = 0; i < p; i++) r *= a;
        return r;
    }

    private static double Sign(int j) => (j & 1) == 0 ? 1.0 : -1.0;

    private static double Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        double r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }
}
