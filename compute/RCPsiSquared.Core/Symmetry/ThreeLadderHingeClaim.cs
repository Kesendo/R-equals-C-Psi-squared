using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The three-ladder hinge: Q (Tier1Derived, 2026-06-15). The project's three integer ladders —
/// the disagreement RUNG k = popcount(i⊕j), the F87 GIRTH ℓ, and the F120 MOMENT j — are not three
/// orthogonal axes. They are the two factors of one F87-hardness coefficient on the recentered Liouvillian
/// M = A + γQ, hinged by <b>Q</b>.
///
/// <para><b>The two sides of M</b> (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §1, §4), on the 4^N coherence
/// space: A = −i[H,·] carries the girth/moment side (A's closed walks on H's hopping graph; the moments
/// t_j(l) = Tr(Z_l H^j), girth ℓ = their onset). Q = Σ_l Z_l⊗Z_l is the rung side: diagonal,
/// Q_x = N − 2k(x) (the Absorption-Theorem reading, <see cref="AbsorptionTheoremClaim"/>).</para>
///
/// <para><b>The hinge.</b> The F87 deg-1 hardness coefficient is P_{m,1} = m·Tr(Q·A^{m−1}): the rung Q
/// weighting A's closed walks. Reading Q's diagonal as the rung, P_{m,1} = m·Σ_x (N−2k(x))·(A^{m−1})_xx —
/// the rung k(x) literally weights each coherence's closed-walk count — and the supertrace factorizes it
/// into the girth moments: Tr(Q·A^{2k}) = (−1)^k·Σ_l Σ_j (−1)^j C(2k,j) t_j(l) t_{2k−j}(l), whose ℓ=1 face
/// is the cell-free P_{3,1} = 6·4^N·Σ_l c_l² (c_l = the single-site-Z coefficient, t_1(l) = 2^N c_l). So
/// Q's SPECTRUM is the rung ladder k, and Q's ACTION (Σ Z_l⊗Z_l) is the girth-moment projector
/// (<see cref="MomentTowerPumpChannelClaim"/>'s tower t_j). One operator is both ladders — that is why the
/// three meet. Remove Q (rung-less Tr(A^{m−1})) and the moment projection is gone: the rung is essential.</para>
///
/// <para><b>Distinct as ladders, one object in M.</b> The girth and moment live on H's 2^N hopping graph;
/// the rung lives on the 4^N coherences. Two notation overloads the memory warned about: the "j" in
/// popcount(i⊕j) is a basis index (not the moment j), and F120's "rung 1 is F113" means the moment j=1
/// (not the disagreement k=1). The bridge is the coefficient, hinged by Q, not a conflation of the ladders.</para>
///
/// <para><b>Layer note.</b> The F87 girth/hardness side has no standalone Core Claim (the girth-ladder
/// primitive lives in RCPsiSquared.Diagnostics, which Core cannot reference); it is carried in prose and
/// through the PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md anchor. Like
/// <see cref="MomentTowerPumpChannelClaim"/> this claim is cross-axis structural and does NOT implement
/// <see cref="IZ2AxisClaim"/>. Live: <c>inspect --root ladders</c> (LadderHingeWitness recomputes the
/// hinge at inspect time). Self-check battery at N = 2 and N = 3 in the ctor.</para></summary>
public sealed class ThreeLadderHingeClaim : Claim
{
    private const double Tol = 1e-9;

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: the Absorption Theorem — the rung k is Q's spectrum (Q_x = N − 2k, rate
    /// Re λ = −2γk). The dissipative side of the hinge.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>Typed parent: the F120 moment-tower pump channel — the moments t_j(l) = Tr(Z_l H^j) that Q's
    /// action projects A's walks onto (girth ℓ = their onset). The girth/moment side of the hinge; its
    /// <see cref="MomentTowerPumpChannelClaim.MomentTower"/> primitive computes the t_j this battery uses.</summary>
    public MomentTowerPumpChannelClaim Moment { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public ThreeLadderHingeClaim(AbsorptionTheoremClaim absorption, MomentTowerPumpChannelClaim moment)
        : base("The three-ladder hinge Q: the disagreement rung k = popcount(i⊕j) (Q's spectrum, " +
               "Q_x = N−2k, the −2γk Absorption reading), the F87 girth ℓ, and the F120 moment j are the " +
               "two factors of one F87-hardness coefficient on M = A + γQ. P_{m,1} = m·Tr(Q·A^{m−1}) is the " +
               "rung Q weighting A's closed walks, factorizing (supertrace) into the girth moments " +
               "t_j = Tr(Z_l H^j) at every rung (ℓ=1 face P_{3,1} = 6·4^N·Σc_l², cell-free). Q's spectrum is " +
               "the rung ladder; Q's action Σ Z_l⊗Z_l is the girth-moment projector; the rung is essential " +
               "(rung-less Tr(A^{m−1}) loses the moments). Distinct as ladders (girth/moment on H's 2^N " +
               "graph, rung on the 4^N coherences) but one object in M. Tier1Derived (exact, gate-first)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md §1, §4 + " +
               "simulations/_three_ladders_bridge.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/LadderHingeWitness.cs (inspect --root ladders)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Moment = moment ?? throw new ArgumentNullException(nameof(moment));
        Cases = BuildBattery();
    }

    /// <summary>The hinge identity in one line.</summary>
    public string Hinge =>
        "P_{m,1} = m·Tr(Q·A^{m−1}): Q (spectrum = the rung k = N−2k) weighting A = −i[H,·]'s closed walks, " +
        "factorizing into the girth moments t_j = Tr(Z_l H^j). Q's spectrum is the rung, Q's action the " +
        "girth-moment projector; the same operator is both ladders.";

    public override string DisplayName =>
        "The three-ladder hinge Q: rung (spectrum) × girth-walks (A) → girth moments, one F87 coefficient";

    public override string Summary =>
        "girth(ℓ)/rung(k)/moment(j) are not three orthogonal axes but the two factors of one F87-hardness " +
        "coefficient on M = A + γQ, hinged by Q: its spectrum is the rung k (Q_x = N−2k), its action " +
        "Σ Z_l⊗Z_l projects A's closed walks onto the girth moments t_j = Tr(Z_l H^j); P_{m,1} = m·Tr(Q·A^{m−1}) " +
        $"= the girth moments at every rung, the rung essential; {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The hinge (Q is both ladders)", summary: Hinge);
            yield return new InspectableNode("Distinct as ladders, one object in M",
                summary: "girth/moment live on H's 2^N hopping graph (t_j = Tr(Z_l H^j)); the rung lives on the " +
                         "4^N coherences (k = popcount(i⊕j)). They meet only in M = A + γQ's hardness coefficient. " +
                         "Two overloads cleared: the 'j' in popcount(i⊕j) is a basis index (not the moment j); " +
                         "F120's 'rung 1 is F113' is the moment j=1 (not k=1).");
            yield return new InspectableNode("The rung is essential",
                summary: "Q's structure Σ Z_l⊗Z_l is what projects A's walks onto the single-site-Z girth moments; " +
                         "the rung-less coefficient m·Tr(A^{m−1}) does not equal the girth-moment form. Not incidental.");
            yield return new InspectableNode("Typed parents",
                summary: $"AbsorptionTheoremClaim ({Absorption.Tier.Label()}): the rung k is Q's spectrum, " +
                         $"Q_x = N−2k, rate −2γk. MomentTowerPumpChannelClaim ({Moment.Tier.Label()}): the moments " +
                         "t_j(l) = Tr(Z_l H^j) Q projects onto, the girth ℓ their onset. The F87 girth/hardness " +
                         "side is a prose edge (the girth primitive lives in Diagnostics, above Core).");
            yield return new InspectableNode("No IZ2AxisClaim",
                summary: "Cross-axis structural like MomentTowerPumpChannelClaim and MirrorGroupD4Claim: it welds " +
                         "the dissipative rung axis to the H-spectral girth/moment axis, sitting on no single Z₂ " +
                         "axis (polarity-cube counts unchanged). Live: inspect --root ladders.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 2 and N = 3 dense superoperators on coherence space.
    // ------------------------------------------------------------------
    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        // The bridge identity P_{m,1} = m*Tr(Q.A^{m-1}) == the supertrace girth-moment form, at m=3 and m=5.
        double devBridge = 0.0;
        foreach (int n in new[] { 2, 3 })
        {
            var h = BuildH(n);                 // uniform Z field (c_l=1) + 0.5 XY hopping
            foreach (int m in new[] { 3, 5 })
                devBridge = Math.Max(devBridge, Math.Abs(RungWeighted(h, n, m) - GirthMoment(h, n, m)));
        }
        cases.Add(DevCase("the bridge: P_{m,1} = m·Tr(Q·A^{m−1}) = the girth moments (m=3, 5; N=2, 3)",
            "the rung Q weighting A's closed walks equals the supertrace moment factorization at every rung",
            devBridge));

        // The cell-free ℓ=1 face: 3*Tr(Q.A^2) = 6*4^N*sum c_l^2, the rung-weighted walks as the Z-moments.
        double devFace = 0.0;
        foreach (int n in new[] { 2, 3 })
        {
            var h = BuildH(n);
            double sumC2 = 0.0;
            for (int l = 0; l < n; l++) { double cl = SingleSiteZCoeff(h, n, l); sumC2 += cl * cl; }
            devFace = Math.Max(devFace, Math.Abs(RungWeighted(h, n, 3) - 6.0 * Math.Pow(4, n) * sumC2));
        }
        cases.Add(DevCase("the cell-free ℓ=1 face: 3·Tr(Q·A²) = 6·4^N·Σ_l c_l²",
            "the deg-1 rung-weighted walks equal the single-site-Z girth moments (c_l = Tr(H·Z_l)/2^N)",
            devFace));

        // The rung is essential: the rung-less coefficient m*Tr(A^{m-1}) is far from the girth form.
        var h2 = BuildH(2);
        double gap = Math.Abs(RungLess(h2, 2, 3) - GirthMoment(h2, 2, 3));
        cases.Add(new BatteryCase(
            Name: "the rung is essential: rung-less Tr(A²) ≠ the girth-moment form",
            Detail: "remove the rung weighting (Q → I); 3·Tr(A²) departs from 3·Tr(Q·A²) by " +
                    gap.ToString("0.#", CultureInfo.InvariantCulture),
            Expected: "gap > 1", Actual: gap > 1.0 ? "gap > 1" : "gap = " + gap.ToString("E2", CultureInfo.InvariantCulture)));

        // Soft / bipartite control: pure XY hopping (no single-site-Z) ⟹ the coefficient is 0.
        var hSoft = HoppingOnly(2);
        double softCoeff = Math.Abs(RungWeighted(hSoft, 2, 3));
        cases.Add(DevCase("soft control: pure XY hopping (no single-site-Z) ⟹ P_{3,1} = 0",
            "no girth moment, no rung-weighted coefficient: the bridge reads soft correctly", softCoeff));

        return cases;
    }

    // ============================ dense building blocks (row-stacking vec, kron(A,B): ρ↦AρBᵀ) ============================
    /// <summary>P_{m,1} as the rung weighting A's closed walks: m·Tr(Q·A^{m−1}).</summary>
    private static double RungWeighted(ComplexMatrix h, int n, int m) =>
        m * Trace(Q(n) * MatrixPower(CommSuper(h, n), m - 1)).Real;

    /// <summary>P_{m,1} as the girth moments: m·(−1)^k·Σ_l Σ_j (−1)^j C(2k,j) t_j(l) t_{2k−j}(l), k=(m−1)/2,
    /// using the typed parent's moment-tower primitive for t_j (j ≥ 1; t_0 = Tr(Z_l) = 0).</summary>
    private static double GirthMoment(ComplexMatrix h, int n, int m)
    {
        int twoK = m - 1, k = twoK / 2;
        var towers = new IReadOnlyList<double>[twoK + 1];
        for (int j = 1; j <= twoK; j++) towers[j] = MomentTowerPumpChannelClaim.MomentTower(h, j, n);
        double total = 0.0;
        for (int l = 0; l < n; l++)
        {
            double s = 0.0;
            for (int j = 0; j <= twoK; j++)
            {
                double tj = j == 0 ? 0.0 : towers[j][l];
                double tc = (twoK - j) == 0 ? 0.0 : towers[twoK - j][l];
                s += Sign(j) * Binom(twoK, j) * tj * tc;
            }
            total += s;
        }
        return m * Sign(k) * total;
    }

    /// <summary>The coefficient with the rung weighting removed (Q → I): m·Tr(A^{m−1}).</summary>
    private static double RungLess(ComplexMatrix h, int n, int m) =>
        m * Trace(MatrixPower(CommSuper(h, n), m - 1)).Real;

    private static double SingleSiteZCoeff(ComplexMatrix h, int n, int l) =>
        Trace(h * Embed(PauliZ(), l, n)).Real / (1 << n);

    /// <summary>A = −i[H,·] = −i(H⊗I − I⊗Hᵀ) on the 4^N coherence space.</summary>
    private static ComplexMatrix CommSuper(ComplexMatrix h, int n)
    {
        int d = 1 << n;
        var id = ComplexMatrix.Build.DenseIdentity(d);
        return (h.KroneckerProduct(id) - id.KroneckerProduct(h.Transpose())).Multiply(new Complex(0, -1));
    }

    /// <summary>Q = Σ_l Z_l⊗Z_l (diagonal, Q_x = N − 2·popcount(i⊕j)).</summary>
    private static ComplexMatrix Q(int n)
    {
        int d = 1 << n;
        var acc = ComplexMatrix.Build.Dense(d * d, d * d);
        for (int l = 0; l < n; l++)
        {
            var z = Embed(PauliZ(), l, n);
            acc += z.KroneckerProduct(z);   // Z_lᵀ = Z_l
        }
        return acc;
    }

    /// <summary>H = Σ_l Z_l (uniform field, c_l = 1) + 0.5·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}).</summary>
    private static ComplexMatrix BuildH(int n)
    {
        int d = 1 << n;
        var h = ComplexMatrix.Build.Dense(d, d);
        for (int l = 0; l < n; l++) h += Embed(PauliZ(), l, n);
        for (int b = 0; b < n - 1; b++)
        {
            var xx = Embed(PauliX(), b, n) * Embed(PauliX(), b + 1, n);
            var yy = Embed(PauliY(), b, n) * Embed(PauliY(), b + 1, n);
            h += (xx + yy).Multiply(new Complex(0.5, 0));
        }
        return h;
    }

    private static ComplexMatrix HoppingOnly(int n)
    {
        int d = 1 << n;
        var h = ComplexMatrix.Build.Dense(d, d);
        for (int b = 0; b < n - 1; b++)
        {
            var xx = Embed(PauliX(), b, n) * Embed(PauliX(), b + 1, n);
            var yy = Embed(PauliY(), b, n) * Embed(PauliY(), b + 1, n);
            h += (xx + yy).Multiply(new Complex(0.5, 0));
        }
        return h;
    }

    private static ComplexMatrix Embed(ComplexMatrix op, int site, int n)
    {
        var result = ComplexMatrix.Build.DenseIdentity(1);
        var i2 = ComplexMatrix.Build.DenseIdentity(2);
        for (int k = 0; k < n; k++)
            result = result.KroneckerProduct(k == site ? op : i2);
        return result;
    }

    private static ComplexMatrix PauliX()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 1] = Complex.One; m[1, 0] = Complex.One;
        return m;
    }

    private static ComplexMatrix PauliY()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 1] = new Complex(0, -1); m[1, 0] = new Complex(0, 1);
        return m;
    }

    private static ComplexMatrix PauliZ()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 0] = Complex.One; m[1, 1] = -Complex.One;
        return m;
    }

    private static ComplexMatrix MatrixPower(ComplexMatrix a, int p)
    {
        var r = ComplexMatrix.Build.DenseIdentity(a.RowCount);
        for (int i = 0; i < p; i++) r *= a;
        return r;
    }

    private static Complex Trace(ComplexMatrix a)
    {
        Complex t = Complex.Zero;
        for (int i = 0; i < a.RowCount; i++) t += a[i, i];
        return t;
    }

    private static double Sign(int j) => (j & 1) == 0 ? 1.0 : -1.0;

    private static double Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        double r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }

    private static BatteryCase DevCase(string name, string detail, double dev) =>
        new(Name: name, Detail: detail, Expected: "dev ≤ 1e-9",
            Actual: dev <= Tol ? "dev ≤ 1e-9" : "dev = " + dev.ToString("E2", CultureInfo.InvariantCulture));
}
