using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The mirror's order-sorting law (F131, Tier1Derived, minted 2026-07-16): a mirror
/// does not throw information away, it sorts it. For a mirror M (a unitary or antiunitary
/// involution) whose conjugation reflects a parameter scan,
///
/// <code>
///   M · G(x₀ + s·δ) · M⁻¹ = G(x₀ + σ_eff·s·δ),   σ_eff = σ_op·χ_M,
/// </code>
///
/// with σ_op the operator parity of the direction δ under M and χ_M = +1 for a linear
/// (unitary) mirror, −1 for an antilinear (antiunitary) one (the antiunitary flips the i in
/// the exponent, which is how an M-FIXED perturbation can still scan reflected), and a
/// readout of definite mirror parity q, the response orders sort by the product q·σ_eff
/// into four cells: <b>generic</b> (+,+), <b>EVEN response</b> (+,−; only even orders
/// enter, the first order vanishes), <b>ODD response</b> (−,−), and <b>IDENTICALLY
/// ZERO</b> (−,+; a static selection rule, not a response-order statement). The law asserts
/// PARITIES only, never magnitudes.
///
/// <para><b>Theorem A (the unitary column, unconditional):</b> R = the F71 site reversal on
/// the open XX chain with local Z dephasing; parameters (J per bond, γ per site, h per
/// site) split into an R-even base plus t times an R-odd direction (σ_op = −1, χ_R = +1,
/// σ_eff = −1). The conjugation identity is entry-wise and exact,
/// (R⊗R)·L(base + t·dir)·(R⊗R) = L(base − t·dir), and for an OPERATOR-R-even preparation
/// (RρR = ρ; expectation-even is provably insufficient) and a readout with ROR = qO,
/// ⟨O⟩(t; time) = q·⟨O⟩(−t; time) at every evolution time. The full Liouville spectrum is
/// even in t (corollary); the owned F91/F92/F93 statements are the STRONGER all-orders
/// pair-sum invariance of a sub-object (the F71-refined diagonal blocks), not the same
/// object.</para>
///
/// <para><b>Theorem B (the antiunitary column):</b> the ζ² anti-protection law
/// (PROOF_ZETA2_ANTI_PROTECTION.md, the F129 derived law): Θ = T·K (class BDI, owner
/// <see cref="ChiralKClaim"/>) commutes with the brickwork Floquet step and FIXES the
/// occupation-diagonal ZZ perturbation (σ_op = +1), but antiunitarity gives χ_Θ = −1, so
/// σ_eff = −1 and Θ·W(ζ)·Θ⁻¹ = W(−ζ): the mirror-pair DIFFERENCE reads even (the ζ² law,
/// exact factor 2), the pair SUM reads odd. This column pays real extra hypotheses (a
/// simple, isolated eigenphase tracked continuously); the multiset eigenphase identity is
/// unconditional, exactly as in Theorem A.</para>
///
/// <para><b>The boundary is the physics:</b> the preparation hypothesis is operator-level;
/// an ε-impure preparation (ρ = ρ_even + ε·ρ_odd) opens the forbidden parity EXACTLY
/// affinely in ε (the master equation is linear in ρ₀), leading order O(ε·t). In both
/// examined linear leaks (the ε-leak; the flown b_qs budget on the Θ side) the failure was
/// a broken parity hypothesis, not a broken mirror. Π-protected observables are a
/// DIFFERENT zero mechanism (degenerate-cluster cancellation, computed per initial state)
/// and are deliberately NOT the zero cell's witness.</para>
///
/// <para><b>Self-check battery (Theorem A, N = 3):</b> built in the constructor on the
/// 64×64 Liouville superoperator, moment-level, with no eigensolver and no matrix
/// exponential. The exact conjugation identity (checked entry-wise, machine zero) gives
/// the per-order parity Tr[O·L(−t)^k·ρ₀] = q·Tr[O·L(+t)^k·ρ₀] for ALL k, and by
/// analyticity of e^{L·T} the full moment tower is the all-times trajectory statement;
/// the battery spot-checks the moments at k = 0..6. Ten cases: the conjugation identity on all three axes (machine zero), the
/// mixed-scan fence (O(1) rejection), the σ_eff = +1 static column, the preparation
/// parities as operator identities, the even/odd/zero cells (odd cell non-vacuous, generic
/// cell responding), and the exactly-affine ε-leak (halving ratio 2). Mirrors the committed
/// gate simulations/mirror_order_sorting.py (which adds the spectral-evenness face by
/// optimal-assignment matching) and the MirrorWorld adoption
/// compute/MirrorWorld/OrderSorting.cs (twin RK4, run mode `sorting N`).</para>
///
/// <para><b>Layer note:</b> like <see cref="AntilinearTriangleClaim"/> this claim is
/// cross-axis structural and deliberately does NOT implement <see cref="IZ2AxisClaim"/>.
/// Typed parents: <see cref="ChiralKClaim"/> (Θ = T·K, the owner of Theorem B's mirror),
/// <see cref="AntilinearTriangleClaim"/> (χ_M: the linear/antilinear character and the
/// transport law that flips the exponent's i), and the F91 family
/// (<see cref="RCPsiSquared.Core.BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance"/> +
/// <see cref="RCPsiSquared.Core.SymmetryFamily.F92BondAntiPalindromicJSpectralInvariance"/> +
/// <see cref="RCPsiSquared.Core.SymmetryFamily.F93DetuningAntiPalindromicSpectralInvariance"/>,
/// the owned witnesses of the σ_eff = −1 pencil identity on the three axes).</para></summary>
public sealed class MirrorOrderSortingClaim : Claim
{
    private const double Tol = 1e-10;
    private const int BatteryN = 3;

    /// <summary>One self-check tying the claim to the order-sorting identities at N = 3.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: Θ = T·K, class BDI (the K sublattice/chiral symmetry). The
    /// antiunitary mirror of Theorem B is built on K; the σ_op = +1 of the ZZ direction
    /// under Θ is K's occupation-diagonal fixing.</summary>
    public ChiralKClaim ChiralK { get; }

    /// <summary>Typed parent: the antilinear triangle. χ_M = −1 for antiunitary mirrors is
    /// the triangle's transport law (the antilinear maps pay the sign of the −i); the
    /// order-sorting principle reads σ_eff = σ_op·χ_M off exactly that grading.</summary>
    public AntilinearTriangleClaim Triangle { get; }

    /// <summary>Typed parent: F91 (γ axis). The anti-palindromic pair-sum invariance is the
    /// STRONGER all-orders sub-object law on the same axis the σ_eff = −1 scan runs.</summary>
    public BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance F91 { get; }

    /// <summary>Typed parent: F92 (J axis), the bond twin.</summary>
    public SymmetryFamily.F92BondAntiPalindromicJSpectralInvariance F92 { get; }

    /// <summary>Typed parent: F93 (h axis), the detuning twin.</summary>
    public SymmetryFamily.F93DetuningAntiPalindromicSpectralInvariance F93 { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public MirrorOrderSortingClaim(
        ChiralKClaim chiralK,
        AntilinearTriangleClaim triangle,
        BlockSpectrum.F71AntiPalindromicGammaSpectralInvariance f91,
        SymmetryFamily.F92BondAntiPalindromicJSpectralInvariance f92,
        SymmetryFamily.F93DetuningAntiPalindromicSpectralInvariance f93)
        : base("F131, the mirror's order-sorting law: mirror conjugation reflects a parameter scan, " +
               "M·G(x₀ + s·δ)·M⁻¹ = G(x₀ + σ_eff·s·δ) with σ_eff = σ_op·χ_M (direction parity × " +
               "linear/antilinear character), and for a readout of definite mirror parity q the response " +
               "orders sort by q·σ_eff into four cells: generic / EVEN response / ODD response / " +
               "IDENTICALLY ZERO. Theorem A: the unitary F71 site reversal on the Lindblad chain, " +
               "unconditional, ⟨O⟩(t) = q·⟨O⟩(−t) for operator-R-even preparation. Theorem B: the ζ² " +
               "anti-protection law (antiunitary Floquet Θ = T·K; branch tracking hypotheses). Parities " +
               "only, never magnitudes; the ε-leak under a broken preparation hypothesis is exactly " +
               "affine. Tier1Derived (assembly of proven Tier-1 results; every cell witness owned)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_MIRROR_ORDER_SORTING.md + " +
               "docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md + " +
               "simulations/mirror_order_sorting.py (gate, PASS N = 4, 5 default; N = 6 via --deep) + " +
               "simulations/zeta2_anti_protection.py (Theorem B gate) + " +
               "compute/MirrorWorld/OrderSorting.cs (the MirrorWorld adoption, run mode sorting N)")
    {
        ChiralK = chiralK ?? throw new ArgumentNullException(nameof(chiralK));
        Triangle = triangle ?? throw new ArgumentNullException(nameof(triangle));
        F91 = f91 ?? throw new ArgumentNullException(nameof(f91));
        F92 = f92 ?? throw new ArgumentNullException(nameof(f92));
        F93 = f93 ?? throw new ArgumentNullException(nameof(f93));
        Cases = BuildBattery();
    }

    /// <summary>The principle in one line.</summary>
    public string Principle =>
        "M·G(x₀ + s·δ)·M⁻¹ = G(x₀ + σ_eff·s·δ) with σ_eff = σ_op·χ_M; a readout of definite mirror " +
        "parity q reads only the orders allowed by q·σ_eff: (+,+) generic, (+,−) even orders, " +
        "(−,−) odd orders, (−,+) identically zero. Parities only, never magnitudes.";

    /// <summary>Theorem A in one line.</summary>
    public string TheoremA =>
        "Unitary column: (R⊗R)·L(base + t·dir)·(R⊗R) = L(base − t·dir) for the F71 site reversal and " +
        "an R-odd direction on any of the three axes (J, γ, h); operator-R-even preparation and " +
        "ROR = qO give ⟨O⟩(t; time) = q·⟨O⟩(−t; time) at every evolution time, unconditionally.";

    /// <summary>Theorem B in one line.</summary>
    public string TheoremB =>
        "Antiunitary column: Θ = T·K fixes the ZZ direction (σ_op = +1) but pays χ_Θ = −1, so " +
        "Θ·W(ζ)·Θ⁻¹ = W(−ζ): the Floquet mirror-pair difference reads even (the ζ² law, exact " +
        "factor 2), the pair sum reads odd; branch-resolved statements pay tracking hypotheses, " +
        "the multiset identity is unconditional.";

    /// <summary>The boundary in one line.</summary>
    public string Boundary =>
        "The hypotheses are the physics: operator-even preparation (expectation-even provably " +
        "insufficient); the ε-leak is exactly affine (master equation linear in ρ₀), leading O(ε·t); " +
        "in both examined linear leaks the failure was a broken parity hypothesis, not a broken " +
        "mirror. Π-protected observables are a different zero mechanism (cluster cancellation), " +
        "deliberately not the zero cell's witness.";

    public override string DisplayName =>
        "F131: the mirror's order-sorting law — response orders sort by q·σ_eff (four cells, two theorems)";

    public override string Summary =>
        "mirror conjugation reflects a parameter scan (σ_eff = σ_op·χ_M) and a definite-parity readout q " +
        "reads only the orders q·σ_eff allows: generic / EVEN / ODD / IDENTICALLY ZERO; Theorem A the " +
        "unitary F71 column (unconditional), Theorem B the antiunitary ζ² column (tracking hypotheses); " +
        $"{PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The principle (σ_eff = σ_op·χ_M)", summary: Principle);
            yield return new InspectableNode("Theorem A (unitary, unconditional)", summary: TheoremA);
            yield return new InspectableNode("Theorem B (antiunitary, the ζ² law)", summary: TheoremB);
            yield return new InspectableNode("The boundary (the hypotheses are the physics)", summary: Boundary);
            yield return new InspectableNode("Spectral face vs trajectory face",
                summary: "The full Liouville spectrum is even in t (corollary of the same identity); the " +
                         "owned F91/F92/F93 pair-sum invariance is the STRONGER all-orders law of a " +
                         "sub-object (the F71-refined diagonal blocks). Two faces of one identity, not " +
                         "one object; the battery here carries the trajectory face (moments), the " +
                         "committed gate adds the spectral face by optimal-assignment matching.");
            yield return new InspectableNode("Hardware, honestly sized",
                summary: "Confirmation 24 (ibm_kingston, the standing fringe) sights the FIRST-order " +
                         "protection only; the even-order anti-protection stayed inside its " +
                         "pre-registered budget and is not independently measured; no other cell has a " +
                         "hardware witness. The ζ²-meter is a candidate future instrument, not a result.");
            yield return new InspectableNode("MirrorWorld adoption",
                summary: "compute/MirrorWorld/OrderSorting.cs (2026-07-16, run mode sorting N): Theorem " +
                         "A's trajectory face by twin RK4 on the γ axis, the pencil face as " +
                         "ParameterKlein.MirrorConjugationResidual; guarded by OrderSortingTests (6). " +
                         "Theorem B and the spectral evenness stay out (the drawn boundary).");
            yield return new InspectableNode("Typed parents",
                summary: $"ChiralKClaim ({ChiralK.Tier.Label()}): Θ = T·K, Theorem B's mirror; " +
                         $"AntilinearTriangleClaim ({Triangle.Tier.Label()}): χ_M, the linear/antilinear " +
                         $"character the σ_eff formula reads; F91/F92/F93 " +
                         $"({F91.Tier.Label()}/{F92.Tier.Label()}/{F93.Tier.Label()}): the three axes' " +
                         "owned pair-sum invariance, the stronger sub-object law beside the σ_eff = −1 scan.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 3, 64×64 Liouville superoperator, moments.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        int n = BatteryN;
        int d = 1 << n;                    // 8
        int d2 = d * d;                    // 64

        // --- single-site operators and the chain builder (Pauli convention, XX+YY bonds) ---
        var sx = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
        var sy = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, new Complex(0, -1) }, { new Complex(0, 1), 0 } });
        var sz = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });
        var s0 = ComplexMatrix.Build.DenseIdentity(2);

        ComplexMatrix Op(ComplexMatrix single, int site)
        {
            var m = site == 0 ? single : s0;
            for (int k = 1; k < n; k++)
                m = m.KroneckerProduct(k == site ? single : s0);
            return m;
        }

        var idD = ComplexMatrix.Build.DenseIdentity(d);
        var idD2 = ComplexMatrix.Build.DenseIdentity(d2);
        var minusI = new Complex(0, -1);
        var zOps = Enumerable.Range(0, n).Select(l => Op(sz, l)).ToArray();

        ComplexMatrix Superop(double[] js, double[] gammas, double[] hs)
        {
            var h = ComplexMatrix.Build.Dense(d, d);
            for (int b = 0; b < n - 1; b++)
                h += (Op(sx, b) * Op(sx, b + 1) + Op(sy, b) * Op(sy, b + 1)).Multiply(new Complex(js[b], 0));
            for (int l = 0; l < n; l++)
                h += zOps[l].Multiply(new Complex(hs[l], 0));
            var m = (h.KroneckerProduct(idD) - idD.KroneckerProduct(h.Transpose())).Multiply(minusI);
            for (int l = 0; l < n; l++)
                m += (zOps[l].KroneckerProduct(zOps[l].Transpose()) - idD2).Multiply(new Complex(gammas[l], 0));
            return m;
        }

        // --- the F71 site reversal R (bit reversal on basis configurations) and R⊗R ---
        var r = ComplexMatrix.Build.Dense(d, d);
        for (int c = 0; c < d; c++)
        {
            int rev = 0;
            for (int s = 0; s < n; s++)
                if (((c >> s) & 1) == 1) rev |= 1 << (n - 1 - s);
            r[rev, c] = Complex.One;
        }
        var rr = r.KroneckerProduct(r);

        // --- profiles: R-even bases, R-odd directions (bitwise exact), scan scale t ---
        const double t = 0.37;
        double[] jBase = { 1.0, 1.0 }, jOdd = { -0.3, 0.3 };
        double[] gBase = { 0.05, 0.08, 0.05 }, gOdd = { -0.1, 0.0, 0.1 };
        double[] hBase = { 0.3, 0.3, 0.3 }, hOdd = { -0.2, 0.0, 0.2 };
        double[] Scan(double[] baseP, double[] dir, double s)
        {
            var res = new double[baseP.Length];
            for (int l = 0; l < res.Length; l++) res[l] = baseP[l] + s * dir[l];
            return res;
        }

        var cases = new List<BatteryCase>();

        // (a) The conjugation identity on the three axes: (R⊗R)·L(+t)·(R⊗R) = L(−t).
        double DevConj(double[] jsP, double[] gsP, double[] hsP, double[] jsM, double[] gsM, double[] hsM)
            => MaxAbsDiff(rr * Superop(jsP, gsP, hsP) * rr, Superop(jsM, gsM, hsM));
        cases.Add(DevCase("conjugation identity, gamma axis",
            "(R⊗R)·L(base + t·dir_γ)·(R⊗R) = L(base − t·dir_γ), entry-wise on the 64×64 superoperator",
            DevConj(jBase, Scan(gBase, gOdd, +t), hBase, jBase, Scan(gBase, gOdd, -t), hBase)));
        cases.Add(DevCase("conjugation identity, J axis",
            "same identity with the R-odd bond direction (bond b mirrors to N−2−b)",
            DevConj(Scan(jBase, jOdd, +t), gBase, hBase, Scan(jBase, jOdd, -t), gBase, hBase)));
        cases.Add(DevCase("conjugation identity, h axis",
            "same identity with the R-odd longitudinal detuning direction",
            DevConj(jBase, gBase, Scan(hBase, hOdd, +t), jBase, gBase, Scan(hBase, hOdd, -t))));

        // (b) The fence: a direction that does not negate EVERY component is rejected at O(1).
        var broken = (double[])gOdd.Clone();
        broken[0] = 0.0;
        double devFence = DevConj(jBase, Scan(gBase, broken, +t), hBase, jBase, Scan(gBase, broken, -t), hBase);
        cases.Add(new BatteryCase(
            Name: "the mixed-scan fence",
            Detail: "zero one component of the odd gamma direction: no longer R-odd, the identity must fail at O(1)",
            Expected: "dev > 1e-6",
            Actual: devFence > 1e-6 ? "dev > 1e-6" : "dev = " + devFence.ToString("E2", CultureInfo.InvariantCulture)));

        // (c) The σ_eff = +1 column is static: an R-even direction commutes with the mirror.
        double[] gEvenDir() => new[] { 0.1, 0.0, 0.1 };
        double devStatic = MaxAbsDiff(
            rr * Superop(jBase, Scan(gBase, gEvenDir(), +t), hBase) * rr,
            Superop(jBase, Scan(gBase, gEvenDir(), +t), hBase));
        cases.Add(DevCase("sigma_eff = +1 column is static",
            "an R-even gamma direction: (R⊗R)·L(t)·(R⊗R) = L(t) at the SAME t (no reflection, a static selection rule)",
            devStatic));

        // (d) The preparations as operator identities: ρ_even is R-even, the leak admixture R-odd.
        int u = 1, v = 1 << (n - 1);
        var rhoEven = ComplexVector.Build.Dense(d2);
        rhoEven[u * d + u] = new Complex(0.5, 0);
        rhoEven[v * d + v] = new Complex(0.5, 0);
        rhoEven[u * d + v] = new Complex(0.5, 0);
        rhoEven[v * d + u] = new Complex(0.5, 0);
        var rhoOdd = ComplexVector.Build.Dense(d2);
        rhoOdd[u * d + u] = new Complex(0.5, 0);
        rhoOdd[v * d + v] = new Complex(-0.5, 0);
        double devPrep = Math.Max(
            (rr * rhoEven - rhoEven).AbsoluteMaximum().Magnitude,
            (rr * rhoOdd + rhoOdd).AbsoluteMaximum().Magnitude);
        cases.Add(DevCase("preparation parities",
            "ρ_even = |ψ⟩⟨ψ| with ψ = (|e₀⟩ + |e_{N−1}⟩)/√2 is operator-R-even; the leak admixture " +
            "(|e₀⟩⟨e₀| − |e_{N−1}⟩⟨e_{N−1}|)/2 is operator-R-odd",
            devPrep));

        // --- moments: Tr[O·L^k·ρ₀], k = 0..6 (all-times statement by analyticity) ---
        var oEven = zOps[0] + zOps[n - 1];
        var oOdd = zOps[0] - zOps[n - 1];
        double[] Moments(ComplexMatrix sup, ComplexVector rho0, ComplexMatrix readout, int kMax)
        {
            var moments = new double[kMax + 1];
            var w = rho0.Clone();
            for (int k = 0; k <= kMax; k++)
            {
                Complex tr = Complex.Zero;
                for (int i = 0; i < d; i++)
                    for (int j = 0; j < d; j++)
                        tr += readout[i, j] * w[j * d + i];
                moments[k] = tr.Real;
                if (k < kMax) w = sup * w;
            }
            return moments;
        }

        const int kMax = 6;
        var supPlus = Superop(jBase, Scan(gBase, gOdd, +t), hBase);
        var supMinus = Superop(jBase, Scan(gBase, gOdd, -t), hBase);
        var mePlus = Moments(supPlus, rhoEven, oEven, kMax);
        var meMinus = Moments(supMinus, rhoEven, oEven, kMax);
        var moPlus = Moments(supPlus, rhoEven, oOdd, kMax);
        var moMinus = Moments(supMinus, rhoEven, oOdd, kMax);

        // (e) The even cell: every moment of the q = +1 readout is EVEN in t.
        double devEven = 0;
        for (int k = 0; k <= kMax; k++)
            devEven = Math.Max(devEven, Math.Abs(mePlus[k] - meMinus[k]) / Math.Max(1.0, Math.Abs(mePlus[k])));
        cases.Add(DevCase("even cell (+,−): moments even in t",
            "Tr[O_even·L(+t)^k·ρ_even] = Tr[O_even·L(−t)^k·ρ_even] for k = 0..6 (relative dev)",
            devEven, tol: 1e-10));

        // (f) The odd cell: every moment of the q = −1 readout is ODD in t, and non-vacuously so.
        double devOdd = 0, oddMag = 0;
        for (int k = 0; k <= kMax; k++)
        {
            devOdd = Math.Max(devOdd, Math.Abs(moPlus[k] + moMinus[k]) / Math.Max(1.0, Math.Abs(moPlus[k])));
            oddMag = Math.Max(oddMag, Math.Abs(moPlus[k]));
        }
        bool oddOk = devOdd <= 1e-10 && oddMag > 1e-3;
        cases.Add(new BatteryCase(
            Name: "odd cell (−,−): moments odd in t, non-vacuous",
            Detail: "Tr[O_odd·L(+t)^k·ρ_even] = −Tr[O_odd·L(−t)^k·ρ_even] for k = 0..6, with at least one moment O(1)",
            Expected: "dev ≤ 1e-10 and max |moment| > 1e-3",
            Actual: oddOk
                ? "dev ≤ 1e-10 and max |moment| > 1e-3"
                : "dev = " + devOdd.ToString("E2", CultureInfo.InvariantCulture) +
                  ", max |moment| = " + oddMag.ToString("E2", CultureInfo.InvariantCulture)));

        // (g) The zero cell: R-even generator, R-even prep, R-odd readout: every moment vanishes;
        //     and the generic cell responds (an R-even direction, q = +1 readout, +t vs −t differ).
        var supBase = Superop(jBase, gBase, hBase);
        var moZero = Moments(supBase, rhoEven, oOdd, kMax);
        double devZero = moZero.Max(Math.Abs);
        var supGenP = Superop(jBase, Scan(gBase, gEvenDir(), +t), hBase);
        var supGenM = Superop(jBase, Scan(gBase, gEvenDir(), -t), hBase);
        var meGenP = Moments(supGenP, rhoEven, oEven, kMax);
        var meGenM = Moments(supGenM, rhoEven, oEven, kMax);
        double genericGap = 0;
        for (int k = 0; k <= kMax; k++)
            genericGap = Math.Max(genericGap, Math.Abs(meGenP[k] - meGenM[k]));
        bool zeroOk = devZero <= Tol && genericGap > 1e-3;
        cases.Add(new BatteryCase(
            Name: "zero cell (−,+) vanishes; generic cell (+,+) responds",
            Detail: "R-even generator + R-even prep + R-odd readout: all moments 0 (a pure selection rule); " +
                    "the same even readout under an R-even scan direction responds at O(1)",
            Expected: "zero dev ≤ 1e-10 and generic gap > 1e-3",
            Actual: zeroOk
                ? "zero dev ≤ 1e-10 and generic gap > 1e-3"
                : "zero dev = " + devZero.ToString("E2", CultureInfo.InvariantCulture) +
                  ", generic gap = " + genericGap.ToString("E2", CultureInfo.InvariantCulture)));

        // (h) The leak: ρ = ρ_even + ε·ρ_odd opens the forbidden odd-in-t component of the even
        //     readout EXACTLY affinely in ε (linearity of the moments in ρ₀): halving ε halves it.
        double LeakOdd(double eps, int k)
        {
            var rho = rhoEven + rhoOdd.Multiply(new Complex(eps, 0));
            var p = Moments(supPlus, rho, oEven, k);
            var m = Moments(supMinus, rho, oEven, k);
            return p[k] - m[k];
        }
        double l1 = LeakOdd(0.02, 3), l2 = LeakOdd(0.01, 3);
        bool leakOk = Math.Abs(l1) > 1e-8 && Math.Abs(l1 / l2 - 2.0) <= 1e-9;
        cases.Add(new BatteryCase(
            Name: "the leak is exactly affine in eps",
            Detail: "the odd-in-t part of the k = 3 even-readout moment under prep ρ_even + ε·ρ_odd " +
                    "(k = 3 is the lowest order carrying the scan to a diagonal readout: the dissipator " +
                    "acts only off-diagonal, so the first reading chain is L_H·D·L_H); " +
                    "moments are linear in ρ₀, so halving ε halves it exactly",
            Expected: "leak nonzero and ratio = 2 to 1e-9",
            Actual: leakOk
                ? "leak nonzero and ratio = 2 to 1e-9"
                : "leak = " + l1.ToString("E2", CultureInfo.InvariantCulture) +
                  ", ratio = " + (l1 / l2).ToString("F9", CultureInfo.InvariantCulture)));

        return cases;
    }

    private static BatteryCase DevCase(string name, string detail, double dev, double tol = 1e-12)
    {
        string expected = "dev ≤ " + tol.ToString("E0", CultureInfo.InvariantCulture);
        return new BatteryCase(
            Name: name,
            Detail: detail,
            Expected: expected,
            Actual: dev <= tol ? expected : "dev = " + dev.ToString("E2", CultureInfo.InvariantCulture));
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
