using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The dead-set law (F132, Tier1Derived, minted 2026-07-16): on the open XY chain
/// with a longitudinal field h_l·Z_l under local Z-dephasing (any watching profile), for
/// each sublattice gauge U_g (Z on the even or on the odd sites) the antiunitary
///
/// <code>
///   V_g(ρ) = W_g · conj(ρ) · W_g†,   W_g = U_g·X^N
/// </code>
///
/// is an EXACT symmetry of the Lindblad flow at every fixed h: W_g flips H(h) to −H(h) as
/// an operator (the gauge flips the hopping, X^N flips the field), conjugation is
/// antiunitary (H real, so conj maps the H-world's trajectories to the −H-world's), and
/// the dephasing mask, real and a function of i⊕j only, rides along. Two h-flipping
/// mirrors (F131's third mirror X^N; the conj∘chiral partner) compose to an h-preserving
/// one. Its kill sign on a Pauli readout O is a pure, N-free function of the left
/// Jordan-Wigner Majorana degree d(O):
///
/// <code>
///   ε_odd-gauge(O) = (−1)^(d(d−1)/2),   ε_even-gauge(O) = (−1)^(d(d+1)/2)
/// </code>
///
/// (the per-site sign of W_g∘conj is uniform on the chain; the Hermitian phase
/// i^(d(d−1)/2) conjugates to a sign). The prep family ρ(0) = (P_s + P_~s)/2 (+ optional
/// cross coherence) splits into V-eigen-sectors, population +1 under both gauges,
/// coherence (−1)^|g|; a sector's contribution to ⟨O⟩ dies unless ε_g matches its sign
/// for BOTH gauges; at d = N the coherence signs hold automatically (⌊N/2⌋, ⌈N/2⌉ are
/// exactly the sublattice sizes). With the popcount blocks and the conserved degree
/// (H quadratic, dissipator diagonal on Pauli strings) the identically-zero-at-every-h
/// set closes into one line:
///
/// <code>
///   alive ⟺ (mask connects a populated diagonal block ∧ d ≡ 0 mod 4)
///           ∨ (coherence on ∧ mask connects the coherence blocks ∧ d = N).
/// </code>
///
/// <para><b>Tiers, honestly split:</b> the necessity direction (forbidden ⟹ identically
/// zero) is derived; the sufficiency direction (allowed ⟹ alive) is the gated
/// observation (twelve full 4^N−1 censuses, N = 3..6, plus the two zz controls in the
/// committed gate). The per-sector form is load-bearing: at N ≡ 2 (mod 4) with coherence
/// no V_g stabilizes ρ(0) globally and the global-stabilizer shortcut goes blind (N = 6:
/// 2047 predicted vs 1055 actual; the collapsed form stays exact). <b>Boundary as
/// content:</b> a ZZ coupling kills BOTH non-kinematic ingredients (the chiral gauge no
/// longer flips H; the degree grading collapses to parity), their kills revive, and the
/// lone fermionic survivor is the conserved parity Z^N: a law of the FREE world, with
/// the interaction knob as its gated discriminator.</para>
///
/// <para><b>Self-check battery (N = 3):</b> built in the constructor, matrix and moment
/// level, no eigensolver and no matrix exponential. Eight cases: the composed mirror
/// flips H (both gauges, machine zero); V_g commutes with the Lindblad right-hand side
/// on a generic dense matrix (the flow symmetry, both gauges, machine zero); the mod-4
/// identity on ALL 64 Pauli strings three ways at once (letter formula == degree formula
/// == direct matrix conjugation); the automatic coherence signs as arithmetic through
/// N = 12; degree conservation of the adjoint flow ([H, O] lands only on strings of the
/// same degree, the dissipator is diagonal on strings); the FULL 63-string moment census
/// against the collapsed rule for BOTH preps (population 15/63, coherence 35/63,
/// first-firing k ≤ 6, relative threshold); and the zz fence (with interaction the gauge
/// no longer flips H, O(1) rejection). Mirrors the committed gate
/// simulations/lattice_dead_set_rule.py (21 checks, which adds N = 4..6, the N = 6
/// divergence, and the binomial anatomy).</para>
///
/// <para><b>Layer note:</b> cross-axis structural, deliberately NOT an
/// <see cref="IZ2AxisClaim"/>. Typed parents: <see cref="ChiralKClaim"/> (the sublattice
/// gauge K·H·K = −H, PROOF_K_PARTNERSHIP; here U_g at the many-body level),
/// <see cref="AntilinearTriangleClaim"/> (conj as an antilinear leg; the transport that
/// flips the exponent's i), and <see cref="MirrorOrderSortingClaim"/> (F131: X^N is the
/// third sighted mirror and its h-axis scan opened the doubly-mirrored zeros this law
/// closes). No MirrorWorld adoption yet (second-consumer stance); no hardware
/// witness.</para></summary>
public sealed class DeadSetLawClaim : Claim
{
    private const int BatteryN = 3;

    /// <summary>One self-check tying the claim to the dead-set identities at N = 3.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: the sublattice gauge (K·H·K = −H for NN hopping,
    /// PROOF_K_PARTNERSHIP). U_g is exactly this gauge at the many-body level; it is the
    /// ingredient that flips the hopping back after conjugation.</summary>
    public ChiralKClaim ChiralK { get; }

    /// <summary>Typed parent: the antilinear triangle. Entrywise conjugation is its
    /// antilinear leg; that it maps the H-world to the −H-world (H real) is the
    /// triangle's transport read on the Lindblad chain.</summary>
    public AntilinearTriangleClaim Triangle { get; }

    /// <summary>Typed parent: F131, the order-sorting law. X^N is its third sighted
    /// mirror (σ_op = −1 on the h axis); the doubly-mirrored zeros of that sighting are
    /// the w = 0 cell of this law, and the composition of the two h-flipping mirrors
    /// into the h-preserving V_g is the move that closes the rest.</summary>
    public MirrorOrderSortingClaim OrderSorting { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public DeadSetLawClaim(
        ChiralKClaim chiralK,
        AntilinearTriangleClaim triangle,
        MirrorOrderSortingClaim orderSorting)
        : base("F132, the dead-set law: the mirror-composition antiunitaries V_g = (U_g·X^N)∘conj are " +
               "exact symmetries of the XY+field Lindblad flow at every fixed h; their kill sign is the " +
               "N-free mod-4 function of the Majorana degree, ε_odd = (−1)^(d(d−1)/2), " +
               "ε_even = (−1)^(d(d+1)/2); the prep splits into V-eigen-sectors (population +1, coherence " +
               "(−1)^|g|, automatic at d = N), and the identically-zero-at-every-h readout set closes into " +
               "one line: alive ⟺ (K_pop ∧ d ≡ 0 mod 4) ∨ (coherence ∧ K_coh ∧ d = N). Necessity derived; " +
               "sufficiency gated (twelve censuses + two zz controls); the zz knob kills V and F and is the " +
               "free-world fence. Tier1Derived (the necessity face; every ingredient owned)",
               Tier.Tier1Derived,
               "experiments/LATTICE_DEAD_SET_RULE.md (owner: the three-layer discovery + the collapse) + " +
               "simulations/lattice_dead_set_rule.py (gate, 21 checks incl. the N = 6 divergence " +
               "2047-vs-1055 and the binomial anatomy) + " +
               "docs/ANALYTICAL_FORMULAS.md (F132, minted 2026-07-16)")
    {
        ChiralK = chiralK ?? throw new ArgumentNullException(nameof(chiralK));
        Triangle = triangle ?? throw new ArgumentNullException(nameof(triangle));
        OrderSorting = orderSorting ?? throw new ArgumentNullException(nameof(orderSorting));
        Cases = BuildBattery();
    }

    /// <summary>The mod-4 kill sign: ε_odd-gauge = (−1)^(d(d−1)/2), ε_even-gauge =
    /// (−1)^(d(d+1)/2). N-free; the odd-sites gauge kills d ≡ 2, 3 (mod 4), the
    /// even-sites gauge kills d ≡ 1, 2 (mod 4).</summary>
    public static int EpsDegree(int d, bool evenGauge)
    {
        long half = evenGauge ? (long)d * (d + 1) / 2 : (long)d * (d - 1) / 2;
        return (half & 1) == 0 ? +1 : -1;
    }

    /// <summary>Left Jordan-Wigner Majorana degree of a Pauli string (a Majorana at site
    /// m puts Z on every site left of m), computed right-to-left via the tail parity:
    /// X/Y contribute one cap; Z is a pair when the tail above is even, absorbed when
    /// odd; I is empty when the tail is even, a pair when odd.</summary>
    public static int MajoranaDegree(string name)
    {
        int tau = 0, d = 0;
        for (int l = name.Length - 1; l >= 0; l--)
        {
            char ch = name[l];
            int n = ch is 'X' or 'Y' ? 1
                : ch == 'Z' ? (tau % 2 == 0 ? 2 : 0)
                : (tau % 2 == 0 ? 0 : 2);
            d += n;
            tau += n;
        }
        return d;
    }

    /// <summary>The law in one line.</summary>
    public string Law =>
        "alive ⟺ (mask connects a populated diagonal block ∧ d ≡ 0 mod 4) ∨ (coherence on ∧ mask " +
        "connects the coherence blocks ∧ d = N); necessity derived, sufficiency the gated observation.";

    /// <summary>The symmetry in one line.</summary>
    public string Symmetry =>
        "V_g = (U_g·X^N)∘conj is an exact flow symmetry at every fixed h: the gauge flips the hopping, " +
        "X^N flips the field, conj (antiunitary, H real) maps the H-world to the −H-world, and the " +
        "dephasing mask (real, a function of i⊕j) rides along; two h-flipping mirrors compose to an " +
        "h-preserving one.";

    /// <summary>The identity in one line.</summary>
    public string Identity =>
        "ε_odd-gauge = (−1)^(d(d−1)/2), ε_even-gauge = (−1)^(d(d+1)/2), N-free: the per-site sign of " +
        "W_g∘conj is uniform on the chain and the Hermitian phase i^(d(d−1)/2) conjugates to a sign; " +
        "at d = N both signs equal (−1)^|g| automatically, so the coherence channel is never V-killed.";

    /// <summary>The boundary in one line.</summary>
    public string Boundary =>
        "A ZZ coupling kills both non-kinematic ingredients (the gauge no longer flips H; the degree " +
        "grading collapses to parity): their kills revive, only the popcount blocks and the conserved " +
        "parity Z^N remain. A law of the FREE world, with the interaction knob as its gated discriminator.";

    public override string DisplayName =>
        "F132: the dead-set law — the mirror-composition kill is the mod-4 face of the conserved fermion degree";

    public override string Summary =>
        "V_g = (U_g·X^N)∘conj is an exact flow symmetry at every fixed h and its kill sign is the N-free " +
        "mod-4 function of the Majorana degree; alive ⟺ (K_pop ∧ d ≡ 0 mod 4) ∨ (coherence ∧ K_coh ∧ " +
        $"d = N); necessity derived, sufficiency gated; {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The law (one line)", summary: Law);
            yield return new InspectableNode("The symmetry (two mirrors compose)", summary: Symmetry);
            yield return new InspectableNode("The mod-4 identity (N-free)", summary: Identity);
            yield return new InspectableNode("The boundary (the free world)", summary: Boundary);
            yield return new InspectableNode("The divergence pin",
                summary: "The per-sector form is load-bearing: at N ≡ 2 (mod 4) with coherence no V_g " +
                         "stabilizes ρ(0) globally and the global-stabilizer shortcut goes blind; the " +
                         "committed gate pins N = 6 (4095-census): 2047 predicted by the global form, " +
                         "1055 actual, the collapsed form exact. Alive-by-degree counts are binomial " +
                         "coefficients cut by kinematics (C(10,5) = 252 fully alive at N = 5, " +
                         "popcount-2 seed).");
            yield return new InspectableNode("Typed parents",
                summary: $"ChiralKClaim ({ChiralK.Tier.Label()}): the sublattice gauge, the hopping " +
                         $"flip; AntilinearTriangleClaim ({Triangle.Tier.Label()}): conj as the " +
                         $"antilinear leg; MirrorOrderSortingClaim ({OrderSorting.Tier.Label()}): F131, " +
                         "whose third-mirror sighting opened the doubly-mirrored zeros this law closes.");
            yield return new InspectableNode("Not yet",
                summary: "No hardware witness; no MirrorWorld adoption (second-consumer stance; the " +
                         "collapsed one-liner is the adoptable face).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 3, matrix + moment level, no eigensolver.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        int n = BatteryN;
        int d = 1 << n;                    // 8

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

        int Bit(int i, int l) => (i >> (n - 1 - l)) & 1;

        var zOps = Enumerable.Range(0, n).Select(l => Op(sz, l)).ToArray();

        ComplexMatrix Hamiltonian(double[] hs, double zz)
        {
            var h = ComplexMatrix.Build.Dense(d, d);
            for (int b = 0; b < n - 1; b++)
            {
                h += Op(sx, b) * Op(sx, b + 1) + Op(sy, b) * Op(sy, b + 1);
                if (zz != 0.0)
                    h += (zOps[b] * zOps[b + 1]).Multiply(new Complex(zz, 0));
            }
            for (int l = 0; l < n; l++)
                h += zOps[l].Multiply(new Complex(hs[l], 0));
            return h;
        }

        var xn = Op(sx, 0);
        for (int site = 1; site < n; site++)
            xn *= Op(sx, site);
        ComplexMatrix Gauge(bool even)
        {
            var u = ComplexMatrix.Build.DenseIdentity(d);
            for (int l = 0; l < n; l++)
                if (l % 2 == (even ? 0 : 1))
                    u *= zOps[l];
            return u;
        }
        var wEven = Gauge(true) * xn;
        var wOdd = Gauge(false) * xn;
        ComplexMatrix Vg(ComplexMatrix w, ComplexMatrix a) => w * a.Conjugate() * w.ConjugateTranspose();

        double[] gammas = { 0.3, 0.1, 0.5 };                    // full non-uniform watching
        var mask = new double[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                double m = 0;
                for (int l = 0; l < n; l++)
                    m += -2.0 * gammas[l] * (Bit(i, l) ^ Bit(j, l));
                mask[i, j] = m;
            }

        double[] hProfile = { 0.4, -0.24, 0.12 };                // 0.4·(1, −0.6, 0.3)
        var hField = Hamiltonian(hProfile, 0.0);

        ComplexMatrix Rhs(ComplexMatrix h, ComplexMatrix rho)
        {
            var outM = (h * rho - rho * h).Multiply(new Complex(0, -1));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    outM[i, j] += rho[i, j] * mask[i, j];
            return outM;
        }

        var cases = new List<BatteryCase>();

        // (1) The composed mirror flips H: W_g·H(h)·W_g† = −H(h), both gauges.
        double devFlip = Math.Max(
            MaxAbs(wEven * hField * wEven.ConjugateTranspose() + hField),
            MaxAbs(wOdd * hField * wOdd.ConjugateTranspose() + hField));
        cases.Add(DevCase("the composed mirror flips H (both gauges)",
            "W_g·H(h)·W_g† = −H(h) for W_g = U_g·X^N: the gauge flips the hopping, X^N flips the field",
            devFlip));

        // (2) V_g is a flow symmetry at fixed h: V_g(RHS(A)) = RHS(V_g(A)) on a generic
        //     dense matrix (not a state; the identity is linear), both gauges.
        var generic = ComplexMatrix.Build.Dense(d, d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                generic[i, j] = new Complex(((i * 17 + j * 29) % 13 - 6) / 10.0,
                                            ((i * 5 + j * 7) % 11 - 5) / 10.0);
        double devFlow = Math.Max(
            MaxAbsDiff(Vg(wEven, Rhs(hField, generic)), Rhs(hField, Vg(wEven, generic))),
            MaxAbsDiff(Vg(wOdd, Rhs(hField, generic)), Rhs(hField, Vg(wOdd, generic))));
        cases.Add(DevCase("V_g is a flow symmetry at fixed h (both gauges, generic matrix)",
            "V_g(−i[H,A] + mask∘A) = −i[H,V_g(A)] + mask∘V_g(A) on a dense generic A, full " +
            "non-uniform watching: prep-independent, the conserved structure itself",
            devFlow));

        // (3) The mod-4 identity, all 64 strings, three ways at once:
        //     letter formula == degree formula == direct matrix conjugation.
        var letters = new[] { 'I', 'X', 'Y', 'Z' };
        var single = new Dictionary<char, ComplexMatrix> { ['I'] = s0, ['X'] = sx, ['Y'] = sy, ['Z'] = sz };
        int mismatches = 0, strings = 0;
        foreach (var name in AllStrings(n, letters))
        {
            strings++;
            var o = Op(single[name[0]], 0);
            for (int site = 1; site < n; site++)
                o *= Op(single[name[site]], site);
            int deg = MajoranaDegree(name);
            int nz = name.Count(ch => ch == 'Z');
            foreach (bool evenGauge in new[] { true, false })
            {
                int xyG = 0;
                for (int l = 0; l < n; l++)
                    if ((name[l] == 'X' || name[l] == 'Y') && l % 2 == (evenGauge ? 0 : 1))
                        xyG++;
                int epsLetters = (nz + xyG) % 2 == 0 ? +1 : -1;
                int epsDeg = EpsDegree(deg, evenGauge);
                var w = evenGauge ? wEven : wOdd;
                double res = MaxAbsDiff(
                    (w * o * w.ConjugateTranspose()).Conjugate(),
                    o.Multiply(new Complex(epsDeg, 0)));
                if (epsLetters != epsDeg || res > 1e-12)
                    mismatches++;
            }
        }
        cases.Add(new BatteryCase(
            Name: "the mod-4 identity, all strings, three ways",
            Detail: $"conj(W_g·O·W_g†) = ε·O with ε from the letter formula (−1)^(n_Z+xy_g) AND the " +
                    $"degree formula, for all {strings} strings × both gauges",
            Expected: "0 mismatches",
            Actual: mismatches == 0 ? "0 mismatches" : mismatches + " mismatches"));

        // (4) The automatic coherence signs: at d = N both gauges' ε equal (−1)^|g|,
        //     because N(N∓1)/2 ≡ ⌊N/2⌋ / ⌈N/2⌉ (mod 2) and those are the sublattice sizes.
        int badAuto = 0;
        for (int nn = 2; nn <= 12; nn++)
        {
            int oddSites = nn / 2, evenSites = (nn + 1) / 2;
            if (EpsDegree(nn, evenGauge: false) != ((oddSites & 1) == 0 ? +1 : -1)) badAuto++;
            if (EpsDegree(nn, evenGauge: true) != ((evenSites & 1) == 0 ? +1 : -1)) badAuto++;
        }
        cases.Add(new BatteryCase(
            Name: "the automatic coherence signs (N = 2..12)",
            Detail: "at d = N: ε_odd = (−1)^⌊N/2⌋ and ε_even = (−1)^⌈N/2⌉, exactly the sublattice " +
                    "sizes, so ε_g = (−1)^|g| and the coherence channel is never V-killed",
            Expected: "0 violations",
            Actual: badAuto == 0 ? "0 violations" : badAuto + " violations"));

        // (5) The adjoint flow preserves the Majorana degree: [H, O] expands only on
        //     strings of the same degree, and every string is a dephasing eigenvector.
        var basis = new List<(string Name, ComplexMatrix O, int Deg)>();
        foreach (var name in AllStrings(n, letters))
        {
            var o = Op(single[name[0]], 0);
            for (int site = 1; site < n; site++)
                o *= Op(single[name[site]], site);
            basis.Add((name, o, MajoranaDegree(name)));
        }
        int degViolations = 0, dissViolations = 0;
        foreach (var (name, o, deg) in basis)
        {
            var comm = hField * o - o * hField;
            foreach (var (name2, p, deg2) in basis)
            {
                var coeff = (p * comm).Trace() / d;
                if (coeff.Magnitude > 1e-10 && deg2 != deg)
                    degViolations++;
            }
            foreach (var z in zOps)
            {
                var conj = z * o * z;
                if (MaxAbsDiff(conj, o) > 1e-12 && MaxAbsDiff(conj, o.Multiply(new Complex(-1, 0))) > 1e-12)
                    dissViolations++;
            }
        }
        cases.Add(new BatteryCase(
            Name: "the adjoint flow preserves the Majorana degree",
            Detail: "[H(h), O] expands only on Pauli strings of the SAME degree (H quadratic under " +
                    "left-JW, field included), and Z_l·O·Z_l = ±O for every string (the dissipator is " +
                    "diagonal on strings): the F face, both premises",
            Expected: "0 violations",
            Actual: degViolations + dissViolations == 0
                ? "0 violations"
                : degViolations + dissViolations + " violations"));

        // (6)+(7) The full 63-string moment census against the collapsed rule, both preps.
        bool BlockReadable(int w, IEnumerable<(int P, int Q)> support)
        {
            foreach (var (p, q) in support)
            {
                if ((p + w - q) % 2 != 0) continue;
                int a = (p + w - q) / 2;
                if (a >= 0 && a <= w && a <= p && p - a <= n - w)
                    return true;
            }
            return false;
        }
        var supPop = new[] { (1, 1), (2, 2) };                   // ps = 1, pt = N − 1 = 2
        var supCoh = new[] { (1, 2), (2, 1) };

        foreach (bool coh in new[] { false, true })
        {
            var rho0 = ComplexMatrix.Build.Dense(d, d);
            rho0[1, 1] = new Complex(0.5, 0);
            rho0[6, 6] = new Complex(0.5, 0);
            if (coh)
            {
                rho0[1, 6] = new Complex(0.3, 0);
                rho0[6, 1] = new Complex(0.3, 0);
            }

            int misPredicted = 0, aliveCount = 0, predCount = 0;
            foreach (var (name, o, deg) in basis)
            {
                if (name == new string('I', n)) continue;
                int w = name.Count(ch => ch is 'X' or 'Y');
                bool pred = (BlockReadable(w, supPop) && deg % 4 == 0)
                            || (coh && BlockReadable(w, supCoh) && deg == n);

                bool fired = false;
                double oNorm = o.FrobeniusNorm();
                var wk = rho0;
                for (int k = 0; k <= 6 && !fired; k++)
                {
                    Complex tr = Complex.Zero;
                    for (int i = 0; i < d; i++)
                        for (int j = 0; j < d; j++)
                            tr += o[i, j] * wk[j, i];
                    if (tr.Magnitude > 1e-8 * Math.Max(1.0, wk.FrobeniusNorm()) * oNorm)
                        fired = true;
                    if (k < 6) wk = Rhs(hField, wk);
                }
                if (fired) aliveCount++;
                if (pred) predCount++;
                if (fired != pred) misPredicted++;
            }
            int expectAlive = coh ? 35 : 15;
            bool ok = misPredicted == 0 && aliveCount == expectAlive && predCount == expectAlive;
            cases.Add(new BatteryCase(
                Name: coh
                    ? "the collapsed rule: full moment census, coherence prep"
                    : "the collapsed rule: full moment census, population prep",
                Detail: $"all 63 readouts, moments k = 0..6 (relative threshold; the scout showed every " +
                        $"alive string fires by k = 4): predicted set == fired set, {expectAlive}/63 alive",
                Expected: $"0 mispredictions, {expectAlive} alive",
                Actual: ok
                    ? $"0 mispredictions, {expectAlive} alive"
                    : $"{misPredicted} mispredictions, {aliveCount} alive (predicted {predCount})"));
        }

        // (8) The zz fence: with an interaction the gauge no longer flips H (O(1) rejection).
        var hZz = Hamiltonian(hProfile, 0.7);
        double devZz = Math.Min(
            MaxAbs(wEven * hZz * wEven.ConjugateTranspose() + hZz),
            MaxAbs(wOdd * hZz * wOdd.ConjugateTranspose() + hZz));
        cases.Add(new BatteryCase(
            Name: "the zz fence: the free-world scope",
            Detail: "with zz = 0.7 the chiral gauge flips the hopping but NOT the interaction, so " +
                    "W_g·H·W_g† ≠ −H at O(1) for both gauges: the law's V and F faces are free-world only",
            Expected: "dev > 0.1",
            Actual: devZz > 0.1 ? "dev > 0.1" : "dev = " + devZz.ToString("E2", CultureInfo.InvariantCulture)));

        return cases;
    }

    private static IEnumerable<string> AllStrings(int n, char[] letters)
    {
        var name = new char[n];
        for (int code = 0; code < 1 << (2 * n); code++)
        {
            int c = code;
            for (int l = 0; l < n; l++)
            {
                name[l] = letters[c & 3];
                c >>= 2;
            }
            yield return new string(name);
        }
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

    private static double MaxAbs(ComplexMatrix a)
    {
        double m = 0.0;
        for (int i = 0; i < a.RowCount; i++)
            for (int j = 0; j < a.ColumnCount; j++)
            {
                double v = a[i, j].Magnitude;
                if (v > m) m = v;
            }
        return m;
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
