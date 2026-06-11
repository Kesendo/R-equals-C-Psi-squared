using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The antilinear triangle (Tier1Derived, 2026-06-11): the three involutions of
/// operator space, the transpose θ(A) = Aᵀ, the entrywise conjugation conj(A) = Ā, and the
/// adjoint †(A) = A†, close with id into one Klein four-group, † = θ∘conj = conj∘θ, graded
/// by two characters: linearity ℓ (θ is ℂ-linear; conj and † are antilinear) and
/// multiplicativity m (conj is an automorphism; θ and † are antiautomorphisms). On a Pauli
/// string everything collapses to one number: θ(σ) = conj(σ) = (−1)^{n_Y}·σ and †(σ) = σ,
/// because Y is the only Pauli that is both antisymmetric and imaginary.
///
/// <para><b>The engine (transport law):</b> for every H, Hermitian or not, and every
/// μ ∈ {id, θ, conj, †},</para>
/// <code>
///   μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)},   L_H = −i[H, ·]
/// </code>
/// <para>the sign being the product of the two characters: one sign from the −i (paid only
/// by antilinear maps), one from the commutator's order (paid only by reversing maps); θ and
/// conj each pay exactly one, † pays both and so pays nothing. The four instances:
/// id → +L_H, θ → −L_{Hᵀ}, conj → −L_{H̄}, † → +L_{H†}.</para>
///
/// <para><b>Five legs, one engine.</b> Five independently derived repository proofs are five
/// faces of the one line: (1) F114's sign law D·L_σ·D = (−1)^{n_Y+1}·L_σ
/// (<see cref="CommutatorDConjugationSign"/>) is θ on a commutator; (2) the F87 girth-ladder
/// reversal kill Tr(reversed word) = (−1)^{n_Y}·Tr(word) is θ at word length j; (3) F112's
/// Lemmas A+B (<see cref="LindbladBitBPiBalance"/>), the HS adjoint (L_H)* = −L_{H†}, skew
/// for Hermitian H, and the dagger conjugating Π-eigenvalues, are the † face at the
/// Hilbert-Schmidt pairing; (4) F113 Lemma C / the F117 Hermitian conjugacy
/// (H = H† ⟺ Hᵀ = H̄, the ket leg is the entrywise conjugate of the bra leg) is the
/// fixed-point collapse, θ = conj on Hermitian operators; (5) the F86/PTF K_b mode mirror
/// T = Σ₁∘conj is the conj face, dressed by the K₁ shadow.</para>
///
/// <para><b>The antilinear double:</b> in the Pauli basis θ IS the D₄ mirror
/// D = diag((−1)^{n_Y}) of <see cref="MirrorGroupD4Claim"/>, † is the antilinear unit 𝒦
/// (pure coefficient conjugation), and conj = D∘𝒦. The closure
/// ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ has order 16 with eight antiunitary members: the antilinear double
/// of the mirror group, every mirror acquiring an antiunitary twin. The two ℓ·m = −1
/// vertices θ and conj invert every dial Ad_{R_z(φ)} ↦ Ad_{R_z(−φ)}; † commutes with every
/// unitary conjugation.</para>
///
/// <para><b>Layer note:</b> the conj-leg's typed claim (ChiralMirrorTrajectoryClaim, the
/// K₁ trajectory identity) lives in RCPsiSquared.Diagnostics, which Core cannot reference;
/// that edge is carried in prose and through the PROOF_ANTILINEAR_TRIANGLE.md §3.4 anchor,
/// not as a typed parent. Like <see cref="MirrorGroupD4Claim"/>, this claim is cross-axis
/// structural and deliberately does NOT implement <see cref="IZ2AxisClaim"/>.</para>
///
/// <para><b>Self-check battery:</b> built in the constructor at N = 2 (4×4 operators and
/// 16×16 superoperators, &lt; 100 ms, machine-exact to 1e-12): the full V₄ composition
/// table on a fixed operator, the Pauli action on all 16 strings, the transport-law sign
/// for all four maps on a fixed non-Hermitian H, the HS adjoint M_H† = −M_{H†} in the
/// row-stacking vec convention (with the convention pinned against the commutator on
/// probes) plus skewness for Hermitian H, the Hermitian collapse Hᵀ = H̄ with a
/// non-Hermitian counterexample, the word-reversal kill on six fixed words covering both
/// signs, and the BFS closure |⟨R, D, 𝒦⟩| = 16 with 8 antilinear members on the 16-dim
/// Pauli coefficient space. Mirrors the blocks of
/// <c>simulations/antilinear_triangle.py</c>; the §6 qudit generalization (the Weyl-Heisenberg
/// lattice) is verified in <c>simulations/qudit_mirror_group_family.py</c>.</para></summary>
public sealed class AntilinearTriangleClaim : Claim
{
    private const double Tol = 1e-12;
    private const int BatteryN = 2;

    /// <summary>One self-check tying the claim to the triangle identities at N = 2.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: the mirror group D₄ = ⟨R, D⟩. The triangle docks onto it in
    /// the Pauli basis (θ = D, the y_par axis as char(θ)); the antilinear double
    /// ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ is its order-16 thickening by the antilinear unit.</summary>
    public MirrorGroupD4Claim MirrorGroup { get; }

    /// <summary>Typed parent: F114, the θ-leg. D·L_σ·D = (−1)^{n_Y+1}·L_σ is the transport
    /// law specialized to μ = θ and H = σ: one sign from commutator reversal (m = −1), the
    /// n_Y from σᵀ = (−1)^{n_Y}·σ.</summary>
    public CommutatorDConjugationSign F114 { get; }

    /// <summary>Typed parent: F112, the †-leg. Lemma B's anti-Hermiticity (L_H)* = −L_H for
    /// Hermitian H is the transport law read through the Hilbert-Schmidt pairing; Lemma A's
    /// antilinear isometry conjugating Π-eigenvalues is the same vertex one level up, with
    /// antilinearity (ℓ = −1) the load-bearing property in both.</summary>
    public LindbladBitBPiBalance F112 { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public AntilinearTriangleClaim(
        MirrorGroupD4Claim mirrorGroup,
        CommutatorDConjugationSign f114,
        LindbladBitBPiBalance f112)
        : base("The antilinear triangle: θ (transpose), conj (entrywise conjugation), † (adjoint) " +
               "close with id into one Klein four-group († = θ∘conj), graded by linearity ℓ and " +
               "multiplicativity m; the transport law μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is the one engine " +
               "behind five proofs (F114 θ-leg, the girth-ladder reversal kill, F112 Lemmas A+B †-leg, " +
               "F113/F117 Hermitian conjugacy, the K_b mode mirror conj-leg); on Pauli strings " +
               "θ(σ) = conj(σ) = (−1)^{n_Y}·σ and †(σ) = σ; in the Pauli basis θ = D, † = 𝒦, and " +
               "⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16), the antilinear double of the mirror group. " +
               "Tier1Derived (one-line algebraic identities, exact)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md + " +
               "simulations/antilinear_triangle.py + " +
               "simulations/qudit_mirror_group_family.py")
    {
        MirrorGroup = mirrorGroup ?? throw new ArgumentNullException(nameof(mirrorGroup));
        F114 = f114 ?? throw new ArgumentNullException(nameof(f114));
        F112 = f112 ?? throw new ArgumentNullException(nameof(f112));
        Cases = BuildBattery();
    }

    /// <summary>The triangle in one line.</summary>
    public string Triangle =>
        "θ(A) = Aᵀ, conj(A) = Ā, †(A) = A†: three involutions, any two compose to the third; " +
        "{id, θ, conj, †} ≅ V₄, graded by linearity ℓ (θ linear; conj, † antilinear) and " +
        "multiplicativity m (conj an automorphism; θ, † antiautomorphisms). On a Pauli string " +
        "everything collapses to one number: θ(σ) = conj(σ) = (−1)^{n_Y}·σ and †(σ) = σ.";

    /// <summary>The transport law (the engine) in one line.</summary>
    public string TransportLaw =>
        "μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)} for every H and every μ ∈ {id, θ, conj, †}, with " +
        "L_H = −i[H, ·]. One sign from the −i (paid only by antilinear maps), one from the " +
        "commutator's order (paid only by reversing maps): id → +L_H, θ → −L_{Hᵀ}, " +
        "conj → −L_{H̄}, † → +L_{H†}.";

    /// <summary>The five legs in one line.</summary>
    public string FiveLegs =>
        "F114 D·L_σ·D = (−1)^{n_Y+1}·L_σ is θ on a commutator; the girth-ladder reversal kill " +
        "Tr(reversed word) = (−1)^{n_Y}·Tr(word) is θ at word length j; F112 Lemmas A+B " +
        "((L_H)* = −L_{H†}, skew for Hermitian H, dagger conjugates Π-eigenvalues) are the † face " +
        "at the Hilbert-Schmidt pairing; F113 Lemma C / F117 Hermitian conjugacy (H = H† ⟺ Hᵀ = H̄, " +
        "ket leg = conj of bra leg) is the fixed-point collapse θ = conj on Hermitian operators; " +
        "the F86/PTF K_b mode mirror T = Σ₁∘conj is the conj face, dressed.";

    /// <summary>The antilinear double in one line.</summary>
    public string AntilinearDouble =>
        "In the Pauli basis θ IS the D₄ mirror D = diag((−1)^{n_Y}), † is the antilinear unit 𝒦 " +
        "(pure coefficient conjugation), conj = D∘𝒦; the closure ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ of order 16 " +
        "with 8 antiunitary members is the antilinear double of the mirror group: every mirror " +
        "acquires an antiunitary twin.";

    /// <summary>The qudit generalization (Weyl-Heisenberg) in one line.</summary>
    public string QuditGeneralization =>
        "The triangle is the d = 2 shadow of a Weyl-Heisenberg lattice action. On the qudit " +
        "operators P_{a,b} = X^a Z^b (a, b ∈ Z_d; clock X, phase Z, ZX = ωXZ, ω = e^{2πi/d}) the " +
        "involutions act with a symplectic phase: θ(P_{a,b}) = ω^{−ab}P_{−a,b}, " +
        "conj(P_{a,b}) = P_{a,−b}, †(P_{a,b}) = ω^{ab}P_{−a,−b} (verified d = 2..5), and the " +
        "transport law is basis-free at every d. The qubit (−1)^{n_Y} is the degeneration: at " +
        "d = 2 the flip a ↦ −a is trivial and ω^{ab} = −1 only on (1,1) = Y. For d > 2 the " +
        "involutions move the labels, so the sign becomes a reflection of the Z_d × Z_d lattice; " +
        "⟨Π_d, D⟩ (F121's Z_d ≀ Z₂) is that lattice's symmetry group. The antilinear double over " +
        "general d and the d → ∞ rotation circle are the open arc. " +
        "Verified in simulations/qudit_mirror_group_family.py.";

    public override string DisplayName =>
        "The antilinear triangle: θ, conj, † as one Klein four-group; five proofs, one engine";

    public override string Summary =>
        "the three involutions θ, conj, † of operator space form with id a Klein four-group graded by " +
        "linearity ℓ and multiplicativity m, and the transport law μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is the " +
        "one engine behind five independently derived proofs (F114, the reversal kill, F112 Lemmas A+B, " +
        "F113/F117 Hermitian conjugacy, the K_b mode mirror); in the Pauli basis θ = D and † = 𝒦, and " +
        "⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16) is the antilinear double of the mirror group; " +
        $"{PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The triangle (V₄ and the two characters)", summary: Triangle);
            yield return new InspectableNode("The transport law (the engine)", summary: TransportLaw);
            yield return new InspectableNode("Five legs, one engine", summary: FiveLegs);
            yield return new InspectableNode("The antilinear double ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂", summary: AntilinearDouble);
            yield return new InspectableNode("The qudit generalization (Weyl-Heisenberg lattice)", summary: QuditGeneralization);
            yield return new InspectableNode("Conj-leg claim (prose edge only)",
                summary: "The conj-leg's typed claim (ChiralMirrorTrajectoryClaim, the K₁ trajectory " +
                         "identity behind the K_b mode mirror T = Σ₁∘conj) lives in " +
                         "RCPsiSquared.Diagnostics; Core cannot reference it, so the edge is carried by " +
                         "the PROOF_ANTILINEAR_TRIANGLE.md §3.4 anchor instead of a ctor parent.");
            yield return new InspectableNode("No IZ2AxisClaim",
                summary: "The triangle is a cross-axis structural identity (the three involutions act " +
                         "on all of operator space; only θ's restriction to Pauli strings reads the " +
                         "y_par axis, and that reading is already typed on MirrorGroupD4Claim). Like " +
                         "MirrorGroupD4Claim, this claim does not sit on a single Z₂ axis.");
            yield return new InspectableNode("Typed parents",
                summary: $"MirrorGroupD4Claim ({MirrorGroup.Tier.Label()}): the D₄ the triangle docks " +
                         $"onto (θ = D) and doubles; CommutatorDConjugationSign ({F114.Tier.Label()}): " +
                         $"F114, the θ-leg; LindbladBitBPiBalance ({F112.Tier.Label()}): F112, whose " +
                         "Lemmas A+B are the †-leg at the Hilbert-Schmidt pairing.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 2 dense operators and superoperators.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        int d = 1 << BatteryN;     // 4
        int d2 = d * d;            // 16 coherence dimensions = 4^N Pauli strings

        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var minusI = new Complex(0.0, -1.0);

        // Fixed deterministic operators (LCG-filled; no RNG dependency, run-to-run stable).
        var a0 = FixedMatrix(seed: 1, dim: d);
        var hNh = FixedMatrix(seed: 2, dim: d);          // non-Hermitian probe Hamiltonian
        var hHermRaw = FixedMatrix(seed: 6, dim: d);
        var hHerm = (hHermRaw + hHermRaw.ConjugateTranspose()).Multiply(new Complex(0.5, 0.0));
        var probes = new[] { FixedMatrix(3, d), FixedMatrix(4, d), FixedMatrix(5, d) };

        ComplexMatrix Theta(ComplexMatrix x) => x.Transpose();
        ComplexMatrix Conj(ComplexMatrix x) => x.Conjugate();
        ComplexMatrix Dag(ComplexMatrix x) => x.ConjugateTranspose();
        ComplexMatrix Lh(ComplexMatrix h, ComplexMatrix rho) => (h * rho - rho * h).Multiply(minusI);

        var cases = new List<BatteryCase>();

        // (a) The full V₄ composition table on a fixed operator.
        double devTriangle = 0.0;
        void Acc(ComplexMatrix x, ComplexMatrix y) => devTriangle = Math.Max(devTriangle, MaxAbsDiff(x, y));
        Acc(Dag(a0), Theta(Conj(a0)));
        Acc(Dag(a0), Conj(Theta(a0)));
        Acc(Theta(Theta(a0)), a0);
        Acc(Conj(Conj(a0)), a0);
        Acc(Dag(Dag(a0)), a0);
        Acc(Conj(Dag(a0)), Theta(a0));
        Acc(Dag(Conj(a0)), Theta(a0));
        Acc(Theta(Dag(a0)), Conj(a0));
        Acc(Dag(Theta(a0)), Conj(a0));
        cases.Add(DevCase("† = θ∘conj (the triangle closes to V₄)",
            "full Klein table on a fixed operator: three involutions, any two compose to the third",
            devTriangle));

        // (b) Pauli action on all 16 strings at N = 2.
        int okPauli = 0;
        for (int k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, BatteryN);
            var sigma = PauliString.Build(letters);
            int nY = letters.Count(l => l == PauliLetter.Y);
            var signed = sigma.Multiply(SignOf(nY));
            bool ok = MaxAbsDiff(Theta(sigma), signed) <= Tol
                   && MaxAbsDiff(Conj(sigma), signed) <= Tol
                   && MaxAbsDiff(Dag(sigma), sigma) <= Tol;
            if (ok) okPauli++;
        }
        cases.Add(new BatteryCase(
            Name: "Pauli action: θ(σ) = conj(σ) = (−1)^{n_Y}·σ, †(σ) = σ",
            Detail: "Y is the only antisymmetric Pauli and the only imaginary one; all 16 strings at N = 2",
            Expected: "16/16",
            Actual: okPauli.ToString(CultureInfo.InvariantCulture) + "/16"));

        // (c) Transport law μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} for all four maps on non-Hermitian H.
        var transports = new (string Symbol, string Rhs, Func<ComplexMatrix, ComplexMatrix> Mu, double Sign, ComplexMatrix MuH)[]
        {
            ("id", "+L_H", x => x, +1.0, hNh),
            ("θ", "−L_{Hᵀ}", Theta, -1.0, hNh.Transpose()),
            ("conj", "−L_{H̄}", Conj, -1.0, hNh.Conjugate()),
            ("†", "+L_{H†}", Dag, +1.0, hNh.ConjugateTranspose()),
        };
        foreach (var t in transports)
        {
            double dev = 0.0;
            foreach (var rho in probes)
                dev = Math.Max(dev, MaxAbsDiff(
                    t.Mu(Lh(hNh, t.Mu(rho))),
                    Lh(t.MuH, rho).Multiply(new Complex(t.Sign, 0.0))));
            cases.Add(DevCase($"transport {t.Symbol}: {t.Symbol}∘L_H∘{t.Symbol} = {t.Rhs}",
                "the engine on a fixed non-Hermitian H over 3 fixed probe ρ; sign = ℓ(μ)·m(μ)",
                dev));
        }

        // (d) HS adjoint at the superoperator level, row-stacking vec convention
        //     (|i⟩⟨j| ↦ e_i ⊗ e_j; kron(A, B): ρ ↦ A·ρ·Bᵀ, matching MirrorGroupD4Claim).
        ComplexMatrix Sup(ComplexMatrix h) =>
            (h.KroneckerProduct(idH) - idH.KroneckerProduct(h.Transpose())).Multiply(minusI);
        double devSup = 0.0;
        foreach (var rho in probes)
        {
            // Pin the convention: vec(L_H(ρ)) = M_H·vec(ρ).
            var vec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d2);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    vec[i * d + j] = rho[i, j];
            var w = Sup(hNh) * vec;
            var lRho = Lh(hNh, rho);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    devSup = Math.Max(devSup, (w[i * d + j] - lRho[i, j]).Magnitude);
        }
        devSup = Math.Max(devSup, MaxAbsDiff(
            Sup(hNh).ConjugateTranspose(),
            Sup(hNh.ConjugateTranspose()).Multiply(-Complex.One)));
        devSup = Math.Max(devSup, MaxAbsDiff(
            Sup(hHerm).ConjugateTranspose(),
            Sup(hHerm).Multiply(-Complex.One)));
        cases.Add(DevCase("HS adjoint: M_H† = −M_{H†}; Hermitian H ⟹ M skew",
            "M_H = −i(H⊗I − I⊗Hᵀ) on the 16×16 coherence space, convention pinned via vec(L_H ρ) = M_H·vec(ρ); " +
            "F112 Lemma B at the matrix level",
            devSup));

        // (e) Hermitian collapse: H = H† ⟺ Hᵀ = H̄, with a non-Hermitian counterexample.
        double devCollapse = MaxAbsDiff(hHerm.Transpose(), hHerm.Conjugate());
        double devCounter = MaxAbsDiff(hNh.Transpose(), hNh.Conjugate());
        bool collapseOk = devCollapse <= Tol && devCounter > 1e-6;
        cases.Add(new BatteryCase(
            Name: "Hermitian collapse: Hᵀ = H̄ exactly on Hermitian H (θ = conj there)",
            Detail: "the fixed-point collapse behind F113 Lemma C / F117; the fixed non-Hermitian H " +
                    "splits the two maps (dev = " + devCounter.ToString("E2", CultureInfo.InvariantCulture) + ")",
            Expected: "collapse dev ≤ 1e-12 and counterexample dev > 1e-6",
            Actual: collapseOk
                ? "collapse dev ≤ 1e-12 and counterexample dev > 1e-6"
                : "collapse dev = " + devCollapse.ToString("E2", CultureInfo.InvariantCulture) +
                  ", counterexample dev = " + devCounter.ToString("E2", CultureInfo.InvariantCulture)));

        // (f) Word reversal: Tr(reversed Pauli word) = (−1)^{n_Y}·Tr(word), both signs covered.
        var words = new[]
        {
            new[] { "XI", "YI", "ZI" },         // n_Y = 1 (odd, sign −1), trace ±4i
            new[] { "XX", "YY", "ZZ" },         // n_Y = 2 (even, sign +1), trace −4
            new[] { "YY", "XY", "ZI" },         // n_Y = 3 (odd, sign −1), trace ±4i
            new[] { "ZY", "XX", "YZ" },         // n_Y = 2 (even, sign +1), trace 4
            new[] { "YI", "ZI", "XI" },         // n_Y = 1 (odd, sign −1), trace ±4i
            new[] { "XI", "YI", "XI", "YI" },   // n_Y = 2 (even, sign +1), length 4, trace −4
        };
        int okWords = 0;
        foreach (var word in words)
        {
            var factors = word
                .Select(w => PauliString.Build(w.Select(PauliLetterExtensions.FromSymbol).ToArray()))
                .ToArray();
            int nY = word.Sum(w => w.Count(c => c == 'Y'));
            var fwd = factors.Aggregate((x, y) => x * y);
            var rev = factors.Reverse().Aggregate((x, y) => x * y);
            if ((rev.Trace() - SignOf(nY) * fwd.Trace()).Magnitude <= Tol) okWords++;
        }
        cases.Add(new BatteryCase(
            Name: "word reversal: Tr(reversed word) = (−1)^{n_Y}·Tr(word)",
            Detail: "the girth-ladder reversal kill = θ inside a trace; 6 fixed words at N = 2, " +
                    "three odd-n_Y (sign −1, nonzero traces) and three even-n_Y (sign +1)",
            Expected: "6/6",
            Actual: okWords.ToString(CultureInfo.InvariantCulture) + "/6"));

        // (g) The antilinear double: BFS closure of {D, R, 𝒦} on the 16-dim Pauli coefficient
        //     space, elements as (matrix, antilinear flag) with composition
        //     (M₁, a₁)∘(M₂, a₂) = (M₁·(a₁ ? conj(M₂) : M₂), a₁ ⊕ a₂).
        var basis = Matrix<Complex>.Build.Dense(d2, d2);   // Pauli-string → coherence basis
        for (int k = 0; k < d2; k++)
        {
            var sigma = PauliString.Build(PauliIndex.FromFlat(k, BatteryN));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    basis[i * d + j, k] = sigma[i, j];
        }
        var basisInv = basis.ConjugateTranspose().Divide(d);
        var f = PauliString.Build(new[] { PauliLetter.X, PauliLetter.X });
        // R = right multiplication by X^⊗2, expressed on Pauli coefficients (a phase
        // permutation with per-string phase i^{n_Z − n_Y}).
        var rp = basisInv * idH.KroneckerProduct(f) * basis;
        var dp = Matrix<Complex>.Build.Dense(d2, d2);      // D = diag((−1)^{n_Y})
        for (int k = 0; k < d2; k++)
        {
            int nY = PauliIndex.FromFlat(k, BatteryN).Count(l => l == PauliLetter.Y);
            dp[k, k] = SignOf(nY);
        }
        var id16 = Matrix<Complex>.Build.DenseIdentity(d2);

        var elems = new List<(ComplexMatrix M, bool Anti)> { (id16, false) };
        var gens = new (ComplexMatrix M, bool Anti)[] { (dp, false), (rp, false), (id16, true) };
        bool changed = true;
        while (changed)
        {
            changed = false;
            foreach (var g in gens)
            {
                foreach (var e in elems.ToList())
                {
                    var cand = (M: g.M * (g.Anti ? e.M.Conjugate() : e.M), Anti: g.Anti ^ e.Anti);
                    if (!elems.Any(x => x.Anti == cand.Anti && MaxAbsDiff(cand.M, x.M) <= Tol))
                    {
                        elems.Add(cand);
                        changed = true;
                    }
                }
            }
        }
        int antiCount = elems.Count(e => e.Anti);
        cases.Add(new BatteryCase(
            Name: "antilinear double: |⟨R, D, 𝒦⟩| = 16 with 8 antilinear (D₄ × Z₂)",
            Detail: "BFS closure of {D = diag((−1)^{n_Y}), R = right-mult by X^⊗2, 𝒦 = pure coefficient " +
                    "conjugation} on the 16-dim Pauli coefficient space",
            Expected: "16 elements, 8 antilinear",
            Actual: elems.Count.ToString(CultureInfo.InvariantCulture) + " elements, " +
                    antiCount.ToString(CultureInfo.InvariantCulture) + " antilinear"));

        return cases;
    }

    /// <summary>Deterministic pseudo-random complex matrix (splitmix-seeded LCG); fixed
    /// across runs so the battery is reproducible without an RNG dependency.</summary>
    private static ComplexMatrix FixedMatrix(ulong seed, int dim)
    {
        var m = Matrix<Complex>.Build.Dense(dim, dim);
        ulong s = seed * 0x9E3779B97F4A7C15UL + 0xD1B54A32D192ED03UL;
        double Next()
        {
            s = s * 6364136223846793005UL + 1442695040888963407UL;
            return (double)((s >> 33) % 2001UL) / 1000.0 - 1.0;   // in [−1, 1], step 1e-3
        }
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                m[i, j] = new Complex(Next(), Next());
        return m;
    }

    private static Complex SignOf(int parity) => (parity & 1) == 0 ? Complex.One : -Complex.One;

    private static BatteryCase DevCase(string name, string detail, double dev) =>
        new(Name: name,
            Detail: detail,
            Expected: "dev ≤ 1e-12",
            Actual: dev <= Tol
                ? "dev ≤ 1e-12"
                : "dev = " + dev.ToString("E2", CultureInfo.InvariantCulture));

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
