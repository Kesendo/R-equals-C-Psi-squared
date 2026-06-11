using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The moment-tower pump channel (F120, Tier1Derived, 2026-06-11): the device's
/// own amplitude damping reads the girth-ladder moment tower linearly.
///
/// <para><b>The engine (pump-slope law):</b> for the standard Lindbladian
/// L = −i[H,·] + Σ_l γ^deph_l·D[Z_l] + Σ_l γ↓_l·D[σ⁻_l] + Σ_l γ↑_l·D[σ⁺_l]
/// (σ⁻ = (X+iY)/2 = [[0,1],[0,0]], D[c](ρ) = cρc† − ½{c†c, ρ}), amplitude damping is the
/// unique non-unital piece and its pump direction is a pure local Z:
/// D[σ⁻_l](I) = +Z_l and D[σ⁺_l](I) = −Z_l, while dephasing and the unitary part
/// annihilate I. Hence for every observable A and every generator Hamiltonian,</para>
/// <code>
///   d/dt Tr(A ρ) |_{ρ = I/d} = (1/d) · Σ_l Δγ_l · Tr(A·Z_l),   Δγ_l = γ↓_l − γ↑_l,
/// </code>
/// <para>exactly. With A = H^j the right side is (1/d)·Σ_l Δγ_l·t_j(l): the slope of the
/// j-th energy moment under nothing but the chip's own damping is the girth-ladder tower
/// t_j(l) = Tr(Z_l·H^j), read linearly, rung by rung.</para>
///
/// <para><b>Three blindnesses</b>, all exact: dephasing-blind (γ^deph never enters, the
/// dephasing channel is unital); evolution-blind (the slope does not contain the
/// generator's Hamiltonian at all, only the measured polynomial A = H_p^j enters); and
/// detailed-balance closure (Δγ ≡ 0 implies slope ≡ 0, the channel is powered exactly by
/// <see cref="F84ThermalAmplitudeDampingPi2Inheritance"/>'s temperature-independent vacuum
/// component).</para>
///
/// <para><b>Rung one is F113:</b> for H = Σ_l (ω_l/2)·Z_l with σ⁻ rate γ_T1,l and σ⁺ rate
/// γ_pump,l, the static F113 polarity asymmetry
/// (<see cref="LindbladBitBPiBreakMagnitude.PredictAsymmetry"/>, a Frobenius norm imbalance
/// of the generator with no time in it) equals −4^N times the dynamic pump slope of ⟨H⟩:
/// asymmetry = −4^N·slope⟨H⟩, exactly. The static spectral-coordinate imbalance IS a
/// measurable pump rate.</para>
///
/// <para><b>The curvature fingerprint:</b> one order up,
/// d²/dt² ⟨A⟩|_{I/d} = (1/d)·Σ_l Δγ_l·Tr(A·L(Z_l)) is exactly affine in the generator (L
/// appears once, no small-δ expansion). A parasite δV in the generator shifts it by
/// (δ/d)·Σ_l Δγ_l·(−i)·Tr(V·[Z_l, A]): linear against the commutator probes [Z_l, H_p^j].
/// Z-flavored parasites are exactly invisible ([Z_l, Z-string] = 0); X/Y-flavored parasites
/// read at first order. This is the complementary channel to F113's Z-drive reader: the two
/// channels partition the single-site parasite algebra (F113's balance reads the Z-flavor,
/// the pump curvature reads the X/Y-flavor), both linearly with closed-form
/// coefficients.</para>
///
/// <para><b>The girth certificate, honestly one-sided:</b> the first rung j whose slope
/// fires is the girth ℓ (the tower below the girth is identically zero), and then
/// m* = 2ℓ+1, deg 1, hard at every γ &gt; 0 by the girth dichotomy. Silence of the deg-1
/// tower is NOT softness: the k = 4 witness IIXY+ZXZY has t_j(l) = 0 for every j and every
/// site, yet is hard at m* = 11 through its deg-5 class.</para>
///
/// <para><b>Layer note (GirthLadder):</b> the C# girth-ladder primitive
/// (<c>compute/RCPsiSquared.Diagnostics/F87/GirthLadder.cs</c>) computes the same per-site
/// tower t_ℓ(l) = Tr(Z_l·H^ℓ) and the m* = 2ℓ+1 forecast this channel reads. It is a
/// compute primitive, not a Claim, and it lives in RCPsiSquared.Diagnostics, which Core
/// cannot reference; the edge is carried here and through the
/// PROOF_MOMENT_TOWER_PUMP_CHANNEL.md anchor, not as a typed parent. Like
/// <see cref="AntilinearTriangleClaim"/>, this claim is cross-axis structural and
/// deliberately does NOT implement <see cref="IZ2AxisClaim"/>.</para>
///
/// <para><b>Self-check battery:</b> built in the constructor at N = 2 and N = 3 (dense
/// complex matrices and superoperators up to 64×64, machine-exact to 1e-12): the pump
/// directions on every site, the slope law against a dense superoperator applied to
/// vec(I/d) for j = 1, 2, 3 with site-dependent rates, the three blindnesses, the F113
/// bridge through the typed parent's <c>PredictAsymmetry</c>, the curvature fingerprint
/// with the Y₀ visible parasite and the Z₀ exact invisibility, and the girth-2 witness
/// H = X₀ + X₀Z₁ (t₁ ≡ 0, t₂ fires). Mirrors the blocks of
/// <c>simulations/moment_tower_pump_channel.py</c>.</para></summary>
public sealed class MomentTowerPumpChannelClaim : Claim
{
    private const double Tol = 1e-12;

    /// <summary>One self-check tying the claim to the pump-channel identities.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: F113, the bridge. The channel's first rung is F113's
    /// closed-form polarity asymmetry: within F113's scope,
    /// <see cref="LindbladBitBPiBreakMagnitude.PredictAsymmetry"/> = −4^N·slope⟨H⟩
    /// exactly. F113's static Frobenius imbalance is this channel's dynamic pump rate.</summary>
    public LindbladBitBPiBreakMagnitude F113 { get; }

    /// <summary>Typed parent: F84, the pump weight. The channel is powered by F84's
    /// temperature-independent vacuum rate Δγ_l = γ↓_l − γ↑_l (thermal photon-number
    /// contributions cancel); F84's detailed-balance regime γ↓ = γ↑ is blindness #3
    /// (slope ≡ 0). The pump direction D[σ⁻](I) = +Z is F82's Π²-antisymmetric entry,
    /// carried transitively through F84's mother-claim edge.</summary>
    public F84ThermalAmplitudeDampingPi2Inheritance F84 { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public MomentTowerPumpChannelClaim(
        LindbladBitBPiBreakMagnitude f113,
        F84ThermalAmplitudeDampingPi2Inheritance f84)
        : base("F120 moment-tower pump channel: amplitude damping is the unique non-unital " +
               "piece of the standard Lindbladian and pumps along pure local Z " +
               "(D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l), so " +
               "d/dt Tr(A ρ)|_{ρ=I/d} = (1/d)·Σ_l Δγ_l·Tr(A·Z_l) with Δγ_l = γ↓_l − γ↑_l; " +
               "with A = H^j the slope reads the girth-ladder tower t_j(l) = Tr(Z_l·H^j) " +
               "linearly. Dephasing-blind, evolution-blind, closes at detailed balance; " +
               "rung 1 is F113 (asymmetry = −4^N·slope⟨H⟩ exactly); the curvature is exactly " +
               "affine in the generator and reads X/Y-flavored parasites against the " +
               "commutator probes [Z_l, H_p^j] while Z-flavored parasites stay exactly " +
               "invisible (the complementary channel to F113's Z-drive reader); the girth " +
               "certificate is one-sided (a firing rung proves m* = 2ℓ+1 hard at all γ; " +
               "silence is not softness, witness IIXY+ZXZY). " +
               "Tier1Derived (one-line identities, exact)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md + " +
               "simulations/moment_tower_pump_channel.py")
    {
        F113 = f113 ?? throw new ArgumentNullException(nameof(f113));
        F84 = f84 ?? throw new ArgumentNullException(nameof(f84));
        Cases = BuildBattery();
    }

    /// <summary>The pump-slope law in one line.</summary>
    public string Theorem =>
        "d/dt Tr(A ρ)|_{ρ=I/d} = (1/d)·Σ_l Δγ_l·Tr(A·Z_l), Δγ_l = γ↓_l − γ↑_l, for every " +
        "observable A and every generator Hamiltonian; with A = H^j the slope is " +
        "(1/d)·Σ_l Δγ_l·t_j(l), the girth-ladder tower t_j(l) = Tr(Z_l·H^j) read linearly.";

    /// <summary>The pump directions in one line.</summary>
    public string PumpDirections =>
        "D[σ⁻_l](I) = +Z_l and D[σ⁺_l](I) = −Z_l (times the rate); dephasing D[Z_l](I) = 0 " +
        "and the unitary part −i[H, I] = 0 are unital. Amplitude damping is the unique " +
        "non-unital piece of the standard noise model, and its pump direction is the same " +
        "(Z_l, I) entry F82 identified as the dissipator's entire Π²-antisymmetric content.";

    /// <summary>The three blindnesses in one line.</summary>
    public string ThreeBlindnesses =>
        "Dephasing-blind: γ^deph never enters the slope (unital). Evolution-blind: the slope " +
        "is independent of the generator's Hamiltonian, only the measured polynomial " +
        "A = H_p^j enters. Detailed-balance closure: Δγ ≡ 0 implies slope ≡ 0, the channel " +
        "is powered exactly by F84's temperature-independent vacuum component.";

    /// <summary>The F113 bridge in one line.</summary>
    public string F113Bridge =>
        "For H = Σ_l (ω_l/2)·Z_l with σ⁻ rate γ_T1,l and σ⁺ rate γ_pump,l: the static F113 " +
        "polarity asymmetry (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l) equals −4^N × the dynamic " +
        "pump slope of ⟨H⟩, exactly. Rung 1 of the tower IS F113; rungs j ≥ 2 are its ladder.";

    /// <summary>The curvature fingerprint in one line.</summary>
    public string CurvatureFingerprint =>
        "d²/dt² ⟨A⟩|_{I/d} = (1/d)·Σ_l Δγ_l·Tr(A·L(Z_l)) is exactly affine in the generator; " +
        "a parasite δV shifts it by (δ/d)·Σ_l Δγ_l·(−i)·Tr(V·[Z_l, A]), linear against the " +
        "commutator probes. Z-flavored parasites are exactly invisible; X/Y-flavored ones " +
        "read at first order: the complementary channel to F113's Z-drive reader, the two " +
        "partition the single-site parasite algebra.";

    /// <summary>The one-sided girth certificate in one line.</summary>
    public string GirthCertificate =>
        "The first rung j whose slope fires is the girth ℓ, and then m* = 2ℓ+1, deg 1, hard " +
        "at every γ > 0 (positive monomial, no positive root). Silence of the deg-1 tower is " +
        "NOT softness: witness IIXY+ZXZY has t_j ≡ 0 at every rung and site yet is hard at " +
        "m* = 11 through its deg-5 class.";

    // ============================================================
    // Static helpers: predict the pump slope from the moment tower
    // ============================================================

    /// <summary>Predict the pump slope of the j-th moment of <paramref name="hMeasured"/>
    /// at the maximally mixed state:
    ///
    /// <para>  slope_j = (1/d)·Σ_l Δγ_l·t_j(l),   t_j(l) = Tr(Z_l·H^j),   d = 2^N,</para>
    ///
    /// <para>with N = <paramref name="deltaGamma"/>.Count and Δγ_l = γ↓_l − γ↑_l the
    /// per-site net damping rate. Site l is the l-th tensor factor from the left (most
    /// significant qubit first). <paramref name="hMeasured"/> is the MEASURED Hamiltonian
    /// (the polynomial read out on the device); by evolution-blindness the generator's
    /// Hamiltonian never enters. Expects Hermitian <paramref name="hMeasured"/> (the
    /// moments t_j(l) are then real; the real part is returned). The moment index
    /// <paramref name="j"/> must be ≥ 1; dimension mismatches throw
    /// <see cref="ArgumentException"/>.</para></summary>
    public static double PredictPumpSlope(
        ComplexMatrix hMeasured,
        int j,
        IReadOnlyList<double> deltaGamma)
    {
        if (hMeasured is null) throw new ArgumentNullException(nameof(hMeasured));
        if (deltaGamma is null) throw new ArgumentNullException(nameof(deltaGamma));
        int n = deltaGamma.Count;
        if (n < 1)
            throw new ArgumentException("deltaGamma must carry at least one per-site rate", nameof(deltaGamma));

        var tower = MomentTower(hMeasured, j, n);
        double sum = 0.0;
        for (int l = 0; l < n; l++)
            sum += deltaGamma[l] * tower[l];
        return sum / (1 << n);
    }

    /// <summary>The per-site moment tower t_j(l) = Tr(Z_l·H^j) for l = 0..N−1 (the same
    /// tower the Diagnostics girth-ladder primitive computes at j = ℓ). Site l is the l-th
    /// tensor factor from the left. Real parts are returned (exact for Hermitian
    /// <paramref name="h"/>). <paramref name="j"/> must be ≥ 1, <paramref name="n"/> ≥ 1,
    /// and <paramref name="h"/> must be 2^n × 2^n.</summary>
    public static IReadOnlyList<double> MomentTower(ComplexMatrix h, int j, int n)
    {
        if (h is null) throw new ArgumentNullException(nameof(h));
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 1; got {n}");
        if (j < 1) throw new ArgumentOutOfRangeException(nameof(j), $"moment index j must be ≥ 1; got {j}");
        int d = 1 << n;
        if (h.RowCount != d || h.ColumnCount != d)
            throw new ArgumentException(
                $"h is {h.RowCount}×{h.ColumnCount}; expected {d}×{d} for N = {n} sites",
                nameof(h));

        var power = MatrixPower(h, j);
        var t = new double[n];
        for (int l = 0; l < n; l++)
        {
            int bit = n - 1 - l;
            double sum = 0.0;
            for (int i = 0; i < d; i++)
                sum += ((i >> bit) & 1) == 0 ? power[i, i].Real : -power[i, i].Real;
            t[l] = sum;
        }
        return t;
    }

    public override string DisplayName =>
        "F120 moment-tower pump channel: the device's own damping reads the girth ladder linearly";

    public override string Summary =>
        "amplitude damping pumps along pure local Z (D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l), so " +
        "d/dt ⟨A⟩|_{I/d} = (1/d)·Σ_l Δγ_l·Tr(A·Z_l) reads the girth-ladder tower t_j(l) = Tr(Z_l·H^j) " +
        "linearly; dephasing-blind, evolution-blind, closed at detailed balance; rung 1 is F113 " +
        "(asymmetry = −4^N·slope⟨H⟩); the curvature is exactly affine and fingerprints X/Y-flavored " +
        "parasites while Z-flavored ones stay invisible; the girth certificate is one-sided; " +
        $"{PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem (the pump-slope law)", summary: Theorem);
            yield return new InspectableNode("Pump directions", summary: PumpDirections);
            yield return new InspectableNode("Three blindnesses", summary: ThreeBlindnesses);
            yield return new InspectableNode("Rung one is F113 (the bridge)", summary: F113Bridge);
            yield return new InspectableNode("Curvature fingerprint (the complementary channel)",
                summary: CurvatureFingerprint);
            yield return new InspectableNode("Girth certificate (one-sided)", summary: GirthCertificate);
            yield return new InspectableNode("GirthLadder primitive (prose edge only)",
                summary: "The C# girth-ladder primitive (compute/RCPsiSquared.Diagnostics/F87/" +
                         "GirthLadder.cs) computes the same tower t_ℓ(l) = Tr(Z_l·H^ℓ) and the " +
                         "m* = 2ℓ+1 forecast this channel reads. It is a compute primitive, not a " +
                         "Claim, and Diagnostics sits above Core, so the edge is carried in prose " +
                         "and through the PROOF_MOMENT_TOWER_PUMP_CHANNEL.md anchor, not as a " +
                         "typed parent.");
            yield return new InspectableNode("Typed parents",
                summary: $"LindbladBitBPiBreakMagnitude ({F113.Tier.Label()}): F113, the bridge; the " +
                         "channel's rung 1 reproduces PredictAsymmetry exactly via " +
                         $"asymmetry = −4^N·slope⟨H⟩. F84ThermalAmplitudeDampingPi2Inheritance " +
                         $"({F84.Tier.Label()}): the pump weight Δγ_l = γ↓_l − γ↑_l is F84's " +
                         "temperature-independent vacuum rate, and F84's detailed-balance regime " +
                         "is blindness #3; F82's pump direction arrives transitively through F84's " +
                         "mother-claim edge.");
            yield return new InspectableNode("No IZ2AxisClaim",
                summary: "The channel is cross-axis structural: the pump direction is a Z-flavored " +
                         "operator-space statement, but the tower it reads is the F87 girth ladder " +
                         "and the curvature reads X/Y flavors. Like AntilinearTriangleClaim and " +
                         "MirrorGroupD4Claim, this claim does not sit on a single Z₂ axis " +
                         "(cube-map counts unchanged).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    // ------------------------------------------------------------------
    // Self-check battery: N = 2 and N = 3 dense operators and superoperators.
    // ------------------------------------------------------------------

    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        // Fixed site-dependent rates at N = 2, shared across the slope-law cases.
        double[] gDeph2 = { 0.05, 0.08 };
        double[] gDown2 = { 0.021, 0.034 };
        double[] gUp2 = { 0.004, 0.011 };
        double[] delta2 = { gDown2[0] - gUp2[0], gDown2[1] - gUp2[1] };

        // Fixed non-trivial measured Hamiltonian at N = 2: diagonal content (Z₀, Z₀Z₁) so
        // the tower fires at j = 1, plus X₀ / X₀X₁ / Y₁ flavor for the curvature probes.
        var hp2 = Combo(
            (0.3, Embed(PauliZ(), 0, 2)),
            (0.5, Embed(PauliX(), 0, 2)),
            (0.45, Embed(PauliZ(), 0, 2) * Embed(PauliZ(), 1, 2)),
            (0.7, Embed(PauliX(), 0, 2) * Embed(PauliX(), 1, 2)),
            (0.2, Embed(PauliY(), 1, 2)));

        // (1) Pump directions: D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l, dephasing and the
        //     unitary part annihilate I. Every site at N = 2 and N = 3.
        double devPump = 0.0;
        foreach (int n in new[] { 2, 3 })
        {
            int d = 1 << n;
            var id = ComplexMatrix.Build.DenseIdentity(d);
            var zero = ComplexMatrix.Build.Dense(d, d);
            for (int l = 0; l < n; l++)
            {
                var zl = Embed(PauliZ(), l, n);
                devPump = Math.Max(devPump, MaxAbsDiff(ApplyDissipator(Embed(SigmaMinus(), l, n), id), zl));
                devPump = Math.Max(devPump, MaxAbsDiff(
                    ApplyDissipator(Embed(SigmaPlus(), l, n), id), zl.Multiply(-Complex.One)));
                devPump = Math.Max(devPump, MaxAbsDiff(ApplyDissipator(zl, id), zero));
            }
        }
        // The unitary part: −i[H_p, I] = 0 on the fixed N = 2 Hamiltonian.
        {
            var id4 = ComplexMatrix.Build.DenseIdentity(4);
            var minusI = new Complex(0.0, -1.0);
            devPump = Math.Max(devPump, MaxAbsDiff(
                (hp2 * id4 - id4 * hp2).Multiply(minusI), ComplexMatrix.Build.Dense(4, 4)));
        }
        cases.Add(DevCase("pump directions: D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l, rest unital",
            "every site at N = 2 and N = 3, plus −i[H, I] = 0; amplitude damping is the unique " +
            "non-unital piece and pumps along pure local Z",
            devPump));

        // (2) The slope law at N = 2, j = 1, 2, 3: dense superoperator applied to vec(I/d)
        //     versus PredictPumpSlope, site-dependent rates.
        var super2 = Lindbladian(hp2, gDeph2, gDown2, gUp2, 2);
        double devSlope2 = 0.0;
        for (int j = 1; j <= 3; j++)
            devSlope2 = Math.Max(devSlope2, Math.Abs(
                MixedStateDerivative(super2, MatrixPower(hp2, j), 4, order: 1)
                - PredictPumpSlope(hp2, j, delta2)));
        cases.Add(DevCase("slope law at N = 2: slope_j = (1/d)·Σ_l Δγ_l·t_j(l), j = 1, 2, 3",
            "dense 16×16 superoperator applied to vec(I/4) vs the closed form, site-dependent " +
            "rates γ↓ = (0.021, 0.034), γ↑ = (0.004, 0.011), γ^deph = (0.05, 0.08)",
            devSlope2));

        // (3) The slope law at N = 3, j = 1, 2, 3: XXZ chain + site fields, site-dependent rates.
        var hp3 = Combo(
            (0.5, Embed(PauliX(), 0, 3) * Embed(PauliX(), 1, 3)),
            (0.5, Embed(PauliY(), 0, 3) * Embed(PauliY(), 1, 3)),
            (0.3, Embed(PauliZ(), 0, 3) * Embed(PauliZ(), 1, 3)),
            (0.5, Embed(PauliX(), 1, 3) * Embed(PauliX(), 2, 3)),
            (0.5, Embed(PauliY(), 1, 3) * Embed(PauliY(), 2, 3)),
            (0.3, Embed(PauliZ(), 1, 3) * Embed(PauliZ(), 2, 3)),
            (0.2, Embed(PauliZ(), 0, 3)),
            (-0.35, Embed(PauliZ(), 1, 3)),
            (0.15, Embed(PauliZ(), 2, 3)),
            (0.25, Embed(PauliX(), 1, 3)));
        double[] gDeph3 = { 0.05, 0.02, 0.07 };
        double[] gDown3 = { 0.012, 0.027, 0.018 };
        double[] gUp3 = { 0.003, 0.009, 0.001 };
        double[] delta3 = { gDown3[0] - gUp3[0], gDown3[1] - gUp3[1], gDown3[2] - gUp3[2] };
        var super3 = Lindbladian(hp3, gDeph3, gDown3, gUp3, 3);
        double devSlope3 = 0.0;
        for (int j = 1; j <= 3; j++)
            devSlope3 = Math.Max(devSlope3, Math.Abs(
                MixedStateDerivative(super3, MatrixPower(hp3, j), 8, order: 1)
                - PredictPumpSlope(hp3, j, delta3)));
        cases.Add(DevCase("slope law at N = 3: XXZ chain + site fields, j = 1, 2, 3",
            "dense 64×64 superoperator vs the closed form, site-dependent rates " +
            "γ↓ = (0.012, 0.027, 0.018), γ↑ = (0.003, 0.009, 0.001)",
            devSlope3));

        // (4) Dephasing-blind: a very different γ^deph profile leaves the dense slope fixed.
        var superDeph = Lindbladian(hp2, new double[] { 0.4, 0.9 }, gDown2, gUp2, 2);
        double devDeph = 0.0;
        for (int j = 1; j <= 3; j++)
            devDeph = Math.Max(devDeph, Math.Abs(
                MixedStateDerivative(superDeph, MatrixPower(hp2, j), 4, order: 1)
                - MixedStateDerivative(super2, MatrixPower(hp2, j), 4, order: 1)));
        cases.Add(DevCase("dephasing-blind: γ^deph (0.05, 0.08) → (0.4, 0.9), slope identical",
            "blindness #1: the dephasing channel is unital and never enters the slope; " +
            "j = 1, 2, 3 at N = 2",
            devDeph));

        // (5) Evolution-blind: replace the generator's Hamiltonian entirely; the slope of the
        //     measured polynomial A = H_p^j is unchanged and still equals the closed form.
        var hGen = Combo(
            (1.3, Embed(PauliY(), 0, 2)),
            (0.8, Embed(PauliX(), 0, 2) * Embed(PauliZ(), 1, 2)),
            (0.6, Embed(PauliX(), 1, 2)));
        var superGen = Lindbladian(hGen, gDeph2, gDown2, gUp2, 2);
        double devEvo = 0.0;
        for (int j = 1; j <= 3; j++)
            devEvo = Math.Max(devEvo, Math.Abs(
                MixedStateDerivative(superGen, MatrixPower(hp2, j), 4, order: 1)
                - PredictPumpSlope(hp2, j, delta2)));
        cases.Add(DevCase("evolution-blind: generator H replaced, slope of ⟨H_p^j⟩ unchanged",
            "blindness #2: only the measured polynomial enters; generator swapped to " +
            "1.3·Y₀ + 0.8·X₀Z₁ + 0.6·X₁, j = 1, 2, 3",
            devEvo));

        // (6) Detailed-balance closure: γ↓ = γ↑ kills every rung.
        var superBal = Lindbladian(hp2, gDeph2, gDown2, gDown2, 2);
        double devBal = 0.0;
        for (int j = 1; j <= 3; j++)
            devBal = Math.Max(devBal, Math.Abs(
                MixedStateDerivative(superBal, MatrixPower(hp2, j), 4, order: 1)));
        cases.Add(DevCase("detailed-balance closure: γ↓ = γ↑ ⟹ slope ≡ 0, j = 1, 2, 3",
            "blindness #3: the channel is powered exactly by F84's vacuum component Δγ; " +
            "at Δγ ≡ 0 every rung reads zero",
            devBal));

        // (7) The F113 bridge: asymmetry = −4^N·slope⟨H⟩ for an F113-scope generator,
        //     against the typed parent's closed form (the typed edge in action).
        double devBridge = 0.0;
        {
            // N = 2: ω = (0.13, 0.27), σ⁻ rates γ_T1, σ⁺ rates γ_pump, plus dephasing
            // (which the bridge never sees).
            double[] omega2 = { 0.13, 0.27 };
            double[] gT1 = { 0.0011, 0.0023 };
            double[] gPump = { 0.0003, 0.0007 };
            var hDrive2 = Combo(
                (omega2[0] / 2.0, Embed(PauliZ(), 0, 2)),
                (omega2[1] / 2.0, Embed(PauliZ(), 1, 2)));
            var superDrive2 = Lindbladian(hDrive2, gDeph2, gT1, gPump, 2);
            double slope2 = MixedStateDerivative(superDrive2, hDrive2, 4, order: 1);
            double asym2 = LindbladBitBPiBreakMagnitude.PredictAsymmetry(omega2, gT1, gPump, 2);
            devBridge = Math.Max(devBridge, Math.Abs(asym2 - (-Math.Pow(4.0, 2) * slope2)));

            // N = 3: pure cooling (γ_pump = 0), the F113 non-uniform-rate table row.
            double[] omega3 = { 0.05, 0.1, 0.2 };
            double[] gT13 = { 0.001, 0.002, 0.003 };
            double[] gPump3 = { 0.0, 0.0, 0.0 };
            var hDrive3 = Combo(
                (omega3[0] / 2.0, Embed(PauliZ(), 0, 3)),
                (omega3[1] / 2.0, Embed(PauliZ(), 1, 3)),
                (omega3[2] / 2.0, Embed(PauliZ(), 2, 3)));
            var superDrive3 = Lindbladian(hDrive3, new double[] { 0.04, 0.04, 0.04 }, gT13, gPump3, 3);
            double slope3 = MixedStateDerivative(superDrive3, hDrive3, 8, order: 1);
            double asym3 = LindbladBitBPiBreakMagnitude.PredictAsymmetry(omega3, gT13, gPump3, 3);
            devBridge = Math.Max(devBridge, Math.Abs(asym3 - (-Math.Pow(4.0, 3) * slope3)));
        }
        cases.Add(DevCase("F113 bridge: asymmetry = −4^N·slope⟨H⟩ (typed parent's PredictAsymmetry)",
            "static Frobenius polarity imbalance = −4^N × dynamic pump rate, exactly; " +
            "N = 2 (two-rate) and N = 3 (pure-cooling non-uniform row)",
            devBridge));

        // (8) The curvature fingerprint: exactly affine in the generator; the Y₀ parasite
        //     reads linearly against the commutator probes [Z_l, H_p^j].
        const double deltaParasite = 0.7;
        double devCurv = 0.0;
        double shiftJ1 = 0.0;
        {
            var parasite = Embed(PauliY(), 0, 2);
            var superFull = Lindbladian(hp2 + parasite.Multiply(new Complex(deltaParasite, 0.0)),
                gDeph2, gDown2, gUp2, 2);
            for (int j = 1; j <= 2; j++)
            {
                var a = MatrixPower(hp2, j);
                double shift = MixedStateDerivative(superFull, a, 4, order: 2)
                             - MixedStateDerivative(super2, a, 4, order: 2);
                double pred = CommutatorProbeShift(parasite, a, deltaParasite, delta2, 2);
                devCurv = Math.Max(devCurv, Math.Abs(shift - pred));
                if (j == 1) shiftJ1 = shift;
            }
        }
        bool curvOk = devCurv <= Tol && Math.Abs(shiftJ1) > 1e-6;
        cases.Add(new BatteryCase(
            Name: "curvature fingerprint: exact affinity, Y₀ parasite reads on the commutator probes",
            Detail: "Δcurvature = (δ/d)·Σ_l Δγ_l·(−i)·Tr(V·[Z_l, A]) with δ = 0.7, V = Y₀, " +
                    "A = H_p^j for j = 1, 2; no small-δ expansion, the curvature is exactly " +
                    "affine in the generator (shift(j=1) = " +
                    shiftJ1.ToString("E2", CultureInfo.InvariantCulture) + ")",
            Expected: "affinity dev ≤ 1e-12 and |shift(j=1)| > 1e-6",
            Actual: curvOk
                ? "affinity dev ≤ 1e-12 and |shift(j=1)| > 1e-6"
                : "affinity dev = " + devCurv.ToString("E2", CultureInfo.InvariantCulture) +
                  ", |shift(j=1)| = " + Math.Abs(shiftJ1).ToString("E2", CultureInfo.InvariantCulture)));

        // (9) The blind spot: a Z₀ parasite is exactly invisible ([Z_l, Z-string] = 0).
        double devBlind = 0.0;
        {
            var parasite = Embed(PauliZ(), 0, 2);
            var superFull = Lindbladian(hp2 + parasite.Multiply(new Complex(deltaParasite, 0.0)),
                gDeph2, gDown2, gUp2, 2);
            for (int j = 1; j <= 2; j++)
            {
                var a = MatrixPower(hp2, j);
                devBlind = Math.Max(devBlind, Math.Abs(
                    MixedStateDerivative(superFull, a, 4, order: 2)
                    - MixedStateDerivative(super2, a, 4, order: 2)));
                devBlind = Math.Max(devBlind, Math.Abs(
                    CommutatorProbeShift(parasite, a, deltaParasite, delta2, 2)));
            }
        }
        cases.Add(DevCase("curvature blind spot: Z₀ parasite exactly invisible",
            "[Z_l, Z-string] = 0, so both the measured curvature shift and the probe " +
            "prediction vanish; the division of labor with F113's Z-drive reader",
            devBlind));

        // (10) The girth-2 witness: H = X₀ + X₀Z₁ at N = 2 has t₁ ≡ 0 (silent rung) and
        //      t₂ = (0, 8) (firing rung); the dense slope agrees rung by rung.
        {
            var hWitness = Embed(PauliX(), 0, 2) + Embed(PauliX(), 0, 2) * Embed(PauliZ(), 1, 2);
            var t1 = MomentTower(hWitness, 1, 2);
            var t2 = MomentTower(hWitness, 2, 2);
            var superW = Lindbladian(hWitness, gDeph2, gDown2, gUp2, 2);
            double slope1 = MixedStateDerivative(superW, hWitness, 4, order: 1);
            double slope2 = MixedStateDerivative(superW, MatrixPower(hWitness, 2), 4, order: 1);
            double devW = Math.Max(
                Math.Max(Math.Abs(t1[0]), Math.Abs(t1[1])),
                Math.Max(Math.Abs(t2[0]), Math.Abs(t2[1] - 8.0)));
            devW = Math.Max(devW, Math.Abs(slope1));
            devW = Math.Max(devW, Math.Abs(slope2 - PredictPumpSlope(hWitness, 2, delta2)));
            bool fires = Math.Abs(slope2) > 1e-6;
            bool witnessOk = devW <= Tol && fires;
            cases.Add(new BatteryCase(
                Name: "girth-2 witness X₀ + X₀Z₁: t₁ ≡ 0, t₂ = (0, 8), slope₂ fires",
                Detail: "the deg-1 tower is silent at j = 1 and fires at j = 2 (the girth); " +
                        "slope₂ = Δγ₁·t₂(1)/d = " + slope2.ToString("E2", CultureInfo.InvariantCulture) +
                        "; a firing rung is a complete hardness certificate (m* = 2ℓ+1)",
                Expected: "t₁ = (0, 0), t₂ = (0, 8), slope₁ = 0, slope₂ fires (dev ≤ 1e-12)",
                Actual: witnessOk
                    ? "t₁ = (0, 0), t₂ = (0, 8), slope₁ = 0, slope₂ fires (dev ≤ 1e-12)"
                    : "dev = " + devW.ToString("E2", CultureInfo.InvariantCulture) +
                      ", |slope₂| = " + Math.Abs(slope2).ToString("E2", CultureInfo.InvariantCulture)));
        }

        return cases;
    }

    // ------------------------------------------------------------------
    // Dense building blocks (row-stacking vec convention: |i⟩⟨j| ↦ e_i ⊗ e_j,
    // kron(A, B): ρ ↦ A·ρ·Bᵀ, matching MirrorGroupD4Claim / AntilinearTriangleClaim).
    // ------------------------------------------------------------------

    private static ComplexMatrix PauliX()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 1] = Complex.One; m[1, 0] = Complex.One;
        return m;
    }

    private static ComplexMatrix PauliY()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 1] = new Complex(0.0, -1.0); m[1, 0] = new Complex(0.0, 1.0);
        return m;
    }

    private static ComplexMatrix PauliZ()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 0] = Complex.One; m[1, 1] = -Complex.One;
        return m;
    }

    /// <summary>σ⁻ = |0⟩⟨1| = [[0, 1], [0, 0]] (standard physics convention, the lowering
    /// operator; same convention as <see cref="LindbladBitBPiBreakMagnitude"/>).</summary>
    private static ComplexMatrix SigmaMinus()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[0, 1] = Complex.One;
        return m;
    }

    /// <summary>σ⁺ = [[0, 0], [1, 0]], the raising operator (pumping channel).</summary>
    private static ComplexMatrix SigmaPlus()
    {
        var m = ComplexMatrix.Build.Dense(2, 2);
        m[1, 0] = Complex.One;
        return m;
    }

    /// <summary>Embed a single-site operator at <paramref name="site"/> (leftmost tensor
    /// factor = site 0 = most significant qubit) into the N-site space.</summary>
    private static ComplexMatrix Embed(ComplexMatrix op, int site, int n)
    {
        var result = ComplexMatrix.Build.DenseIdentity(1);
        var i2 = ComplexMatrix.Build.DenseIdentity(2);
        for (int k = 0; k < n; k++)
            result = result.KroneckerProduct(k == site ? op : i2);
        return result;
    }

    /// <summary>Σ_i c_i·M_i for a fixed list of real-weighted dense terms.</summary>
    private static ComplexMatrix Combo(params (double C, ComplexMatrix M)[] terms)
    {
        var r = terms[0].M.Multiply(new Complex(terms[0].C, 0.0));
        for (int i = 1; i < terms.Length; i++)
            r += terms[i].M.Multiply(new Complex(terms[i].C, 0.0));
        return r;
    }

    private static ComplexMatrix MatrixPower(ComplexMatrix h, int j)
    {
        var p = h;
        for (int k = 1; k < j; k++)
            p *= h;
        return p;
    }

    /// <summary>D[c](ρ) = cρc† − ½{c†c, ρ}, applied directly to a dense ρ.</summary>
    private static ComplexMatrix ApplyDissipator(ComplexMatrix c, ComplexMatrix rho)
    {
        var cdc = c.ConjugateTranspose() * c;
        return c * rho * c.ConjugateTranspose()
             - (cdc * rho + rho * cdc).Multiply(new Complex(0.5, 0.0));
    }

    /// <summary>The dense d²×d² Lindblad superoperator −i[H,·] + Σ_l γ^deph_l·D[Z_l]
    /// + γ↓_l·D[σ⁻_l] + γ↑_l·D[σ⁺_l] in the row-stacking vec convention.</summary>
    private static ComplexMatrix Lindbladian(
        ComplexMatrix h,
        IReadOnlyList<double> gammaDeph,
        IReadOnlyList<double> gammaDown,
        IReadOnlyList<double> gammaUp,
        int n)
    {
        int d = 1 << n;
        var id = ComplexMatrix.Build.DenseIdentity(d);
        var minusI = new Complex(0.0, -1.0);
        var m = (h.KroneckerProduct(id) - id.KroneckerProduct(h.Transpose())).Multiply(minusI);
        for (int l = 0; l < n; l++)
        {
            m += DissipatorSuper(Embed(PauliZ(), l, n), id).Multiply(new Complex(gammaDeph[l], 0.0));
            m += DissipatorSuper(Embed(SigmaMinus(), l, n), id).Multiply(new Complex(gammaDown[l], 0.0));
            m += DissipatorSuper(Embed(SigmaPlus(), l, n), id).Multiply(new Complex(gammaUp[l], 0.0));
        }
        return m;
    }

    private static ComplexMatrix DissipatorSuper(ComplexMatrix c, ComplexMatrix id)
    {
        var cdc = c.ConjugateTranspose() * c;
        return c.KroneckerProduct(c.Conjugate())
             - (cdc.KroneckerProduct(id) + id.KroneckerProduct(cdc.Transpose()))
                   .Multiply(new Complex(0.5, 0.0));
    }

    /// <summary>Re Tr(A·L^order(I/d)): the slope (order 1) or curvature (order 2) of ⟨A⟩
    /// at the maximally mixed state, via the dense superoperator applied to vec(I/d).</summary>
    private static double MixedStateDerivative(ComplexMatrix superL, ComplexMatrix aObs, int d, int order)
    {
        var vec = ComplexVector.Build.Dense(d * d);
        for (int i = 0; i < d; i++)
            vec[i * d + i] = new Complex(1.0 / d, 0.0);
        for (int k = 0; k < order; k++)
            vec = superL * vec;
        Complex tr = Complex.Zero;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                tr += aObs[i, j] * vec[j * d + i];
        return tr.Real;
    }

    /// <summary>The predicted curvature shift of a parasite δV against the commutator
    /// probes: (δ/d)·Σ_l Δγ_l·(−i)·Tr(V·[Z_l, A]).</summary>
    private static double CommutatorProbeShift(
        ComplexMatrix v, ComplexMatrix a, double delta, IReadOnlyList<double> deltaGamma, int n)
    {
        int d = 1 << n;
        var minusI = new Complex(0.0, -1.0);
        Complex sum = Complex.Zero;
        for (int l = 0; l < n; l++)
        {
            var zl = Embed(PauliZ(), l, n);
            sum += deltaGamma[l] * (minusI * (v * (zl * a - a * zl)).Trace());
        }
        return delta / d * sum.Real;
    }

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
