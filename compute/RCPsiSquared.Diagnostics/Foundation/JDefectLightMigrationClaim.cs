using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The jdefect axis's light-migration mechanism, the first typed claim living ON a
/// navigator axis (<see cref="JDefectField"/>). Along the J-defect axis (uniform XY chain +
/// Z-dephasing, bond defect δJ; a Π-invariant NON-similarity: the palindrome holds at every δJ
/// but the spectrum moves) the per-mode absorption identity tracks the motion exactly:
///
/// <code>
///   Re λ(δJ) = −2·Σ_l γ_l·light_l(v(δJ))   for EVERY eigenmode v(δJ), at EVERY δJ,
///   light(v) = Σ_x w(x)·|v_x|² / Σ_x |v_x|²,  w(x) = popcount(i ⊕ j),  x = i·d + j
/// </code>
///
/// The identity is δJ-pointwise: the Re-drift of every mark IS light migration of its
/// eigenvector, and the Hamiltonian contributes zero throughout (H(δJ) is Hermitian at every
/// δJ, so its Liouvillian part is anti-Hermitian and drops out of the Rayleigh quotient).
/// Two constraints ride on it: (a) the N+1 kernel modes stay at light ≡ 0 at every δJ (the
/// identity read backwards at Re λ = 0; U(1) protection, no migration); (b) palindrome
/// partners migrate oppositely, light_s(δJ) + light_f(δJ) = N at every δJ (complementarity is
/// δJ-pointwise because the palindrome is Π-invariant along the whole axis). So the axis's
/// "marks move O(δJ)" reading gains a mechanism: the marks move exactly as much as the
/// eigenvectors' light migrates, and partner migrations cancel.
///
/// <para><b>Tier1Derived honestly stated as a composition</b>: this claim is two proven
/// identities composed along the axis, plus the composition itself. (1) The per-mode
/// absorption identity (<see cref="Core.Symmetry.AbsorptionTheoremClaim"/>, vector form):
/// exact for ANY Hermitian H because Herm(L) is the dephasing dissipator alone; H(δJ) is
/// Hermitian at every δJ, so the identity holds pointwise along the axis. (2) The F1
/// palindrome (<see cref="Core.F1.F1PalindromeIdentity"/>, Π·L·Π⁻¹ = −L − 2Σγ·I): valid for
/// XY on any graph with any bond weights, so the detuned bond keeps the complete pairing
/// λ_s + λ_f = −2Σγ at every δJ; substituting (1) into (2) gives the pointwise light
/// complementarity. What is new is the composition read ON the axis: the Π-invariant
/// non-similarity's spectral motion is exactly eigenvector light migration.</para>
///
/// <para><b>Verification</b> (<c>simulations/jdefect_light_migration.py</c>, N=4 and N=5,
/// γ=0.05 uniform, δJ ∈ {0, 0.05, 0.1} on bond (0,1)): identity residual ≤ 8.4·10⁻¹⁵ (N=4)
/// and ≤ 1.6·10⁻¹⁴ (N=5) over every mode at every δJ; kernel exactly N+1 modes with light
/// ≤ 5.4·10⁻²⁸; complete pairing with |light_s + light_f − N| ≤ 7.9·10⁻¹⁴ per pair at every
/// δJ; max per-mode light migration δJ 0 → 0.1: 0.0268 at N=4 (88/256 modes move),
/// 0.0449 at N=5 (672/1024 modes move). The migration is the content of the axis; the
/// identity is its bookkeeping.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> +
/// <c>simulations/jdefect_light_migration.py</c>.</para></summary>
public sealed class JDefectLightMigrationClaim : Claim
{
    /// <summary>One dense self-check of the composed identity on the axis at N = 3.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public JDefectLightMigrationClaim()
        : base("JDefect light migration: along the jdefect axis Re λ(δJ) = −2·Σ_l γ_l·light_l(v(δJ)) for every eigenmode at every δJ (the Re-drift IS eigenvector light migration, H contributes zero); the N+1 kernel modes stay at light ≡ 0 (U(1)); palindrome partners migrate oppositely, light_s + light_f = N δJ-pointwise",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/proofs/MIRROR_SYMMETRY_PROOF.md + " +
               "simulations/jdefect_light_migration.py")
    {
        Cases = BuildBattery();
    }

    public string Identity =>
        "Re λ = ⟨v|Herm(L)|v⟩/‖v‖² for any right eigenmode of any L; here Herm(L) is the " +
        "Z-dephasing dissipator alone (diagonal, −2γ·popcount(i⊕j)) because L_H = −i[H, ·] is " +
        "anti-Hermitian for every Hermitian H. H(δJ) is Hermitian at every δJ, so the identity " +
        "is δJ-pointwise along the whole axis: the marks move exactly as the light migrates.";

    public string KernelConstraint =>
        "Re λ = 0 forces light(v) = 0 (the identity read backwards): the N+1 U(1) steady states " +
        "are diagonal at every δJ and carry no light. The kernel cannot migrate; it is the dark " +
        "anchor of the moving spectrum.";

    public string Complementarity =>
        "The F1 palindrome Π·L·Π⁻¹ = −L − 2Σγ·I holds for XY on any graph with any bond weights, " +
        "so the defect keeps the complete pairing λ_s + λ_f = −2Σγ at every δJ; substituting the " +
        "absorption identity gives light_s(δJ) + light_f(δJ) = N pointwise. Partner migrations " +
        "cancel exactly: what one bank gains in light, the mirror bank loses.";

    public string Composition =>
        "Two proven identities composed along the axis: the per-mode absorption identity " +
        "(AbsorptionTheoremClaim, vector form) + the F1 palindrome (F1PalindromeIdentity), each " +
        "applied δJ-pointwise. The new content is the reading ON the axis: JDefectField's " +
        "Π-invariant non-similarity (palindrome kept, spectrum moves) gains its mechanism, the " +
        "Re-motion of the marks is eigenvector light migration, nothing else.";

    public override string DisplayName =>
        "JDefect light migration (the absorption identity δJ-pointwise on the jdefect axis, Tier1Derived)";

    public override string Summary =>
        "Re λ(δJ) = −2γ·light(v(δJ)) for every mode at every δJ; kernel dark (light ≡ 0), " +
        "partners complementary (light_s + light_f = N) pointwise; the marks move exactly as " +
        $"the eigenvectors' light migrates; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The identity (δJ-pointwise, H contributes zero)", summary: Identity);
            yield return new InspectableNode("The kernel constraint (U(1), no migration)", summary: KernelConstraint);
            yield return new InspectableNode("The complementarity (partner migrations cancel)", summary: Complementarity);
            yield return new InspectableNode("The composition (what is proven, what is new)", summary: Composition);
            yield return new InspectableNode("Verification (simulations/jdefect_light_migration.py)",
                summary: "N=4: identity ≤ 8.4e-15, kernel 5 modes light ≤ 1.1e-28, complementarity ≤ 2.0e-14, " +
                         "max migration 0.0268 (88/256 modes); N=5: identity ≤ 1.6e-14, kernel 6 modes light " +
                         "≤ 5.3e-28, complementarity ≤ 7.9e-14, max migration 0.0449 (672/1024 modes); " +
                         "all at δJ ∈ {0, 0.05, 0.1}, γ=0.05 uniform, defect bond (0,1).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    /// <summary>light(v) = Σ_x w(x)·|v_x|² / Σ_x |v_x|² with w(x) = popcount(i ⊕ j) at the
    /// row-major coherence index x = i·d + j, matching <see cref="PauliDephasingDissipator"/>'s
    /// vec convention.</summary>
    private static double Light(ComplexVector v, int dim)
    {
        double num = 0.0, den = 0.0;
        for (int x = 0; x < v.Count; x++)
        {
            var c = v[x];
            double p = c.Real * c.Real + c.Imaginary * c.Imaginary;
            num += BitOperations.PopCount((uint)((x / dim) ^ (x % dim))) * p;
            den += p;
        }
        return num / den;
    }

    /// <summary>Dense battery at N = 3 (64 modes per point) on the axis's own Hamiltonian
    /// builder (<see cref="DimensionAxis.JDefect"/>, defect bond (0,1), γ = 0.05 uniform):
    /// the identity at δJ ∈ {0, 0.1}, the kernel zeros, one partner pair's complementarity,
    /// and the nonzero light movement (the content).</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        const int N = 3;
        const double Gamma = 0.05;
        const double tol = 1e-9;
        int dim = 1 << N;
        string Zero(double v) => v < tol ? "0" : v.ToString("E3", System.Globalization.CultureInfo.InvariantCulture);

        var axis = DimensionAxis.JDefect(N, gamma: Gamma, defectBond: 0, deltaJMax: 0.1, points: 2);
        var cases = new List<BatteryCase>();
        var lightsPerPoint = new double[2][];

        for (int p = 0; p < 2; p++)
        {
            double deltaJ = p == 0 ? 0.0 : 0.1;
            string deltaJLabel = p == 0 ? "0" : "0.1";
            var L = PauliDephasingDissipator.BuildZ(axis.Hamiltonian(deltaJ), axis.GammaPerSite);
            var evd = L.Evd();
            var vals = evd.EigenValues.ToArray();
            var right = evd.EigenVectors;

            var light = new double[vals.Length];
            for (int k = 0; k < vals.Length; k++) light[k] = Light(right.Column(k), dim);
            lightsPerPoint[p] = light;

            // (a) the absorption identity, every mode
            double worstIdentity = 0.0;
            for (int k = 0; k < vals.Length; k++)
                worstIdentity = Math.Max(worstIdentity, Math.Abs(vals[k].Real + 2.0 * Gamma * light[k]));
            cases.Add(new BatteryCase(
                Name: $"absorption identity at δJ={deltaJLabel}",
                Detail: $"max |Re λ + 2γ·light(v)| over all {vals.Length} modes (N=3)",
                Expected: "0",
                Actual: Zero(worstIdentity)));

            // (b) the kernel: exactly N+1 modes at λ ≈ 0, all dark
            int kernelCount = 0;
            double kernelLight = 0.0;
            for (int k = 0; k < vals.Length; k++)
                if (vals[k].Magnitude < 1e-8) // far below the smallest nonzero |λ| ≈ 0.086 at N=3
                {
                    kernelCount++;
                    kernelLight = Math.Max(kernelLight, light[k]);
                }
            cases.Add(new BatteryCase(
                Name: $"kernel dark at δJ={deltaJLabel}",
                Detail: $"modes with |λ| < 1e-8 (expect N+1 = {N + 1}), max light among them",
                Expected: $"{N + 1} modes, light 0",
                Actual: $"{kernelCount} modes, light {Zero(kernelLight)}"));
        }

        // (c) one partner pair's complementarity at δJ = 0.1 (the moved spectrum):
        // the slowest strictly decaying mode and its palindrome partner.
        {
            var L = PauliDephasingDissipator.BuildZ(axis.Hamiltonian(0.1), axis.GammaPerSite);
            var evd = L.Evd();
            var vals = evd.EigenValues.ToArray();
            var right = evd.EigenVectors;
            var target = new Complex(-2.0 * N * Gamma, 0.0);

            int s = -1;
            for (int k = 0; k < vals.Length; k++)
                if (vals[k].Real < -1e-8 && (s < 0 || vals[k].Real > vals[s].Real)) s = k;
            int f = -1;
            double bestPair = double.MaxValue;
            for (int k = 0; k < vals.Length; k++)
            {
                if (k == s) continue;
                double r = (vals[s] + vals[k] - target).Magnitude;
                if (r < bestPair) { bestPair = r; f = k; }
            }
            double comp = Math.Abs(Light(right.Column(s), dim) + Light(right.Column(f), dim) - N);
            cases.Add(new BatteryCase(
                Name: "partner complementarity at δJ=0.1",
                Detail: $"slowest decaying mode λ_s = {vals[s].ToString("0.####", System.Globalization.CultureInfo.InvariantCulture)} and its palindrome partner: |λ_s + λ_f + 2Nγ| and |light_s + light_f − N|",
                Expected: "0 and 0",
                Actual: $"{Zero(bestPair)} and {Zero(comp)}"));
        }

        // (d) the content: the light distribution genuinely moves δJ 0 → 0.1 (sorted-vector
        // comparison, a matching-free lower bound on per-mode migration; ≈ 0.0136 at N=3).
        {
            var l0 = (double[])lightsPerPoint[0].Clone();
            var l1 = (double[])lightsPerPoint[1].Clone();
            Array.Sort(l0);
            Array.Sort(l1);
            double move = 0.0;
            for (int k = 0; k < l0.Length; k++) move = Math.Max(move, Math.Abs(l1[k] - l0[k]));
            cases.Add(new BatteryCase(
                Name: "the light migrates (the content of the axis)",
                Detail: $"max sorted-light movement δJ 0 → 0.1 = {move.ToString("0.######", System.Globalization.CultureInfo.InvariantCulture)} (matching-free lower bound; expect > 1e-3)",
                Expected: "nonzero",
                Actual: move > 1e-3 ? "nonzero" : Zero(move)));
        }

        return cases;
    }
}
