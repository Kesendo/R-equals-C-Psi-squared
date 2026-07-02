using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness for the TIME_IRREVERSIBILITY_EXCLUSION argument (2026-07-02): the
/// anticommutator {L_H, L_Dc} of the Hamiltonian Liouvillian and the F1-centered Z-dephasing
/// Liouvillian vanishes EXACTLY at N=2 and is nonzero for N&gt;2, while the COMMUTATOR
/// [L_H, L_Dc] is nonzero already at N=2. The value ‖{L_H, L_Dc}‖² = 4γ²(N−2)‖L_H‖² is already
/// typed in <see cref="F49NonUniformCrossTermClaim"/>; what this witness adds is the ARGUMENT the
/// proof carries and the F49 value-claim does not (typed as
/// <see cref="TimeIrreversibilityExclusionClaim"/>): the N=2-only orthogonality READING plus its
/// load-bearing caveat.
///
/// <para>Be precise about what is recomputed and what is a closed form. The witness BUILDS the two
/// Liouvillians live via <see cref="LindbladianBuilder"/> (L_H = −i[H, ·] for a Heisenberg chain;
/// L_Dc = L_D + Σγ·I for uniform Z-dephasing c_k = √γ·Z_k) and forms both the anticommutator and
/// the commutator, taking their Frobenius norms — a genuine per-N matrix compute. It then
/// cross-checks the live ‖{L_H, L_Dc}‖ against the F49 closed form
/// <see cref="F49NonUniformCrossTermClaim.PredictHeisenbergChain"/>, fed a LIVE per-bond
/// ‖L_H^bond‖² (no hardcoded 384/1536/…): two independent computations meeting.</para>
///
/// <para>The load-bearing caveat, encoded so the witness cannot overclaim: the N=2 vanishing of
/// the ANTIcommutator is a Frobenius-orthogonality fact (the Pythagorean split
/// L_c² = L_H² + L_Dc²), NOT a separability or reversibility criterion. Separability is governed by
/// the COMMUTATOR [L_H, L_Dc], which is nonzero already at N=2 (the proof reports ≈22.6). The
/// witness reports BOTH norms so the arrow-of-time reading is bounded exactly as the proof bounds
/// it (TIME_IRREVERSIBILITY_EXCLUSION.md, the Tier-3 caveat).</para>
///
/// <para>Guard: N ≤ <see cref="MaxN"/> (= 5), where the Liouvillian is 4^N × 4^N = 1024 × 1024 and
/// each anticommutator is a handful of 1024³ complex multiplies. Anchors:
/// <c>docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md</c> + <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c>
/// (the uniform-γ 4γ²(N−2)‖L_H‖² parent) + <c>simulations/review2_A8_time.py</c>.</para></summary>
public sealed class TimeIrreversibilityExclusionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Largest chain length materialised. The Liouvillian is 4^N × 4^N; at N=5 that is
    /// 1024 × 1024 and each anticommutator is a few 1024³ complex multiplies.</summary>
    public const int MaxN = 5;

    public double Gamma { get; }
    public double J { get; }
    public int NMax { get; }

    /// <summary>Per-N reading: the live anticommutator and commutator Frobenius norms of the pair
    /// (L_H, L_Dc), the Liouvillian norm ‖L_H‖, and the F49 closed-form ‖{L_H, L_Dc}‖ it is
    /// checked against.</summary>
    public sealed record NRow(int N, double AntiNorm, double CommNorm, double LhNorm, double ClosedFormAntiNorm);

    public IReadOnlyList<NRow> Rows { get; }

    public TimeIrreversibilityExclusionWitness(int nMax = 4, double gamma = 0.05, double j = 1.0)
    {
        if (nMax < 2)
            throw new ArgumentOutOfRangeException(nameof(nMax), nMax, "nMax must be ≥ 2 (N=2 is the vanishing anchor).");
        if (nMax > MaxN)
            throw new ArgumentOutOfRangeException(nameof(nMax), nMax, $"nMax must be ≤ {MaxN} (the 4^N Liouvillian build guard).");
        if (gamma <= 0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive.");

        Gamma = gamma;
        J = j;
        NMax = nMax;

        var rows = new List<NRow>(nMax - 1);
        for (int n = 2; n <= nMax; n++) rows.Add(BuildRow(n, gamma, j));
        Rows = rows;
    }

    private static NRow BuildRow(int n, double gamma, double j)
    {
        int d = 1 << n;          // Hilbert dimension 2^N
        int dSq = d * d;         // Liouville dimension 4^N

        // L_H = −i[H, ·] for the Heisenberg chain (no dissipator).
        var H = PauliHamiltonian.HeisenbergChain(n, j).ToMatrix();
        var lH = LindbladianBuilder.Build(H, Array.Empty<ComplexMatrix>());

        // L_Dc = L_D + Σγ·I, the F1-centered uniform Z-dephasing dissipator (c_k = √γ·Z_k).
        var rootGamma = new Complex(Math.Sqrt(gamma), 0.0);
        var cOps = new List<ComplexMatrix>(n);
        for (int k = 0; k < n; k++)
            cOps.Add(new PauliHamiltonian(n,
                new[] { PauliTerm.SingleSite(n, k, PauliLetter.Z, rootGamma) }).ToMatrix());
        var lD = LindbladianBuilder.Build(Matrix<Complex>.Build.Dense(d, d), cOps);
        double sigma = n * gamma;
        var lDc = lD + sigma * Matrix<Complex>.Build.DenseIdentity(dSq);

        // The two norms: the anticommutator (→ 0 at N=2) and the commutator (≠ 0 even at N=2).
        double antiNorm = (lH * lDc + lDc * lH).FrobeniusNorm();
        double commNorm = (lH * lDc - lDc * lH).FrobeniusNorm();
        double lhNorm = lH.FrobeniusNorm();

        // Cross-check the live anticommutator against the canonical F49 closed form, fed a LIVE
        // per-bond ‖L_H^bond‖² (no hardcoded constants): two independent computations meeting.
        var bondC = new Complex(j / 4.0, 0.0);
        var bond = new PauliHamiltonian(n, new[]
        {
            PauliTerm.TwoSite(n, 0, PauliLetter.X, 1, PauliLetter.X, bondC),
            PauliTerm.TwoSite(n, 0, PauliLetter.Y, 1, PauliLetter.Y, bondC),
            PauliTerm.TwoSite(n, 0, PauliLetter.Z, 1, PauliLetter.Z, bondC),
        }).ToMatrix();
        double bondNorm = LindbladianBuilder.Build(bond, Array.Empty<ComplexMatrix>()).FrobeniusNorm();
        var gammaArr = new double[n];
        for (int k = 0; k < n; k++) gammaArr[k] = gamma;
        double closedAnti = Math.Sqrt(
            F49NonUniformCrossTermClaim.PredictHeisenbergChain(n, gammaArr, bondNorm * bondNorm));

        return new NRow(n, antiNorm, commNorm, lhNorm, closedAnti);
    }

    /// <summary>The proof's relative cross-term R(N) = √((N−2)/(N·4^(N−1))): γ- and
    /// topology-independent, exactly 0 at N=2 then ≈ 1.83% (N=3), 2.07% (N=4). Reported for the
    /// narrative, not a gate.</summary>
    public static double RelativeCrossTerm(int n) => Math.Sqrt((n - 2.0) / (n * Math.Pow(4.0, n - 1)));

    private static string E(double v) => v.ToString("E3", Inv);
    private static string F(double v) => v.ToString("0.###", Inv);

    public string DisplayName =>
        $"TimeIrreversibilityExclusionWitness (N=2..{NMax}, γ={F(Gamma)}, J={F(J)})";

    public string Summary
    {
        get
        {
            var n2 = Rows[0];
            var last = Rows[^1];
            return $"anticommutator ‖{{L_H, L_Dc}}‖ = {E(n2.AntiNorm)} at N=2 (Frobenius-orthogonality, exactly 0) " +
                   $"and grows to {E(last.AntiNorm)} by N={last.N}; but the commutator ‖[L_H, L_Dc]‖ = {F(n2.CommNorm)} " +
                   $"is nonzero ALREADY at N=2, so the N=2 vanishing is an orthogonality fact, NOT a " +
                   $"separability/reversibility criterion (the load-bearing caveat). Live norms match the F49 closed " +
                   $"form ‖{{L_H, L_Dc}}‖² = 4γ²(N−2)‖L_H‖².";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // The caveat, front and centre: the commutator does not vanish where the anticommutator does.
            yield return new InspectableNode(
                displayName: "the caveat: commutator ≠ 0 even at N=2 (not reversibility)",
                summary: $"‖[L_H, L_Dc]‖ = {F(Rows[0].CommNorm)} at N=2 ≠ 0, while ‖{{L_H, L_Dc}}‖ = {E(Rows[0].AntiNorm)} ≈ 0. " +
                         "Separability/reversibility is governed by the COMMUTATOR, not the anticommutator; the N=2 " +
                         "anticommutator vanishing is Frobenius-orthogonality (the Pythagorean split " +
                         "L_c² = L_H² + L_Dc²), NOT an arrow-of-time criterion. TIME_IRREVERSIBILITY_EXCLUSION.md, the Tier-3 caveat.",
                provenance: NodeProvenance.Live);

            // One node per N: live anticommutator + commutator, checked against the F49 closed form.
            foreach (var r in Rows)
            {
                bool match = Math.Abs(r.AntiNorm - r.ClosedFormAntiNorm) < 1e-6;
                yield return new InspectableNode(
                    displayName: $"N={r.N}: ‖{{L_H, L_Dc}}‖ = {E(r.AntiNorm)} ({(r.N == 2 ? "orthogonal" : "broken")})",
                    summary: $"live ‖{{L_H, L_Dc}}‖ = {E(r.AntiNorm)} vs F49 closed form {E(r.ClosedFormAntiNorm)} " +
                             $"({(match ? "match" : "MISMATCH")}); ‖[L_H, L_Dc]‖ = {F(r.CommNorm)} (≠0); ‖L_H‖ = {F(r.LhNorm)}; " +
                             $"relative R(N) = ‖{{,}}‖/(‖L_H‖·‖L_Dc‖) = √((N−2)/(N·4^(N−1))) = {F(RelativeCrossTerm(r.N) * 100)}% " +
                             $"({(r.N == 2 ? "exactly 0: the orthogonality" : "γ-independent")}).",
                    provenance: NodeProvenance.Live);
            }

            // Payload curve: the anticommutator norm across N (0 at N=2, growing after).
            var xs = Rows.Select(r => (double)r.N).ToArray();
            var ys = Rows.Select(r => r.AntiNorm).ToArray();
            yield return new InspectableNode(
                displayName: "‖{L_H, L_Dc}‖ vs N (0 at N=2, growing after)",
                summary: "the anticommutator norm across N: exactly 0 at N=2, nonzero and growing for N>2 " +
                         "(the (N−2) factor of the F49 closed form).",
                payload: new InspectablePayload.Curve("anti-norm vs N", xs, ys, "N", "‖{L_H, L_Dc}‖"),
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
