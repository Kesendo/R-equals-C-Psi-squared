using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness: the F89 path-3 octic branch locus is a palindrome. The EP/diabolic collisions,
/// as q = J/γ varies in the complex plane, are mirror-symmetric about Re λ = −σ = −4 (the midpoint of the
/// absorption rungs −2γ, −6γ, the AT-midpoint). This is FORCED, not coincidental: the F1 palindrome is
/// carried on the (SE,DE) block as an ANTIUNITARY symmetry T = P·K (P the weight-complement rung swap,
/// K conjugation), T L(q) T⁻¹ = −L(q̄) − 2σ (a same-q fold at real q), so the merged-eigenvalue locus is invariant under
/// the antilinear mirror λ ↦ −λ̄ − 2σ (reflect Re about −σ, preserve Im). Every EP lies on the line or in
/// a mirror pair across it, no orphan. The line is the palindrome's gift; the diabolic's silence
/// (semisimple) is free-fermion integrability's, a separate gift (DIABOLIC_BY_INTEGRABILITY: XXZ Δ≠0
/// stays on the line yet defects). Typed home: F89BranchLocusPalindromeClaim; reading:
/// reflections/ON_WHO_WATCHES_WHOM.md; the q-direction sibling of --root galoismonodromy.</summary>
public sealed class BranchLocusPalindromeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Sigma = 4.0;   // Σγ = N·γ = 4 (path-3 octic: N=4, γ=1); palindrome centre = −σ = −4

    /// <summary>The two-sided gate. Max residual of closing the octic roots at q0 under the ANTIUNITARY
    /// mirror λ ↦ −λ̄ − 2σ (Re ↦ −Re−2σ, Im preserved) — should be ≈ 0 — and under the bare LINEAR
    /// palindrome λ ↦ −λ − 2σ (Im flipped) — should be large, so a sign/construction error fails it.</summary>
    public static (double antiunitary, double linear) MirrorClosureResiduals(double q0)
    {
        var roots = GaloisMonodromyWitness.OcticRootsAt(new Complex(q0, 0));
        double MaxResidual(Func<Complex, Complex> mirror)
        {
            double max = 0;
            foreach (var r in roots)
            {
                var image = mirror(r);
                max = Math.Max(max, roots.Min(s => (s - image).Magnitude));
            }
            return max;
        }
        double anti = MaxResidual(z => new Complex(-z.Real - 2 * Sigma, z.Imaginary));
        double lin = MaxResidual(z => new Complex(-z.Real - 2 * Sigma, -z.Imaginary));
        return (anti, lin);
    }

    private static Complex Collision(Complex q)
    {
        var r = GaloisMonodromyWitness.OcticRootsAt(q);
        int bi = 0, bj = 1; double bg = double.PositiveInfinity;
        for (int i = 0; i < r.Length; i++)
            for (int j = i + 1; j < r.Length; j++)
            {
                double g = (r[i] - r[j]).Magnitude;
                if (g < bg) { bg = g; bi = i; bj = j; }
            }
        return (r[bi] + r[bj]) / 2;
    }

    /// <summary>Each branch point's collision λ classified: on the line (|Re+σ| small) or with a mirror
    /// partner (another collision at −Re−2σ, same Im). Returns (total, onLine, paired, orphans). The
    /// cell-scan q are resolution-limited, so the tolerances are generous; the exact closure is the
    /// algebraic gate above and the discriminant factorisation in the experiment.</summary>
    public static (int total, int onLine, int paired, int orphans) BranchLocusStructure(
        double reLo = -2.0, double reHi = 2.0, double imLo = -0.3, double imHi = 0.3, double cell = 0.04)
    {
        var lams = GaloisMonodromyWitness.FindBranchPoints(reLo, reHi, imLo, imHi, cell)
            .Select(ep => Collision(ep.Q)).ToList();
        const double lineTol = 0.08, pairTol = 0.18;
        int onLine = 0, paired = 0, orphans = 0;
        for (int i = 0; i < lams.Count; i++)
        {
            if (Math.Abs(lams[i].Real + Sigma) < lineTol) { onLine++; continue; }
            var mirror = new Complex(-lams[i].Real - 2 * Sigma, lams[i].Imaginary);
            bool partner = lams.Where((_, j) => j != i).Any(o => (o - mirror).Magnitude < pairTol);
            if (partner) paired++; else orphans++;
        }
        return (lams.Count, onLine, paired, orphans);
    }

    public string DisplayName => "The F89 octic branch locus is a palindrome (mirror about Re λ = −4, forced by Π)";

    public string Summary
    {
        get
        {
            var (anti, lin) = MirrorClosureResiduals(2.0);
            return $"the EP/diabolic collisions are mirror-symmetric about Re λ = −σ = −4: the octic roots close " +
                   $"under the antiunitary palindrome λ→−λ̄−2σ (residual {anti.ToString("E1", Inv)}) but NOT under the " +
                   $"linear λ→−λ−2σ (residual {lin.ToString("0.0", Inv)}); so the branch locus is a palindrome, forced.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (anti, lin) = MirrorClosureResiduals(2.0);
            yield return new InspectableNode("the gate: the octic obeys the antiunitary mirror, not the linear one",
                summary: $"at q=2 the eight octic roots close under λ↦−λ̄−2σ (reflect Re about −4, keep Im) to " +
                         $"{anti.ToString("E2", Inv)}, and do NOT close under the bare linear λ↦−λ−2σ (Im flipped), " +
                         $"residual {lin.ToString("0.0", Inv)}. Two-sided: a sign or construction error fails it.");

            yield return new InspectableNode("the centre is the palindrome (Re λ = −σ = −4)",
                summary: "−4 is the midpoint of the absorption rungs −2γ (overlap) and −6γ (no-overlap), F1 " +
                         "weight-complement partners n_diff = 1 ↔ 3 (the AbsorptionTheorem). The centre is not a " +
                         "separate fact, it IS the palindrome.");

            var lamEp = new Complex(-4, 2 * GaloisMonodromyWitness.QEp);
            yield return new InspectableNode("the diabolic sits on the line (Re λ = −4), but for a different reason",
                summary: $"the diabolic collision λ_EP = −4γ + 2iJ ({lamEp.Real.ToString("0.0", Inv)}{lamEp.Imaginary.ToString("+0.000;-0.000", Inv)}i " +
                         $"at q_EP={GaloisMonodromyWitness.QEp.ToString("0.000", Inv)}, γ=J=1) is on the line because its pair is overlap-balanced " +
                         "(p=½, the AT-midpoint), an integrability fact, not Π. On-line does NOT imply silent: XXZ Δ≠0 stays on the line yet " +
                         "defects (DIABOLIC_BY_INTEGRABILITY). Line = palindrome's gift; silence = free-fermion's.");

            var (total, onLine, paired, orphans) = BranchLocusStructure();
            yield return new InspectableNode($"every branch point on the line or in a mirror pair ({orphans} orphans)",
                summary: $"over a symmetric box, {total} branch points: {onLine} on the line (self-mirror, Re=−4), " +
                         $"{paired} in mirror pairs across it, {orphans} orphans. The map of the seams is itself a " +
                         "palindrome (reflections/ON_WHO_WATCHES_WHOM.md; the picture is visualizations/f89_octic_branch_locus.png).");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
