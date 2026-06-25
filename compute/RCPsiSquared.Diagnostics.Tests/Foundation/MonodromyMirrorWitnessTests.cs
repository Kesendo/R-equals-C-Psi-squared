using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The mirror meets the monodromy: does the q ↦ −q̄ reflection intertwine the F89 path-3 octic
/// braiding? Exploratory gate. This file asserts only the two spectral sanities (which the palindrome
/// guarantees) and PRINTS the per-EP intertwining verdicts, so the precise true statement can be read
/// before the F89MonodromyMirrorClaim text is fixed.</summary>
public class MonodromyMirrorWitnessTests
{
    private readonly ITestOutputHelper _out;
    public MonodromyMirrorWitnessTests(ITestOutputHelper o) => _out = o;

    [Fact]
    public void MirrorMonodromy_SpectralSanitiesHold_AndReportIntertwining()
    {
        var (specK, specT, reports, all) = GaloisMonodromyWitness.MirrorMonodromy();
        _out.WriteLine($"=== C-K (q -> -conj(q)) ===");
        _out.WriteLine($"specK(conj@+2 vs @-2)={specK:E3}  specT(lam->-conj(lam)-8 @+2)={specT:E3}");
        foreach (var r in reports)
        {
            var sk = r.TauPlus.Select(s => r.SigmaK[s]).OrderBy(x => x).ToArray();
            _out.WriteLine($"q={r.Q}  ReLamEP={r.LambdaMidRe:F3}  tau+={Fmt(r.TauPlus)} tau-={Fmt(r.TauMinus)} " +
                           $"sigmaK(tau+)={Fmt(sk)} intertwines={r.Intertwines}");
        }
        _out.WriteLine($"sigmaK={Fmt(GetSigmaK(reports))}  reports={reports.Count} allIntertwine={all}");

        var (sigmaT, fixedPts, twoCyc, braidInv, nBraids) = GaloisMonodromyWitness.PalindromeStrandPairing();
        _out.WriteLine($"=== C-T (fixed-q palindrome lam -> -conj(lam)-8) ===");
        _out.WriteLine($"sigmaT={Fmt(sigmaT)}  fixedPoints={fixedPts} (strands on Re=-4)  twoCycles={twoCyc}  " +
                       $"braidSetInvariantUnderSigmaT={braidInv}  nBraids={nBraids}");

        Assert.True(specK < 1e-9, $"L(2)*=L(-2) family symmetry (got {specK:E3})");
        Assert.True(specT < 1e-9, $"fixed-q palindrome lam->-conj(lam)-8 (got {specT:E3})");
    }

    private static int[] GetSigmaK(System.Collections.Generic.IReadOnlyList<GaloisMonodromyWitness.MirrorReport> reports)
        => reports.Count > 0 ? reports[0].SigmaK : System.Array.Empty<int>();

    private static string Fmt(int[] a) => "(" + string.Join(" ", a) + ")";
}
