using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Head-1 "function thread": track the defective eigenvalue branch of the (1,2) block as a smooth
/// function of real q and check whether ⟨Ô⟩(q) (= 3/2 + Re λ/4, the AT identity) has a nicer form than the
/// non-radical per-locus values. The AT identity ⟨Ô⟩_fromλ = ⟨Ô⟩_fromVec is a genuine regression gate along
/// the whole branch; the CSV dump is the exploration artifact.</summary>
public class SectorBraidBranchSweepTests
{
    private readonly ITestOutputHelper _out;
    public SectorBraidBranchSweepTests(ITestOutputHelper o) => _out = o;

    [Fact(DisplayName = "N=5 defective branch: AT identity holds along the branch + CSV dump")]
    [Trait("Category", "SLOW_MSM")]
    public void Sweep_N5_DefectiveBranch_ATIdentity_AndDump()
    {
        // Seed at the known real defective locus (locus 1) and its defective λ.
        var samples = SectorBraidModeGeometry.SweepDefectiveBranch(
            N: 5, qStart: 0.620878, lambdaStart: new Complex(-4.6189, 0.0),
            qMin: 0.25, qMax: 1.6, steps: 270);

        Assert.NotEmpty(samples);

        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder("q,ReLambda,ImLambda,Ohat_fromLambda,Ohat_fromVec,ndiff\n");
        double maxAtResidual = 0;
        foreach (var s in samples)
        {
            sb.Append($"{s.Q.ToString("0.######", inv)},{s.Lambda.Real.ToString("0.########", inv)},");
            sb.Append($"{s.Lambda.Imaginary.ToString("0.########", inv)},{s.OhatFromLambda.ToString("0.########", inv)},");
            sb.Append($"{s.OhatFromVec.ToString("0.########", inv)},{s.NDiff.ToString("0.########", inv)}\n");
            // AT gate only where the branch is real (Im λ ~ 0); on a complex pair ⟨Ô⟩_fromVec is the
            // single-vector overlap and Re λ = ⟨Ô⟩-from-λ still holds, so compare unconditionally.
            maxAtResidual = Math.Max(maxAtResidual, Math.Abs(s.OhatFromLambda - s.OhatFromVec));
        }

        string path = Path.Combine(Path.GetTempPath(), "sectorbraid_branch_sweep_N5.csv");
        File.WriteAllText(path, sb.ToString());
        _out.WriteLine($"CSV: {path}  ({samples.Count} points)  maxATresidual={maxAtResidual:E3}");

        // The AT identity Re λ = -2⟨n_diff⟩ is exact at real q (C anti-Hermitian): ⟨Ô⟩ from λ == from vec.
        Assert.True(maxAtResidual < 1e-9, $"AT identity broken along branch: max residual {maxAtResidual:E3}");
    }
}
