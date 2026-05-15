using System.Diagnostics;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>One-shot stress probe: push the native F89 path polynomial pipeline
/// well beyond the int-typed tabulation boundary and emit table-ready Markdown
/// rows for the F89 D_k verification table in PROOF_F89_PATH_D_CLOSED_FORM.md.
/// Run via <c>dotnet test --filter "PipelineStress"</c>.</summary>
public class F89PathPolynomialPipelineStressTest
{
    private readonly ITestOutputHelper _out;
    public F89PathPolynomialPipelineStressTest(ITestOutputHelper output) => _out = output;

    [Fact]
    public void EmitMarkdownTableRows()
    {
        var targets = new[] { 25, 30, 40, 46, 50, 60, 75, 100, 150, 200, 300 };
        _out.WriteLine("| k | v2(k) | odd(k) | FA | deg | E(k) | D_pred | bits | match | factored |");
        _out.WriteLine("|---|-------|--------|-----|-----|------|--------|------|-------|----------|");
        foreach (var k in targets)
        {
            int v2 = V2(k);
            int oddPart = k >> v2;
            int FA = (k + 1) / 2;
            int deg = FA - 1;
            int E = Math.Max(0, (k - 5) / 2) + v2 + Math.Max(0, v2 - 2);

            var formulaD = F89UnifiedFaClosedFormClaim.PredictDenominatorBig(k);
            var (_, pipelineD) = F89PathPolynomialPipeline.Compute(k);
            bool match = formulaD == pipelineD;

            string factored = $"{oddPart}^2·2^{E}";
            int bits = (int)formulaD.GetBitLength();
            string matchStr = match ? "OK" : "FAIL";
            _out.WriteLine($"| {k} | {v2} | {oddPart} | {FA} | {deg} | {E} | {formulaD} | {bits} | {matchStr} | {factored} |");
            Assert.True(match);
        }
    }

    [Fact]
    public void ReportLargeKDenominators()
    {
        var targets = new[] { 50, 60, 80, 100, 150, 200, 300 };
        _out.WriteLine($"{"k",6}  {"FA",4}  {"deg(P_k)",10}  {"#bits(D_k)",12}  {"#digits(D_k)",14}  {"time(ms)",10}");
        _out.WriteLine(new string('-', 80));
        foreach (var k in targets)
        {
            var sw = Stopwatch.StartNew();
            var (coefs, d) = F89PathPolynomialPipeline.Compute(k);
            sw.Stop();
            int FA = (k + 1) / 2;
            int degP = coefs.Length - 1;
            int bits = (int)d.GetBitLength();
            int digits = d.ToString().Length;
            _out.WriteLine($"{k,6}  {FA,4}  {degP,10}  {bits,12}  {digits,14}  {sw.ElapsedMilliseconds,10}");
        }
    }

    private static int V2(int n)
    {
        if (n <= 0) return 0;
        int v = 0;
        while ((n & 1) == 0) { n >>= 1; v++; }
        return v;
    }
}
