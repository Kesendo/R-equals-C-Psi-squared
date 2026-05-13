using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Symmetry;

var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
var anatomy = C2FullBlockSigmaAnatomy.Build(block);

double faRateTarget = -2.0 * block.GammaZero;
double faRateTolerance = 1e-5 * block.GammaZero;
double q = 1.5;
double j = q * block.GammaZero;

Console.WriteLine($"N={block.N}, γ₀={block.GammaZero}, Q={q}");
Console.WriteLine($"F_a rate target: {faRateTarget:G6}, tolerance: {faRateTolerance:G6}");
Console.WriteLine($"J = {j:G6}");

var rateMatches = anatomy.SigmaSpectrum.Where(w => Math.Abs(w.EigenvalueReal - faRateTarget) <= faRateTolerance).ToList();
Console.WriteLine($"\nModes with Re(λ) ≈ {faRateTarget:G6}: {rateMatches.Count}");
foreach(var w in rateMatches)
{
    Console.WriteLine($"  λ = {w.EigenvalueReal:G6} {w.EigenvalueImag:+G6;-G6}i");
}

var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(7);
Console.WriteLine($"\nBloch orbit for N=7: {string.Join(", ", orbit)}");
foreach(int n in orbit)
{
    double yn = F89PathKAtLockMechanismClaim.BlochEigenvalueY(7, n);
    Console.WriteLine($"  n={n}: y_n = {yn:G6}");
}

Console.WriteLine("\nAll eigenvalues (real parts only):");
var realParts = anatomy.SigmaSpectrum.Select(w => w.EigenvalueReal).OrderBy(x => x).ToList();
foreach(var re in realParts)
{
    Console.WriteLine($"  {re:G6}");
}
