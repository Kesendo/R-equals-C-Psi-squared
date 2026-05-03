using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Probe: dump the actual Klein-cell aggregate masses for each f83 fingerprint
/// case, to surface whether 1/2 (the framework's recurring half-anchor: F81 50/50 split,
/// F83 anti-fraction at r=0, balanced Π² partition) appears as a structural quantity in
/// the slow-mode × observable-cell overlap.</summary>
public class Pi2KleinHalfProbe
{
    private readonly ITestOutputHelper _out;
    public Pi2KleinHalfProbe(ITestOutputHelper output) => _out = output;

    private static ChainSystem Chain3() => new(N: 3, J: 1.0, GammaZero: 0.05);

    private static ComplexMatrix BuildL(IReadOnlyList<PauliPairBondTerm> terms, ChainSystem chain)
    {
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, terms.ToBilinearSpec(chain.J)).ToMatrix();
        var gammaList = Enumerable.Repeat(chain.GammaZero, chain.N).ToArray();
        return PauliDephasingDissipator.BuildZ(H, gammaList);
    }

    private static double CellMass(KleinSpectralMode m, int z, int x) =>
        (z, x) switch
        {
            (+1, +1) => m.MassPp,
            (+1, -1) => m.MassPm,
            (-1, +1) => m.MassMp,
            (-1, -1) => m.MassMm,
            _ => 0,
        };

    [Fact(DisplayName = "Probe: f83 slow-mode aggregate in observable cell, all 4 H-classes")]
    public void DumpAggregates()
    {
        var cases = new (string label, PauliPairBondTerm[] terms, PauliLetter[] obs)[]
        {
            ("Truly XX+YY → Y₀IZ₂",
                new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.Y, PauliLetter.Y) },
                new[] { PauliLetter.Y, PauliLetter.I, PauliLetter.Z }),
            ("Pi2EvenNonTruly YZ+ZY → X₀IX₂",
                new[] { Term(PauliLetter.Y, PauliLetter.Z), Term(PauliLetter.Z, PauliLetter.Y) },
                new[] { PauliLetter.X, PauliLetter.I, PauliLetter.X }),
            ("Pi2OddPure XY+YX → X₀IZ₂",
                new[] { Term(PauliLetter.X, PauliLetter.Y), Term(PauliLetter.Y, PauliLetter.X) },
                new[] { PauliLetter.X, PauliLetter.I, PauliLetter.Z }),
            ("Mixed XX+XY → Z₀IX₂",
                new[] { Term(PauliLetter.X, PauliLetter.X), Term(PauliLetter.X, PauliLetter.Y) },
                new[] { PauliLetter.Z, PauliLetter.I, PauliLetter.X }),
        };

        var chain = Chain3();
        _out.WriteLine($"# Klein-cell mass per slow-mode budget at N={chain.N}, J={chain.J}, γ_Z={chain.GammaZero}");
        _out.WriteLine($"# 4^N = 64 modes total; aggregating over slowest 16 = quarter, then full 64.");
        _out.WriteLine($"# Each mode normalised to total Klein-mass = 1.");
        _out.WriteLine($"# Threshold reference: 1/2 = 0.5 (framework's recurring half-anchor).");
        _out.WriteLine("");

        foreach (var (label, terms, obs) in cases)
        {
            var L = BuildL(terms, chain);
            var modes = Pi2KleinSpectralView.ComputeFor(L, chain.N)
                .OrderBy(m => Math.Abs(m.Eigenvalue.Real)).ToList();
            int eigZ = PiOperator.SquaredEigenvalue(obs, PauliLetter.Z);
            int eigX = PiOperator.SquaredEigenvalue(obs, PauliLetter.X);
            string cellLabel = (eigZ, eigX) switch
            {
                (+1, +1) => "Pp", (+1, -1) => "Pm",
                (-1, +1) => "Mp", (-1, -1) => "Mm",
                _ => "??",
            };

            double slow16 = modes.Take(16).Sum(m => CellMass(m, eigZ, eigX));
            double full64 = modes.Sum(m => CellMass(m, eigZ, eigX));
            double slow16AsHalfMultiple = slow16 / 0.5;
            double full64AsHalfMultiple = full64 / 0.5;

            _out.WriteLine($"## {label}");
            _out.WriteLine($"  observable Klein-cell:           {cellLabel}");
            _out.WriteLine($"  slow-16 aggregate in cell:       {slow16:F6}     (= {slow16AsHalfMultiple:F4} × 1/2)");
            _out.WriteLine($"  full-64 aggregate in cell:       {full64:F6}     (= {full64AsHalfMultiple:F4} × 1/2)");
            // Also dump per-cell distribution of slow-16:
            double slowPp = modes.Take(16).Sum(m => m.MassPp);
            double slowPm = modes.Take(16).Sum(m => m.MassPm);
            double slowMp = modes.Take(16).Sum(m => m.MassMp);
            double slowMm = modes.Take(16).Sum(m => m.MassMm);
            _out.WriteLine($"  slow-16 cell distribution:       Pp={slowPp:F4} Pm={slowPm:F4} Mp={slowMp:F4} Mm={slowMm:F4}");
            _out.WriteLine($"  slow-16 sum (sanity):            {slowPp + slowPm + slowMp + slowMm:F4} (should be 16)");
            _out.WriteLine("");
        }
    }

    private static PauliPairBondTerm Term(PauliLetter a, PauliLetter b) => new(a, b);
}
