using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Theorem A's D-half (remainder (c) of PROOF_CODIM1_BY_ADDITIVITY): is a residual diabolic
/// coalescence "twin-scalar", i.e. is the dephasing D scalar on the coalescing 2-plane?
///
/// <para>The 2026-07-02 survey reframed this: given the H-half (single-multiplet descent ⟹ qC₂ scalar) and
/// a coalescence (M₂ = A₂ + qC₂ has a double eigenvalue), the D-half FOLLOWS, because A₂ = M₂ − qC₂ is
/// Hermitian with a repeated eigenvalue, hence scalar. So at REAL q a diabolic is automatically twin-scalar
/// (the Hermitian/anti-Hermitian split of λI), while at COMPLEX q a diabolic could carry A₂ and qC₂ both
/// non-scalar summing to λI (a coincidence not from additivity). This gate measures the DEPHASE (D-half) and
/// HOP (H-half) scalar-departures on the coalescing plane at two controls and at the N=5 complex-q diabolics.</para>
///
/// <para>Run: <c>dotnet test compute/RCPsiSquared.Diagnostics.Tests --filter "Category=TWINSCALAR"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class TwinScalarDHalfTests
{
    private readonly ITestOutputHelper _out;
    public TwinScalarDHalfTests(ITestOutputHelper o) => _out = o;

    private void Report(string tag, Complex q, Complex lambda, EpCharacter.PencilReading r)
    {
        _out.WriteLine($"{tag}: q={q}  λ={lambda}");
        _out.WriteLine($"    kind={r.Character.Kind} alg={r.Character.Algebraic} geo={r.Character.Geometric} " +
                       $"dep={r.Character.Departure:E2}");
        _out.WriteLine($"    D-half(dephase) scalarDep={r.DephaseScalarDeparture:E3}  " +
                       $"H-half(hop) scalarDep={r.HopScalarDeparture:E3}  antiHerm(A₂)={r.DephaseAntiHermiticity:E3}");
    }

    // Control 1 (positive) + noise-floor calibration: the N=4 (k=3) octic diabolic, at REAL q. The theory
    // says a real-q diabolic is EXACTLY twin-scalar (D-half = 0). We find it by the same scan the 11 came
    // from (so the locus precision matches), then measure the D-half at the refined locus AND at loci
    // perturbed in λ by 1e-3 / 1e-2, to see how much D-half a locus error induces (the noise floor against
    // which the complex-q departures must be judged).
    [Fact]
    [Trait("Category", "TWINSCALAR")]
    public void N4_OcticDiabolic_Calibration()
    {
        var pts = PathKMonodromyScout.FindDiabolicsExact(
            k: 3, reLo: 0.30, reHi: 1.10, imLo: -0.03, imHi: 0.03, cell: 0.02);
        _out.WriteLine($"N=4 (k=3) scan found {pts.Count} coalescence points on the real-q strip");
        var diab = pts.Where(p => p.IsSemisimple).OrderBy(p => p.Gap).FirstOrDefault();
        Assert.NotNull(diab);
        _out.WriteLine($"picked diabolic: q={diab!.QValue} λ={diab.MergeLambda} gap={diab.Gap:E3} " +
                       $"loopId={diab.LoopIsIdentity} expo={diab.GapScalingExponent:F3}");

        var r = PathKMonodromyScout.CharacterizeCoherencePencilAt(4, 2, 1, diab.QValue, diab.MergeLambda);
        Report("N=4 octic diabolic (refined)", diab.QValue, diab.MergeLambda, r);

        Assert.True(r.DephaseAntiHermiticity < 1e-9,
            $"the dephase restriction must be Hermitian; got {r.DephaseAntiHermiticity:E3}");
        Assert.Equal(EpCharacter.EpKind.Diabolic, r.Character.Kind);
        // a real-q diabolic is twin-scalar (D-half ~ gap): the calibration "zero" against which the
        // complex-q loci are judged. Computed in the HS-orthonormal coherence basis.
        Assert.True(r.DephaseScalarDeparture < 1e-6,
            $"N=4 real-q diabolic must be twin-scalar; got D-half={r.DephaseScalarDeparture:E3}");
        _out.WriteLine($"N=4 TWIN-SCALAR CALIBRATION: D-half = {r.DephaseScalarDeparture:E3} (gap {diab.Gap:E2})");
    }

    // Control 2 (negative): the N=5 (k=4) REAL DEFECTIVE seed q* = 0.620878, λ = −4.6189. A defective EP is
    // NOT twin-scalar (the dephase restriction carries the Jordan off-diagonal). Validates the discriminator.
    [Fact]
    [Trait("Category", "TWINSCALAR")]
    public void N5_RealDefectiveSeed_IsNotTwinScalar()
    {
        var q = new Complex(0.620878, 0);
        var lambda = new Complex(-4.6189, 0);
        var r = PathKMonodromyScout.CharacterizeCoherencePencilAt(5, 2, 1, q, lambda);
        Report("N=5 real defective seed", q, lambda, r);
        Assert.Equal(EpCharacter.EpKind.Defective, r.Character.Kind);
        Assert.True(r.DephaseScalarDeparture > 1e-3,
            $"a defective EP must NOT be twin-scalar; got {r.DephaseScalarDeparture:E3}");
    }

    // the mean of the closest residual-root pair at q, and their gap (the coalescence precision).
    private static (Complex Lambda, double Gap) CoalescingPair(int k, Complex q)
    {
        var roots = PathKMonodromyScout.ResidualRootsExact(k, q);
        int ai = -1, bi = -1; double best = double.PositiveInfinity;
        for (int i = 0; i < roots.Length; i++)
            for (int j = i + 1; j < roots.Length; j++)
            {
                double g = (roots[i] - roots[j]).Magnitude;
                if (g < best) { best = g; ai = i; bi = j; }
            }
        return ((roots[ai] + roots[bi]) / 2, best);
    }

    // Control 3: at PURE-IMAGINARY q the coherence block is Hermitian (the anti-Hermitian coherent part
    // q·C becomes real-symmetric when q = iy), so its degeneracies are semisimple AUTOMATICALLY, for a
    // reason distinct from twin-scalar. This is why the imaginary-q diabolics (λ real) are non-twin-scalar
    // yet diabolic: Hermiticity, not additivity. Confirms the mechanism that explains that family.
    [Fact]
    [Trait("Category", "TWINSCALAR")]
    public void ImaginaryQ_CoherenceBlockIsHermitian()
    {
        foreach (double y in new[] { 0.381, 0.879, 1.695 })
        {
            var l = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(5, 2, 1, new Complex(0, y)));
            double antiHerm = (l - l.ConjugateTranspose()).FrobeniusNorm() / l.FrobeniusNorm();
            _out.WriteLine($"q={y}i: ‖L − Lᴴ‖/‖L‖ = {antiHerm:E3}");
            Assert.True(antiHerm < 1e-14, $"at pure-imaginary q the block must be Hermitian; got {antiHerm:E3}");
        }
    }

    // The full survey with a NOISE-FLOOR discriminant. Enumerate the k=4 coalescences over the tabulated
    // region, DEEPLY REFINE each locus (drive the gap far below the scan's 1e-3 acceptance), then measure the
    // D-half at the refined locus. The verdict is D-half vs gap: twin-scalar ⟹ D-half ~ O(gap) → 0 as the
    // locus sharpens (calibrated by N=4: gap 3e-9 → D-half 4e-10); a genuine complex-q coincidence ⟹ D-half
    // stays ≫ gap (independent of the locus precision). Report-only; the split is the answer.
    [Fact]
    [Trait("Category", "TWINSCALAR")]
    public void N5_ComplexQDiabolics_TwinScalarSurvey()
    {
        var points = PathKMonodromyScout.FindDiabolicsExact(
            k: 4, reLo: -0.05, reHi: 2.2, imLo: -0.05, imHi: 1.6, cell: 0.05);
        _out.WriteLine($"FindDiabolicsExact(k=4) returned {points.Count} coalescence points; refining each deeper");
        int diabTwin = 0, diabCoin = 0, other = 0;
        foreach (var d in points.OrderBy(p => p.QValue.Magnitude))
        {
            var qref = PathKMonodromyScout.GapRefineResidualLocalAt(4, d.QValue, cell: 0.02);
            var (lam, gap) = CoalescingPair(4, qref);
            var r = PathKMonodromyScout.CharacterizeCoherencePencilAt(5, 2, 1, qref, lam);
            double ratio = r.DephaseScalarDeparture / Math.Max(gap, 1e-300);
            string verdict;
            if (r.Character.Kind == EpCharacter.EpKind.Diabolic)
            {
                // genuine coincidence iff D-half is far above the gap (not a finite-split artifact)
                bool twin = ratio < 10.0;
                verdict = twin ? "twin-scalar (D-half ~ gap)" : "COINCIDENCE (D-half >> gap)";
                if (twin) diabTwin++; else diabCoin++;
            }
            else { verdict = r.Character.Kind.ToString(); other++; }
            var de = r.DephaseEigenvalues.Select(e => e.Real).OrderBy(x => x).ToArray();
            string rates = de.Length == 2 ? $"[{de[0]:F3}, {de[1]:F3}] Δ={de[1] - de[0]:F3}" : string.Join(",", de.Select(x => x.ToString("F3")));
            _out.WriteLine($"q={qref}  λ={lam}  gap={gap:E2}  D-half={r.DephaseScalarDeparture:E3}  " +
                           $"D/gap={ratio:E1}  dephaseRates={rates}  {verdict}");
        }
        _out.WriteLine($"SUMMARY (diabolics): {diabTwin} twin-scalar (D-half~gap), " +
                       $"{diabCoin} genuine COINCIDENCE (D-half>>gap); ({other} non-diabolic)");
        Assert.NotEmpty(points);
    }
}
