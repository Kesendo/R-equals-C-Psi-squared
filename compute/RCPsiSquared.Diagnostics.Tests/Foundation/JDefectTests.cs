using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Ptf;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class JDefectTests
{
    private readonly ITestOutputHelper _out;
    public JDefectTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void JDefect_Axis_StaysPalindromic_AcrossDeltaJ()
    {
        // The J-defect is Π-invariant (Π H Π⁻¹ = −H), so it keeps the spectrum palindromic at every
        // δJ: the contract holds even as the spectrum moves. This is the property that makes it the
        // contrast to the crossover (which freezes the eigenvalues by an exact similarity instead).
        var axis = DimensionAxis.JDefect(N: 4, gamma: 0.5, defectBond: 0, deltaJMax: 0.1, points: 3);

        Assert.Equal("jdefect", axis.Name);
        Assert.Equal(3, axis.Theta.Count);
        Assert.Equal(0.0, axis.Theta[0], 1e-12);
        Assert.Equal(0.1, axis.Theta[^1], 1e-12);
        Assert.Empty(axis.LitSites);

        double sigmaGamma = axis.N * 0.5;
        foreach (double deltaJ in new[] { 0.0, 0.05, 0.1 })
        {
            var H = axis.Hamiltonian(deltaJ);
            var L = PauliDephasingDissipator.BuildZ(H, axis.GammaPerSite);
            double residual = PalindromeResidual.Build(L, axis.N, sigmaGamma).FrobeniusNorm();
            Assert.True(residual < 1e-9,
                $"J-defect must stay palindromic at δJ={deltaJ}; residual {residual:E2}");
        }
    }

    [Fact]
    public void SlowModeMixing_KernelProtected_CoherencesShift_EigenvectorsMix()
    {
        // The in-between of the J-defect is eigenvector mixing, not a pure rotation. Build the
        // unperturbed Liouvillian L_A and the bond perturbation V_L = ∂L/∂J on the defect bond, then
        // read ⟨W_s|V_L|M_s'⟩ on the slow modes. The shift profile tells the story: the kernel (the
        // steady states, Re λ ≈ 0) is protected to machine zero by U(1) conservation; the slow
        // coherences (Re λ ≈ −2γ) shift at first order (the spectrum moves, though it stays
        // palindromic); and the off-diagonal mixing is alive throughout (the eigenvectors mix, unlike
        // the crossover's rigid rotation).
        int N = 4;
        var axis = DimensionAxis.JDefect(N, gamma: 0.5, defectBond: 0);
        var lA = PauliDephasingDissipator.BuildZ(axis.Hamiltonian(0.0), axis.GammaPerSite);
        var vL = BondPerturbation.Build(N, 0, 1, BondPerturbation.Kind.XY);

        var reading = SlowModeMixing.Compute(lA, vL, slowCount: 16);

        _out.WriteLine($"off-diagonal mass = {reading.OffDiagonalMass:E3}");
        for (int s = 0; s < Math.Min(10, reading.DiagonalShiftMagnitudes.Count); s++)
            _out.WriteLine($"  shift[{s}] = {reading.DiagonalShiftMagnitudes[s]:E3}   (Re λ={reading.SlowEigenvalues[s].Real:0.000})");

        // The kernel (Re λ ≈ 0) is protected to machine zero: its first-order eigenvalue shift is ≈ 0.
        double maxKernelShift = 0.0;
        int kernelCount = 0;
        for (int s = 0; s < reading.DiagonalShiftMagnitudes.Count; s++)
            if (Math.Abs(reading.SlowEigenvalues[s].Real) < 1e-6)
            {
                kernelCount++;
                maxKernelShift = Math.Max(maxKernelShift, reading.DiagonalShiftMagnitudes[s]);
            }
        Assert.Equal(N + 1, kernelCount); // N+1 steady states (F4)
        Assert.True(maxKernelShift < 1e-9,
            $"the kernel should be first-order protected (U(1)); max shift {maxKernelShift:E3}");

        // At least one slow coherence shifts at first order: the spectrum moves (no similarity).
        Assert.True(reading.DiagonalShiftMagnitudes.Max() > 0.1,
            $"a slow coherence should shift (the spectrum moves); max shift {reading.DiagonalShiftMagnitudes.Max():E3}");

        // The eigenvectors mix: the off-diagonal coupling is substantial (not a pure rotation).
        Assert.True(reading.OffDiagonalMass > 1e-3,
            $"eigenvectors should mix (off-diagonal alive); got {reading.OffDiagonalMass:E3}");
    }

    [Fact]
    public void JDefectField_Surfaces_Contract_SpectrumMoves_Mixing_Fan()
    {
        var field = new JDefectField(N: 4, gamma: 0.5, defectBond: 0, deltaJMax: 0.1, points: 5, slowCount: 16);
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)field).Children.ToList();
        var labels = children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("contract"));
        Assert.Contains(labels, l => l.Contains("spectrum moves"));
        Assert.Contains(labels, l => l.Contains("mixing"));
        Assert.Contains(labels, l => l.Contains("fan"));

        // The contract holds: the palindrome residual curve is flat at the machine floor across δJ.
        var contract = children.First(c => c.DisplayName.Contains("contract"));
        var curve = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(contract.Payload);
        Assert.All(curve.Y, r => Assert.True(r < 1e-9, $"palindrome residual should be ~0; got {r:E2}"));

        // The spectrum moves, but the movement is reorder-robust: O(δJ), not an index-alignment
        // artifact (sort swaps across moving/crossing eigenvalues would otherwise inflate it to O(Σγ)).
        var moves = children.First(c => c.DisplayName.Contains("spectrum moves"));
        var moveCurve = Assert.IsType<RCPsiSquared.Core.Inspection.InspectablePayload.Curve>(moves.Payload);
        Assert.True(moveCurve.Y.Max() > 1e-6, "the spectrum should move (not frozen like the crossover)");
        Assert.True(moveCurve.Y.Max() < 1.0,
            $"spectrum movement should be O(δJ), reorder-robust, not a sort artifact; got {moveCurve.Y.Max():E2}");

        // Smoke: renders to JSON without throwing and carries the J-defect story.
        var json = RCPsiSquared.Core.Inspection.InspectionJsonExporter.ToJson(field);
        Assert.Contains("palindrome", json);
        Assert.Contains("mixing", json);
    }
}
