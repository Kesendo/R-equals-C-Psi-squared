using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the sideways spin ladder (typed: <see cref="SidewaysSpinLadderClaim"/>;
/// open arc <c>sideways_spin_ladder</c>). Where the claim carries the closed forms (ℓ = (N−3)/2, the CG
/// norms √(ℓ(ℓ+1) − m(m+1)), orbit 4N−8 = 4(N−2)), this witness RUNS the chain walk at inspect time
/// through the shared construction <see cref="SidewaysSpinLadderChain"/> (the same one the gate asserts
/// on, so the gate's numbers and the inspect-time numbers cannot drift apart).
///
/// <para>Be precise about what is derived and what is recomputed. The seed (q*, λ_A) is READ FROM the
/// recorded census table (<c>RealDefectiveSeeds</c>, not rediscovered); the CG coefficients are WRITTEN
/// FROM the closed form. What is genuinely recomputed live: both ladders' intertwining residuals on every
/// rung of the two chains (compared to 0.0 exactly), the block spectra and the derived tolerance, the
/// fold/band sector sets, and the transport norms of the fold eigenvector, rung by rung, terminal ratio
/// and survivors control included.</para>
///
/// <para>Guard: N ∈ {5, 7}, the two N with a recorded seed. N=5 is sub-second; N=7 builds and
/// eigendecomposes blocks up to 1225² and takes ~half a minute at inspect time.</para></summary>
public sealed class SidewaysSpinLadderWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public int N { get; }
    private readonly SidewaysSpinLadderChain.LadderRun _run;

    public SidewaysSpinLadderWitness(int n = 5)
    {
        if (n != 5 && n != 7)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"N must be 5 or 7 (the dense chain walk); got {n}. " +
                "N=9 lives in the sparse walk (SidewaysSpinLadderSparse, LU shift-invert past the dense " +
                "wall) and is gated in SidewaysSpinLadderSparseTests (SLOW_SIDEWAYS9), confirmed 2026-08-09.");
        N = n;
        _run = SidewaysSpinLadderChain.Run(n);
    }

    public string DisplayName =>
        $"SidewaysSpinLadderWitness (S⁺ chain walk live, N={N}, q*={_run.QStar.ToString("0.######", Inv)})";

    public string Summary
    {
        get
        {
            bool shape = _run.FoldSet.SequenceEqual(_run.PredictedSet);
            int count = _run.FoldSet.Count + _run.BandSet.Count;
            var norms = _run.Chains.SelectMany(c => c.Rungs.Select(r => r.Norm)).ToList();
            bool cg = _run.Chains.All(c =>
                c.Rungs.Count == _run.CgPredicted.Count
                && c.Rungs.Select(r => r.Norm).Zip(_run.CgPredicted).All(t => Math.Abs(t.First - t.Second) < 1e-6));
            return $"N={N}: CONTROL worst residual {_run.WorstControlResidual.ToString("E1", Inv)} " +
                   $"({(_run.WorstControlResidual == 0.0 ? "exactly 0.0" : "NOT exact")}), " +
                   $"fold sectors {(shape ? "= the two chain interiors" : "≠ prediction")}, " +
                   $"fold+band {count} vs 4N−8 = {4 * N - 8}, " +
                   $"norms {(cg ? "= CG coefficients" : "≠ CG")} " +
                   $"[{string.Join(", ", norms.Select(x => x.ToString("0.000000", Inv)).Distinct())}]; " +
                   "seed from RealDefectiveSeeds, CG from the closed form, everything else recomputed live.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "CONTROL: both ladders intertwine, residual vs 0.0 exactly",
                summary: $"worst residual of L_target·M − M·L_source over Φ and S⁺ on every rung of both " +
                         $"chains p+q̃ = N∓1: {_run.WorstControlResidual.ToString("E1", Inv)} " +
                         $"({(_run.WorstControlResidual == 0.0 ? "EXACTLY 0.0, the operator identity cancels pairwise" : "NOT exactly zero: a construction finding, not a tolerance case")}). " +
                         "Φ's half is what both proofs derive, so it is the control that makes the S⁺ half readable.",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "SHAPE + COUNT: the fold family is the two chain interiors, fold+band = 4N−8",
                summary: $"interior blocks carrying the cross-fold partner −λ_A − 2N = {_run.Fold.ToString("0.####", Inv)}: " +
                         $"[{string.Join(", ", _run.FoldSet)}] vs the two S⁺ chain interiors at p+q̃ = N∓1: " +
                         $"[{string.Join(", ", _run.PredictedSet)}] " +
                         $"({(_run.FoldSet.SequenceEqual(_run.PredictedSet) ? "MATCH" : "MISMATCH")}); " +
                         $"band sectors {_run.BandSet.Count}, fold+band = {_run.FoldSet.Count + _run.BandSet.Count} " +
                         $"vs the F125 orbit size 4N−8 = {4 * N - 8}; tolerance {_run.Tolerance.ToString("E1", Inv)}, " +
                         "DERIVED from the run (λ_A is recorded to 4 decimals).",
                provenance: NodeProvenance.Live);

            foreach (var chain in _run.Chains)
            {
                var norms = chain.Rungs.Select(r => r.Norm).ToList();
                yield return new InspectableNode(
                    displayName: $"chain p+q̃ = {chain.Total}: norms vs CG, terminal death, survivors",
                    summary: $"interior transport norms [{string.Join(", ", norms.Select(x => x.ToString("0.000000", Inv)))}] " +
                             $"vs CG √(ℓ(ℓ+1)−m(m+1)) [{string.Join(", ", _run.CgPredicted.Select(x => x.ToString("0.000000", Inv)))}] " +
                             $"for ℓ = {SidewaysSpinLadderClaim.ChainSpin(N).ToString("0.#", Inv)}; " +
                             $"terminal ‖S⁺v‖ = {chain.Terminal.ToString("E1", Inv)}, ratio to interior " +
                             $"{chain.TerminalRatio.ToString("E1", Inv)} (the eigensolver floor, NOT an exact zero: " +
                             $"v comes from an eigensolver, the highest-weight death S⁺|ℓ,ℓ⟩ = 0 is read as this ratio); " +
                             $"survivors {chain.Survivors} of {chain.TerminalSourceDim} " +
                             $"{(chain.Survivors >= chain.TerminalTargetDim ? "≥" : "BELOW")} target dim {chain.TerminalTargetDim} " +
                             "(negative control: the terminal map is not the zero map).",
                    provenance: NodeProvenance.Live);
            }

            yield return new InspectableNode(
                displayName: "N=9 confirmed; even N stays open",
                summary: "N=9, ℓ=3 was confirmed 2026-08-09 on both chains to six decimals " +
                         "(√6, √10, √12, √12, √10, √6; the 10584²/15876² middle blocks walked by LU " +
                         "shift-invert, SidewaysSpinLadderSparse, gate SLOW_SIDEWAYS9). Even N stays open: no " +
                         "even-N defective seed is recorded (RealDefectiveSeeds lists odd N, lower bounds), so " +
                         "the even-N chain reading is untested, not predicted.",
                provenance: NodeProvenance.Stored);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
