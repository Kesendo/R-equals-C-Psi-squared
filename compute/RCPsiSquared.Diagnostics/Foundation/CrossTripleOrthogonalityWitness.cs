using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live F127 witness (<c>inspect --root crosstriple</c>): recomputes, at inspect time,
/// an INDEPENDENT C# slice of the cross-triple orthogonality along the residue-collapse chain of
/// <c>docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md</c>. Three recomputations, each with its own
/// discriminating control:
/// <list type="bullet">
/// <item>the full form 𝔉 vanishes on the two-constraint variety over split primes, with an
/// off-variety nonzero control (<see cref="CrossFormCertificate.CertifySlice"/>);</item>
/// <item>the §3 core function T vanishes on the three-constraint variety, with an off-variety
/// nonzero control (<see cref="CrossFormCertificate.CertifyCoreIdentitySlice"/>);</item>
/// <item>the §2 sheet lattice: 72 atoms → 288 pole events → 32 sheets, all coefficients ±1, nine
/// events per sheet, exact integer combinatorics (<see cref="SheetLattice"/>).</item>
/// </list>
/// This is the second implementation the F127 code-trust caveat asks for, at sample scale; the proof
/// objects remain the 527/527 grid+CRT wall and the six committed <c>f127_*.py</c> gates.
/// Sub-second; fully deterministic (xorshift streams seeded from the prime).</summary>
public sealed class CrossTripleOrthogonalityWitness : IInspectable
{
    // the wall's own first two 30-bit split primes (core_grid.gen_primes(17)[0..1]; both = 1 mod 4)
    private static readonly long[] Primes = { 1073741833L, 1073741857L };
    private const int SamplesPerPrime = 6;
    private const int CoreSamplesPerPrime = 12;

    public string DisplayName =>
        "F127 cross-triple orthogonality (live: 𝔉-slice + core-identity T + sheet lattice)";

    public string Summary
    {
        get
        {
            int onVariety = 0, bad = 0, controls = 0, samples = 0;
            int coreEvals = 0, coreBad = 0, coreCtrl = 0, coreCtrlEval = 0;
            foreach (long p in Primes)
            {
                var (ev, b, c, s) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
                onVariety += ev; bad += b; controls += c; samples += s;
                var (cev, cb, cc, cce) = CrossFormCertificate.CertifyCoreIdentitySlice(p, CoreSamplesPerPrime);
                coreEvals += cev; coreBad += cb; coreCtrl += cc; coreCtrlEval += cce;
            }
            var lattice = SheetLattice.Analyze();
            bool latticeOk = lattice.Atoms == SheetLattice.AtomCount &&
                             lattice.Events == SheetLattice.EventCount &&
                             lattice.AllCoefficientsPmOne &&
                             lattice.DistinctSheets == SheetLattice.SheetCount &&
                             lattice.MinEventsPerSheet == SheetLattice.EventsPerSheet &&
                             lattice.MaxEventsPerSheet == SheetLattice.EventsPerSheet &&
                             lattice.EachSheetOneEventPerBlock;
            bool pass = bad == 0 && controls >= (int)(0.9 * samples) &&
                        coreBad == 0 && coreCtrl >= (int)(0.8 * coreCtrlEval) && latticeOk;
            return $"{(pass ? "PASS" : "FAIL")}: 𝔉-slice {onVariety} zeros ({bad} bad), " +
                   $"{controls}/{samples} controls nonzero; core-T {coreEvals} zeros ({coreBad} bad), " +
                   $"{coreCtrl}/{coreCtrlEval} controls nonzero; lattice {lattice.DistinctSheets} sheets × " +
                   $"{lattice.MaxEventsPerSheet} events {(latticeOk ? "OK" : "BAD")}";
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return SheetLatticeNode();
            yield return CoreIdentityNode();
            yield return CrossFormNode();
            yield return GatesNode();
            yield return new InspectableNode("What this is",
                summary: "the residue-collapse chain recomputed live in C#: the §2 lattice (exact integer), " +
                         "the §3 core identity T (GF(p), independent of 𝔉), and the 𝔉-slice; a SLICE of F127, " +
                         "not the wall (that is the 527/527 sweep and the six committed f127_* gates)");
        }
    }

    private static InspectableNode SheetLatticeNode()
    {
        var r = SheetLattice.Analyze();
        return new InspectableNode("§2 sheet lattice (exact integer)",
            summary: $"{r.Atoms} atoms → {r.Events} pole events → {r.DistinctSheets} sheets; " +
                     $"all coefficients ±1: {r.AllCoefficientsPmOne}; " +
                     $"{r.MinEventsPerSheet}..{r.MaxEventsPerSheet} events per sheet; " +
                     $"one event per (i,j) block per sheet: {r.EachSheetOneEventPerBlock}",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode CoreIdentityNode()
    {
        int evals = 0, bad = 0, ctrl = 0, ctrlEval = 0;
        var perPrime = new List<IInspectable>();
        foreach (long p in Primes)
        {
            var (ev, b, c, ce) = CrossFormCertificate.CertifyCoreIdentitySlice(p, CoreSamplesPerPrime);
            evals += ev; bad += b; ctrl += c; ctrlEval += ce;
            perPrime.Add(new InspectableNode($"p = {p}",
                summary: $"{ev} on-variety points, T = 0 (nonzero: {b}, must be 0); " +
                         $"{c}/{ce} off-variety controls nonzero (must be ≥ 80%)",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("§3 core identity T (live GF(p))",
            summary: $"{evals} three-constraint variety points, {bad} nonzero (must be 0); " +
                     $"{ctrl}/{ctrlEval} off-variety controls nonzero; T independent of 𝔉",
            children: perPrime, provenance: NodeProvenance.Live);
    }

    private static InspectableNode CrossFormNode()
    {
        int onVariety = 0, bad = 0, ctrl = 0, samples = 0;
        var perPrime = new List<IInspectable>();
        foreach (long p in Primes)
        {
            var (ev, b, c, s) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
            onVariety += ev; bad += b; ctrl += c; samples += s;
            long iRoot = CrossFormCertificate.SqrtMinusOne(p);
            perPrime.Add(new InspectableNode($"p = {p}",
                summary: $"i = {iRoot}; {ev} on-variety evaluations, {b} nonzero (must be 0); " +
                         $"{c}/{s} off-variety controls nonzero (must be ≥ 90%)",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("𝔉-slice (live GF(p))",
            summary: $"{onVariety} on-variety evaluations, {bad} nonzero (must be 0); " +
                     $"{ctrl}/{samples} off-variety controls nonzero",
            children: perPrime, provenance: NodeProvenance.Live);
    }

    private static InspectableNode GatesNode() =>
        new("The proof chain (committed gates)",
            summary: "docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md + six f127_* gate scripts",
            children: new IInspectable[]
            {
                new InspectableNode("PROOF_F127_RESIDUE_COLLAPSE.md",
                    summary: "the derived chain: sheet lattice (§2), core identity by resultant divisibility (§3), " +
                             "transport (§4), the window + oddness (§5), the mirror anchor (§6)"),
                new InspectableNode("f127_core_identity.py",
                    summary: "§3: naked-T falsification, N assembly (540 terms), Res | E exact, branch check"),
                new InspectableNode("f127_core_locus_patch.py",
                    summary: "§3.2: forcing chain, gcd(Res, S_w−m₀S_z) = 1, degenerate sublocus r₀ = r₁ = 0"),
                new InspectableNode("f127_exact_gates.py",
                    summary: "S1 (72 atoms), S2 (288 events), S3/S5/S6 (oddness bijections), S4 (mirror gate); LEM/tanLEM"),
                new InspectableNode("f127_transport_lemma.py",
                    summary: "§4: incidence table, 31 discriminants, norm identity (62), step-2 divisibility (384)"),
                new InspectableNode("f127_relabel_pin.py",
                    summary: "§5: the literal w₂-distinguished rerun (inventory, all poles simple, window ⊆ {−1,0,1})"),
                new InspectableNode("f127_sheet_lattice.py",
                    summary: "§2 numerics: the 32×9 lattice survey and the residue-bookkeeping gates"),
            });
}
