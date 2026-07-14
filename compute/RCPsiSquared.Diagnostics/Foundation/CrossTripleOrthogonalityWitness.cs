using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live F127 witness (<c>inspect --root crosstriple</c>): recomputes, at inspect
/// time, an INDEPENDENT C# slice of the cross-triple orthogonality via
/// <see cref="CrossFormCertificate"/>: deterministic variety points over split primes, all four
/// (z₃, w₃) root combinations each, must evaluate 𝔉 = 0 in GF(p); one off-variety control per
/// point must be NONZERO (the discriminating half: a zero function would pass the first check).
/// This is the second implementation the F127 code-trust caveat asks for, at sample scale; the
/// proof object remains the 527/527 grid+CRT wall (<c>grid_proof_sweep.py --assert</c>).
/// Sub-second; fully deterministic (xorshift stream seeded from the prime).</summary>
public sealed class CrossTripleOrthogonalityWitness : IInspectable
{
    // the wall's own first two 30-bit split primes (core_grid.gen_primes(17)[0..1]; both = 1 mod 4)
    private static readonly long[] Primes = { 1073741833L, 1073741857L };
    private const int SamplesPerPrime = 6;

    public string DisplayName =>
        "F127 cross-triple orthogonality (live: GF(p) variety slice + off-variety control)";

    public string Summary
    {
        get
        {
            int onVariety = 0, bad = 0, controls = 0, samples = 0;
            foreach (long p in Primes)
            {
                var (ev, b, c, s) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
                onVariety += ev; bad += b; controls += c; samples += s;
            }
            string verdict = bad == 0 && controls >= (int)(0.9 * samples) ? "PASS" : "FAIL";
            return $"{verdict}: {onVariety} on-variety evaluations all zero ({bad} bad), " +
                   $"{controls}/{samples} controls nonzero, {Primes.Length} split primes";
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (long p in Primes)
            {
                var (ev, bad, ctrl, samples) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
                long iRoot = CrossFormCertificate.SqrtMinusOne(p);
                yield return new InspectableNode($"p = {p}",
                    summary: $"i = {iRoot}; {ev} on-variety evaluations, {bad} nonzero (must be 0); " +
                             $"{ctrl}/{samples} off-variety controls nonzero (must be ≥ 90%)");
            }
            yield return new InspectableNode("What this is",
                summary: "an independent C# transcription of 𝔉 over GF(p) (Tonelli roots of z₃²+S·z₃+1, " +
                         "both roots each side); a SLICE of F127, not the wall (that is the 527/527 sweep)");
        }
    }
}
