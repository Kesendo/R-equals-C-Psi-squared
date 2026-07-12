using System.Numerics;

namespace MirrorWorld;

// The local page (adopted 2026-07-12, the k=2 doorway): the reduced state of a chosen set of sites,
// obtained by tracing out the rest -- m[a,b] = sum_r cloud[assemble(a,r), assemble(b,r)], summing over
// every configuration r of the traced-out sites (the full indices agree outside the keep-set). This is
// the one computation no whole-cloud read produces: what a watcher standing on the kept sites alone can
// see. The selection it obeys is F70 (docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md, Tier 1): a |S|-site
// page carries only |Delta popcount| <= |S| content -- one line, Tr_{~S}(|x><y|) = <x_~S|y_~S> |x_S><y_S|.
// Lifted from the two inline kernels the F72 and F73 pins already carried (MirrorWorld.Tests/SmokeTests).
// Reads are LIVE: the page is recomputed from the cloud at read time, so it follows every Step.
// Bit convention (little-endian, the world's): page bit s <-> site keep[s]; the order GIVEN in keep
// defines the page's bit order. The read delegate must return the Hermitian-correct (i,j) entry:
// Restless stores the full complex matrix (safe); Field stores only the upper triangle and mirrors the
// read UN-conjugated, which is correct only because Field's weights are real (the empty world).
public sealed class Marginal : GameObject
{
    public int N { get; }                    // the cloud's site count
    readonly int[] keep;                     // the sites the page stands on, in page-bit order
    readonly Func<int, int, Complex> read;   // the cloud entry (i,j), Hermitian-correct
    readonly int[] traced;                   // the sites traced out
    readonly int pageDim;

    public Marginal(Field cloud, int[] keep) : this((GameObject)cloud, cloud.N, (i, j) => cloud[i, j], keep) { }
    public Marginal(Restless cloud, int[] keep) : this((GameObject)cloud, cloud.N, (i, j) => cloud[i, j], keep) { }

    Marginal(GameObject cloud, int n, Func<int, int, Complex> read, int[] keep) : base(cloud)
    {
        if (keep.Length == 0 || keep.Distinct().Count() != keep.Length || keep.Any(s => s < 0 || s >= n))
            throw new ArgumentException("keep must be nonempty, distinct sites in [0, N)", nameof(keep));
        N = n;
        this.keep = (int[])keep.Clone();
        this.read = read;
        pageDim = 1 << keep.Length;
        traced = Enumerable.Range(0, n).Where(s => !keep.Contains(s)).ToArray();
    }

    public int PageDim => pageDim;

    // assemble the full index: page bits a at the kept sites, configuration r at the traced sites.
    int Assemble(int a, int r)
    {
        int u = 0;
        for (int s = 0; s < keep.Length; s++) u |= ((a >> s) & 1) << keep[s];
        for (int s = 0; s < traced.Length; s++) u |= ((r >> s) & 1) << traced[s];
        return u;
    }

    // the page entry (a, b): the partial trace, live. keep = all sites is the identity page (the
    // r-sum has exactly one term, r = 0, and the page IS the cloud).
    public Complex this[int a, int b]
    {
        get
        {
            var m = Complex.Zero;
            int rDim = 1 << traced.Length;
            for (int r = 0; r < rDim; r++)
                m += read(Assemble(a, r), Assemble(b, r));
            return m;
        }
    }

    // left: the page's own split, the same cut the cloud makes -- structure = the diagonal (real: the
    // partial trace of a Hermitian cloud has a real diagonal), novelty = the off-diagonal magnitude,
    // each upper cell standing for its mirror twin (the Field/Restless aggregate convention).
    public double Structure
    {
        get { double s = 0; for (int a = 0; a < pageDim; a++) s += this[a, a].Real; return s; }
    }

    public double Novelty
    {
        get
        {
            double s = 0;
            for (int a = 0; a < pageDim; a++)
                for (int b = a + 1; b < pageDim; b++)
                    s += 2.0 * this[a, b].Magnitude;
            return s;
        }
    }

    public override IReadOnlyList<string> Own => new[] { "page", "structure", "novelty" };
}
