using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Head-1 derivation instrument (multi-sector monodromy): decompose the DEFECTIVE exceptional-point
/// eigenvector of the (1,2)=(SE,DE) coherence block into the free-fermion mode-product basis |k⟩⟨k₁,k₂| and
/// read the coincidence mixture ⟨Ô⟩ = (3 − ⟨n_diff⟩)/2 FROM THE MODE GEOMETRY, to expose the closed form of
/// ⟨n_diff⟩(q).
///
/// <para><b>Reuses only built tools</b> (no reimplementation): the L(q) block via
/// <see cref="CoherenceBlock"/> + <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> (H = XYChain(N, 2q) so
/// the spectrum matches <c>WeightCoherenceBlock.Build(N,1,2,q)</c> element-for-element, γ=1); the mode-product
/// transform via <see cref="JwBlockBasis"/> (Uinv·v gives the (k,k₁,k₂) amplitudes); the sine modes ε_k, ψ_k
/// via <see cref="XyJordanWignerModes"/>.</para>
///
/// <para><b>The closed-form skeleton.</b> The AT diagonal of a pure mode-product |k⟩⟨k₁,k₂| is
/// ⟨Ô⟩_diag(k;k₁,k₂) = I(k,k₁)+I(k,k₂), the ket-mode/bra-mode density overlap I(a,b)=Σ_l ψ_a(l)²ψ_b(l)², which
/// is EXACTLY quantized: I = 1/(N+1) generic, 3/(2(N+1)) if a=b or a+b=N+1 (chiral partners). The EP
/// eigenvector is a MIXTURE within a δ-degenerate multiplet (δ = ε_k − ε_{k₁} − ε_{k₂}); the dephasing is
/// off-diagonal in the mode basis and couples equal-δ triples. So ⟨Ô⟩ = Σ_αβ c*_α c_β ⟨α|Ô|β⟩ over that small
/// reduced eigenvector: the diagonal part is Σ|c_α|²·⟨Ô⟩_diag(α), the remainder is the off-diagonal mode
/// mixing. This method reports both, plus the δ-multiplet support, to drive the analytic reduction.</para></summary>
public static class SectorBraidModeGeometry
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Build the (1,2) block L(q) (native CoherenceBlock/JwBlockBasis basis) at real q and the
    /// per-basis-index n_diff (γ=1, H=XYChain(N,2q) to match WeightCoherenceBlock).</summary>
    public static (Matrix<Complex> L, double[] Ndiff) BuildBlock(int N, double q)
    {
        var block = new CoherenceBlock(N, 1, 1.0);
        BlockBasis basis = block.Basis;
        int d = 1 << N;
        var H = PauliHamiltonian.XYChain(N, 2.0 * q).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(1.0, N).ToList();
        var flatIndices = new int[basis.MTotal];
        for (int pIdx = 0; pIdx < basis.Mp; pIdx++)
            for (int qIdx = 0; qIdx < basis.Mq; qIdx++)
            {
                long p = basis.StatesP[pIdx], b = basis.StatesQ[qIdx];
                flatIndices[basis.FlatIndex(p, b)] = (int)(p * d + b);
            }
        var L = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
        var ndiff = new double[basis.MTotal];
        for (int i = 0; i < basis.MTotal; i++)
            ndiff[i] = System.Numerics.BitOperations.PopCount((uint)((flatIndices[i] / d) ^ (flatIndices[i] % d)));
        return (L, ndiff);
    }

    /// <summary>One sample of the defective branch: coupling, the tracked eigenvalue, ⟨Ô⟩ from the AT identity
    /// (3/2 + Re λ/4) and ⟨Ô⟩ from the eigenvector contraction (v†Ôv/v†v), and ⟨n_diff⟩ = −Re λ/2.</summary>
    public readonly record struct BranchSample(double Q, Complex Lambda, double OhatFromLambda, double OhatFromVec, double NDiff);

    /// <summary>Track the defective eigenvalue branch (the one that coalesces at the EP) as a smooth function of
    /// REAL q by continuity: seed at (qStart, lambdaStart), step outward to [qMin,qMax] picking the eigenvalue
    /// nearest the running value. Records ⟨Ô⟩(q) = 3/2 + Re λ/4 (AT identity) and, as a gate, ⟨Ô⟩ from the
    /// eigenvector. This is the "function thread": does ⟨Ô⟩(q) along the branch have a nicer form than the
    /// non-radical per-locus values?</summary>
    public static List<BranchSample> SweepDefectiveBranch(int N, double qStart, Complex lambdaStart,
        double qMin, double qMax, int steps)
    {
        double dq = (qMax - qMin) / steps;
        int startIdx = (int)Math.Round((qStart - qMin) / dq);
        startIdx = Math.Clamp(startIdx, 0, steps);
        var samples = new BranchSample?[steps + 1];

        void Walk(int from, int dir, Complex seed)
        {
            Complex prev = seed;
            for (int i = from; i >= 0 && i <= steps; i += dir)
            {
                double q = qMin + i * dq;
                var (L, ndiff) = BuildBlock(N, q);
                var evd = L.Evd();
                var eigs = evd.EigenValues;
                int best = 0; double bd = double.MaxValue;
                for (int j = 0; j < eigs.Count; j++)
                {
                    double d = (eigs[j] - prev).Magnitude;
                    if (d < bd) { bd = d; best = j; }
                }
                var lam = eigs[best];
                var v = evd.EigenVectors.Column(best);
                double num = 0, den = 0;
                for (int r = 0; r < v.Count; r++)
                {
                    double w = v[r].Real * v[r].Real + v[r].Imaginary * v[r].Imaginary;
                    num += w * ndiff[r]; den += w;
                }
                double ndv = num / den;
                samples[i] = new BranchSample(q, lam, 1.5 + lam.Real / 4.0, (3.0 - ndv) / 2.0, ndv);
                prev = lam;
            }
        }

        Walk(startIdx, -1, lambdaStart);
        Walk(Math.Min(startIdx + 1, steps), +1, samples[startIdx]?.Lambda ?? lambdaStart);
        return samples.Where(s => s.HasValue).Select(s => s!.Value).OrderBy(s => s.Q).ToList();
    }

    /// <summary>I(a,b) = Σ_l ψ_a(l)²·ψ_b(l)² (a,b ∈ [1,N]); the closed-form ket/bra mode-density overlap.</summary>
    public static double ModeOverlap(XyJordanWignerModes modes, int a, int b)
    {
        double s = 0;
        for (int l = 0; l < modes.N; l++)
        {
            double pa = modes.SineMode(a, l), pb = modes.SineMode(b, l);
            s += pa * pa * pb * pb;
        }
        return s;
    }

    /// <summary>The mode-geometry reading at (N, real q, defective λ): mixture ⟨Ô⟩, the diagonal-approximation
    /// from the quantized overlaps, the off-diagonal remainder, and the dominant δ-multiplet support.</summary>
    public static string Analyze(int N, double q, Complex lambdaDefective, int topTriples = 10)
    {
        var block = new CoherenceBlock(N, 1, 1.0);
        BlockBasis basis = block.Basis;
        int d = 1 << N;

        var H = PauliHamiltonian.XYChain(N, 2.0 * q).ToMatrix();     // J = 2q matches WeightCoherenceBlock
        var gammaPerSite = Enumerable.Repeat(1.0, N).ToList();
        var flatIndices = new int[basis.MTotal];
        for (int pIdx = 0; pIdx < basis.Mp; pIdx++)
            for (int qIdx = 0; qIdx < basis.Mq; qIdx++)
            {
                long p = basis.StatesP[pIdx], b = basis.StatesQ[qIdx];
                flatIndices[basis.FlatIndex(p, b)] = (int)(p * d + b);
            }
        var L = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);

        // n_diff per flat basis index (ket = f/d, bra = f%d).
        var ndiff = new double[basis.MTotal];
        for (int i = 0; i < basis.MTotal; i++)
        {
            int f = flatIndices[i];
            ndiff[i] = System.Numerics.BitOperations.PopCount((uint)((f / d) ^ (f % d)));
        }

        var evd = L.Evd();
        var eigs = evd.EigenValues;
        var Vv = evd.EigenVectors;
        int best = 0; double bd = double.MaxValue;
        for (int i = 0; i < eigs.Count; i++)
        {
            double dd = (eigs[i] - lambdaDefective).Magnitude;
            if (dd < bd) { bd = dd; best = i; }
        }
        var v = Vv.Column(best);

        // ⟨n_diff⟩ = v† diag(ndiff) v / v†v  (the trusted AT rate, computational basis; the gate).
        double num = 0, den = 0;
        for (int i = 0; i < v.Count; i++)
        {
            double w = v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
            num += w * ndiff[i]; den += w;
        }
        double nd = num / den, overlap = (3.0 - nd) / 2.0;

        // Mode-product amplitudes c = Uinv·v (q-independent U).
        var jw = JwBlockBasis.Build(block);
        var c = jw.Uinv * v;
        var modes = jw.Modes;
        double cnorm = 0;
        for (int a = 0; a < c.Count; a++) cnorm += c[a].Real * c[a].Real + c[a].Imaginary * c[a].Imaginary;

        double Delta(int k, int k1, int k2) => modes.Dispersion[k - 1] - modes.Dispersion[k1 - 1] - modes.Dispersion[k2 - 1];
        double ODiag(int k, int k1, int k2) => ModeOverlap(modes, k, k1) + ModeOverlap(modes, k, k2);

        // diagonal-approximation ⟨Ô⟩ = Σ |c_α|²·⟨Ô⟩_diag(α); the off-diagonal remainder is the mode mixing.
        double diagApprox = 0;
        for (int a = 0; a < c.Count; a++)
        {
            double w = (c[a].Real * c[a].Real + c[a].Imaginary * c[a].Imaginary) / cnorm;
            var t = jw.Triples[a];
            diagApprox += w * ODiag(t.K, t.K1, t.K2);
        }

        var sb = new StringBuilder();
        sb.Append($"q={q.ToString("0.######", Inv)}: λ={eigs[best].Real.ToString("0.####", Inv)}");
        sb.Append($"{(eigs[best].Imaginary >= 0 ? "+" : "")}{eigs[best].Imaginary.ToString("0.####", Inv)}i, ");
        sb.Append($"⟨n_diff⟩={nd.ToString("F6", Inv)}, ⟨Ô⟩={overlap.ToString("F6", Inv)}; ");
        sb.Append($"diag-approx ⟨Ô⟩={diagApprox.ToString("F6", Inv)} (off-diag mixing {(overlap - diagApprox).ToString("F6", Inv)}). ");

        // dominant support, grouped by δ.
        var order = Enumerable.Range(0, c.Count)
            .OrderByDescending(a => c[a].Real * c[a].Real + c[a].Imaginary * c[a].Imaginary).Take(topTriples).ToList();
        sb.Append("top |c|² triples (k,k₁,k₂ | δ | weight | I-overlap): ");
        foreach (var a in order)
        {
            var t = jw.Triples[a];
            double w = (c[a].Real * c[a].Real + c[a].Imaginary * c[a].Imaginary) / cnorm;
            sb.Append($"({t.K},{t.K1},{t.K2}|δ={Delta(t.K, t.K1, t.K2).ToString("0.###", Inv)}|");
            sb.Append($"{w.ToString("F4", Inv)}|{ODiag(t.K, t.K1, t.K2).ToString("F4", Inv)}) ");
        }
        return sb.ToString();
    }
}
