using System;
using System.Collections.Generic;
using RCPsiSquared.Core.Numerics;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>Spin 3/2 as two qubits in Souza's labelling, |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩,
/// |10⟩, |11⟩ (the first letter on the high bit): I_z = ZI + ½IZ, √6·T₂₀ = 3·ZZ, the jumps ZZ at
/// 3J₀/2, ZX and ZY at 3J₁/2, XI and YI at 3J₂/2, each γ(PρP − ρ), H = q·ZZ + ω·(ZI + ½IZ)
/// − 3K₂·ZI − 3K₁·IZ. Data for string instruments, and a second construction of the generator.</summary>
public static class Spin32Strings
{
    public static IReadOnlyList<(string Label, BigRational Rate)> Jumps(QuadrupolarRates r) => new[]
    {
        ("ZZ", 3 * r.J0 / 2), ("ZX", 3 * r.J1 / 2), ("ZY", 3 * r.J1 / 2), ("XI", 3 * r.J2 / 2), ("YI", 3 * r.J2 / 2),
    };

    public static IReadOnlyList<(string Label, BigRational Coefficient)> Hamiltonian(SecularHamiltonian h)
    {
        if (!h.Nu.IsZero || h.ExtraPhysical is not null)
            throw new NotSupportedException("the lock's I_x = (√3/2)·IX + ½(XX + YY) is not rational in strings.");
        return new[] { ("ZZ", h.Q), ("ZI", h.Omega - 3 * h.K2), ("IZ", h.Omega / 2 - 3 * h.K1) };
    }

    /// <summary>The 4 × 4 matrix of a two-letter Pauli string, first letter on the high bit.</summary>
    public static GaussianRational[,] Matrix(string label)
    {
        if (label is null || label.Length != 2) throw new ArgumentException("two letters.", nameof(label));
        var a = Letter(label[0]);
        var b = Letter(label[1]);
        var r = M.Zero(4);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++) r[2 * i + k, 2 * j + l] = a[i, j] * b[k, l];
        return r;
    }

    public static GaussianRational[,] Generator(QuadrupolarRates rates, SecularHamiltonian h)
    {
        var hm = M.Zero(4);
        foreach (var (label, c) in Hamiltonian(h)) hm = M.Add(hm, M.Scale(Matrix(label), c));
        var l = QuadrupolarGenerator.CommutatorPart(hm);
        var one = M.Identity(16);
        foreach (var (label, rate) in Jumps(rates))
        {
            if (rate.IsZero) continue;
            var p = Matrix(label);
            l = M.Add(l, M.Scale(M.Subtract(M.Sandwich(p, p), one), rate));
        }
        return l;
    }

    private static GaussianRational[,] Letter(char c)
    {
        var m = M.Zero(2);
        switch (c)
        {
            case 'I': m[0, 0] = 1; m[1, 1] = 1; break;
            case 'X': m[0, 1] = 1; m[1, 0] = 1; break;
            case 'Y': m[0, 1] = -GaussianRational.I; m[1, 0] = GaussianRational.I; break;
            case 'Z': m[0, 0] = 1; m[1, 1] = -1; break;
            default: throw new ArgumentException($"not a Pauli letter: {c}");
        }
        return m;
    }
}
