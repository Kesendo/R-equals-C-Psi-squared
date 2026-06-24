using System;
using ScottPlot;

namespace RCPsiSquared.Visualization;

/// <summary>A cyberpunk + Matrix colormap: deep black through dark and neon "Matrix" green to a
/// cyan-white core. Position 0 = low (far from a branch point), 1 = high (at a branch point). Paired
/// with magenta EP markers, it reads as a Matrix rain field lit by neon at the exceptional points.</summary>
public sealed class CyberpunkMatrixColormap : IColormap
{
    public string Name => "CyberpunkMatrix";

    private static readonly (double Pos, Color Col)[] Stops =
    {
        (0.00, Color.FromHex("#000500")),   // near-black, faint green
        (0.22, Color.FromHex("#011d0c")),   // very dark green
        (0.42, Color.FromHex("#004d16")),   // dark matrix green
        (0.60, Color.FromHex("#00b32d")),   // matrix green
        (0.74, Color.FromHex("#39ff77")),   // bright matrix green
        (0.85, Color.FromHex("#39ffd0")),   // green-cyan
        (0.93, Color.FromHex("#7df9ff")),   // cyan
        (1.00, Color.FromHex("#eaffff")),   // cyan-white core
    };

    public Color GetColor(double position)
    {
        double p = Math.Clamp(position, 0.0, 1.0);
        for (int i = 1; i < Stops.Length; i++)
        {
            if (p <= Stops[i].Pos)
            {
                var (p0, c0) = Stops[i - 1];
                var (p1, c1) = Stops[i];
                double t = p1 > p0 ? (p - p0) / (p1 - p0) : 0.0;
                return new Color(
                    (byte)Math.Round(c0.R + (c1.R - c0.R) * t),
                    (byte)Math.Round(c0.G + (c1.G - c0.G) * t),
                    (byte)Math.Round(c0.B + (c1.B - c0.B) * t));
            }
        }
        return Stops[^1].Col;
    }
}
