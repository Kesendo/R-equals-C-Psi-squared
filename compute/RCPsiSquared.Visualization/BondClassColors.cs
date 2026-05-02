using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Visualization;

public static class BondClassColors
{
    public static readonly ScottPlot.Color Interior = ScottPlot.Colors.RoyalBlue;
    public static readonly ScottPlot.Color Endpoint = ScottPlot.Colors.Crimson;

    public static ScottPlot.Color Of(BondClass cls) => cls switch
    {
        BondClass.Interior => Interior,
        BondClass.Endpoint => Endpoint,
        _ => ScottPlot.Colors.DarkSlateGray,
    };
}
