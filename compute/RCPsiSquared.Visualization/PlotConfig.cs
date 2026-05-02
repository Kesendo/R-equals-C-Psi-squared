namespace RCPsiSquared.Visualization;

public sealed record PlotConfig(int Width = 1024, int Height = 640)
{
    public static readonly PlotConfig Default = new();
}
