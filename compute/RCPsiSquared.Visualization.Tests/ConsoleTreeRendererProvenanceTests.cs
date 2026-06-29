using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Visualization.Inspection;

namespace RCPsiSquared.Visualization.Tests;

/// <summary>The [live]/[stored] provenance badge as the ConsoleTreeRenderer prints it
/// (stranger_door's fifth door): its own bracketed token, never folded into the Summary
/// string, so the substring-content tests elsewhere stay green.</summary>
public class ConsoleTreeRendererProvenanceTests
{
    [Fact]
    public void Render_StoredNode_ShowsStoredBadge()
    {
        var node = new InspectableNode("glossary", "the house language");   // bare carrier -> Stored
        var text = ConsoleTreeRenderer.Render(node);
        Assert.Contains("[stored]", text);
    }

    [Fact]
    public void Render_LiveNode_ShowsLiveBadge()
    {
        var node = new InspectableNode("census", "60/81 pairs", provenance: NodeProvenance.Live);
        var text = ConsoleTreeRenderer.Render(node);
        Assert.Contains("[live]", text);
    }

    [Fact]
    public void Render_BadgeIsItsOwnColumn_ShowsEvenWhenSummaryEmpty()
    {
        // The badge is never concatenated into Summary: a node with no summary still shows it,
        // and the summary dash is absent (proving the two are independent columns).
        var node = new InspectableNode("bare", "", provenance: NodeProvenance.Live);
        var text = ConsoleTreeRenderer.Render(node);
        Assert.Contains("[live]", text);
        Assert.DoesNotContain("—", text);
    }
}
