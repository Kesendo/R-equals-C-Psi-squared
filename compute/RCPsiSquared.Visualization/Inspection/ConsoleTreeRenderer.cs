using System.Text;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Visualization.Inspection;

/// <summary>Walks an <see cref="IInspectable"/> tree and emits an indented text view —
/// the simplest possible Object Manager. Each line is "<c>display name — summary</c>",
/// with payload type appended for leaves.
///
/// <para>Use <see cref="Render"/> for a string, or <see cref="WriteTo"/> to stream to a
/// <see cref="TextWriter"/> (Console, file, xUnit's <c>ITestOutputHelper</c>).</para>
/// </summary>
public static class ConsoleTreeRenderer
{
    /// <summary>Render the tree rooted at <paramref name="root"/> as a multi-line string.</summary>
    public static string Render(IInspectable root, int maxDepth = 4, bool drawPayloads = false, int maxWidth = 64)
    {
        var sb = new StringBuilder();
        WriteNode(root, sb, prefix: "", isLast: true, depth: 0, maxDepth, drawPayloads, maxWidth);
        return sb.ToString();
    }

    /// <summary>Stream the tree to <paramref name="writer"/>.</summary>
    public static void WriteTo(IInspectable root, TextWriter writer, int maxDepth = 4,
        bool drawPayloads = false, int maxWidth = 64)
    {
        writer.Write(Render(root, maxDepth, drawPayloads, maxWidth));
    }

    private static void WriteNode(IInspectable node, StringBuilder sb, string prefix, bool isLast,
        int depth, int maxDepth, bool drawPayloads, int maxWidth)
    {
        string connector = depth == 0 ? "" : (isLast ? "└── " : "├── ");
        string label = node.DisplayName;
        string summary = string.IsNullOrEmpty(node.Summary) ? "" : $"  —  {node.Summary}";
        string payloadTag = PayloadTag(node.Payload);

        sb.Append(prefix);
        sb.Append(connector);
        sb.Append(label);
        sb.Append(summary);
        if (!string.IsNullOrEmpty(payloadTag))
        {
            sb.Append("  [");
            sb.Append(payloadTag);
            sb.Append(']');
        }
        sb.AppendLine();

        // The payload's picture, hung under the node and indented into its child region. Drawn for
        // leaves too (before the maxDepth cut), so a vector/matrix/curve at the depth limit still shows.
        if (drawPayloads)
        {
            var plot = AsciiPayloadPlot.Render(node.Payload, maxWidth);
            if (plot.Count > 0)
            {
                string plotPrefix = prefix + (depth == 0 ? "" : isLast ? "    " : "│   ") + "    ";
                foreach (var line in plot)
                {
                    sb.Append(plotPrefix);
                    sb.Append(line);
                    sb.AppendLine();
                }
            }
        }

        if (depth >= maxDepth) return;

        var childList = node.Children.ToList();
        for (int i = 0; i < childList.Count; i++)
        {
            string childPrefix = depth == 0 ? "" : prefix + (isLast ? "    " : "│   ");
            WriteNode(childList[i], sb, childPrefix, isLast: i == childList.Count - 1,
                depth: depth + 1, maxDepth, drawPayloads, maxWidth);
        }
    }

    private static string PayloadTag(InspectablePayload payload) => payload switch
    {
        InspectablePayload.None => "",
        InspectablePayload.Real r => $"real {r.Value:G4}",
        InspectablePayload.Scalar s => $"complex |{s.Value.Magnitude:G4}|",
        InspectablePayload.Vector v => $"vector[{v.Values.Count}]",
        InspectablePayload.MatrixView m => $"matrix[{m.Values.RowCount}×{m.Values.ColumnCount}]",
        InspectablePayload.Curve c => $"curve[{c.X.Count}]",
        _ => "",
    };
}
