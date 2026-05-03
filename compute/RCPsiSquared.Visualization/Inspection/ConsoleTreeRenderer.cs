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
    public static string Render(IInspectable root, int maxDepth = 4)
    {
        var sb = new StringBuilder();
        WriteNode(root, sb, prefix: "", isLast: true, depth: 0, maxDepth);
        return sb.ToString();
    }

    /// <summary>Stream the tree to <paramref name="writer"/>.</summary>
    public static void WriteTo(IInspectable root, TextWriter writer, int maxDepth = 4)
    {
        writer.Write(Render(root, maxDepth));
    }

    private static void WriteNode(IInspectable node, StringBuilder sb, string prefix, bool isLast,
        int depth, int maxDepth)
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

        if (depth >= maxDepth) return;

        var childList = node.Children.ToList();
        for (int i = 0; i < childList.Count; i++)
        {
            string childPrefix = depth == 0 ? "" : prefix + (isLast ? "    " : "│   ");
            WriteNode(childList[i], sb, childPrefix, isLast: i == childList.Count - 1,
                depth: depth + 1, maxDepth);
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
