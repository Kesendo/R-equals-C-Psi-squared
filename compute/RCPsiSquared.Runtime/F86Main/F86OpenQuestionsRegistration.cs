using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Marker type for <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the F86
/// family.</summary>
public sealed class F86Marker { }

/// <summary>Registers <see cref="F86OpenQuestions.Standard"/> as a single
/// <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the OpenQuestion Tier. The
/// collection inherits from <see cref="ChiralAiiiClassification"/> as a representative
/// F86 anchor; any other F86 root claim would do as the "what is open downstream of F86"
/// edge target.</summary>
public static class F86OpenQuestionsRegistration
{
    public static ClaimRegistryBuilder RegisterF86OpenQuestions(this ClaimRegistryBuilder builder) =>
        builder.Register<OpenQuestionCollection<F86Marker>>(b =>
        {
            _ = b.Get<ChiralAiiiClassification>();
            return new OpenQuestionCollection<F86Marker>(
                familyName: "F86",
                items: F86OpenQuestions.Standard,
                anchor: "compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs");
        });
}
