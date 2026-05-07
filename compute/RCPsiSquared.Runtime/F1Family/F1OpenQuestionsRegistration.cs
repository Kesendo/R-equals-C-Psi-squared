using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F1Family;

/// <summary>Marker type used to instantiate
/// <see cref="OpenQuestionCollection{TFamilyMarker}"/> for the F1 family. Each F-family
/// declares its own empty marker so its open-question collection is a distinct Type in
/// the Runtime registry.</summary>
public sealed class F1Marker { }

/// <summary>Registers <see cref="F1OpenQuestions.Standard"/> as a single
/// <see cref="OpenQuestionCollection{TFamilyMarker}"/> at the OpenQuestion Tier. The
/// collection inherits from <see cref="F1PalindromeIdentity"/>, declaring the
/// "what is open downstream of F1" relationship.</summary>
public static class F1OpenQuestionsRegistration
{
    public static ClaimRegistryBuilder RegisterF1OpenQuestions(this ClaimRegistryBuilder builder) =>
        builder.Register<OpenQuestionCollection<F1Marker>>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();
            return new OpenQuestionCollection<F1Marker>(
                familyName: "F1",
                items: F1OpenQuestions.Standard,
                anchor: "compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs");
        });
}
