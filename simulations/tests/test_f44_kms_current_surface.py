from pathlib import Path


ROOT = Path(__file__).parents[2]


def test_oq168_assigns_no_temperature_to_pure_z_dephasing():
    text = (ROOT / "review" / "OPEN_QUESTIONS_INDEX.md").read_text(encoding="utf-8")
    section = text.split("### OQ-168", 1)[1].split("### OQ-170", 1)[0]
    assert "pure Z-dephasing assigns no bath temperature" in " ".join(section.split())
    assert "pure dephasing (infinite-temperature bath)" not in section


def test_f44_typed_sources_use_only_algebraic_rate_language():
    paths = (
        ROOT / "compute" / "RCPsiSquared.Core" / "Symmetry" / "F44AlgebraicPairRateLogIdentityPi2Inheritance.cs",
        ROOT / "compute" / "RCPsiSquared.Runtime" / "PolarityArchitecture" / "F44AlgebraicPairRateLogIdentityPi2InheritanceRegistration.cs",
    )
    forbidden = ("β_eff", "EffectiveInverseTemperature", "EmpiricalJarzynskiMean",
                 "IsCrooksFluctuationTheorem", "Crooks-like")
    for path in paths:
        text = path.read_text(encoding="utf-8")
        assert "inverse-rate linear coefficient" in text
        for phrase in forbidden:
            assert phrase not in text, (path, phrase)


def test_kms_doc_scopes_alhambra_woods_petz_statement_to_its_hypotheses():
    text = (ROOT / "docs" / "KMS_DETAILED_BALANCE.md").read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    required = (
        "For a QDB dissipative semigroup with no unitary part",
        "equal to its Petz recovery map",
        "With an additional commuting unitary part",
        "reverses the unitary sign while retaining the same dissipative evolution",
    )
    forbidden = (
        "showed that quantum detailed balance is equivalent to the Petz recovery map",
        "QDB = Petz recovery map being exact channel reversal",
    )
    for phrase in required:
        assert phrase in normalized
    for phrase in forbidden:
        assert phrase not in text


def test_shifted_sls_outbound_adapter_makes_no_priority_or_ownership_claim():
    text = (ROOT / "docs" / "outbound" / "SHIFTED_ORDER4_CHIRAL_SYMMETRY.md").read_text(
        encoding="utf-8")
    normalized = " ".join(text.split())
    required = (
        "Prior-art coverage and equivalence to published constructions remain OPEN",
        "No priority or novelty claim is made",
    )
    forbidden = (
        "Priority bookkeeping",
        "does not state the operator or the spectral reflection",
        "The repository-specific object here",
    )
    for phrase in required:
        assert phrase in normalized
    for phrase in forbidden:
        assert phrase not in text
