from pathlib import Path


ROOT = Path(__file__).parents[2]


def test_legacy_delta_locator_carries_no_positive_delta_character_verdict():
    source = (ROOT / "simulations" / "f89_zz_break_gate.py").read_text(encoding="utf-8")
    normalized = " ".join(source.split())

    for required in (
        "historical locator without a positive-Delta character verdict",
        "strict full-block coincidence/correspondence tolerance",
        "positive-Delta proposals are Uncertified",
        "neither persistence, defectiveness, nor lifting",
    ):
        assert required in normalized

    for forbidden in (
        "Delta>0 off-axis/defective",
        "BECOMES DEFECTIVE: breaking free-fermion",
        "FREE-FERMION INTEGRABILITY WAS THE PROTECTION",
        "the EP has MOVED OFF the real axis",
        "genuine defective Jordan EP",
    ):
        assert forbidden not in source


def test_current_hypothesis_calls_the_script_a_locator_not_a_character_gate():
    text = (ROOT / "hypotheses" / "DIABOLIC_BY_INTEGRABILITY.md").read_text(encoding="utf-8")
    assert "historical locator without a positive-Delta character verdict" in " ".join(text.split())
    assert "gate-first; Stage 0 reproduces" not in text
