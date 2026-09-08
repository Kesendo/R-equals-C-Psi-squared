import copy
import json
from pathlib import Path

import pytest

from simulations.route_b_a2_n6_atlas import render_atlas, validate_manifest


ROOT = Path(__file__).parents[2]
MANIFEST = ROOT / "simulations/results/route_b_a2_n6_atlas.json"


def test_manifest_validation_and_render_contract(tmp_path):
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    validate_manifest(data, ROOT)
    outputs = render_atlas(MANIFEST, tmp_path)
    assert {path.suffix for path in outputs} == {".png", ".svg"}
    assert len(outputs) == 6
    assert all(path.stat().st_size > 1_000 for path in outputs)
    assert all("<svg" in path.read_text(encoding="utf-8")[:500]
               for path in outputs if path.suffix == ".svg")
    repeated = render_atlas(MANIFEST, tmp_path / "repeated")
    by_name = {path.name: path.read_bytes() for path in outputs}
    assert {path.name: path.read_bytes() for path in repeated} == by_name


def test_render_can_resolve_source_when_manifest_is_outside_the_repository(tmp_path, monkeypatch):
    external_manifest = tmp_path / "atlas.json"
    external_manifest.write_bytes(MANIFEST.read_bytes())
    monkeypatch.chdir(ROOT)
    outputs = render_atlas(external_manifest, tmp_path / "rendered")
    assert len(outputs) == 6


def test_manifest_validation_rejects_an_orbit_that_reuses_a_locus():
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    broken = copy.deepcopy(data)
    broken["orbits"][1]["memberIds"][0] = broken["orbits"][0]["memberIds"][0]
    with pytest.raises(ValueError, match="orbit partition"):
        validate_manifest(broken, ROOT)


@pytest.mark.parametrize("mutation, message", [
    ("coordinates", "coordinate mismatch"),
    ("partner", "partner mismatch"),
    ("evidence", "evidence"),
    ("orbit-closure", "orbit closure"),
])
def test_manifest_validation_rejects_local_drift_hidden_by_correct_aggregates(mutation, message):
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    broken = copy.deepcopy(data)
    if mutation == "coordinates":
        for locus in broken["loci"]:
            locus["tReal"] = 0.0
            locus["tImag"] = 0.0
    elif mutation == "partner":
        broken["loci"][0]["parityPartnerId"] = broken["loci"][1]["id"]
    elif mutation == "evidence":
        locus = next(item for item in broken["loci"]
                     if item["characterSource"] == "EpCharacterStable")
        locus["contours"][0]["isolationMargin"] = -1.0
        locus["contours"][0]["relativeDeparture"] = 1.0
    else:
        same_size = {}
        left = right = None
        for orbit in broken["orbits"]:
            size = len(orbit["memberIds"])
            if size in same_size:
                left, right = same_size[size], orbit
                break
            same_size[size] = orbit
        assert left is not None and right is not None
        left["memberIds"][0], right["memberIds"][0] = (
            right["memberIds"][0], left["memberIds"][0])
    with pytest.raises(ValueError, match=message):
        validate_manifest(broken, ROOT)


def test_manifest_validation_rejects_source_residual_and_reconciliation_mutations():
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))

    swapped = copy.deepcopy(data)
    axis = next(locus for locus in swapped["loci"]
                if locus["characterSource"] == "HermitianAxis")
    off_axis = next(locus for locus in swapped["loci"]
                    if locus["characterSource"] == "EpCharacterStable")
    axis["characterSource"] = "EpCharacterStable"
    axis["fullBlockHermiticityResidual"] = None
    for contour in axis["contours"]:
        contour["relativeDeparture"] = 0.0
    off_axis["characterSource"] = "HermitianAxis"
    off_axis["fullBlockHermiticityResidual"] = 0.0
    for contour in off_axis["contours"]:
        contour["relativeDeparture"] = None
    with pytest.raises(ValueError, match="character source does not match carrier geometry"):
        validate_manifest(swapped, ROOT)

    negative_residual = copy.deepcopy(data)
    axis = next(locus for locus in negative_residual["loci"]
                if locus["characterSource"] == "HermitianAxis")
    axis["fullBlockHermiticityResidual"] = -1.0
    with pytest.raises(ValueError, match="invalid HermitianAxis evidence"):
        validate_manifest(negative_residual, ROOT)

    incomplete = copy.deepcopy(data)
    incomplete["reconciliation"]["totalLoci"] = 0
    with pytest.raises(ValueError, match="reconciliation is incomplete"):
        validate_manifest(incomplete, ROOT)
