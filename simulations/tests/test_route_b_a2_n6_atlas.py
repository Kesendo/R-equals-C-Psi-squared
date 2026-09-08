import copy
import json
from pathlib import Path

import pytest
from matplotlib.colors import to_rgba

import simulations.route_b_a2_n6_atlas as atlas
from simulations.route_b_a2_n6_atlas import (GOLD, MAGENTA, _legend, render_atlas,
                                              validate_manifest)


ROOT = Path(__file__).parents[2]
MANIFEST = ROOT / "simulations/results/route_b_a2_n6_atlas.json"


def _artifact_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    return data.replace(b"\r\n", b"\n") if path.suffix == ".svg" else data


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
    assert all(_artifact_bytes(ROOT / "visualizations" / path.name) == _artifact_bytes(path)
               for path in outputs)


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


def test_verdict_legend_uses_marker_outlines_instead_of_nonexistent_lines():
    verdict_handles = _legend()[-2:]

    assert [handle.get_linestyle() for handle in verdict_handles] == ["None", "None"]
    assert [handle.get_marker() for handle in verdict_handles] == ["o", "o"]
    assert [handle.get_markeredgecolor() for handle in verdict_handles] == [GOLD, MAGENTA]
    assert [handle.get_label() for handle in verdict_handles] == [
        "Outline: Diabolic (266 loci)",
        "Outline: Defective (none in atlas)",
    ]


def test_constellation_facets_parity_and_binds_actual_outlines_to_legend(monkeypatch, tmp_path):
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    captured = []

    def capture_figure(fig, _output_dir, _stem):
        captured.append(fig)
        return []

    monkeypatch.setattr(atlas, "_save", capture_figure)
    atlas._constellation(data, tmp_path)
    fig = captured[0]
    assert len(fig.axes) == 2
    assert [axis.get_title() for axis in fig.axes] == [
        "R-even (133 loci)", "R-odd (133 loci)"]
    assert not fig.axes[0].get_shared_x_axes().joined(fig.axes[0], fig.axes[1])
    assert atlas.SOURCE_STYLE == {
        "HermitianAxis": ("D", 18), "EpCharacterStable": ("o", 12)}

    assert len(fig.legends) == 1
    legend = fig.legends[0]
    assert [text.get_text() for text in legend.get_texts()] == [
        "Fill: R-even (R = chain reflection)", "Fill: R-odd",
        "Shape: diamond = real-t Hermitian classifier", "Shape: circle = off-axis EP check (stable at 3 radii)",
        "Outline: Diabolic (266 loci)", "Outline: Defective (none in atlas)"]
    legend_x = legend.get_bbox_to_anchor().transformed(fig.transFigure.inverted()).x0
    assert max(axis.get_position().x1 for axis in fig.axes) < legend_x
    point_collections = [collection for axis in fig.axes for collection in axis.collections
                         if len(collection.get_offsets()) > 0]
    assert sum(len(collection.get_offsets()) for collection in point_collections) == 266
    assert all(tuple(edge) == to_rgba(GOLD) for collection in point_collections
               for edge in collection.get_edgecolors())
    assert [[len(collection.get_offsets()) for collection in axis.collections
             if len(collection.get_offsets()) > 0] for axis in fig.axes] == [[59, 74], [59, 74]]
    for axis, fill in zip(fig.axes, (atlas.CYAN, atlas.VIOLET), strict=True):
        assert all(tuple(face) == to_rgba(fill) for collection in axis.collections
                   for face in collection.get_facecolors())
    footer = " ".join(text.get_text() for text in fig.texts)
    assert "t = i qCSharp" in footer
    assert "qCSharp = J/γ" in footer
    assert "Non-real qCSharp is analytic continuation" in footer
    assert "real t is a classification axis" in footer
    assert "t ↦ -t" in footer
    assert "Dense clusters overlap" in footer
    assert "visible dots are not the census" in footer
    assert "simulations/results/route_b_a2_n6_atlas.json" in footer

    mutated = copy.deepcopy(data)
    mutated["loci"][0]["verdict"] = "Defective"
    atlas._constellation(mutated, tmp_path)
    mutated_fig = captured[1]
    mutated_edges = [tuple(edge) for axis in mutated_fig.axes
                     for collection in axis.collections
                     for edge in collection.get_edgecolors()]
    assert mutated_edges.count(to_rgba(MAGENTA)) == 1
    assert mutated_edges.count(to_rgba(GOLD)) == 265
    atlas.plt.close(fig)
    atlas.plt.close(mutated_fig)
