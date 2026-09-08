"""Render the N=6 Route-B A2 discovery atlas from its gated manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "route-b-a2-n6-atlas-v1"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SOURCE_SHA256 = "c4b19015a750e56a55b0c8fdba3d05534c0c133bc5355a74eb89b77020532d47"
CYAN = "#28d7e5"
VIOLET = "#9b7bea"
GOLD = "#e2b83b"
MAGENTA = "#ed4ba4"
INK = "#e8edf6"
MUTED = "#8f9aad"
GRID = "#303846"
BACKGROUND = "#11151c"
VERDICT_OUTLINE = {"Diabolic": GOLD, "Defective": MAGENTA}
SOURCE_STYLE = {"HermitianAxis": ("D", 18), "EpCharacterStable": ("o", 12)}


def _rational(value: dict) -> Fraction:
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def _box_midpoint(box: dict, component: str) -> Fraction:
    interval = box[component]
    return (_rational(interval["lower"]) + _rational(interval["upper"])) / 2


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _display_double(value: Fraction) -> float:
    # Match ExactRational.ToDouble: convert numerator and denominator separately.
    return float(value.numerator) / float(value.denominator)


def validate_manifest(data: dict, repo_root: Path) -> None:
    """Reject any drift between the manifest, its source and the approved atlas semantics."""
    require = lambda condition, message: condition or (_raise(message))
    require(data.get("schemaVersion") == 1 and data.get("n") == 6,
            "unexpected atlas schema or N")
    require(data.get("sourceSha256") == SOURCE_SHA256, "unexpected source digest")
    source = repo_root / data.get("sourceArtifact", "")
    require(source.is_file(), "source artifact is missing")
    require(hashlib.sha256(source.read_bytes()).hexdigest() == SOURCE_SHA256,
            "source artifact digest mismatch")
    carrier = json.loads(source.read_text(encoding="utf-8"))
    require(data.get("coordinateMeaning") ==
            "exact t-box midpoint converted to display double; not an algebraic root",
            "coordinate meaning drift")

    loci = data.get("loci", [])
    orbits = data.get("orbits", [])
    require(len(loci) == 266 and len({locus["id"] for locus in loci}) == 266,
            "expected 266 unique loci")
    source_by_id = {locus["id"]: locus for locus in carrier["loci"]}
    require(set(source_by_id) == {locus["id"] for locus in loci}, "locus ID mismatch")
    for locus in loci:
        source_locus = source_by_id[locus["id"]]
        re_mid = _box_midpoint(source_locus["tBox"], "real")
        im_mid = _box_midpoint(source_locus["tBox"], "imag")
        require(locus["tRealExactMidpoint"] == _fraction_text(re_mid) and
                locus["tImagExactMidpoint"] == _fraction_text(im_mid) and
                math.isclose(locus["tReal"], _display_double(re_mid), rel_tol=0, abs_tol=1e-15) and
                math.isclose(locus["tImag"], _display_double(im_mid), rel_tol=0, abs_tol=1e-15),
                f"coordinate mismatch at {locus['id']}")
        require(locus["parity"] == source_locus["parity"] and
                locus["algebraicMultiplicity"] == source_locus["algebraicMultiplicity"] and
                locus["conjugationPartnerId"] == source_locus["conjugationPartnerId"] and
                locus["parityPartnerId"] == source_locus["parityPartnerId"],
                f"partner mismatch at {locus['id']}")
        contours = locus["contours"]
        require(len(contours) == 3 and
                all(math.isfinite(c["radius"]) and c["radius"] > 0 and
                    math.isfinite(c["isolationMargin"]) and c["isolationMargin"] > 0
                    for c in contours) and
                all(contours[index]["radius"] < contours[index + 1]["radius"]
                    for index in range(2)) and
                all(contours[index]["isolationMargin"] > contours[index + 1]["isolationMargin"]
                    for index in range(2)),
                f"invalid contour evidence at {locus['id']}")
        if locus["characterSource"] == "HermitianAxis":
            require(im_mid == 0,
                    f"character source does not match carrier geometry at {locus['id']}")
            residual = locus["fullBlockHermiticityResidual"]
            require(math.isfinite(residual) and 0 <= residual < 1e-12 and
                    all(c["relativeDeparture"] is None for c in contours),
                    f"invalid HermitianAxis evidence at {locus['id']}")
        elif locus["characterSource"] == "EpCharacterStable":
            require(im_mid != 0,
                    f"character source does not match carrier geometry at {locus['id']}")
            require(locus["fullBlockHermiticityResidual"] is None and
                    all(c["relativeDeparture"] is not None and
                        math.isfinite(c["relativeDeparture"]) and
                        0 <= c["relativeDeparture"] < 1e-6 for c in contours),
                    f"invalid EpCharacter evidence at {locus['id']}")
        else:
            _raise(f"invalid evidence source at {locus['id']}")
    require(sum(locus["parity"] == "E" for locus in loci) == 133 and
            sum(locus["parity"] == "O" for locus in loci) == 133,
            "parity census drift")
    require(sum(locus["characterSource"] == "HermitianAxis" for locus in loci) == 118 and
            sum(locus["characterSource"] == "EpCharacterStable" for locus in loci) == 148,
            "character-source census drift")
    require(all(locus["verdict"] == "Diabolic" and
                locus["algebraicMultiplicity"] == locus["geometricMultiplicity"] == 2 and
                len(locus["contours"]) == 3 for locus in loci),
            "character or contour drift")
    member_ids = [member for orbit in orbits for member in orbit["memberIds"]]
    require(len(orbits) == 96 and
            sum(len(orbit["memberIds"]) == 2 for orbit in orbits) == 59 and
            sum(len(orbit["memberIds"]) == 4 for orbit in orbits) == 37,
            "orbit census drift")
    require(sorted(member_ids) == sorted(locus["id"] for locus in loci),
            "orbit partition does not cover every locus exactly once")
    manifest_by_id = {locus["id"]: locus for locus in loci}
    for orbit in orbits:
        members = set(orbit["memberIds"])
        require(all(manifest_by_id[member]["conjugationPartnerId"] in members and
                    manifest_by_id[member]["parityPartnerId"] in members
                    for member in members),
                f"orbit closure fails at {orbit['id']}")
    reconciliation = data.get("reconciliation", {})
    require(reconciliation.get("totalLoci") == 266 and
            reconciliation.get("consumedLoci") == 266 and
            reconciliation.get("unresolvedLoci") == 0,
            "reconciliation is incomplete")
    require(reconciliation.get("byAlgebraicSource") == {"ExactAlgebraic": 266} and
            reconciliation.get("byCharacterSource") == {
                "HermitianAxis": 118, "EpCharacterStable": 148, "ExactRankExecuted": 0} and
            reconciliation.get("byVerdict") == {"Diabolic": 266},
            "reconciliation census mismatch")
    middle_margins = [locus["contours"][1]["isolationMargin"] for locus in loci
                      if locus["characterSource"] == "EpCharacterStable"]
    require(math.isclose(reconciliation.get("minimumIsolationMargin", math.nan),
                         min(middle_margins), rel_tol=0, abs_tol=1e-15),
            "reconciliation margin mismatch")


def _raise(message: str) -> bool:
    raise ValueError(message)


def _style_axis(ax, xlabel: str, ylabel: str) -> None:
    ax.set_facecolor(BACKGROUND)
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.grid(color=GRID, linewidth=0.55, alpha=0.55)


def _save(fig, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.svg"]
    metadata = {"Creator": "RCPsiSquared route_b_a2_n6_atlas.py", "Date": "2026-09-08"}
    fig.savefig(paths[0], dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor(),
                metadata=metadata)
    fig.savefig(paths[1], bbox_inches="tight", facecolor=fig.get_facecolor(), metadata=metadata)
    svg = paths[1].read_text(encoding="utf-8")
    normalized_svg = "\n".join(line.rstrip() for line in svg.splitlines()) + "\n"
    paths[1].write_text(normalized_svg, encoding="utf-8", newline="\n")
    plt.close(fig)
    return paths


def _legend() -> list[Line2D]:
    return [
        Line2D([], [], marker="s", linestyle="none", markerfacecolor=CYAN,
               markeredgecolor=CYAN, label="Fill: R-even (R = chain reflection)"),
        Line2D([], [], marker="s", linestyle="none", markerfacecolor=VIOLET,
               markeredgecolor=VIOLET, label="Fill: R-odd"),
        Line2D([], [], marker="D", linestyle="none", markerfacecolor=INK,
               markeredgecolor=INK, label="Shape: diamond = real-t Hermitian classifier"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=INK,
               markeredgecolor=INK, label="Shape: circle = off-axis EP check (stable at 3 radii)"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor=VERDICT_OUTLINE["Diabolic"], markeredgewidth=1.8,
               label="Outline: Diabolic (266 loci)"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor=VERDICT_OUTLINE["Defective"], markeredgewidth=1.8,
               label="Outline: Defective (none in atlas)"),
    ]


def _constellation(data: dict, output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 7.2), facecolor=BACKGROUND,
                             sharey=True)
    for ax, parity, color, title in (
            (axes[0], "E", CYAN, "R-even (133 loci)"),
            (axes[1], "O", VIOLET, "R-odd (133 loci)")):
        for source, (marker, size) in SOURCE_STYLE.items():
            points = [locus for locus in data["loci"]
                      if locus["characterSource"] == source and locus["parity"] == parity]
            ax.scatter([p["tReal"] for p in points], [p["tImag"] for p in points],
                       s=size, marker=marker, facecolors=color,
                       edgecolors=[VERDICT_OUTLINE[p["verdict"]] for p in points],
                       linewidths=0.75, alpha=1.0, zorder=3)
        ax.axhline(0, color=MUTED, linewidth=0.8, alpha=0.8)
        ax.axvline(0, color=MUTED, linewidth=0.8, alpha=0.8)
        _style_axis(ax, r"Re $t$  (box-midpoint display coordinate)",
                    r"Im $t$  (box-midpoint display coordinate)")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, color=color, pad=10)
    fig.suptitle(r"N=6 Route-B $A_2$ midpoint overview: 266 certified boxes",
                 color=INK, y=0.94)
    fig.subplots_adjust(left=0.07, right=0.77, bottom=0.18, top=0.88, wspace=0.08)
    fig.legend(handles=_legend(), loc="upper left", bbox_to_anchor=(0.79, 0.88),
               frameon=False, labelcolor=INK, fontsize=9, ncol=1)
    fig.text(0.5, 0.078,
             "Displayed points are exact rational t-box midpoints converted to float; they are not algebraic roots.",
             color=MUTED, ha="center", fontsize=8)
    fig.text(0.5, 0.046,
             "t = i qCSharp, qCSharp = J/γ. Non-real qCSharp is analytic continuation; real t is a classification axis.",
             color=MUTED, ha="center", fontsize=8)
    fig.text(0.5, 0.014,
             "Panels obey t ↦ -t. Dense clusters overlap; visible dots are not the census. Exact table: simulations/results/route_b_a2_n6_atlas.json",
             color=MUTED, ha="center", fontsize=8)
    return _save(fig, output_dir, "route_b_a2_n6_constellation")


def _orbit_map(data: dict, output_dir: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(9.4, 6.8), facecolor=BACKGROUND)
    by_id = {locus["id"]: locus for locus in data["loci"]}
    drawn_conjugation: set[tuple[str, str]] = set()
    drawn_parity: set[tuple[str, str]] = set()
    for locus in data["loci"]:
        for field, drawn, color, style in (
                ("conjugationPartnerId", drawn_conjugation, MUTED, "-"),
                ("parityPartnerId", drawn_parity, GOLD, (0, (2, 2)))):
            edge = tuple(sorted((locus["id"], locus[field])))
            if edge in drawn or edge[0] == edge[1]:
                continue
            drawn.add(edge)
            partner = by_id[locus[field]]
            ax.plot([locus["tReal"], partner["tReal"]],
                    [locus["tImag"], partner["tImag"]], color=color,
                    linestyle=style, linewidth=0.55, alpha=0.32, zorder=1)
    for parity, color in (("E", CYAN), ("O", VIOLET)):
        points = [locus for locus in data["loci"] if locus["parity"] == parity]
        ax.scatter([p["tReal"] for p in points], [p["tImag"] for p in points],
                   s=18, c=color, edgecolors=GOLD, linewidths=0.55, zorder=3)
    ax.axhline(0, color=MUTED, linewidth=0.8, alpha=0.8)
    ax.axvline(0, color=MUTED, linewidth=0.8, alpha=0.8)
    _style_axis(ax, r"Re $t$", r"Im $t$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Partner-orbit map — 59 axis pairs + 37 off-axis quartets",
                 color=INK, pad=14)
    handles = [
        Line2D([], [], color=MUTED, linewidth=1.2, label="complex conjugation"),
        Line2D([], [], color=GOLD, linewidth=1.2, linestyle=(0, (2, 2)),
               label=r"parity transport $t\mapsto-t$"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=CYAN,
               markeredgecolor=GOLD, label="R-even"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=VIOLET,
               markeredgecolor=GOLD, label="R-odd"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False,
              labelcolor=INK, fontsize=8, ncol=2)
    fig.text(0.5, 0.015,
             "Only the two partner maps present in the certified carrier are drawn; no crossfold edge is inferred.",
             color=MUTED, ha="center", fontsize=8)
    return _save(fig, output_dir, "route_b_a2_n6_orbit_map")


def _evidence_profile(data: dict, output_dir: Path) -> list[Path]:
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.2, 4.8), facecolor=BACKGROUND)
    loci = data["loci"]
    for source, color, label in (("HermitianAxis", CYAN, "HermitianAxis (118)"),
                                  ("EpCharacterStable", VIOLET, "EpCharacterStable (148)")):
        values = sorted(min(contour["isolationMargin"] for contour in locus["contours"])
                        for locus in loci if locus["characterSource"] == source)
        left.scatter(range(1, len(values) + 1), values, s=12, c=color,
                     edgecolors=GOLD, linewidths=0.35, label=label)
    _style_axis(left, "rank within evidence source", "minimum isolation margin over 3 contours")
    left.set_title("Isolation from the third eigenvalue", color=INK)
    left.legend(frameon=False, labelcolor=INK, fontsize=8)

    departures = sorted(max(contour["relativeDeparture"] for contour in locus["contours"]
                            if contour["relativeDeparture"] is not None)
                        for locus in loci if locus["characterSource"] == "EpCharacterStable")
    floor = 1e-18
    exact_zero_count = sum(value == 0 for value in departures)
    right.scatter(range(1, len(departures) + 1), [max(value, floor) for value in departures],
                  s=13, c=VIOLET, edgecolors=GOLD, linewidths=0.35)
    right.axhline(1e-6, color=MAGENTA, linewidth=1.1, linestyle=(0, (4, 3)),
                  label=r"Diabolic gate: departure $<10^{-6}$")
    right.set_yscale("log")
    _style_axis(right, "rank among EpCharacterStable loci",
                "maximum relative departure over 3 contours")
    right.set_title("Three-contour semisimplicity stability", color=INK)
    right.legend(frameon=False, labelcolor=INK, fontsize=8)
    right.text(0.03, 0.04,
               f"{exact_zero_count} exact zeros are shown at the display floor {floor:.0e}",
               transform=right.transAxes, color=MUTED, fontsize=7.5)
    fig.suptitle("Executed character evidence — 266 diabolic, 0 defective",
                 color=INK, y=1.02)
    fig.tight_layout()
    fig.text(0.5, -0.015,
             "Hermitian-axis character is structural (full 45-dimensional block residual = 0); the log display floor is not a measured departure.",
             color=MUTED, ha="center", fontsize=8)
    return _save(fig, output_dir, "route_b_a2_n6_evidence_profile")


def _find_repo_root(manifest_path: Path) -> Path:
    cwd = Path.cwd().resolve()
    candidates = list(manifest_path.parents) + [cwd] + list(cwd.parents)
    for candidate in candidates:
        if (candidate / "simulations/results/route_b_a2_n6.json").is_file():
            return candidate
    raise FileNotFoundError(
        "Cannot locate simulations/results/route_b_a2_n6.json from the manifest path or working directory.")


def render_atlas(manifest_path: Path, output_dir: Path) -> list[Path]:
    manifest_path = manifest_path.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    repo_root = _find_repo_root(manifest_path)
    validate_manifest(data, repo_root)
    outputs: list[Path] = []
    outputs.extend(_constellation(data, output_dir))
    outputs.extend(_orbit_map(data, output_dir))
    outputs.extend(_evidence_profile(data, output_dir))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path,
                        default=Path("simulations/results/route_b_a2_n6_atlas.json"))
    parser.add_argument("--out-dir", type=Path, default=Path("visualizations"))
    args = parser.parse_args()
    for output in render_atlas(args.manifest, args.out_dir):
        print(output)


if __name__ == "__main__":
    main()
