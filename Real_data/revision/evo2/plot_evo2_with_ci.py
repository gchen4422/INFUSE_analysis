#!/usr/bin/env python3
"""Preserve the original INFUSE-only bar layout and add saved bootstrap CIs."""
from pathlib import Path
import csv
import hashlib
import json
import os

OUT = Path(os.environ.get("INFUSE_REVISION", "/path/to/revision_output")) / "INFUSE_variant_availability/evo2_ci_20261020"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REVISION = OUT.parent.parent
SOURCE = REVISION / "evo2_raw_score_external_validation" / "results"
STEM = OUT / "Evo2_external_validation_INFUSE_only_two_panel_bar_word_with_CI"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8, "axes.labelweight": "bold",
    "axes.titlesize": 9, "axes.titleweight": "bold",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .7, "pdf.fonttype": 42, "hatch.linewidth": .6,
})
fig, axes = plt.subplots(1, 2, figsize=(170 / 25.4, 78 / 25.4))
fig.subplots_adjust(left=.14, right=.975, top=.80, bottom=.18, wspace=.55)

export = []
inputs = {}
for col, (rule, cutoff, title) in enumerate([
    ("pip_threshold", .5, "PIP threshold: 0.50"),
    ("top_n", 100, "Top 100 SNPs"),
]):
    values, lower, upper = [], [], []
    for trait in ["BMI", "BP"]:
        source = SOURCE / f"raw_evo2_enrichment_by_{rule}_{trait.lower()}.tsv"
        inputs[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        with source.open() as handle:
            row, = [r for r in csv.DictReader(handle, delimiter="\t")
                    if float(r["cutoff"]) == cutoff
                    and r["method"] == "Classic_Annotations"]
        value = float(row["deleterious_shift"])
        low = float(row["deleterious_shift_ci_low"])
        high = float(row["deleterious_shift_ci_high"])
        assert low <= value <= high
        assert value == -float(row["mean_difference"])
        assert low == -float(row["difference_ci_high"])
        assert high == -float(row["difference_ci_low"])
        values.append(value)
        lower.append(low)
        upper.append(high)
        export.append(row)

    ax = axes[col]
    bars = ax.bar([0, 1], values, width=.55, color="#71B8ED",
                  edgecolor="#3779A6", linewidth=.4)
    bars[1].set_hatch("///")
    errors = [[v - lo for v, lo in zip(values, lower)],
              [hi - v for v, hi in zip(values, upper)]]
    ax.errorbar([0, 1], values, yerr=errors, fmt="none", ecolor=".25",
                elinewidth=.8, capsize=3, capthick=.8, zorder=4)
    ax.set_xlim(-.6, 1.6)
    low, high = min(0, *lower), max(0, *upper)
    span = high - low
    ax.set_ylim(low - .10 * span, high + .10 * span)
    ax.axhline(0, color=".6", linestyle="--", linewidth=.7)
    ax.set_xticks([0, 1], ["BMI", "BP"])
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-4, -4), useMathText=True)
    ax.set_title(title, pad=9)
    ax.text(-.02, 1.16, "ab"[col], transform=ax.transAxes, fontsize=11)
    if col == 0:
        ax.set_ylabel("Deleterious Evo 2 shift")
    assert ax.get_ylim()[0] < min(lower) and ax.get_ylim()[1] > max(upper)

fig.savefig(str(STEM) + ".png", dpi=600, facecolor="white")
fig.savefig(str(STEM) + ".pdf", facecolor="white")
plt.close(fig)
with Path(str(STEM) + "_data.tsv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(export[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(export)
with (OUT / "provenance.json").open("w") as handle:
    json.dump({"source_sha256": inputs, "figure_size_mm": [170, 78],
               "png_dpi": 600, "matplotlib_version": matplotlib.__version__,
               "interval": "95% locus-bootstrap percentile CI; 1000 replicates",
               "source": "raw bootstrap summary tables"}, handle, indent=2)
print(str(STEM) + ".pdf")
print(str(STEM) + ".png")
for row in export:
    print(row["Trait"], row["selection_rule"], row["cutoff"],
          row["deleterious_shift"], row["deleterious_shift_ci_low"],
          row["deleterious_shift_ci_high"])
