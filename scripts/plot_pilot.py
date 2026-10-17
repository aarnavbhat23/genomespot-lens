"""Render the archived pilot; requires matplotlib==3.7.4."""

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from genomespot_lens.core import quality_flags

root = Path(__file__).resolve().parents[1]
cohorts = [
    ("ATB_SAMD00000344_longest20_stress_test", "AllTheBacteria: SAMD00000344", "Longest 20% of proteins"),
    ("M_jannaschii_current_reviewed_exact_matches", "Methanocaldococcus jannaschii", "Current reviewed proteins: exact matches"),
]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
rng = np.random.default_rng(23)
for ax, (folder, title, subtitle) in zip(axes, cohorts):
    records = json.load(gzip.open(root / "examples" / folder / "predictions.json.gz", "rt"))
    target = "temperature_optimum"
    baseline = records[0]["predictions"][target]["value"]
    observed = records[1]["predictions"][target]["value"] - baseline
    for x, kind in enumerate(["count_matched", "length_matched"]):
        controls = [r["predictions"][target] for r in records if r["kind"] == kind]
        pairs = [(p["value"] - baseline, bool(quality_flags(target, p))) for p in controls if p["value"] is not None]
        values = [v for v, _ in pairs]
        ax.boxplot([values], positions=[x], widths=.38, whis=(2.5, 97.5), showfliers=False, medianprops={"color": "#17324d", "linewidth": 2}, boxprops={"color": "#536779"}, whiskerprops={"color": "#536779"}, capprops={"color": "#536779"})
        for flagged, color in [(False, "#21867a"), (True, "#cc7a25")]:
            points = [v for v, flag in pairs if flag == flagged]
            ax.scatter(x + rng.uniform(-.16, .16, len(points)), points, s=12, alpha=.60, color=color, zorder=3)
    ax.axhline(0, color="#9ba8b4", linewidth=.8)
    ax.axhline(observed, color="#a52655", linewidth=1.7, linestyle="--", label="Selected subset")
    ax.set_title(title + "\n" + subtitle + f"\nFull: {baseline:.2f} °C  |  Selected change: {observed:+.2f} °C", loc="left", fontsize=10, pad=12)
    ax.set_xticks([0, 1], ["Count-matched\n(n = 100)", "Length-matched\n(n = 100)"])
    ax.set_xlim(-.5, 1.5)
    ax.grid(axis="y", alpha=.15)
axes[0].set_ylabel("Change from full-input temperature prediction (°C)")
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([], [], color="#a52655", linestyle="--", label="Selected subset"), Line2D([], [], color="#21867a", marker="o", linestyle="", label="Unflagged control"), Line2D([], [], color="#cc7a25", marker="o", linestyle="", label="Flagged control")], loc="lower center", ncol=3, bbox_to_anchor=(.5, .075), frameon=False)
fig.suptitle("Protein selection changes predictions — controls determine the context", x=.07, ha="left", fontsize=14, fontweight="bold")
fig.text(.07, .035, "Two-genome exploratory pilot • 100-aa length bins • Box: quartiles; whiskers: 2.5–97.5%\nModel sensitivity only: full-input predictions are not experimental ground truth.", color="#536779", fontsize=9)
fig.subplots_adjust(left=.075, right=.98, top=.79, bottom=.25, wspace=.17)
fig.savefig(root / "docs/pilot-temperature.png", dpi=180)
fig.savefig(root / "docs/pilot-temperature.svg")
svg = root / "docs/pilot-temperature.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
