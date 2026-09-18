"""
plot_results.py
Generates summary result figures for the NA12878 chr20 variant-calling project.
Run from the project root: python scripts/plot_results.py
Outputs PNG files into results/
"""

import matplotlib.pyplot as plt
import os

OUTDIR = "results"
os.makedirs(OUTDIR, exist_ok=True)

# --- Real numbers from this project's pipeline output ---
# (Sourced from: high_impact_variants.txt counts, SNP/indel filtered VCFs,
#  and benchmark_results.summary.csv — update if you re-run the pipeline.)

impact_labels = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
impact_counts = [134206, 604, 325, 41]

snp_indel_labels = ["SNP", "Indel"]
total_counts = [113174, 22002]
pass_counts = [106536, 21646]

metrics = ["Recall", "Precision", "F1"]
snp_pass = [0.992801, 0.991353, 0.992076]
indel_pass = [0.989673, 0.993195, 0.991431]


def plot_impact_distribution():
    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(impact_labels, impact_counts, color="#4C72B0")
    ax.set_yscale("log")
    ax.set_ylabel("Number of variants (log scale)")
    ax.set_title("Variant Impact Distribution — Chromosome 20 (NA12878)")
    for bar, count in zip(bars, impact_counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.1,
                 f"{count:,}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/impact_distribution.png", dpi=150)
    plt.close(fig)


def plot_total_vs_pass():
    x = range(len(snp_indel_labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar([i - width / 2 for i in x], total_counts, width, label="Total calls", color="#DD8452")
    ax.bar([i + width / 2 for i in x], pass_counts, width, label="PASS only", color="#4C72B0")
    ax.set_xticks(list(x))
    ax.set_xticklabels(snp_indel_labels)
    ax.set_ylabel("Number of variants")
    ax.set_title("Total vs. PASS-Filtered Variant Counts")
    ax.legend()
    for i in x:
        ax.text(i - width / 2, total_counts[i] + 1000, f"{total_counts[i]:,}", ha="center", fontsize=9)
        ax.text(i + width / 2, pass_counts[i] + 1000, f"{pass_counts[i]:,}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/total_vs_pass.png", dpi=150)
    plt.close(fig)


def plot_benchmark_scores():
    x = range(len(metrics))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar([i - width / 2 for i in x], snp_pass, width, label="SNP (PASS)", color="#4C72B0")
    ax.bar([i + width / 2 for i in x], indel_pass, width, label="Indel (PASS)", color="#55A868")
    ax.set_xticks(list(x))
    ax.set_xticklabels(metrics)
    ax.set_ylim(0.95, 1.0)
    ax.set_ylabel("Score")
    ax.set_title("Benchmark Accuracy vs. GIAB Truth Set (Confident Regions)")
    ax.legend()
    for i in x:
        ax.text(i - width / 2, snp_pass[i] + 0.001, f"{snp_pass[i]:.3f}", ha="center", fontsize=9)
        ax.text(i + width / 2, indel_pass[i] + 0.001, f"{indel_pass[i]:.3f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/benchmark_scores.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot_impact_distribution()
    plot_total_vs_pass()
    plot_benchmark_scores()
    print(f"Saved 3 charts to {OUTDIR}/")
