import argparse
import os
from typing import List

import matplotlib

# Force a non-GUI backend to avoid Qt conflicts in headless or mixed-Qt environments
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

STAGE_ORDER: List[str] = ["Early", "Mid", "Late"]
PROPERTY_COLUMNS = ["MW_R2", "Ring_R2", "Atoms_R2", "LogP_R2"]


def load_data(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "Stage" not in df.columns:
        raise ValueError("CSV must include a Stage column")
    missing = [c for c in PROPERTY_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"CSV is missing columns: {', '.join(missing)}")
    df = df[df["Stage"].isin(STAGE_ORDER)].copy()
    df["Stage"] = pd.Categorical(df["Stage"], categories=STAGE_ORDER, ordered=True)
    return df.sort_values("Stage")


def plot_line(df: pd.DataFrame, out_path: str) -> None:
    melted = df.melt(id_vars=["Stage"], value_vars=PROPERTY_COLUMNS, var_name="Property", value_name="R2")
    plt.figure(figsize=(6, 4))
    sns.lineplot(data=melted, x="Stage", y="R2", hue="Property", marker="o")
    plt.ylim(0.0, 1.05)
    plt.ylabel("R² score")
    plt.title("R² vs Decoder Stage")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_grouped_bar(df: pd.DataFrame, out_path: str) -> None:
    melted = df.melt(id_vars=["Stage"], value_vars=PROPERTY_COLUMNS, var_name="Property", value_name="R2")
    plt.figure(figsize=(6, 4))
    sns.barplot(data=melted, x="Stage", y="R2", hue="Property")
    plt.ylim(0.0, 1.05)
    plt.ylabel("R² score")
    plt.title("Per-Stage R² Comparison")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_heatmap(df: pd.DataFrame, out_path: str) -> None:
    heatmap_data = df.set_index("Stage")[PROPERTY_COLUMNS]
    plt.figure(figsize=(5, 3))
    sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="Blues", vmin=0.0, vmax=1.0)
    plt.title("R² by Stage and Property")
    plt.ylabel("Stage")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize probe_experiment R² scores")
    parser.add_argument("--csv", help="Path to probe_results.csv")
    parser.add_argument("--outdir", default=None, help="Directory to save figures (default: same as CSV)")
    args = parser.parse_args()

    csv_path = os.path.abspath(args.csv)
    outdir = os.path.abspath(args.outdir) if args.outdir else os.path.dirname(csv_path)
    os.makedirs(outdir, exist_ok=True)

    df = load_data(csv_path)

    plot_line(df, os.path.join(outdir, "probe_r2_line.png"))
    plot_grouped_bar(df, os.path.join(outdir, "probe_r2_grouped_bar.png"))
    plot_heatmap(df, os.path.join(outdir, "probe_r2_heatmap.png"))


if __name__ == "__main__":
    main()
