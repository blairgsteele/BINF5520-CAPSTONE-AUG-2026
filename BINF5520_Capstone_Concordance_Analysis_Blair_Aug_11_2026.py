# Capstone Project — Concordance Analysis
# Blair Steele
# August 11, 2026

#!/usr/bin/env python3
"""
concordance_analysis.py
------------------------------------------------------------------------
Capstone Project - Viral Lineage Surveillance
STUDENT STARTER TEMPLATE

Goal: turn your pipeline's per-sample outputs into:
    1. A concordance table  (known lineage vs. each tool's call)
    2. A figure             (concordance vs. sequencing depth/coverage)
    3. Case study notes      (why did the mismatches happen?)

The function signatures, docstrings, and file layout below are provided
so your code stays organized and matches what run_pipeline.sh produces.
The logic inside each function is up to you, that's the assignment.

Expected inputs (matching run_pipeline.sh's output layout):
    metadata.csv                       accession,reported_lineage,primer_scheme
    lineage/<acc>.nextclade.tsv        Nextclade output (tab-separated)
    lineage/<acc>.pangolin.csv         Pangolin output (comma-separated)
    qc/<acc>.depth.txt                 samtools depth -a output (pos, depth)

Expected outputs:
    results/concordance_table.csv
    results/concordance_vs_depth.png
    results/case_studies.txt

Reminders from lecture / the slide deck:
    - Use a colorblind-safe palette for the figure (e.g. Wong 2011).
    - Keep code modular (one function per stage) so you can test and
      debug each piece independently.
------------------------------------------------------------------------
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# ---- Paths ----------------------------------------------------------------
METADATA_PATH = Path("metadata.csv")
LINEAGE_DIR = Path("lineage")
QC_DIR = Path("qc")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

# TODO: decide on and document any thresholds you use to flag "low coverage"
# or "low breadth" in your case studies. Keep them as named constants here
# (not magic numbers buried in a function) so they're easy to justify/tune.

Min_Depth_Threshold = 10

def load_metadata() -> pd.DataFrame:
    """ Load metadata.csv into a DataFrame with one row per sample. """
    df = pd.read_csv(METADATA_PATH)
    df.columns = df.columns.str.strip().str.lower() #strip whitespace & normalize col names
    return df

def load_nextclade_call(accession: str) -> str | None:
    """ Read Nextclade TSV and return it Nextclade lineage/clade call."""
    nextclade_file = LINEAGE_DIR / f"{accession}.nextclade.tsv"
    if not nextclade_file.exists() or nextclade_file.stat().st_size == 0:
        return None

    df = pd.read_csv(nextclade_file, sep="\t")
    for col in ["Nextclade_pango", "clade", "lineage"]:
        if col in df.columns and pd.notna(df[col].iloc[0]):
            return str(df[col].iloc[0]).strip()
    return None

def load_pangolin_call(accession: str) -> str | None:
    """ Read Pangolin CSV and return called Pango lineage."""
    pangolin_file = LINEAGE_DIR / f"{accession}.pangolin.csv"
    if not pangolin_file.exists() or pangolin_file.stat().st_size == 0:
        return None

    df = pd.read_csv(pangolin_file)
    if "lineage" in df.columns and pd.notna(df["lineage"].iloc[0]):
        return str(df["lineage"].iloc[0]).strip()
    return None

def load_depth_stats(accession: str) -> tuple[float | None, float | None]:
    """ Calcuate mean read depth and breadth of coverage (>0x) from samtools depth. """
    depth_file = QC_DIR / f"{accession}.depth.txt"
    if not depth_file.exists():
        return None, None 

    total_positions = 0
    total_depth = 0
    covered_positions = 0

    # Read line-by-line for memory efficiency
    with open(depth_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3: 
                depth = int(parts[2])
                total_positions += 1
                total_depth += depth
                if depth > 0:
                    covered_positions += 1

    if total_positions == 0: 
        return 0.0, 0.0

    mean_depth = round(total_depth / total_positions, 2)
    breadth = round(covered_positions / total_positions, 4)
    return mean_depth, breadth

# Evaluate match between reported and called lineage
def is_lineage_match(reported: str | None, called: str | None) -> bool:
    if not reported or not called or called in ["None", "Unassigned", ""]:
        return False

    rep = reported.strip().upper()
    cal = called.strip().upper()

    # Case-sensitive matching
    if rep == cal:
        return True
    
    # Match sub-lineage & parent hierarchy
    if cal.startswith(rep + ".") or rep.startswith(cal + "."):
        return True

    return False

def build_concordance_table(metadata: pd.DataFrame) -> pd.DataFrame:
    """ 
    Assembles one row per sample comparing the known
    lineage against each tool's independent call, plus coverage context.
    """

    records = []

    for _, row in metadata.iterrows():
        acc = str(row["accession"]).strip()
        reported = str(row["reported_lineage"]).strip() if pd.notna(row["reported_lineage"]) else None

        # Fetch CLI outputs
        nextclade_call = load_nextclade_call(acc)
        pangolin_call = load_pangolin_call(acc)
        mean_depth, breadth = load_depth_stats(acc)

        # Match evaluation logic
        nextclade_match = is_lineage_match(reported, nextclade_call)
        pangolin_match = is_lineage_match(reported, pangolin_call)

        # Store row data matching requested format
        records.append({
            "accession": acc,
            "reported_lineage": reported,
            "nextclade_call": nextclade_call,
            "pangolin_call": pangolin_call,
            "nextclade_match": nextclade_match,
            "pangolin_match": pangolin_match,
            "mean_depth": mean_depth,
            "breadth_covered": breadth
        })

    return pd.DataFrame(records)

def plot_concordance_vs_depth(table: pd.DataFrame) -> None:
    """
    Plot of mean depth v. breadth of coverage for Nextclade & Pangolin calls.
    The colour-blind-safe pallet (Wong, 2011) was used to differentiate matches and mismatches
    """
    COLOR_MATCH = "#009E73"
    COLOR_MISMATCH = "#D55E00"

    plt.figure(figsize=(10,6))

    # Iterate through each row
    for _, row in table.iterrows():
        depth = row["mean_depth"]
        breadth = row["breadth_covered"] * 100 #convert fraction to percentage

        # Skip rows with missing depth data
        if pd.isna(depth) or pd.isna(breadth):
            continue

        # Nextclade point
        nc_color = COLOR_MATCH if row["nextclade_match"] else COLOR_MISMATCH
        nc_marker = "o" if row["nextclade_match"] else "X"
        plt.scatter(
            depth, breadth, color=nc_color, marker=nc_marker,
            s=80, alpha=0.8, edgecolors="k", linewidths=0.5
        )       

        # Pangolin point
        pango_color = COLOR_MATCH if row["pangolin_match"] else COLOR_MISMATCH 
        pango_marker = "o" if row["pangolin_match"] else "X" 
        plt.scatter(
            depth, breadth, color=pango_color, marker=pango_marker,
            s=80, alpha=0.8, edgecolors="k", linewidths=0.5
        )               

    # Reference threshold
    plt.axvline(x=10, color="gray", linestyle="--", alpha=0.7, label="10x Depth Threshold")
    plt.axhline(y=90, color="darkgray", linestyle=":", alpha=0.7, label="90% Coverage Breadth") 

    plt.title("SARS-CoV-2 Lineage Call Concordance v. Sequencing Coverage", fontsize=14, fontweight="bold") 
    plt.xlabel("Mean Sequencing Depth (x)", fontsize=12) 
    plt.ylabel("Breadth of Coverage (%)", fontsize=12) 
    plt.grid(True, linestyle="--", alpha=0.5)

    # Colorblind-safe legend
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="Concordant Call (Match)",
                markerfacecolor=COLOR_MATCH, markeredgecolor="k", markersize=10), 
        Line2D([0], [0], marker="X", color="w", label="Discrepant Call (Mismatch)",
                markerfacecolor=COLOR_MISMATCH, markeredgecolor="k", markersize=10), 
        Line2D([0], [0], color="gray", linestyle="--", label="10x Depth Cutoff"), 
        Line2D([0], [0], color="darkgray", linestyle=":", label="90% Breadth Cutoff"), 
    ]
    plt.legend(handles=legend_elements, loc="lower right")

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "concordance_vs_depth.png", dpi=300)
    plt.close()

LOW_DEPTH_THRESHOLD = 10.0
LOW_BREADTH_THRESHOLD = 0.90 


def draft_case_studies(table: pd.DataFrame) -> None:
    """
    Write results/case_studies.txt: for each sample where at least one
    tool missed the reported lineage, note the relevant coverage stats
    and your explanation for the mismatch.

    You need at least 2 case studies in your final write-up, but this
    function can draft notes for every mismatch so you have material to
    choose from. Consider: what pattern points to low coverage? What
    points to primer dropout vs. a real disagreement between tools?
    """

    output_path = RESULTS_DIR / "case_studies.txt"

    # Filter for samples where at least one tool failed to match
    mismatches = table[(table["nextclade_match"] == False) | (table["pangolin_match"] == False)]

    with open(output_path, "w") as f:
        f.write("Capstone Project: Lineage Discrepancy & Mismatch Case Studies\n\n")

        if mismatches.empty:
            f.write("No lineage mismatches detected across analyzed samples. \n")

        for idx, row in mismatches.iterrows():
            acc = row["accession"]
            rep = row["reported_lineage"]
            nc_call = row["nextclade_call"]
            pg_call = row["pangolin_call"]
            depth = row["mean_depth"]
            breadth = row["breadth_covered"]

            f.write(f"SAMPLE ACCESSION: {acc}\n") 
            f.write(f" - Reported Lineage : {rep}\n") 
            f.write(f" - Nextclade Call : {nc_call} (Match: {row['nextclade_match']})\n")
            f.write(f" - Pangolin Call : {pg_call} (Match: {row['pangolin_match']})\n") 
            f.write(f" - Mean Depth : {depth}x\n") 
            f.write(f" - Coverage Breadth : {breadth * 100:.2f}%\n") 

            f.write(" - DIAGNOSTIC ANALYSIS: \n")

            # Diagnostic Checks
            reasons = []
            if pd.isna(depth) or depth < LOW_DEPTH_THRESHOLD:
                reasons.append(
                    f"Low mean depth ({depth}x < {LOW_DEPTH_THRESHOLD}x). Insufficient read depth forces iVar to insert ambiguous 'N' bases into the consensus FASTA, breaking lineage assignment."
                )
            if pd.isna(breadth) or breadth < LOW_BREADTH_THRESHOLD:
                reasons.append(
                    f"Low genome coverage breadth ({breadth*100:.2f}% < {LOW_BREADTH_THRESHOLD*100}%). Likely points to amplicon dropout or primer binding site mutations preventing amplification."
                )
            if row["nextclade_match"] != row["pangolin_match"]:
                reasons.append(
                    "Tool Disagreement: Nextclade uses mutation distance on a reference tree, whereas Pangolin uses pUSHER phylogenetic placement. Resolution differences can cause nominal designation mismatches."
                )
            if not reasons: 
                reasons.append("High coverage sample mismatch. Investigate potential outdated ARTIC primer scheme or lineage database stale definitions.")

            for r in reasons:
                f.write(f" * {r}\n")
            f.write("\n" + "-"*72 + "\n\n")


def main() -> None:
    metadata = load_metadata()
    table = build_concordance_table(metadata)

    table.to_csv(RESULTS_DIR / "concordance_table.csv", index=False)
    plot_concordance_vs_depth(table)
    draft_case_studies(table)

    # Summary
    total_samples = len(table)
    nc_acc = (table["nextclade_match"].sum() / total_samples) * 100 if total_samples > 0 else 0
    pg_acc = (table["pangolin_match"].sum() / total_samples) * 100 if total_samples > 0 else 0

    print("Pipeline Concordance Analysis Complete")
    print(f"Total Samples Processed : {total_samples}")
    print(f"Nextclade Accuracy : {nc_acc:.2f}%")
    print(f"Pangolin Accuracy : {pg_acc:.2f}%")
    print(f"Results written to : {RESULTS_DIR.resolve()}")

if __name__ == "__main__":
    main()
