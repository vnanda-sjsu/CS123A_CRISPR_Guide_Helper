"""
Person 2 Script: Computes GC content, PAM distance, self-complementarity,
and merges off-target results into a unified TSV output.
Author: Jumana
"""

import csv
import argparse


def gc_content(seq: str) -> float:
    """Return GC content as a percent (0–100)."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return 100.0 * gc / len(seq)


def calculate_self_complementarity(seq: str, k: int = 4) -> int:
    """
    Counts the number of self-binding k-mer substrings (simple heuristic).
    If this returns > 0, you can treat that as "yes" for self-complementarity.
    """
    seq = seq.upper()
    count = 0
    for i in range(len(seq) - 2 * k + 1):
        sub = seq[i:i + k]
        if sub in seq[i + k:]:
            count += 1
    return count


def parse_genomic_location(loc: str):
    """
    Parse 'chr2:60462265' -> ('chr2', 60462265)
    """
    loc = loc.strip()
    if not loc:
        return None, None
    chrom, pos = loc.split(":")
    return chrom, int(pos)


def pam_distance_from_gene_start(seq: str, genomic_loc: str, gene_start: int) -> int:
    """
    Compute distance of PAM from the beginning of the BCL11A gene.

    Assumptions:
    - 'Genomic location' in results.tsv is the start coordinate of the full
      23-bp sequence (20-mer + NGG).
    - PAM is at the 3' end of the sequence (the last 3 bases, NGG).

    So:
        pam_start = genomic_start + (len(seq) - 3)
        distance = abs(pam_start - gene_start)
    """
    chrom, seq_start = parse_genomic_location(genomic_loc)
    if chrom is None or seq_start is None:
        return -1

    seq = seq.upper()
    if len(seq) < 3:
        return -1

    # Assume PAM is the last 3 bases (NGG) of the sequence
    pam_offset = len(seq) - 3
    pam_genomic_pos = seq_start + pam_offset

    return abs(pam_genomic_pos - gene_start)


def load_results_tsv(path: str):
    """
    Load metrics from results.tsv into a dict keyed by guide sequence.
    Also keep genomic location so we can compute PAM distance.
    """
    metrics = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq = row["Target sequence"].strip().upper()
            genomic_loc = row.get("Genomic location", "").strip()

            try:
                gc = float(row["GC content (%)"])
            except (ValueError, KeyError):
                gc = gc_content(seq)

            try:
                self_comp = int(row["Self-complementarity"])
            except (ValueError, KeyError):
                self_comp = calculate_self_complementarity(seq)

            metrics[seq] = {
                "GC_content": gc,
                "self_complementarity": self_comp,
                "genomic_loc": genomic_loc,
            }
    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="Combine off-target hits with GC content, PAM distance from gene start, and self-complementarity into a TSV."
    )
    parser.add_argument(
        "--offtargets",
        default="offtarget_results_bcl11a.csv",
        help="CSV file with columns: guide,total_hits,on_target_hits,off_target_hits,...",
    )
    parser.add_argument(
        "--results",
        default="results.tsv",
        help="TSV file with guide metrics from CHOPCHOP (results.tsv).",
    )
    parser.add_argument(
        "--out",
        default="person2_guide_metrics.tsv",
        help="Output TSV file.",
    )
    parser.add_argument(
        "--gene-start",
        type=int,
        required=True,
        help="Genomic coordinate of the beginning of the BCL11A gene.",
    )

    args = parser.parse_args()

    results_metrics = load_results_tsv(args.results)

    with open(args.offtargets, newline="") as f_in, open(args.out, "w", newline="") as f_out:
        reader = csv.DictReader(f_in)
        writer = csv.writer(f_out, delimiter="\t")

        header = [
            "guide",
            "total_hits",
            "off_target_hits",
            "GC_content_percent",
            "pam_distance_from_BCL11A_start",
            "self_complementarity_score",
        ]
        writer.writerow(header)

        for row in reader:
            guide = row["guide"].strip().upper()
            total_hits = row.get("total_hits", "")
            off_target_hits = row.get("off_target_hits", "")

            m = results_metrics.get(guide)
            if m is not None:
                gc = m["GC_content"]
                self_comp = m["self_complementarity"]
                genomic_loc = m["genomic_loc"]
            else:
                gc = gc_content(guide)
                self_comp = calculate_self_complementarity(guide)
                genomic_loc = ""

            pam_dist = pam_distance_from_gene_start(
                guide, genomic_loc, args.gene_start
            )

            writer.writerow([
                guide,
                total_hits,
                off_target_hits,
                f"{gc:.2f}",
                pam_dist,
                self_comp,
            ])

    print(f"Done. Wrote TSV to {args.out}")


if __name__ == "__main__":
    main()
