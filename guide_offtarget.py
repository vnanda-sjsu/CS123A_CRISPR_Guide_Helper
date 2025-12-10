import argparse
import csv

def load_fasta_sequence(path):
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)

def find_guide_hits(chrom_seq, guide, max_mismatches=3):
    guide = guide.strip().upper()
    g_len = len(guide)
    hits = []
    chrom_len = len(chrom_seq)
    for i in range(chrom_len - g_len + 1):
        window = chrom_seq[i:i + g_len]
        mismatches = 0
        for a, b in zip(window, guide):
            if a != b:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            hits.append((i, mismatches))
    return hits

def classify_region(hit_start, guide_len, gene_start, gene_end):
    hit_end = hit_start + guide_len
    if hit_end <= gene_start:
        return "upstream"
    if hit_start >= gene_end:
        return "downstream"
    return "inside"

def analyze_guide_for_offtargets(chrom_seq, guide, gene_start, gene_end, max_mismatches=3):
    hits = find_guide_hits(chrom_seq, guide, max_mismatches)
    guide_len = len(guide.strip())
    on_target_hits = 0
    off_target_hits = 0
    best_mismatches = None
    for pos, mismatches in hits:
        region = classify_region(pos, guide_len, gene_start, gene_end)
        if region == "inside":
            on_target_hits += 1
        else:
            off_target_hits += 1
        if best_mismatches is None or mismatches < best_mismatches:
            best_mismatches = mismatches
    total_hits = len(hits)
    if best_mismatches is None:
        best_mismatches = -1
    return {
        "guide": guide.strip().upper(),
        "total_hits": total_hits,
        "on_target_hits": on_target_hits,
        "off_target_hits": off_target_hits,
        "best_mismatches": best_mismatches,
    }

def process_guides(fasta_path, guides_path, output_path, gene_start, gene_end, max_mismatches=3):
    chrom_seq = load_fasta_sequence(fasta_path)
    results = []
    with open(guides_path) as f:
        for line in f:
            guide = line.strip()
            if not guide:
                continue
            result = analyze_guide_for_offtargets(chrom_seq, guide, gene_start, gene_end, max_mismatches)
            results.append(result)
    fieldnames = ["guide", "total_hits", "on_target_hits", "off_target_hits", "best_mismatches"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in fieldnames})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--guides", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--gene-start", type=int, required=True)
    parser.add_argument("--gene-end", type=int, required=True)
    parser.add_argument("--max-mismatches", type=int, default=3)
    args = parser.parse_args()
    process_guides(
        fasta_path=args.fasta,
        guides_path=args.guides,
        output_path=args.out,
        gene_start=args.gene_start,
        gene_end=args.gene_end,
        max_mismatches=args.max_mismatches,
    )

if __name__ == "__main__":
    main()
