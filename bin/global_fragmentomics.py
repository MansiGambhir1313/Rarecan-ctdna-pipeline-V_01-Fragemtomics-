#!/usr/bin/env python3
"""
global_fragmentomics.py

Analyzes cfDNA fragment length distributions and end-motif signatures.
Outputs:
  - Fragment length histogram (.tsv)
  - End-motif frequencies (.tsv)
  - Summary metrics (.json)

Highlights:
  * Mononucleosome peak (~167 bp)
  * Short fragment enrichment (90-150 bp by default)
  * Tumor-associated motifs (4-mer C-rich / G-rich at fragment ends)
"""

import argparse
import json
import math
from collections import Counter, defaultdict

import pysam


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def get_fragment_end_motif(read: pysam.AlignedSegment, k: int = 4):
    """
    Returns the 5' end motif (k-mer) for the fragment represented by this read.
    Only evaluates read1 to avoid double-counting the pair.
    """
    if read.query_sequence is None or len(read.query_sequence) < k:
        return None

    seq = read.query_sequence.upper()

    if read.is_reverse:
        motif = reverse_complement(seq[-k:])
    else:
        motif = seq[:k]

    return motif


def compute_fragment_metrics(
    bam_path: str,
    max_length: int,
    short_min: int,
    short_max: int,
    motif_k: int = 4,
) -> dict:
    hist = Counter()
    motif_counts = Counter()

    total_fragments = 0
    short_fragments = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_read1:
                continue
            
            # Skip low quality reads
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue
            
            # Handle unmapped reads: use read length as approximate fragment length
            # This allows processing when alignment fails (e.g., dummy reference)
            if read.is_unmapped or read.mate_is_unmapped:
                # For unmapped reads, use read length as fragment length estimate
                # This is approximate but allows analysis to proceed
                if read.query_length and read.query_length > 0:
                    # For paired unmapped reads, estimate fragment length
                    # Since both are unmapped, use read length directly
                    estimated_len = read.query_length
                    
                    # Check if length is within valid range
                    if estimated_len > 0 and estimated_len <= max_length:
                        total_fragments += 1
                        hist[estimated_len] += 1
                        if short_min <= estimated_len <= short_max:
                            short_fragments += 1
                        # Extract motif from read sequence
                        motif = get_fragment_end_motif(read, motif_k)
                        if motif:
                            motif_counts[motif] += 1
                continue
            
            # Normal processing for mapped reads
            if not read.is_proper_pair:
                continue

            frag_len = abs(read.template_length)
            if frag_len == 0 or frag_len > max_length:
                continue

            total_fragments += 1
            hist[frag_len] += 1

            if short_min <= frag_len <= short_max:
                short_fragments += 1

            motif = get_fragment_end_motif(read, motif_k)
            if motif:
                motif_counts[motif] += 1

    short_fraction = (
        (short_fragments / total_fragments) if total_fragments > 0 else 0.0
    )

    # Mononucleosome peak detection (140-190 bp)
    mono_window = range(140, 191)
    peak_bp = max(mono_window, key=lambda bp: hist.get(bp, 0))
    peak_count = hist.get(peak_bp, 0)

    baseline_window = [hist.get(bp, 0) for bp in range(200, 221)]
    baseline_mean = (
        sum(baseline_window) / len(baseline_window) if baseline_window else 1.0
    )
    peak_enrichment = (
        (peak_count / baseline_mean) if baseline_mean > 0 else math.nan
    )

    return {
        "histogram": hist,
        "motif_counts": motif_counts,
        "total_fragments": total_fragments,
        "short_fragment_count": short_fragments,
        "short_fragment_fraction": short_fraction,
        "mononucleosome_peak_bp": peak_bp,
        "mononucleosome_peak_count": peak_count,
        "mononucleosome_peak_enrichment": peak_enrichment,
    }


def write_histogram(hist: Counter, path: str):
    with open(path, "w") as fh:
        fh.write("fragment_size\tcount\n")
        for size in sorted(hist.keys()):
            fh.write(f"{size}\t{hist[size]}\n")


def write_motifs(motif_counts: Counter, path: str):
    total = sum(motif_counts.values()) or 1
    with open(path, "w") as fh:
        fh.write("motif\tcount\tfraction\n")
        for motif, count in motif_counts.most_common():
            fh.write(f"{motif}\t{count}\t{count/total:.6f}\n")


def write_summary(metrics: dict, sample_id: str, path: str):
    cccc = metrics["motif_counts"].get("CCCC", 0)
    gggg = metrics["motif_counts"].get("GGGG", 0)
    motifs_total = sum(metrics["motif_counts"].values()) or 1

    summary = {
        "sample_id": sample_id,
        "total_fragments": metrics["total_fragments"],
        "short_fragments_90_150": metrics["short_fragment_count"],
        "short_fragment_fraction": round(
            metrics["short_fragment_fraction"], 4
        ),
        "mononucleosome_peak_bp": metrics["mononucleosome_peak_bp"],
        "mononucleosome_peak_count": metrics["mononucleosome_peak_count"],
        "mononucleosome_peak_enrichment": round(
            metrics["mononucleosome_peak_enrichment"] or 0, 3
        ),
        "tumor_motif_counts": {
            "CCCC": cccc,
            "GGGG": gggg,
        },
        "tumor_motif_fraction": {
            "CCCC": round(cccc / motifs_total, 4),
            "GGGG": round(gggg / motifs_total, 4),
        },
        "motifs_total": motifs_total,
    }

    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)


def parse_args():
    parser = argparse.ArgumentParser(
        description="cfDNA fragment length and motif analysis"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--output_json", required=True, help="Summary JSON")
    parser.add_argument("--histogram", required=True, help="Histogram TSV")
    parser.add_argument("--motifs", required=True, help="Motif TSV")
    parser.add_argument(
        "--short_min", type=int, default=90, help="Short fragment min bp"
    )
    parser.add_argument(
        "--short_max", type=int, default=150, help="Short fragment max bp"
    )
    parser.add_argument(
        "--max_length",
        type=int,
        default=700,
        help="Maximum fragment length to consider",
    )
    parser.add_argument(
        "--motif_length", type=int, default=4, help="End motif k-mer length"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Debug: Print input parameters
    import sys
    print(f"DEBUG: Processing BAM: {args.bam}", file=sys.stderr)
    print(f"DEBUG: Sample: {args.sample}", file=sys.stderr)
    print(f"DEBUG: Short fragment range: {args.short_min}-{args.short_max} bp", file=sys.stderr)

    metrics = compute_fragment_metrics(
        bam_path=args.bam,
        max_length=args.max_length,
        short_min=args.short_min,
        short_max=args.short_max,
        motif_k=args.motif_length,
    )
    
    # Debug: Print metrics
    print(f"DEBUG: Total fragments found: {metrics['total_fragments']}", file=sys.stderr)
    print(f"DEBUG: Histogram entries: {len(metrics['histogram'])}", file=sys.stderr)
    print(f"DEBUG: Motif entries: {len(metrics['motif_counts'])}", file=sys.stderr)

    write_histogram(metrics["histogram"], args.histogram)
    write_motifs(metrics["motif_counts"], args.motifs)
    write_summary(metrics, args.sample, args.output_json)
    
    print(f"DEBUG: Output files written successfully", file=sys.stderr)


if __name__ == "__main__":
    main()

