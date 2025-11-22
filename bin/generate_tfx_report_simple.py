#!/usr/bin/env python3
"""
Simplified TFX report generator that works without pysam
Generates report from available fragmentomics data
"""

import json
import os
import sys
import pandas as pd
from datetime import datetime

def parse_ichor_tfx(tfx_path):
    """Reads the ichorCNA tumor fraction file."""
    if not tfx_path or not os.path.exists(tfx_path) or tfx_path == "":
        return 0.0, "ichorCNA file not available"
    try:
        df = pd.read_csv(tfx_path, sep="\t")
        if "tumorFraction" in df.columns and not df.empty:
            return float(df["tumorFraction"].iloc[0]), "Success"
    except Exception as e:
        return 0.0, f"Error parsing ichorCNA file: {e}"
    return 0.0, "No tumor fraction found"

def parse_global_histogram(hist_path, min_bp, max_bp):
    """Parses the global fragment histogram."""
    if not hist_path or not os.path.exists(hist_path):
        return 0, 0.0, "Global Fragmentomics file not found"
    try:
        df = pd.read_csv(hist_path, sep="\t")
        if df.empty:
            return 0, 0.0, "Histogram file is empty (placeholder data)"
        
        # Handle different column formats
        if "fragment_size" in df.columns and "count" in df.columns:
            # Numeric fragment_size
            total_reads = df["count"].sum()
            short_reads = df[
                (df["fragment_size"] >= min_bp) & 
                (df["fragment_size"] <= max_bp)
            ]["count"].sum()
        elif "fragment_length" in df.columns and "count" in df.columns:
            # Check if fragment_length is numeric or range format (e.g., "50-100")
            try:
                # Try numeric first
                df["fragment_length_num"] = pd.to_numeric(df["fragment_length"], errors='coerce')
                if df["fragment_length_num"].notna().any():
                    total_reads = df["count"].sum()
                    short_reads = df[
                        (df["fragment_length_num"] >= min_bp) & 
                        (df["fragment_length_num"] <= max_bp)
                    ]["count"].sum()
                else:
                    # Handle range format like "50-100", "100-150"
                    total_reads = df["count"].sum()
                    short_reads = 0
                    for idx, row in df.iterrows():
                        frag_range = str(row["fragment_length"])
                        if "-" in frag_range:
                            try:
                                start, end = map(int, frag_range.split("-"))
                                if start <= max_bp and end >= min_bp:
                                    # Range overlaps with short fragment range
                                    overlap_start = max(start, min_bp)
                                    overlap_end = min(end, max_bp)
                                    if overlap_start <= overlap_end:
                                        short_reads += row["count"]
                            except:
                                pass
            except Exception as e:
                return 0, 0.0, f"Error parsing fragment_length: {e}"
        else:
            return 0, 0.0, "Histogram file format not recognized"
        
        if total_reads == 0:
            return 0, 0.0, "Histogram contains no reads (placeholder data)"
        
        short_percent = (short_reads / total_reads) * 100
        return int(total_reads), round(short_percent, 2), "Success"
    except Exception as e:
        return 0, 0.0, f"Error parsing histogram: {e}"

def parse_annotated_vcf_simple(vcf_path):
    """Simplified VCF parser that doesn't require pysam."""
    if not vcf_path or not os.path.exists(vcf_path) or vcf_path == "":
        return 0.0, 0, 0, "Variant Fragmentomics file not available"
    
    # For now, return placeholder values
    # In a real scenario, this would parse the VCF
    return 0.0, 0, 0, "VCF parsing requires pysam (not available in current environment)"

def load_template(template_path):
    """Loads an HTML template file."""
    try:
        with open(template_path, 'r') as f:
            return f.read()
    except FileNotFoundError:
        return "<h1>{sample_id} Report</h1><pre>{json_data}</pre>"

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate TFX report from fragmentomics data")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--ichor_tfx", default="")
    parser.add_argument("--annot_vcf", default="")
    parser.add_argument("--global_hist", required=True)
    parser.add_argument("--min_bp", type=int, default=90)
    parser.add_argument("--max_bp", type=int, default=150)
    parser.add_argument("--ks_pval", type=float, default=0.05)
    parser.add_argument("--cna_lod", type=float, default=0.03)
    parser.add_argument("--template_html", required=True)
    parser.add_argument("--out_json", required=True)
    parser.add_argument("--out_html", required=True)
    
    args = parser.parse_args()
    
    # Parse inputs
    cna_purity, cna_status = parse_ichor_tfx(args.ichor_tfx if args.ichor_tfx else None)
    total_reads, short_frag_pct, hist_status = parse_global_histogram(
        args.global_hist, args.min_bp, args.max_bp
    )
    snv_purity, pass_vars, somatic_vars, vcf_status = parse_annotated_vcf_simple(
        args.annot_vcf if args.annot_vcf else None
    )
    
    # Consensus logic
    consensus_purity = 0.0
    consensus_method = "None"
    
    if cna_purity >= args.cna_lod:
        consensus_purity = cna_purity
        consensus_method = "CNA-based (ichorCNA)"
    else:
        consensus_purity = snv_purity
        consensus_method = "SNV/INDEL-based (Fragmentomics)"
    
    # Compile JSON report
    report_data = {
        "sample_id": args.sample_id,
        "ctdna_purity_consensus": round(consensus_purity, 4),
        "ctdna_purity_percent": f"{consensus_purity * 100:.2f}%",
        "consensus_method": consensus_method,
        "methods": {
            "cna_based_purity": {
                "value": round(cna_purity, 4),
                "percent": f"{cna_purity * 100:.2f}%",
                "tool": "ichorCNA",
                "status": cna_status
            },
            "snv_indel_based_purity": {
                "value": round(snv_purity, 4),
                "percent": f"{snv_purity * 100:.2f}%",
                "tool": "GATK + Variant Fragmentomics (SNV/INDEL)",
                "total_pass_variants": pass_vars,
                "somatic_variants_used": somatic_vars,
                "status": vcf_status
            }
        },
        "qc_metrics": {
            "global_fragmentomics": {
                "total_reads": int(total_reads),
                "short_fragment_percent": short_frag_pct,
                "definition": f"% of fragments {args.min_bp}-{args.max_bp} bp",
                "status": hist_status
            },
            "variant_fragmentomics": {
                "ks_pval_cutoff": args.ks_pval
            },
            "cna_lod_cutoff": args.cna_lod
        }
    }
    
    # Write JSON
    with open(args.out_json, 'w') as f:
        json.dump(report_data, f, indent=4)
    
    # Generate HTML
    template = load_template(args.template_html)
    html_out = template.replace("{{sample_id}}", report_data["sample_id"])
    html_out = html_out.replace("{{ctdna_purity_percent}}", report_data["ctdna_purity_percent"])
    html_out = html_out.replace("{{consensus_method}}", report_data["consensus_method"])
    html_out = html_out.replace("{{cna_based_percent}}", report_data["methods"]["cna_based_purity"]["percent"])
    html_out = html_out.replace("{{snv_based_percent}}", report_data["methods"]["snv_indel_based_purity"]["percent"])
    html_out = html_out.replace("{{somatic_variants_used}}", str(report_data["methods"]["snv_indel_based_purity"]["somatic_variants_used"]))
    html_out = html_out.replace("{{total_pass_variants}}", str(report_data["methods"]["snv_indel_based_purity"]["total_pass_variants"]))
    html_out = html_out.replace("{{short_fragment_percent}}", str(report_data["qc_metrics"]["global_fragmentomics"]["short_fragment_percent"]))
    html_out = html_out.replace("{{fragment_qc_definition}}", report_data["qc_metrics"]["global_fragmentomics"]["definition"])
    html_out = html_out.replace("{{json_data}}", json.dumps(report_data, indent=4))
    
    with open(args.out_html, 'w') as f:
        f.write(html_out)
    
    print(f"TFX report generated for {args.sample_id}")
    print(f"  Consensus Purity: {report_data['ctdna_purity_percent']}")
    print(f"  Method: {consensus_method}")
    print(f"  Short Fragment %: {short_frag_pct}%")

if __name__ == "__main__":
    main()

