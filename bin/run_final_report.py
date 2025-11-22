#!/usr/bin/env python3

"""
run_final_report.py (ADVANCED VERSION)
(Module 3.5 - Final Reporting)

Aggregates all data streams to create a final JSON and HTML report.
Gracefully handles missing input files if a pipeline branch fails.

1.  ichorCNA output (*.tfx.txt) -> CNA-based Purity
2.  Global Fragmentomics (*.global_hist.tsv) -> Global QC metric
3.  Annotated VCF (*.frag.vcf.gz) -> SNV/INDEL-based Purity
    
This script applies the clinical logic from the plan:
-   Global Fragment % = % of fragments between `min_bp` and `max_bp`
-   SNV/INDEL Purity is calculated *only* from variants that are:
    1.  `PASS`
    2.  `FRAG_KS_PVAL < ks_pval` (i.e., statistically shorter, "Somatic")
    3.  VAF > 1% (to avoid noise)
-   **Consensus Purity**:
    -   If CNA Purity > `cna_lod` (e.g., 3%), trust CNA Purity.
    -   Else, trust SNV/INDEL Purity (more sensitive at low TFX).
"""

import argparse
import json
import pysam
import statistics
import pandas as pd
import os

def parse_ichor_tfx(tfx_path):
    """Reads the ichorCNA tumor fraction file."""
    if not tfx_path or not os.path.exists(tfx_path):
        return 0.0, "ichorCNA failed or file not found."
    try:
        df = pd.read_csv(tfx_path, sep="\t")
        if "tumorFraction" in df.columns and not df.empty:
            return float(df["tumorFraction"].iloc[0]), "Success"
    except Exception as e:
        return 0.0, f"Error parsing ichorCNA file: {e}"
    return 0.0, "No tumor fraction found in ichorCNA file."

def parse_global_histogram(hist_path, min_bp, max_bp):
    """
    Parses the global fragment histogram.
    Returns: total_reads (int), short_fragment_percent (float), status (str)
    """
    if not hist_path or not os.path.exists(hist_path):
        return 0, 0.0, "Global Fragmentomics failed or file not found."
    try:
        df = pd.read_csv(hist_path, sep="\t")
        if df.empty or "fragment_size" not in df.columns or "count" not in df.columns:
            return 0, 0.0, "Histogram file was empty or malformed."
            
        total_reads = df["count"].sum()
        if total_reads == 0:
            return 0, 0.0, "Histogram file contained no reads."
            
        short_reads = df[
            (df["fragment_size"] >= min_bp) & 
            (df["fragment_size"] <= max_bp)
        ]["count"].sum()
        
        short_percent = (short_reads / total_reads) * 100
        return int(total_reads), round(short_percent, 2), "Success"
        
    except Exception as e:
        return 0, 0.0, f"Error parsing histogram: {e}"

def parse_annotated_vcf(vcf_path, ks_pval_cutoff, min_vaf):
    """
    Parses the fragment-annotated VCF.
    Calculates purity (2 * VAF) from high-confidence "Somatic" (short) variants.
    Returns: snv_purity (float), total_pass (int), somatic_found (int), status (str)
    """
    if not vcf_path or not os.path.exists(vcf_path):
        return 0.0, 0, 0, "Variant Fragmentomics failed or file not found."

    somatic_vafs = []
    total_pass_variants = 0
    somatic_variants_found = 0
    
    try:
        vcf = pysam.VariantFile(vcf_path, "r")
        for rec in vcf:
            if "PASS" not in rec.filter:
                continue
            total_pass_variants += 1
            
            # Check for our fragmentomics annotations
            if "FRAG_KS_PVAL" not in rec.info or rec.info["FRAG_KS_PVAL"] is None:
                continue

            # --- Clinical Logic ---
            is_shorter = rec.info["FRAG_KS_PVAL"] < ks_pval_cutoff
            
            smp = rec.samples[0]
            if "AD" in smp and smp["AD"] and sum(smp["AD"]) > 1: # Need at least 2 reads
                refc, altc = smp["AD"][0], smp["AD"][1]
                vaf = altc / (refc + altc)
                
                if is_shorter and vaf >= min_vaf:
                    somatic_vafs.append(vaf)
                    somatic_variants_found += 1
            
        if not somatic_vafs:
            return 0.0, total_pass_variants, 0, "Success (No high-confidence somatic variants found)"
            
        # Purity = Median(2 * VAF)
        median_somatic_vaf = statistics.median(somatic_vafs)
        snv_purity = min(1.0, 2 * median_somatic_vaf) # Cap at 100%
        
        return round(snv_purity, 4), total_pass_variants, somatic_variants_found, "Success"
        
    except Exception as e:
        return 0.0, 0, 0, f"Error parsing annotated VCF: {e}"

def load_template(template_path):
    """Loads an HTML template file."""
    try:
        with open(template_path, 'r') as f:
            return f.read()
    except FileNotFoundError:
        print(f"Warning: HTML template not found at {template_path}. Using basic report.", file=sys.stderr)
        return "<h1>{sample_id} Report</h1><pre>{json_data}</pre>"

def main(args):
    
    # --- 1. Run All Analyses (Gracefully handle missing files) ---
    cna_purity, cna_status = parse_ichor_tfx(args.ichor_tfx)
    total_reads, short_frag_pct, hist_status = parse_global_histogram(args.global_hist, args.min_bp, args.max_bp)
    snv_purity, pass_vars, somatic_vars, vcf_status = parse_annotated_vcf(args.annot_vcf, args.ks_pval, 0.01)

    # --- 2. Final Consensus Logic (Module 3.5 - ADVANCED) ---
    consensus_purity = 0.0
    consensus_method = "None"
    
    if cna_purity >= args.cna_lod:
        consensus_purity = cna_purity
        consensus_method = "CNA-based (ichorCNA)"
    else:
        consensus_purity = snv_purity
        consensus_method = "SNV/INDEL-based (Fragmentomics)"
    
    # --- 3. Compile JSON Report ---
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

    # Write JSON output
    with open(args.out_json, 'w') as f:
        json.dump(report_data, f, indent=4)

    # --- 4. Compile HTML Report ---
    template = load_template(args.template_html)
    
    # Populate template
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

    # Write HTML output
    with open(args.out_html, 'w') as f:
        f.write(html_out)

    print(f"Report generated for {args.sample_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate all TFX/Frag data into a final report.")
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    # Optional file args (robust to pipeline failures)
    parser.add_argument("--ichor_tfx", required=False, default=None, help="ichorCNA *.tfx.txt file")
    parser.add_argument("--annot_vcf", required=False, default=None, help="Fragment-annotated *.frag.vcf.gz file")
    parser.add_argument("--global_hist", required=False, default=None, help="Global fragment histogram *.global_hist.tsv file")
    # Param args
    parser.add_argument("--min_bp", type=int, required=True, help="Min BP for global frag %")
    parser.add_argument("--max_bp", type=int, required=True, help="Max BP for global frag %")
    parser.add_argument("--ks_pval", type=float, required=True, help="P-value cutoff for somatic variant fragment test")
    parser.add_argument("--cna_lod", type=float, required=True, help="TFX cutoff to trust CNA over SNV")
    # Output args
    parser.add_argument("--template_html", required=True, help="Path to HTML report template")
    parser.add_argument("--out_json", required=True, help="Output JSON report file")
    parser.add_argument("--out_html", required=True, help="Output HTML report file")
    
    args = parser.parse_args()
    main(args)
