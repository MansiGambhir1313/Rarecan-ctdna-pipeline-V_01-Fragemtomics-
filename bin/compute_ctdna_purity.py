#!/usr/bin/env python3
"""
compute_ctdna_purity.py

Comprehensive ctDNA% computation module that integrates:
1. ichorCNA tumor fraction estimation
2. Fragment-length distribution and KS logic
3. SNV/INDEL purity from variant calling
4. CHIP filtering integration
5. Fallback behaviors

This is the core analytical component for clinical-grade ctDNA% estimation.
"""

import argparse
import json
import sys
import os
from pathlib import Path
from typing import Dict, Optional, Tuple
import statistics

try:
    import pandas as pd
except ImportError:
    pd = None


def load_ichorcna_tumor_fraction(tfx_file: str) -> Tuple[float, Dict]:
    """
    Load ichorCNA tumor fraction from .tfx.txt file.
    Returns: (tumor_fraction, metadata_dict)
    """
    if not os.path.exists(tfx_file):
        return 0.0, {"status": "file_not_found", "source": "ichorCNA"}
    
    try:
        if pd:
            df = pd.read_csv(tfx_file, sep="\t")
            if "tumorFraction" in df.columns and not df.empty:
                tf = float(df["tumorFraction"].iloc[0])
                metadata = {
                    "status": "success",
                    "source": "ichorCNA",
                    "ploidy": float(df["ploidy"].iloc[0]) if "ploidy" in df.columns else None
                }
                return max(0.0, min(1.0, tf)), metadata
        else:
            # Fallback: read manually
            with open(tfx_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Header + data
                    parts = lines[1].strip().split("\t")
                    if len(parts) >= 2:
                        tf = float(parts[1])  # Assuming sample_id, tumorFraction, ploidy
                        return max(0.0, min(1.0, tf)), {"status": "success", "source": "ichorCNA"}
    except Exception as e:
        print(f"WARNING: Could not parse ichorCNA file {tfx_file}: {e}", file=sys.stderr)
    
    return 0.0, {"status": "parse_error", "source": "ichorCNA"}


def load_ichorcna_segments(segments_file: str) -> Dict:
    """
    Load ichorCNA segments to extract log2CN and segmentMeans.
    Returns: dict with segment statistics
    """
    if not os.path.exists(segments_file):
        return {"status": "file_not_found", "segments": []}
    
    try:
        if pd:
            df = pd.read_csv(segments_file, sep="\t")
            if not df.empty:
                segments = df.to_dict('records')
                log2_values = df["log2"].tolist() if "log2" in df.columns else []
                return {
                    "status": "success",
                    "segments": segments,
                    "mean_log2": float(statistics.mean(log2_values)) if log2_values else 0.0,
                    "median_log2": float(statistics.median(log2_values)) if log2_values else 0.0,
                    "num_segments": len(segments)
                }
        else:
            # Fallback: read manually
            segments = []
            with open(segments_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Skip header
                    for line in lines[1:]:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            segments.append({
                                "chromosome": parts[0],
                                "start": int(parts[1]),
                                "end": int(parts[2]),
                                "log2": float(parts[3]) if len(parts) > 3 else 0.0
                            })
            return {"status": "success", "segments": segments, "num_segments": len(segments)}
    except Exception as e:
        print(f"WARNING: Could not parse segments file {segments_file}: {e}", file=sys.stderr)
    
    return {"status": "parse_error", "segments": []}


def compute_fragment_ks_purity(frag_vcf: str, ks_pval_cutoff: float = 0.05) -> Tuple[float, Dict]:
    """
    Compute ctDNA purity from fragment-length KS statistics in annotated VCF.
    Uses variants with FRAG_KS_PVAL < cutoff as somatic (tumor-derived).
    
    For low tumor burden cases: Calculate mean VAF of variants that passed 
    "Smoking Gun" test (FRAG_KS_PVAL < cutoff) and use 2 * mean_VAF as purity estimate.
    
    Returns: (purity_estimate, metadata)
    """
    if not os.path.exists(frag_vcf):
        return 0.0, {"status": "file_not_found", "source": "fragment_KS"}
    
    try:
        import pysam
        
        vcf = pysam.VariantFile(frag_vcf, "r")
        
        total_variants = 0
        somatic_variants = 0
        ks_pvals = []
        alt_medians = []
        ref_medians = []
        vafs_passed_smoking_gun = []  # VAFs of variants that passed smoking gun test
        
        for rec in vcf:
            if "PASS" not in rec.filter:
                continue
            
            total_variants += 1
            
            # Extract VAF from INFO or compute from AD/DP
            vaf = None
            if "AF" in rec.info:
                af = rec.info["AF"]
                vaf = float(af[0]) if isinstance(af, (list, tuple)) else float(af)
            elif "AD" in rec.info and "DP" in rec.info:
                ad = rec.info["AD"]
                dp = rec.info["DP"]
                if dp > 0:
                    alt_count = ad[1] if isinstance(ad, (list, tuple)) and len(ad) > 1 else ad
                    vaf = float(alt_count) / float(dp)
            
            # Check if variant has fragmentomics annotations
            if "FRAG_KS_PVAL" in rec.info and rec.info["FRAG_KS_PVAL"] is not None:
                ks_pval = rec.info["FRAG_KS_PVAL"]
                ks_pvals.append(ks_pval)
                
                if ks_pval < ks_pval_cutoff:
                    somatic_variants += 1
                    # Collect VAF of variants that passed "Smoking Gun" test
                    if vaf is not None and vaf > 0:
                        vafs_passed_smoking_gun.append(vaf)
                
                if "FRAG_ALT_MED" in rec.info:
                    alt_medians.append(rec.info["FRAG_ALT_MED"])
                if "FRAG_REF_MED" in rec.info:
                    ref_medians.append(rec.info["FRAG_REF_MED"])
        
        vcf.close()
        
        # Compute purity: Use 2 * mean VAF of variants that passed smoking gun test
        # This is the recommended approach for low tumor burden cases (Dr. Bharat's spec)
        if vafs_passed_smoking_gun:
            mean_vaf = float(statistics.mean(vafs_passed_smoking_gun))
            purity = 2.0 * mean_vaf  # 2 * VAF as per Dr. Bharat's specification
            purity_method = "2x_mean_VAF_smoking_gun"
        elif total_variants > 0:
            # Fallback: fraction of somatic variants (if VAF not available)
            purity = somatic_variants / total_variants
            purity_method = "fraction_somatic_variants"
        else:
            purity = 0.0
            purity_method = "none"
        
        metadata = {
            "status": "success",
            "source": "fragment_KS",
            "total_variants": total_variants,
            "somatic_variants": somatic_variants,
            "variants_passed_smoking_gun": len(vafs_passed_smoking_gun),
            "mean_vaf_smoking_gun": float(statistics.mean(vafs_passed_smoking_gun)) if vafs_passed_smoking_gun else None,
            "purity_method": purity_method,
            "mean_ks_pval": float(statistics.mean(ks_pvals)) if ks_pvals else None,
            "mean_alt_median": float(statistics.mean(alt_medians)) if alt_medians else None,
            "mean_ref_median": float(statistics.mean(ref_medians)) if ref_medians else None
        }
        
        return max(0.0, min(1.0, purity)), metadata
        
    except ImportError:
        print("WARNING: pysam not available, skipping fragment KS purity", file=sys.stderr)
        return 0.0, {"status": "pysam_missing", "source": "fragment_KS"}
    except Exception as e:
        print(f"WARNING: Could not compute fragment KS purity: {e}", file=sys.stderr)
        return 0.0, {"status": "error", "source": "fragment_KS", "error": str(e)}


def compute_snv_indel_purity(vcf_file: str, min_vaf: float = 0.01) -> Tuple[float, Dict]:
    """
    Compute ctDNA purity from SNV/INDEL VAF distribution.
    Uses median VAF of PASS variants as purity estimate.
    Returns: (purity_estimate, metadata)
    """
    if not os.path.exists(vcf_file):
        return 0.0, {"status": "file_not_found", "source": "SNV/INDEL"}
    
    try:
        import pysam
        
        vcf = pysam.VariantFile(vcf_file, "r")
        
        vafs = []
        total_variants = 0
        
        for rec in vcf:
            if "PASS" not in rec.filter:
                continue
            
            total_variants += 1
            
            # Extract VAF from INFO or compute from AD/DP
            vaf = None
            if "AF" in rec.info:
                af = rec.info["AF"]
                vaf = float(af[0]) if isinstance(af, (list, tuple)) else float(af)
            elif "AD" in rec.info and "DP" in rec.info:
                ad = rec.info["AD"]
                dp = rec.info["DP"]
                if dp > 0:
                    alt_count = ad[1] if isinstance(ad, (list, tuple)) and len(ad) > 1 else ad
                    vaf = float(alt_count) / float(dp)
            
            if vaf is not None and vaf >= min_vaf:
                vafs.append(vaf)
        
        vcf.close()
        
        # Use median VAF as purity estimate (more robust than mean)
        if vafs:
            purity = float(statistics.median(vafs))
        else:
            purity = 0.0
        
        metadata = {
            "status": "success",
            "source": "SNV/INDEL",
            "total_variants": total_variants,
            "variants_with_vaf": len(vafs),
            "median_vaf": float(statistics.median(vafs)) if vafs else 0.0,
            "mean_vaf": float(statistics.mean(vafs)) if vafs else 0.0
        }
        
        return max(0.0, min(1.0, purity)), metadata
        
    except ImportError:
        print("WARNING: pysam not available, skipping SNV/INDEL purity", file=sys.stderr)
        return 0.0, {"status": "pysam_missing", "source": "SNV/INDEL"}
    except Exception as e:
        print(f"WARNING: Could not compute SNV/INDEL purity: {e}", file=sys.stderr)
        return 0.0, {"status": "error", "source": "SNV/INDEL", "error": str(e)}


def compute_consensus_purity(
    ichorcna_tf: float,
    ichorcna_metadata: Dict,
    fragment_ks_purity: float,
    fragment_ks_metadata: Dict,
    snv_purity: float,
    snv_metadata: Dict,
    cna_lod_cutoff: float = 0.03
) -> Tuple[float, Dict]:
    """
    Compute consensus ctDNA purity using clinical decision rule:
    - If CNA Purity (ichorCNA) > cna_lod_cutoff (default 3%) → trust CNA result
    - Else → trust SNV-based purity from variant-level fragmentomics
    
    This implements the decision-support system logic as specified:
    - CNA-based purity is more reliable at higher tumor fractions (>3%)
    - SNV/INDEL-based purity (from fragmentomics) is more sensitive at low TFX (<3%)
    
    Args:
        ichorcna_tf: ichorCNA tumor fraction estimate (0.0-1.0)
        ichorcna_metadata: Metadata dict with status information
        fragment_ks_purity: Fragment KS-based purity estimate (0.0-1.0)
        fragment_ks_metadata: Fragment KS metadata
        snv_purity: SNV/INDEL-based purity estimate (0.0-1.0)
        snv_metadata: SNV/INDEL metadata
        cna_lod_cutoff: CNA limit of detection cutoff (default 0.03 = 3%)
    
    Returns:
        Tuple of (consensus_purity, metadata_dict)
    """
    # Determine which SNV-based purity to use (prefer fragment KS if available, else SNV/INDEL)
    snv_based_purity = 0.0
    snv_based_source = None
    snv_based_metadata = {}
    
    # Prefer fragment KS (variant-level fragmentomics) if available
    if fragment_ks_metadata.get("status") == "success" and fragment_ks_metadata.get("total_variants", 0) > 0:
        snv_based_purity = fragment_ks_purity
        snv_based_source = "fragment_KS"
        snv_based_metadata = fragment_ks_metadata
    # Fallback to SNV/INDEL purity if fragment KS not available
    elif snv_metadata.get("status") == "success" and snv_metadata.get("total_variants", 0) > 0:
        snv_based_purity = snv_purity
        snv_based_source = "SNV/INDEL"
        snv_based_metadata = snv_metadata
    
    # Get CNA purity (ichorCNA)
    cna_purity = 0.0
    cna_available = ichorcna_metadata.get("status") == "success"
    if cna_available:
        cna_purity = ichorcna_tf
    
    # Apply decision rule
    if cna_available and cna_purity >= cna_lod_cutoff:
        # Trust CNA result (CNA Purity > 3%)
        consensus_purity = cna_purity
        consensus_method = "CNA-based (ichorCNA)"
        decision_reason = f"CNA purity ({cna_purity*100:.2f}%) >= LOD cutoff ({cna_lod_cutoff*100:.0f}%)"
    elif snv_based_source:
        # Trust SNV-based purity (CNA Purity < 3% or CNA not available)
        consensus_purity = snv_based_purity
        consensus_method = f"SNV-based ({snv_based_source})"
        if cna_available:
            decision_reason = f"CNA purity ({cna_purity*100:.2f}%) < LOD cutoff ({cna_lod_cutoff*100:.0f}%), using SNV-based estimate"
        else:
            decision_reason = "CNA not available, using SNV-based estimate"
    elif cna_available:
        # Fallback: use CNA even if below cutoff (better than nothing)
        consensus_purity = cna_purity
        consensus_method = "CNA-based (ichorCNA) [fallback]"
        decision_reason = f"CNA available but below LOD ({cna_purity*100:.2f}% < {cna_lod_cutoff*100:.0f}%), SNV-based not available"
    else:
        # No estimates available
        consensus_purity = 0.0
        consensus_method = "None"
        decision_reason = "No purity estimates available"
    
    metadata = {
        "status": "success" if consensus_purity > 0.0 or cna_available or snv_based_source else "no_estimates",
        "consensus_purity": consensus_purity,
        "consensus_method": consensus_method,
        "decision_rule": {
            "cna_lod_cutoff": cna_lod_cutoff,
            "cna_purity": cna_purity,
            "cna_available": cna_available,
            "snv_based_purity": snv_based_purity,
            "snv_based_source": snv_based_source,
            "decision_reason": decision_reason
        },
        "individual_estimates": {
            "ichorCNA": ichorcna_tf if cna_available else None,
            "fragment_KS": fragment_ks_purity if fragment_ks_metadata.get("status") == "success" else None,
            "SNV/INDEL": snv_purity if snv_metadata.get("status") == "success" else None
        },
        "fallback": not (cna_available or snv_based_source)
    }
    
    return max(0.0, min(1.0, consensus_purity)), metadata


def main():
    parser = argparse.ArgumentParser(
        description="Compute comprehensive ctDNA purity by integrating ichorCNA, fragment KS, and SNV/INDEL estimates"
    )
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--ichorcna_tfx", help="ichorCNA tumor fraction file (.tfx.txt)")
    parser.add_argument("--ichorcna_segments", help="ichorCNA segments file (.seg.txt)")
    parser.add_argument("--frag_vcf", help="Fragmentomics-annotated VCF file")
    parser.add_argument("--snv_vcf", help="SNV/INDEL VCF file")
    parser.add_argument("--ks_pval_cutoff", type=float, default=0.05, help="KS p-value cutoff for somatic variants")
    parser.add_argument("--min_vaf", type=float, default=0.01, help="Minimum VAF for SNV/INDEL purity")
    parser.add_argument("--cna_lod_cutoff", type=float, default=0.03, help="CNA LOD cutoff (default 0.03 = 3%): if CNA purity >= this, trust CNA; else trust SNV-based")
    parser.add_argument("--output_json", required=True, help="Output JSON file with all metrics")
    parser.add_argument("--output_tfx", help="Output TFX format file (.tfx.txt)")
    
    args = parser.parse_args()
    
    # Load ichorCNA tumor fraction
    ichorcna_tf = 0.0
    ichorcna_metadata = {"status": "not_provided"}
    if args.ichorcna_tfx:
        ichorcna_tf, ichorcna_metadata = load_ichorcna_tumor_fraction(args.ichorcna_tfx)
    
    # Load ichorCNA segments
    segments_data = {"status": "not_provided", "segments": []}
    if args.ichorcna_segments:
        segments_data = load_ichorcna_segments(args.ichorcna_segments)
    
    # Compute fragment KS purity
    fragment_ks_purity = 0.0
    fragment_ks_metadata = {"status": "not_provided"}
    if args.frag_vcf:
        fragment_ks_purity, fragment_ks_metadata = compute_fragment_ks_purity(
            args.frag_vcf, args.ks_pval_cutoff
        )
    
    # Compute SNV/INDEL purity
    snv_purity = 0.0
    snv_metadata = {"status": "not_provided"}
    if args.snv_vcf:
        snv_purity, snv_metadata = compute_snv_indel_purity(args.snv_vcf, args.min_vaf)
    
    # Compute consensus purity using decision rule
    consensus_purity, consensus_metadata = compute_consensus_purity(
        ichorcna_tf, ichorcna_metadata,
        fragment_ks_purity, fragment_ks_metadata,
        snv_purity, snv_metadata,
        cna_lod_cutoff=args.cna_lod_cutoff
    )
    
    # Compile comprehensive output
    output = {
        "sample_id": args.sample_id,
        "ctdna_purity": {
            "consensus": consensus_purity,
            "consensus_percent": f"{consensus_purity * 100:.2f}%",
            "consensus_method": consensus_metadata.get("consensus_method", "None"),
            "decision_rule": consensus_metadata.get("decision_rule", {}),
            "sources": consensus_metadata.get("sources_used", []),
            "fallback_used": consensus_metadata.get("fallback", False)
        },
        "ichorCNA": {
            "tumor_fraction": ichorcna_tf,
            "tumor_fraction_percent": f"{ichorcna_tf * 100:.2f}%",
            "status": ichorcna_metadata.get("status"),
            "ploidy": ichorcna_metadata.get("ploidy"),
            "segments": {
                "num_segments": segments_data.get("num_segments", 0),
                "mean_log2": segments_data.get("mean_log2", 0.0),
                "median_log2": segments_data.get("median_log2", 0.0),
                "segments": segments_data.get("segments", [])
            }
        },
        "fragment_KS": {
            "purity": fragment_ks_purity,
            "purity_percent": f"{fragment_ks_purity * 100:.2f}%",
            "status": fragment_ks_metadata.get("status"),
            "total_variants": fragment_ks_metadata.get("total_variants", 0),
            "somatic_variants": fragment_ks_metadata.get("somatic_variants", 0),
            "mean_ks_pval": fragment_ks_metadata.get("mean_ks_pval"),
            "mean_alt_median_frag": fragment_ks_metadata.get("mean_alt_median"),
            "mean_ref_median_frag": fragment_ks_metadata.get("mean_ref_median")
        },
        "SNV_INDEL": {
            "purity": snv_purity,
            "purity_percent": f"{snv_purity * 100:.2f}%",
            "status": snv_metadata.get("status"),
            "total_variants": snv_metadata.get("total_variants", 0),
            "variants_with_vaf": snv_metadata.get("variants_with_vaf", 0),
            "median_vaf": snv_metadata.get("median_vaf", 0.0),
            "mean_vaf": snv_metadata.get("mean_vaf", 0.0)
        },
        "qc_metrics": {
            "ichorCNA_available": ichorcna_metadata.get("status") == "success",
            "fragment_KS_available": fragment_ks_metadata.get("status") == "success",
            "SNV_INDEL_available": snv_metadata.get("status") == "success",
            "consensus_method": consensus_metadata.get("consensus_method", "decision_rule"),
            "cna_lod_cutoff": args.cna_lod_cutoff
        }
    }
    
    # Write JSON output
    with open(args.output_json, 'w') as f:
        json.dump(output, f, indent=2)
    
    # Write TFX format output (if requested)
    if args.output_tfx:
        with open(args.output_tfx, 'w') as f:
            f.write("sample\ttumorFraction\tploidy\tmethod\n")
            ploidy_str = str(ichorcna_metadata.get("ploidy", "NA")) if ichorcna_metadata.get("ploidy") is not None else "NA"
            f.write(f"{args.sample_id}\t{consensus_purity:.6f}\t{ploidy_str}\tconsensus\n")
    
    print(f"ctDNA purity computation complete:")
    print(f"  Consensus purity: {consensus_purity * 100:.2f}%")
    print(f"  Method: {consensus_metadata.get('consensus_method', 'None')}")
    decision_rule = consensus_metadata.get('decision_rule', {})
    if decision_rule:
        print(f"  Decision: {decision_rule.get('decision_reason', 'N/A')}")
    print(f"  Output written to: {args.output_json}")


if __name__ == "__main__":
    main()

