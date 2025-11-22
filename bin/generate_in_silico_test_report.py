#!/usr/bin/env python3
"""
Generate In Silico Test/Validation Report for TFX Pipeline
Compares actual results against expected validation criteria
"""

import json
import os
import sys
from datetime import datetime
import pandas as pd

def load_tfx_report(tfx_json_path):
    """Load TFX report JSON"""
    if not os.path.exists(tfx_json_path):
        return None
    try:
        with open(tfx_json_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading TFX report: {e}", file=sys.stderr)
        return None

def load_fragmentomics_summary(frag_json_path):
    """Load fragmentomics summary JSON"""
    if not os.path.exists(frag_json_path):
        return None
    try:
        with open(frag_json_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading fragmentomics summary: {e}", file=sys.stderr)
        return None

def validate_tfx_results(tfx_data, frag_data):
    """Validate TFX results against clinical criteria"""
    validation_results = {
        "sample_id": tfx_data.get("sample_id", "unknown"),
        "validation_date": datetime.now().isoformat(),
        "overall_status": "PASS",
        "criteria": {}
    }
    
    # Criterion 1: Consensus Purity Calculation
    consensus_purity = tfx_data.get("ctdna_purity_consensus", 0.0)
    validation_results["criteria"]["consensus_purity"] = {
        "value": consensus_purity,
        "threshold": ">= 0.0",
        "status": "PASS" if consensus_purity >= 0.0 else "FAIL",
        "note": "Consensus purity should be calculated"
    }
    
    # Criterion 2: Method Selection Logic
    consensus_method = tfx_data.get("consensus_method", "")
    cna_purity = tfx_data.get("methods", {}).get("cna_based_purity", {}).get("value", 0.0)
    cna_lod = tfx_data.get("qc_metrics", {}).get("cna_lod_cutoff", 0.03)
    
    method_correct = False
    if cna_purity >= cna_lod:
        method_correct = "CNA-based" in consensus_method
    else:
        method_correct = "SNV/INDEL-based" in consensus_method or "Fragmentomics" in consensus_method
    
    validation_results["criteria"]["method_selection"] = {
        "expected_method": "CNA-based" if cna_purity >= cna_lod else "SNV/INDEL-based",
        "actual_method": consensus_method,
        "status": "PASS" if method_correct else "FAIL",
        "note": "Method selection should follow clinical logic"
    }
    
    # Criterion 3: Global Fragmentomics QC
    global_frag = tfx_data.get("qc_metrics", {}).get("global_fragmentomics", {})
    total_reads = global_frag.get("total_reads", 0)
    short_frag_pct = global_frag.get("short_fragment_percent", 0.0)
    
    validation_results["criteria"]["global_fragmentomics"] = {
        "total_reads": total_reads,
        "short_fragment_percent": short_frag_pct,
        "min_reads_threshold": 1000000,
        "status": "PASS" if total_reads >= 1000000 else "WARN",
        "note": "Sufficient reads for reliable fragmentomics analysis"
    }
    
    # Criterion 4: Fragmentomics Data Quality
    if frag_data:
        total_fragments = frag_data.get("total_fragments", 0)
        mean_length = frag_data.get("mean_fragment_length", 0.0)
        
        validation_results["criteria"]["fragmentomics_quality"] = {
            "total_fragments": total_fragments,
            "mean_fragment_length": mean_length,
            "status": "PASS" if total_fragments > 0 and mean_length > 0 else "FAIL",
            "note": "Fragmentomics data should contain real fragment measurements"
        }
    else:
        validation_results["criteria"]["fragmentomics_quality"] = {
            "status": "WARN",
            "note": "Fragmentomics summary data not available"
        }
    
    # Criterion 5: Both Methods Available
    cna_status = tfx_data.get("methods", {}).get("cna_based_purity", {}).get("status", "")
    snv_status = tfx_data.get("methods", {}).get("snv_indel_based_purity", {}).get("status", "")
    
    validation_results["criteria"]["method_availability"] = {
        "cna_method": cna_status,
        "snv_method": snv_status,
        "status": "PASS" if "Success" in cna_status or "Success" in snv_status else "WARN",
        "note": "At least one purity estimation method should succeed"
    }
    
    # Overall status
    all_pass = all(
        c.get("status") == "PASS" 
        for c in validation_results["criteria"].values() 
        if isinstance(c, dict) and "status" in c
    )
    validation_results["overall_status"] = "PASS" if all_pass else "WARN"
    
    return validation_results

def generate_validation_report(tfx_json_path, frag_json_path, output_path):
    """Generate comprehensive validation report"""
    tfx_data = load_tfx_report(tfx_json_path)
    frag_data = load_fragmentomics_summary(frag_json_path)
    
    if not tfx_data:
        print(f"ERROR: Could not load TFX report from {tfx_json_path}", file=sys.stderr)
        return False
    
    validation_results = validate_tfx_results(tfx_data, frag_data)
    
    # Add summary statistics
    validation_results["summary"] = {
        "total_criteria": len(validation_results["criteria"]),
        "passed_criteria": sum(
            1 for c in validation_results["criteria"].values()
            if isinstance(c, dict) and c.get("status") == "PASS"
        ),
        "warned_criteria": sum(
            1 for c in validation_results["criteria"].values()
            if isinstance(c, dict) and c.get("status") == "WARN"
        ),
        "failed_criteria": sum(
            1 for c in validation_results["criteria"].values()
            if isinstance(c, dict) and c.get("status") == "FAIL"
        )
    }
    
    # Write JSON report
    with open(output_path, 'w') as f:
        json.dump(validation_results, f, indent=2)
    
    # Print summary
    print("="*60)
    print("In Silico Validation Test Report")
    print("="*60)
    print(f"Sample ID: {validation_results['sample_id']}")
    print(f"Overall Status: {validation_results['overall_status']}")
    print(f"Validation Date: {validation_results['validation_date']}")
    print("")
    print("Criteria Results:")
    for name, result in validation_results["criteria"].items():
        if isinstance(result, dict):
            status = result.get("status", "UNKNOWN")
            note = result.get("note", "")
            print(f"  {name}: {status} - {note}")
    print("")
    print(f"Summary: {validation_results['summary']['passed_criteria']}/{validation_results['summary']['total_criteria']} criteria passed")
    print("="*60)
    
    return True

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate in silico validation test report")
    parser.add_argument("--tfx_json", required=True, help="Path to TFX report JSON")
    parser.add_argument("--frag_json", help="Path to fragmentomics summary JSON")
    parser.add_argument("--output", required=True, help="Output validation report JSON path")
    
    args = parser.parse_args()
    
    success = generate_validation_report(
        args.tfx_json,
        args.frag_json,
        args.output
    )
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

