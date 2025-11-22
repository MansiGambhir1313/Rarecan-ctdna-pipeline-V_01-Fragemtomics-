#!/usr/bin/env python3
"""
Generate Final Comprehensive Fragmentomics Report
Combines all fragmentomics outputs into a single comprehensive report
"""

import json
import os
import sys
import pandas as pd
from datetime import datetime

def load_fragmentomics_data(sample_id, results_dir):
    """Load all fragmentomics data files"""
    data = {
        "sample_id": sample_id,
        "analysis_date": datetime.now().isoformat(),
        "files_loaded": [],
        "files_missing": []
    }
    
    # Load summary JSON
    summary_path = f"{results_dir}/fragmentomics/{sample_id}.fragmentomics.summary.json"
    if os.path.exists(summary_path):
        try:
            with open(summary_path, 'r') as f:
                data["summary"] = json.load(f)
            data["files_loaded"].append("summary")
        except Exception as e:
            data["files_missing"].append(f"summary (error: {e})")
    else:
        data["files_missing"].append("summary")
    
    # Load global histogram
    hist_path = f"{results_dir}/fragmentomics/{sample_id}.global_hist.tsv"
    if os.path.exists(hist_path):
        try:
            data["histogram"] = pd.read_csv(hist_path, sep="\t").to_dict('records')
            data["files_loaded"].append("histogram")
        except Exception as e:
            data["files_missing"].append(f"histogram (error: {e})")
    else:
        data["files_missing"].append("histogram")
    
    # Load fragment motifs
    motifs_path = f"{results_dir}/fragmentomics/{sample_id}.fragment_motifs.tsv"
    if os.path.exists(motifs_path):
        try:
            data["motifs"] = pd.read_csv(motifs_path, sep="\t").to_dict('records')
            data["files_loaded"].append("motifs")
        except Exception as e:
            data["files_missing"].append(f"motifs (error: {e})")
    else:
        data["files_missing"].append("motifs")
    
    return data

def generate_comprehensive_report(fragmentomics_data, output_path):
    """Generate comprehensive fragmentomics report"""
    report = {
        "report_type": "Comprehensive Fragmentomics Analysis Report",
        "sample_id": fragmentomics_data["sample_id"],
        "analysis_date": fragmentomics_data["analysis_date"],
        "data_availability": {
            "files_loaded": fragmentomics_data["files_loaded"],
            "files_missing": fragmentomics_data["files_missing"],
            "completeness": len(fragmentomics_data["files_loaded"]) / 3.0 * 100
        }
    }
    
    # Add summary data
    if "summary" in fragmentomics_data:
        summary = fragmentomics_data["summary"]
        report["fragment_metrics"] = {
            "total_fragments": summary.get("total_fragments", 0),
            "short_fragments_90_150bp": summary.get("short_fragments_90_150bp", 0),
            "short_fragment_percentage": summary.get("short_fragment_percentage", 0.0),
            "mean_fragment_length": summary.get("mean_fragment_length", 0.0),
            "median_fragment_length": summary.get("median_fragment_length", 0.0),
            "fragment_length_distribution": summary.get("fragment_length_distribution", {})
        }
        report["analysis_status"] = summary.get("status", "unknown")
        report["notes"] = summary.get("note", "")
    else:
        report["fragment_metrics"] = {
            "status": "No data available",
            "note": "Fragmentomics summary file not found or could not be loaded"
        }
    
    # Add histogram data
    if "histogram" in fragmentomics_data:
        report["histogram_data"] = fragmentomics_data["histogram"]
    else:
        report["histogram_data"] = []
    
    # Add motif data
    if "motifs" in fragmentomics_data:
        report["motif_analysis"] = fragmentomics_data["motifs"]
    else:
        report["motif_analysis"] = []
    
    # Add interpretation
    if "summary" in fragmentomics_data:
        short_pct = fragmentomics_data["summary"].get("short_fragment_percentage", 0.0)
        total_frags = fragmentomics_data["summary"].get("total_fragments", 0)
        
        report["interpretation"] = {
            "data_quality": "High" if total_frags > 1000000 else "Moderate" if total_frags > 100000 else "Low",
            "short_fragment_enrichment": "High" if short_pct > 20 else "Moderate" if short_pct > 10 else "Low",
            "clinical_relevance": "Short fragment enrichment may indicate ctDNA presence" if short_pct > 10 else "Limited short fragment enrichment observed"
        }
    else:
        report["interpretation"] = {
            "status": "Cannot interpret - no data available"
        }
    
    # Write report
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Print summary
    print("="*60)
    print("Final Fragmentomics Report Generated")
    print("="*60)
    print(f"Sample ID: {report['sample_id']}")
    print(f"Analysis Date: {report['analysis_date']}")
    print(f"Data Completeness: {report['data_availability']['completeness']:.1f}%")
    print(f"Files Loaded: {', '.join(report['data_availability']['files_loaded'])}")
    if report['data_availability']['files_missing']:
        print(f"Files Missing: {', '.join(report['data_availability']['files_missing'])}")
    print("="*60)
    
    return report

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate final comprehensive fragmentomics report")
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--results_dir", required=True, help="Results directory path")
    parser.add_argument("--output", required=True, help="Output report JSON path")
    
    args = parser.parse_args()
    
    fragmentomics_data = load_fragmentomics_data(args.sample_id, args.results_dir)
    report = generate_comprehensive_report(fragmentomics_data, args.output)
    
    print(f"\nReport saved to: {args.output}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

