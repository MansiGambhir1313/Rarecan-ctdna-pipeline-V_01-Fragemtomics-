// modules/report.nf - Clinical Report Generation Module

nextflow.enable.dsl=2

process COLLECT_RESULTS {
    tag "$sample_id"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(clinical_variants), path(cnv_calls), path(sv_calls), path(msi_result), path(tmb_result), path(qc_metrics)
    
    output:
    tuple val(sample_id), path("${sample_id}_collected_results.json"), emit: collected_results
    path "versions.yml", emit: versions
    
    script:
    """
    cat > collect_results.py << 'EOF'
import json
import csv
import os
from datetime import datetime

def load_json_file(filepath):
    \"\"\"Load JSON file if it exists\"\"\"
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        try:
            with open(filepath, 'r') as f:
                return json.load(f)
        except:
            return {}
    return {}

def load_tsv_file(filepath):
    \"\"\"Load TSV file as list of dictionaries\"\"\"
    if not os.path.exists(filepath):
        return []
    
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            return list(reader)
    except:
        return []

def parse_qc_metrics(qc_file):
    \"\"\"Parse QC metrics from various formats\"\"\"
    qc_data = {}
    
    if not os.path.exists(qc_file):
        return qc_data
    
    # Try to parse as JSON first
    try:
        with open(qc_file, 'r') as f:
            qc_data = json.load(f)
    except:
        # Try to parse as text file (Picard metrics format)
        try:
            with open(qc_file, 'r') as f:
                lines = f.readlines()
                
            # Look for metrics in Picard format
            for i, line in enumerate(lines):
                if line.startswith('BAIT_SET') or line.startswith('SAMPLE'):
                    # Found header line, next line should have values
                    if i + 1 < len(lines):
                        headers = line.strip().split('\\t')
                        values = lines[i + 1].strip().split('\\t')
                        
                        if len(headers) == len(values):
                            for h, v in zip(headers, values):
                                try:
                                    qc_data[h] = float(v) if '.' in v else int(v)
                                except:
                                    qc_data[h] = v
                        break
        except:
            pass
    
    return qc_data

# Load all result files
clinical_variants = load_tsv_file("${clinical_variants}")
cnv_calls = load_tsv_file("${cnv_calls}")
sv_calls = load_tsv_file("${sv_calls}")
msi_result = load_json_file("${msi_result}")
tmb_result = load_json_file("${tmb_result}")
qc_metrics = parse_qc_metrics("${qc_metrics}")

# Get current timestamp
current_time = datetime.now().isoformat()

# Collect comprehensive results
collected_results = {
    "sample_id": "${sample_id}",
    "analysis_date": current_time,
    "pipeline_version": "1.0.0",
    "genome_build": "GRCh38",
    
    "quality_control": {
        "on_target_depth": qc_metrics.get("MEAN_TARGET_COVERAGE", 0),
        "on_target_rate": qc_metrics.get("PCT_SELECTED_BASES", 0),
        "coverage_uniformity": qc_metrics.get("PCT_TARGET_BASES_20X", 0) / 100 if qc_metrics.get("PCT_TARGET_BASES_20X") else 0,
        "total_reads": qc_metrics.get("TOTAL_READS", 0),
        "pf_reads": qc_metrics.get("PF_READS", 0),
        "pf_unique_reads": qc_metrics.get("PF_UNIQUE_READS", 0),
        "mean_insert_size": qc_metrics.get("MEAN_INSERT_SIZE", 0),
        "contamination": 0.0  # Would come from contamination analysis
    },
    
    "small_variants": {
        "total_variants": len(clinical_variants),
        "pathogenic_variants": len([v for v in clinical_variants if 'pathogenic' in v.get('clinical_significance', '').lower()]),
        "vus_variants": len([v for v in clinical_variants if 'vus' in v.get('clinical_significance', '').lower()]),
        "variants": clinical_variants[:50]  # Limit to top 50 for report
    },
    
    "copy_number_variants": {
        "total_cnvs": len(cnv_calls),
        "amplifications": len([c for c in cnv_calls if c.get('interpretation', '').lower() in ['amplification', 'gain']]),
        "deletions": len([c for c in cnv_calls if c.get('interpretation', '').lower() in ['deletion', 'deep deletion']]),
        "cnvs": cnv_calls
    },
    
    "structural_variants": {
        "total_svs": len(sv_calls),
        "fusion_candidates": len([s for s in sv_calls if 'fusion' in s.get('clinical_significance', '').lower()]),
        "svs": sv_calls
    },
    
    "biomarkers": {
        "msi": msi_result,
        "tmb": tmb_result
    },
    
    "summary": {
        "actionable_variants": len([v for v in clinical_variants if v.get('evidence_level', '5') in ['1', '2']]),
        "total_alterations": len(clinical_variants) + len(cnv_calls) + len(sv_calls),
        "analysis_quality": "High" if qc_metrics.get("MEAN_TARGET_COVERAGE", 0) >= 300 else "Moderate" if qc_metrics.get("MEAN_TARGET_COVERAGE", 0) >= 100 else "Low"
    }
}

# Write collected results
with open("${sample_id}_collected_results.json", 'w') as f:
    json.dump(collected_results, f, indent=2)

print(f"Results Collection Complete:")
print(f"  Sample: ${sample_id}")
print(f"  Small Variants: {len(clinical_variants)}")
print(f"  CNVs: {len(cnv_calls)}")
print(f"  SVs: {len(sv_calls)}")
print(f"  MSI Status: {msi_result.get('msi_status', 'Unknown')}")
print(f"  TMB: {tmb_result.get('tmb_per_mb', 0)} mut/Mb")
EOF

    python collect_results.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

process GENERATE_CLINICAL_REPORT {
    tag "$sample_id"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(collected_results)
    
    output:
    tuple val(sample_id), path("${sample_id}_clinical_report.json"), emit: json_report
    tuple val(sample_id), path("${sample_id}_clinical_report.html"), emit: html_report
    path "versions.yml", emit: versions
    
    script:
    """
    # Install required packages
    pip install jinja2 --quiet
    
    cat > generate_report.py << 'EOF'
import json
import os
from datetime import datetime
from jinja2 import Template

# Load collected results
with open("${collected_results}", 'r') as f:
    data = json.load(f)

# Create clinical report JSON (structured for clinical systems)
clinical_report = {
    "report_header": {
        "sample_id": data["sample_id"],
        "report_date": datetime.now().isoformat(),
        "pipeline_version": data["pipeline_version"],
        "genome_build": data["genome_build"],
        "laboratory": "Clinical Bioinformatics Laboratory",
        "report_type": "ctDNA Analysis Report"
    },
    
    "quality_metrics": {
        "sequencing_quality": {
            "mean_target_coverage": data["quality_control"]["on_target_depth"],
            "on_target_percentage": data["quality_control"]["on_target_rate"] * 100,
            "coverage_uniformity": data["quality_control"]["coverage_uniformity"] * 100,
            "total_reads": data["quality_control"]["total_reads"]
        },
        "quality_flags": {
            "sufficient_coverage": data["quality_control"]["on_target_depth"] >= 300,
            "adequate_on_target": data["quality_control"]["on_target_rate"] >= 0.7,
            "uniform_coverage": data["quality_control"]["coverage_uniformity"] >= 0.8,
            "overall_quality": data["summary"]["analysis_quality"]
        }
    },
    
    "genomic_findings": {
        "small_variants": [
            {
                "gene": v["gene"],
                "variant_type": "SNV/INDEL",
                "genomic_change": f"{v['chromosome']}:g.{v['position']}{v['ref']}>{v['alt']}",
                "protein_change": v["hgvs_protein"],
                "vaf": v["vaf"],
                "clinical_significance": v["clinical_significance"],
                "evidence_level": v["evidence_level"]
            }
            for v in data["small_variants"]["variants"]
            if v.get("evidence_level", "5") in ["1", "2", "3"]  # Only clinically relevant
        ],
        
        "copy_number_variants": [
            {
                "gene": c["gene"],
                "variant_type": "CNV",
                "alteration": c["interpretation"],
                "copy_number": c.get("copy_number", "Unknown"),
                "clinical_significance": "Significant" if abs(float(c.get("log2", 0))) > 1.0 else "Moderate"
            }
            for c in data["copy_number_variants"]["cnvs"]
        ],
        
        "structural_variants": [
            {
                "variant_type": "SV",
                "sv_type": s["sv_type"],
                "genomic_location": f"{s['chrom']}:{s['pos']}",
                "clinical_significance": s["clinical_significance"]
            }
            for s in data["structural_variants"]["svs"]
        ]
    },
    
    "biomarker_analysis": {
        "microsatellite_instability": {
            "msi_status": data["biomarkers"]["msi"].get("msi_status", "Unknown"),
            "msi_score": data["biomarkers"]["msi"].get("msi_score", 0),
            "interpretation": data["biomarkers"]["msi"].get("interpretation", "")
        },
        
        "tumor_mutational_burden": {
            "tmb_score": data["biomarkers"]["tmb"].get("tmb_per_mb", 0),
            "tmb_interpretation": data["biomarkers"]["tmb"].get("tmb_interpretation", "Unknown"),
            "confidence_interval": data["biomarkers"]["tmb"].get("confidence_interval_95", [0, 0]),
            "immunotherapy_relevant": data["biomarkers"]["tmb"].get("tmb_per_mb", 0) >= 10
        }
    },
    
    "clinical_summary": {
        "total_alterations": data["summary"]["total_alterations"],
        "actionable_findings": data["summary"]["actionable_variants"],
        "key_findings": [],  # Would be populated based on clinical rules
        "recommendations": [],  # Would be populated based on findings
        "limitations": "This analysis is limited to the targeted gene panel. Variants outside the panel regions are not detected. Clinical correlation is recommended."
    }
}

# Add key findings based on results
key_findings = []
if clinical_report["biomarker_analysis"]["microsatellite_instability"]["msi_status"] == "MSI-H":
    key_findings.append("High microsatellite instability detected - consider immunotherapy")

if clinical_report["biomarker_analysis"]["tumor_mutational_burden"]["tmb_score"] >= 10:
    key_findings.append(f"High tumor mutational burden ({clinical_report['biomarker_analysis']['tumor_mutational_burden']['tmb_score']} mut/Mb) - consider immunotherapy")

pathogenic_variants = [v for v in clinical_report["genomic_findings"]["small_variants"] if "pathogenic" in v["clinical_significance"].lower()]
if pathogenic_variants:
    key_findings.append(f"{len(pathogenic_variants)} pathogenic/likely pathogenic variant(s) identified")

clinical_report["clinical_summary"]["key_findings"] = key_findings

# Write JSON report
with open("${sample_id}_clinical_report.json", 'w') as f:
    json.dump(clinical_report, f, indent=2)

# Generate HTML report
html_template = '''
<!DOCTYPE html>
<html>
<head>
    <title>ctDNA Analysis Report - {{ report.report_header.sample_id }}</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        .quality-pass { color: green; font-weight: bold; }
        .quality-fail { color: red; font-weight: bold; }
        .variant-table { border-collapse: collapse; width: 100%; }
        .variant-table th, .variant-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .variant-table th { background-color: #f2f2f2; }
        .pathogenic { background-color: #ffebee; }
        .vus { background-color: #fff3e0; }
        .biomarker { background-color: #e8f5e8; padding: 10px; margin: 10px 0; border-radius: 5px; }
    </style>
</head>
<body>
    <div class="header">
        <h1>ctDNA Analysis Report</h1>
        <p><strong>Sample ID:</strong> {{ report.report_header.sample_id }}</p>
        <p><strong>Report Date:</strong> {{ report.report_header.report_date }}</p>
        <p><strong>Pipeline Version:</strong> {{ report.report_header.pipeline_version }}</p>
        <p><strong>Genome Build:</strong> {{ report.report_header.genome_build }}</p>
    </div>
    
    <div class="section">
        <h2>Quality Metrics</h2>
        <p><strong>Mean Target Coverage:</strong> 
            <span class="{{ 'quality-pass' if report.quality_metrics.quality_flags.sufficient_coverage else 'quality-fail' }}">
                {{ "%.0f"|format(report.quality_metrics.sequencing_quality.mean_target_coverage) }}x
            </span>
        </p>
        <p><strong>On-Target Rate:</strong> 
            <span class="{{ 'quality-pass' if report.quality_metrics.quality_flags.adequate_on_target else 'quality-fail' }}">
                {{ "%.1f"|format(report.quality_metrics.sequencing_quality.on_target_percentage) }}%
            </span>
        </p>
        <p><strong>Coverage Uniformity:</strong> 
            <span class="{{ 'quality-pass' if report.quality_metrics.quality_flags.uniform_coverage else 'quality-fail' }}">
                {{ "%.1f"|format(report.quality_metrics.sequencing_quality.coverage_uniformity) }}%
            </span>
        </p>
        <p><strong>Overall Quality:</strong> {{ report.quality_metrics.quality_flags.overall_quality }}</p>
    </div>
    
    <div class="section">
        <h2>Key Findings</h2>
        {% if report.clinical_summary.key_findings %}
            <ul>
            {% for finding in report.clinical_summary.key_findings %}
                <li>{{ finding }}</li>
            {% endfor %}
            </ul>
        {% else %}
            <p>No significant findings identified.</p>
        {% endif %}
    </div>
    
    <div class="section">
        <h2>Biomarker Analysis</h2>
        <div class="biomarker">
            <h3>Microsatellite Instability (MSI)</h3>
            <p><strong>Status:</strong> {{ report.biomarker_analysis.microsatellite_instability.msi_status }}</p>
            <p><strong>Score:</strong> {{ "%.3f"|format(report.biomarker_analysis.microsatellite_instability.msi_score) }}</p>
        </div>
        
        <div class="biomarker">
            <h3>Tumor Mutational Burden (TMB)</h3>
            <p><strong>TMB Score:</strong> {{ "%.1f"|format(report.biomarker_analysis.tumor_mutational_burden.tmb_score) }} mutations/Mb</p>
            <p><strong>Interpretation:</strong> {{ report.biomarker_analysis.tumor_mutational_burden.tmb_interpretation }}</p>
            {% if report.biomarker_analysis.tumor_mutational_burden.immunotherapy_relevant %}
                <p style="color: green; font-weight: bold;">May be relevant for immunotherapy consideration</p>
            {% endif %}
        </div>
    </div>
    
    {% if report.genomic_findings.small_variants %}
    <div class="section">
        <h2>Small Variants (SNVs/INDELs)</h2>
        <table class="variant-table">
            <tr>
                <th>Gene</th>
                <th>Genomic Change</th>
                <th>Protein Change</th>
                <th>VAF</th>
                <th>Clinical Significance</th>
            </tr>
            {% for variant in report.genomic_findings.small_variants %}
            <tr class="{{ 'pathogenic' if 'pathogenic' in variant.clinical_significance.lower() else 'vus' if 'vus' in variant.clinical_significance.lower() else '' }}">
                <td>{{ variant.gene }}</td>
                <td>{{ variant.genomic_change }}</td>
                <td>{{ variant.protein_change }}</td>
                <td>{{ "%.3f"|format(variant.vaf) }}</td>
                <td>{{ variant.clinical_significance }}</td>
            </tr>
            {% endfor %}
        </table>
    </div>
    {% endif %}
    
    {% if report.genomic_findings.copy_number_variants %}
    <div class="section">
        <h2>Copy Number Variants</h2>
        <table class="variant-table">
            <tr>
                <th>Gene</th>
                <th>Alteration</th>
                <th>Copy Number</th>
                <th>Clinical Significance</th>
            </tr>
            {% for cnv in report.genomic_findings.copy_number_variants %}
            <tr>
                <td>{{ cnv.gene }}</td>
                <td>{{ cnv.alteration }}</td>
                <td>{{ cnv.copy_number }}</td>
                <td>{{ cnv.clinical_significance }}</td>
            </tr>
            {% endfor %}
        </table>
    </div>
    {% endif %}
    
    <div class="section">
        <h2>Limitations</h2>
        <p>{{ report.clinical_summary.limitations }}</p>
    </div>
    
    <div class="section">
        <p><em>This report was generated automatically by the ctDNA analysis pipeline. Clinical correlation and interpretation by a qualified healthcare professional is recommended.</em></p>
    </div>
</body>
</html>
'''

template = Template(html_template)
html_content = template.render(report=clinical_report)

with open("${sample_id}_clinical_report.html", 'w') as f:
    f.write(html_content)

print(f"Clinical Report Generated:")
print(f"  Sample: {clinical_report['report_header']['sample_id']}")
print(f"  Total Alterations: {clinical_report['clinical_summary']['total_alterations']}")
print(f"  Actionable Findings: {clinical_report['clinical_summary']['actionable_findings']}")
print(f"  MSI Status: {clinical_report['biomarker_analysis']['microsatellite_instability']['msi_status']}")
print(f"  TMB Score: {clinical_report['biomarker_analysis']['tumor_mutational_burden']['tmb_score']} mut/Mb")
EOF

    python generate_report.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        jinja2: \$(python -c "import jinja2; print(jinja2.__version__)")
    END_VERSIONS
    """
}

// Report Generation Workflow
workflow REPORT {
    take:
    clinical_variants    // channel: [sample_id, clinical_variants_tsv]
    cnv_calls           // channel: [sample_id, cnv_calls_tsv]
    sv_calls            // channel: [sample_id, sv_calls_tsv]
    msi_result          // channel: [sample_id, msi_result_json]
    tmb_result          // channel: [sample_id, tmb_result_json]
    qc_metrics          // channel: [sample_id, qc_metrics_file]
    
    main:
    // Combine all results by sample_id
    combined_results = clinical_variants
        .join(cnv_calls, by: 0)
        .join(sv_calls, by: 0)
        .join(msi_result, by: 0)
        .join(tmb_result, by: 0)
        .join(qc_metrics, by: 0)
    
    // Collect all results into a single JSON
    COLLECT_RESULTS(combined_results)
    
    // Generate clinical report
    GENERATE_CLINICAL_REPORT(COLLECT_RESULTS.out.collected_results)
    
    // Combine versions
    ch_versions = COLLECT_RESULTS.out.versions
        .mix(GENERATE_CLINICAL_REPORT.out.versions)
    
    emit:
    json_report = GENERATE_CLINICAL_REPORT.out.json_report
    html_report = GENERATE_CLINICAL_REPORT.out.html_report
    collected_results = COLLECT_RESULTS.out.collected_results
    versions = ch_versions
}