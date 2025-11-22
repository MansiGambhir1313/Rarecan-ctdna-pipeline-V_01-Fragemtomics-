// modules/tmb.nf - Tumor Mutational Burden Analysis Module
// ECR-compatible version using deployed containers

nextflow.enable.dsl=2

process CALCULATE_TMB {
    tag "$sample_id"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(vcf)
    path targets_bed
    path filter_stats
    
    output:
    tuple val(sample_id), path("${sample_id}_tmb_calculation.json"), emit: calculation
    path "versions.yml", emit: versions
    
    script:
    """
    cat > calculate_tmb.py << 'EOF'
import gzip
import json
import os
from datetime import datetime

def parse_bed_file(bed_file):
    \"\"\"Calculate total targetable regions from BED file\"\"\"
    total_bases = 0
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\\t')
                if len(fields) >= 3:
                    start = int(fields[1])
                    end = int(fields[2])
                    total_bases += (end - start)
        
        return total_bases
    except Exception as e:
        print(f"Error parsing BED file: {e}")
        return 0

def count_coding_mutations(vcf_file):
    \"\"\"Count coding mutations from VCF file\"\"\"
    coding_mutations = 0
    total_variants = 0
    synonymous_mutations = 0
    
    try:
        # Handle both gzipped and regular VCF files
        if vcf_file.endswith('.gz'):
            file_handle = gzip.open(vcf_file, 'rt')
        else:
            file_handle = open(vcf_file, 'r')
        
        with file_handle as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                total_variants += 1
                
                # Parse VCF line
                fields = line.strip().split('\\t')
                if len(fields) < 8:
                    continue
                
                # Check INFO field for coding consequences
                info_field = fields[7]
                
                # Look for coding consequences (simplified approach)
                coding_consequences = [
                    'missense_variant', 'stop_gained', 'stop_lost',
                    'start_lost', 'frameshift_variant', 'inframe_insertion',
                    'inframe_deletion', 'protein_altering_variant'
                ]
                
                synonymous_consequences = ['synonymous_variant']
                
                # Check if variant has coding consequence
                is_coding = any(cons in info_field.lower() for cons in coding_consequences)
                is_synonymous = any(cons in info_field.lower() for cons in synonymous_consequences)
                
                if is_coding:
                    coding_mutations += 1
                elif is_synonymous:
                    synonymous_mutations += 1
                elif 'CSQ=' in info_field:
                    # If CSQ annotation is present, check for coding consequences
                    csq_part = info_field.split('CSQ=')[1].split(';')[0]
                    if any(cons in csq_part.lower() for cons in coding_consequences):
                        coding_mutations += 1
                    elif any(cons in csq_part.lower() for cons in synonymous_consequences):
                        synonymous_mutations += 1
    
    except Exception as e:
        print(f"Error processing VCF file: {e}")
    
    return {
        'total_variants': total_variants,
        'coding_mutations': coding_mutations,
        'synonymous_mutations': synonymous_mutations
    }

def load_filter_stats(stats_file):
    \"\"\"Load filtering statistics if available\"\"\"
    try:
        if os.path.exists(stats_file):
            with open(stats_file, 'r') as f:
                return json.load(f)
    except:
        pass
    
    return {}

# Main TMB calculation
print("Starting TMB calculation...")

# Parse target regions
target_bases = parse_bed_file("${targets_bed}")
target_mb = target_bases / 1_000_000  # Convert to megabases

print(f"Target region size: {target_bases:,} bases ({target_mb:.2f} Mb)")

# Count mutations
mutation_counts = count_coding_mutations("${vcf}")
coding_mutations = mutation_counts['coding_mutations']
total_variants = mutation_counts['total_variants']
synonymous_mutations = mutation_counts['synonymous_mutations']

print(f"Total variants: {total_variants}")
print(f"Coding mutations: {coding_mutations}")
print(f"Synonymous mutations: {synonymous_mutations}")

# Calculate TMB
if target_mb > 0:
    tmb_per_mb = coding_mutations / target_mb
    tmb_total = coding_mutations / target_mb  # Legacy compatibility
else:
    tmb_per_mb = 0
    tmb_total = 0

# Load filter statistics
filter_stats = load_filter_stats("${filter_stats}")

# Determine TMB interpretation
def interpret_tmb(tmb_score):
    if tmb_score >= 20:
        return "Very High"
    elif tmb_score >= 10:
        return "High"
    elif tmb_score >= 6:
        return "Intermediate"
    elif tmb_score >= 1:
        return "Low"
    else:
        return "Very Low"

tmb_interpretation = interpret_tmb(tmb_per_mb)

# Calculate confidence intervals (simplified approach)
# Using Poisson confidence intervals for mutation counts
import math

def poisson_confidence_interval(count, alpha=0.05):
    \"\"\"Calculate Poisson confidence interval\"\"\"
    if count == 0:
        lower = 0
        upper = -math.log(alpha/2)
    else:
        # Approximation for Poisson CI
        lower = max(0, count - 1.96 * math.sqrt(count))
        upper = count + 1.96 * math.sqrt(count)
    
    return lower, upper

lower_count, upper_count = poisson_confidence_interval(coding_mutations)
lower_tmb = lower_count / target_mb if target_mb > 0 else 0
upper_tmb = upper_count / target_mb if target_mb > 0 else 0

# Create comprehensive TMB result
tmb_result = {
    "sample_id": "${sample_id}",
    "analysis_date": datetime.now().isoformat(),
    "tmb_per_mb": round(tmb_per_mb, 2),
    "tmb_total": round(tmb_total, 2),  # Legacy field
    "coding_mutations": coding_mutations,
    "total_variants": total_variants,
    "synonymous_mutations": synonymous_mutations,
    "target_size_mb": round(target_mb, 3),
    "target_size_bases": target_bases,
    "tmb_interpretation": tmb_interpretation,
    "confidence_interval_95": [round(lower_tmb, 2), round(upper_tmb, 2)],
    "quality_metrics": {
        "sufficient_target_size": target_mb >= 1.0,  # At least 1 Mb
        "adequate_mutation_count": coding_mutations >= 5,
        "analysis_quality": "High" if target_mb >= 1.0 and coding_mutations >= 5 else "Moderate"
    },
    "clinical_significance": {
        "immunotherapy_relevant": tmb_per_mb >= 10,
        "hypermutated": tmb_per_mb >= 100,
        "msi_associated": tmb_per_mb >= 20  # Potentially MSI-associated
    },
    "methodology": {
        "mutation_types_counted": ["missense_variant", "nonsense_variant", "frameshift_variant", "inframe_indel"],
        "exclusions": ["synonymous_variant", "intronic_variant", "intergenic_variant"],
        "normalization": "per_megabase_coding_sequence"
    },
    "filter_statistics": filter_stats
}

# Write results
with open("${sample_id}_tmb_calculation.json", 'w') as f:
    json.dump(tmb_result, f, indent=2)

# Console output
print("="*60)
print("TMB Calculation Complete")
print("="*60)
print(f"Sample: ${sample_id}")
print(f"TMB Score: {tmb_per_mb:.1f} mutations/Mb")
print(f"Interpretation: {tmb_interpretation}")
print(f"Coding Mutations: {coding_mutations}")
print(f"Target Size: {target_mb:.2f} Mb")
print(f"95% CI: {lower_tmb:.1f} - {upper_tmb:.1f}")
print(f"Immunotherapy Relevant: {'Yes' if tmb_per_mb >= 10 else 'No'}")
print("="*60)
EOF

    python calculate_tmb.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

process TMB_SUMMARY {
    tag "$sample_id"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(tmb_calculation)
    
    output:
    tuple val(sample_id), path("${sample_id}_tmb_summary.json"), emit: summary
    tuple val(sample_id), path("${sample_id}_tmb_report.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    cat > tmb_summary.py << 'EOF'
import json
from datetime import datetime

# Load TMB calculation
with open("${tmb_calculation}", 'r') as f:
    tmb_data = json.load(f)

# Create clinical summary
clinical_summary = {
    "sample_id": tmb_data["sample_id"],
    "report_date": datetime.now().isoformat(),
    "tmb_score": tmb_data["tmb_per_mb"],
    "tmb_interpretation": tmb_data["tmb_interpretation"],
    "confidence_interval": tmb_data["confidence_interval_95"],
    "clinical_recommendations": [],
    "quality_assessment": tmb_data["quality_metrics"]["analysis_quality"],
    "methodology_note": "TMB calculated as coding mutations per megabase of sequenced target region"
}

# Add clinical recommendations based on TMB score
tmb_score = tmb_data["tmb_per_mb"]

if tmb_score >= 20:
    clinical_summary["clinical_recommendations"].append(
        "Very high TMB (≥20 mut/Mb) - Strong consideration for immunotherapy"
    )
elif tmb_score >= 10:
    clinical_summary["clinical_recommendations"].append(
        "High TMB (≥10 mut/Mb) - Consider immunotherapy (pembrolizumab, nivolumab)"
    )
elif tmb_score >= 6:
    clinical_summary["clinical_recommendations"].append(
        "Intermediate TMB (6-9.9 mut/Mb) - Clinical trial eligibility may vary"
    )
else:
    clinical_summary["clinical_recommendations"].append(
        "Low TMB (<6 mut/Mb) - Less likely to benefit from TMB-directed immunotherapy"
    )

# Quality considerations
if not tmb_data["quality_metrics"]["sufficient_target_size"]:
    clinical_summary["clinical_recommendations"].append(
        "CAUTION: Limited target size may affect TMB accuracy"
    )

if not tmb_data["quality_metrics"]["adequate_mutation_count"]:
    clinical_summary["clinical_recommendations"].append(
        "CAUTION: Low mutation count may affect statistical confidence"
    )

# Write summary JSON
with open("${sample_id}_tmb_summary.json", 'w') as f:
    json.dump(clinical_summary, f, indent=2)

# Generate human-readable report
with open("${sample_id}_tmb_report.txt", 'w') as f:
    f.write("Tumor Mutational Burden (TMB) Analysis Report\\n")
    f.write("=" * 50 + "\\n")
    f.write(f"Sample ID: {clinical_summary['sample_id']}\\n")
    f.write(f"Analysis Date: {clinical_summary['report_date']}\\n")
    f.write("\\n")
    
    f.write("TMB Results:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"TMB Score: {clinical_summary['tmb_score']:.1f} mutations/Mb\\n")
    f.write(f"Interpretation: {clinical_summary['tmb_interpretation']}\\n")
    f.write(f"95% Confidence Interval: {clinical_summary['confidence_interval'][0]:.1f} - {clinical_summary['confidence_interval'][1]:.1f}\\n")
    f.write(f"Quality Assessment: {clinical_summary['quality_assessment']}\\n")
    f.write("\\n")
    
    f.write("Technical Details:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"Coding Mutations: {tmb_data['coding_mutations']}\\n")
    f.write(f"Total Variants: {tmb_data['total_variants']}\\n")
    f.write(f"Target Size: {tmb_data['target_size_mb']:.2f} Mb\\n")
    f.write("\\n")
    
    f.write("Clinical Recommendations:\\n")
    f.write("-" * 30 + "\\n")
    for i, recommendation in enumerate(clinical_summary['clinical_recommendations'], 1):
        f.write(f"{i}. {recommendation}\\n")
    f.write("\\n")
    
    f.write("Methodology:\\n")
    f.write("-" * 15 + "\\n")
    f.write(f"{clinical_summary['methodology_note']}\\n")
    f.write("\\n")
    
    f.write("Clinical Context:\\n")
    f.write("-" * 20 + "\\n")
    f.write("TMB is a biomarker that may predict response to immune checkpoint inhibitors.\\n")
    f.write("FDA-approved TMB thresholds vary by tumor type and assay.\\n")
    f.write("This analysis should be interpreted in conjunction with other clinical factors.\\n")

print(f"TMB Summary Generated:")
print(f"  Sample: {clinical_summary['sample_id']}")
print(f"  TMB Score: {clinical_summary['tmb_score']} mut/Mb")
print(f"  Interpretation: {clinical_summary['tmb_interpretation']}")
print(f"  Quality: {clinical_summary['quality_assessment']}")
EOF

    python tmb_summary.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

// TMB Workflow
workflow TMB {
    take:
    vcf              // channel: [sample_id, vcf]
    targets_bed      // path: target regions BED
    filter_stats     // channel: [sample_id, filter_stats_json]
    
    main:
    // Combine VCF with filter stats
    vcf_with_stats = vcf.join(filter_stats, by: 0)
    
    // Calculate TMB
    CALCULATE_TMB(
        vcf_with_stats.map { sample_id, vcf_file, stats_file -> [sample_id, vcf_file] },
        targets_bed,
        vcf_with_stats.map { sample_id, vcf_file, stats_file -> stats_file }
    )
    
    // Generate TMB summary
    TMB_SUMMARY(CALCULATE_TMB.out.calculation)
    
    // Combine versions
    ch_versions = CALCULATE_TMB.out.versions
        .mix(TMB_SUMMARY.out.versions)
    
    emit:
    calculation = CALCULATE_TMB.out.calculation
    summary = TMB_SUMMARY.out.summary
    report = TMB_SUMMARY.out.report
    versions = ch_versions
}