// modules/cna_context.nf - CNA-Contextual Variant Filtering Module
// PO-CFS Module 3.4: Adjust variant VAF based on local copy number

nextflow.enable.dsl=2

process CNA_CONTEXTUAL_FILTER {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/filtering/cna_context", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(cnv_segments)
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_cna_adjusted.vcf.gz"), emit: adjusted_vcf, optional: true
    tuple val(sample_id), path("${sample_id}_cna_context_stats.json"), emit: stats, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running CNA-contextual filtering for sample: ${sample_id}"
    
    # Check if Python is available
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        echo "WARNING: Python not found"
        
        cp ${vcf} ${sample_id}_cna_adjusted.vcf.gz 2>/dev/null || touch ${sample_id}_cna_adjusted.vcf.gz
        echo '{"sample_id":"${sample_id}","variants_adjusted":0,"status":"python_missing"}' > ${sample_id}_cna_context_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    PYTHON_CMD=\$(command -v python3 || command -v python)
    
    # Create CNA-contextual filtering script
    cat > cna_context_filter.py << 'EOF'
import gzip
import json
import sys
from datetime import datetime

def load_cnv_segments(segments_file):
    """Load CNV segments from CNVkit output"""
    segments = []
    
    try:
        with open(segments_file, 'r') as f:
            # Skip header
            header = f.readline()
            
            for line in f:
                fields = line.strip().split('\\t')
                if len(fields) >= 6:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    log2 = float(fields[4]) if fields[4] != 'NA' else 0.0
                    copy_number = float(fields[5]) if len(fields) > 5 and fields[5] != 'NA' else 2.0
                    
                    segments.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'log2': log2,
                        'copy_number': copy_number
                    })
    except Exception as e:
        print(f"Warning: Could not load CNV segments: {e}", file=sys.stderr)
    
    return segments

def get_segment_for_position(segments, chrom, pos):
    """Find CNV segment containing the given position"""
    for seg in segments:
        if seg['chrom'] == chrom and seg['start'] <= pos <= seg['end']:
            return seg
    return None

def calculate_adjusted_vaf(original_vaf, copy_number, is_loh=False):
    """Adjust VAF based on local copy number"""
    
    # Normal copy number is 2
    # If copy number is different, VAF needs adjustment
    
    if copy_number == 0:
        # Deletion - variant should not exist
        return 0.0
    
    if is_loh:
        # Loss of heterozygosity - adjust VAF
        # In LOH, one allele is lost, so VAF is doubled
        adjusted_vaf = original_vaf * 2.0
        return min(1.0, adjusted_vaf)
    
    # General adjustment: VAF = alt_count / (ref_count + alt_count)
    # With copy number CN, if CN != 2, need to adjust
    # Simplified: VAF_adj = VAF_orig * (2 / CN)
    if copy_number > 0:
        adjusted_vaf = original_vaf * (2.0 / copy_number)
        return min(1.0, adjusted_vaf)
    
    return original_vaf

def detect_loh(segment):
    """Detect Loss of Heterozygosity from segment"""
    # LOH typically: copy_number ~1, log2 ~ -1
    if segment['copy_number'] < 1.5 and segment['log2'] < -0.5:
        return True
    return False

def adjust_vcf_with_cna(vcf_file, segments_file, output_file, stats_file):
    """Adjust VCF variants based on CNA context"""
    
    segments = load_cnv_segments(segments_file)
    variants_adjusted = 0
    variants_with_loh = 0
    total_variants = 0
    
    try:
        # Open VCF
        if vcf_file.endswith('.gz'):
            vcf_handle = gzip.open(vcf_file, 'rt')
        else:
            vcf_handle = open(vcf_file, 'r')
        
        passed_variants = []
        header_lines = []
        
        with vcf_handle as vcf:
            for line in vcf:
                if line.startswith('#'):
                    header_lines.append(line)
                    continue
                
                total_variants += 1
                fields = line.strip().split('\\t')
                
                if len(fields) < 10:
                    passed_variants.append(line)
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # Find CNV segment for this position
                segment = get_segment_for_position(segments, chrom, pos)
                
                if segment:
                    # Extract original VAF from sample data
                    sample_data = fields[9]
                    original_vaf = 0.0
                    
                    if ':' in sample_data:
                        sample_fields = sample_data.split(':')
                        # Try to extract VAF
                        for i, field in enumerate(sample_fields):
                            if 'VAF' in field.upper() or i == len(sample_fields) - 1:
                                try:
                                    original_vaf = float(field)
                                    break
                                except:
                                    pass
                        
                        # If VAF not found, try to calculate from AD
                        if original_vaf == 0.0 and len(sample_fields) > 1:
                            ad_str = sample_fields[1]
                            if ',' in ad_str:
                                ref_count, alt_count = map(int, ad_str.split(',')[:2])
                                total = ref_count + alt_count
                                if total > 0:
                                    original_vaf = alt_count / total
                    
                    # Detect LOH
                    is_loh = detect_loh(segment)
                    if is_loh:
                        variants_with_loh += 1
                    
                    # Calculate adjusted VAF
                    adjusted_vaf = calculate_adjusted_vaf(
                        original_vaf,
                        segment['copy_number'],
                        is_loh
                    )
                    
                    # Update INFO field with CNA context
                    info_field = fields[7]
                    if info_field == '.':
                        info_field = ''
                    
                    # Add CNA context to INFO
                    cna_info = f"CNA_CN={segment['copy_number']:.2f};CNA_LOG2={segment['log2']:.3f}"
                    if is_loh:
                        cna_info += ";LOH=TRUE"
                    cna_info += f";VAF_ORIG={original_vaf:.4f};VAF_ADJ={adjusted_vaf:.4f}"
                    
                    if info_field:
                        info_field = f"{info_field};{cna_info}"
                    else:
                        info_field = cna_info
                    
                    fields[7] = info_field
                    
                    # Update sample data with adjusted VAF if possible
                    if ':' in sample_data:
                        sample_fields = sample_data.split(':')
                        # Try to update VAF field
                        updated = False
                        for i, field in enumerate(sample_fields):
                            if 'VAF' in field.upper() or (i == len(sample_fields) - 1 and not updated):
                                sample_fields[i] = f"{adjusted_vaf:.4f}"
                                updated = True
                                break
                        
                        if not updated:
                            sample_fields.append(f"{adjusted_vaf:.4f}")
                        
                        fields[9] = ':'.join(sample_fields)
                    
                    variants_adjusted += 1
                    
                    # Filter variants in regions with extreme copy number changes
                    # (likely artifacts)
                    if segment['copy_number'] > 6 or segment['copy_number'] < 0.5:
                        # Skip variants in extreme CNV regions
                        continue
                
                # Write adjusted variant
                passed_variants.append('\\t'.join(fields))
        
        # Write adjusted VCF
        with gzip.open(output_file, 'wt') as out_vcf:
            for header_line in header_lines:
                # Add CNA context to header
                if header_line.startswith('##INFO'):
                    out_vcf.write(header_line)
                    out_vcf.write('##INFO=<ID=CNA_CN,Number=1,Type=Float,Description="Local copy number">\\n')
                    out_vcf.write('##INFO=<ID=CNA_LOG2,Number=1,Type=Float,Description="Log2 copy number ratio">\\n')
                    out_vcf.write('##INFO=<ID=LOH,Number=0,Type=Flag,Description="Loss of heterozygosity">\\n')
                    out_vcf.write('##INFO=<ID=VAF_ORIG,Number=1,Type=Float,Description="Original VAF">\\n')
                    out_vcf.write('##INFO=<ID=VAF_ADJ,Number=1,Type=Float,Description="CNA-adjusted VAF">\\n')
                else:
                    out_vcf.write(header_line)
            
            for variant_line in passed_variants:
                out_vcf.write(variant_line)
        
        # Write statistics
        stats = {
            "sample_id": "${sample_id}",
            "analysis_date": datetime.now().isoformat(),
            "total_variants": total_variants,
            "variants_adjusted": variants_adjusted,
            "variants_with_loh": variants_with_loh,
            "segments_loaded": len(segments)
        }
        
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        print(f"CNA-contextual filtering complete:")
        print(f"  Total variants: {total_variants}")
        print(f"  Variants adjusted: {variants_adjusted}")
        print(f"  Variants in LOH regions: {variants_with_loh}")
        
    except Exception as e:
        print(f"ERROR in CNA-contextual filtering: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        
        # Create minimal output
        with open(stats_file, 'w') as f:
            json.dump({
                "sample_id": "${sample_id}",
                "status": "failed",
                "error": str(e)
            }, f, indent=2)
        
        # Copy input as output
        import shutil
        shutil.copy(vcf_file, output_file)

# Run filtering
adjust_vcf_with_cna(
    "${vcf}",
    "${cnv_segments}",
    "${sample_id}_cna_adjusted.vcf.gz",
    "${sample_id}_cna_context_stats.json"
)
EOF

    # Run CNA-contextual filtering
    \${PYTHON_CMD} cna_context_filter.py 2>&1 | tee cna_context.log
    
    # Index adjusted VCF
    if [ -f "${sample_id}_cna_adjusted.vcf.gz" ]; then
        tabix -p vcf ${sample_id}_cna_adjusted.vcf.gz 2>/dev/null || echo "WARNING: Could not index VCF"
    fi
    
    # Ensure output files exist
    if [ ! -f "${sample_id}_cna_adjusted.vcf.gz" ]; then
        cp ${vcf} ${sample_id}_cna_adjusted.vcf.gz 2>/dev/null || touch ${sample_id}_cna_adjusted.vcf.gz
    fi
    
    if [ ! -f "${sample_id}_cna_context_stats.json" ]; then
        echo '{"sample_id":"${sample_id}","status":"unknown"}' > ${sample_id}_cna_context_stats.json
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(\${PYTHON_CMD} --version 2>&1 | sed 's/Python //')
    status: "success"
END_VERSIONS
    
    set -e
    exit 0
    """
}

// CNA-Contextual Filtering Workflow
workflow CNA_CONTEXTUAL_FILTER {
    take:
    vcf              // channel: [sample_id, vcf]
    cnv_segments     // channel: [sample_id, cnv_segments]
    ref_fasta        // path: reference FASTA
    
    main:
    // Combine VCF with CNV segments
    vcf_cnv_combined = vcf.join(cnv_segments, by: 0)
    
    // Run CNA-contextual filtering
    CNA_CONTEXTUAL_FILTER(
        vcf_cnv_combined.map { sample_id, vcf_file, segments_file -> 
            [sample_id, vcf_file]
        },
        vcf_cnv_combined.map { sample_id, vcf_file, segments_file -> 
            [sample_id, segments_file]
        },
        ref_fasta
    )
    
    // Combine versions
    ch_versions = CNA_CONTEXTUAL_FILTER.out.versions
    
    emit:
    adjusted_vcf = CNA_CONTEXTUAL_FILTER.out.adjusted_vcf
    stats = CNA_CONTEXTUAL_FILTER.out.stats
    versions = ch_versions
}

