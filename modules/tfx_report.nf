// modules/tfx_report.nf - TFX Clinical Report Generation
// Wraps run_final_report.py from the RareCan TFX pipeline

nextflow.enable.dsl=2

process TFX_FINAL_REPORT {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/tfx", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple val(sample_id), 
          path(ichor_tfx), 
          path(global_hist), 
          path(frag_vcf), 
          path(frag_vcf_index)
    path template_html
    
    output:
    tuple val(sample_id), path("${sample_id}.tfx.report.json"), emit: json
    tuple val(sample_id), path("${sample_id}.tfx.report.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    def min_bp = params.frag_min_bp ?: 90
    def max_bp = params.frag_max_bp ?: 150
    def ks_pval = params.frag_ks_pval ?: 0.05
    def cna_lod = params.cna_lod_cutoff ?: 0.03
    """
    set -euo pipefail
    
    # Validate required files exist
    if [ ! -f "${template_html}" ]; then
        echo "ERROR: Template not found: ${template_html}" >&2
        exit 1
    fi
    
    if [ ! -f "${projectDir}/bin/run_final_report.py" ]; then
        echo "ERROR: Script not found: ${projectDir}/bin/run_final_report.py" >&2
        exit 1
    fi
    
    # Handle optional input files - check existence in bash
    ICHOR_FILE="${ichor_tfx}"
    HIST_FILE="${global_hist}"
    VCF_FILE="${frag_vcf}"
    
    if [ "$ICHOR_FILE" != "NO_FILE" ] && [ ! -f "$ICHOR_FILE" ]; then
        echo "WARNING: ichorCNA file not found: $ICHOR_FILE, using NO_FILE" >&2
        ICHOR_FILE="NO_FILE"
    fi
    
    if [ "$HIST_FILE" != "NO_FILE" ] && [ ! -f "$HIST_FILE" ]; then
        echo "WARNING: Histogram file not found: $HIST_FILE, using NO_FILE" >&2
        HIST_FILE="NO_FILE"
    fi
    
    if [ "$VCF_FILE" != "NO_FILE" ] && [ ! -f "$VCF_FILE" ]; then
        echo "WARNING: VCF file not found: $VCF_FILE, using NO_FILE" >&2
        VCF_FILE="NO_FILE"
    fi
    
    python3 ${projectDir}/bin/run_final_report.py \\
        --sample_id ${sample_id} \\
        --ichor_tfx "$ICHOR_FILE" \\
        --annot_vcf "$VCF_FILE" \\
        --global_hist "$HIST_FILE" \\
        --min_bp ${min_bp} \\
        --max_bp ${max_bp} \\
        --ks_pval ${ks_pval} \\
        --cna_lod ${cna_lod} \\
        --template_html ${template_html} \\
        --out_json ${sample_id}.tfx.report.json \\
        --out_html ${sample_id}.tfx.report.html
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """
}

workflow TFX_REPORT {
    take:
    ichor_tfx
    global_hist
    frag_vcf
    frag_vcf_index
    template_html
    
    main:
    // Join VCF with its index first
    frag_vcf_combined = frag_vcf.join(frag_vcf_index, by: 0)
        .map { sid, vcf, idx -> [sid, vcf, idx] }
    
    // Combine all inputs - join by sample_id
    // Use map to handle empty channels gracefully
    combined = ichor_tfx
        .join(global_hist, by: 0)
        .join(frag_vcf_combined, by: 0)
        .map { sid, tfx, hist, vcf, idx -> 
            // Ensure all paths are valid or use NO_FILE placeholder
            def tfx_file = tfx ?: file("NO_FILE")
            def hist_file = hist ?: file("NO_FILE")
            def vcf_file = vcf ?: file("NO_FILE")
            def idx_file = idx ?: file("NO_FILE")
            [sid, tfx_file, hist_file, vcf_file, idx_file]
        }
    
    TFX_FINAL_REPORT(combined, template_html)
    
    emit:
    json = TFX_FINAL_REPORT.out.json
    html = TFX_FINAL_REPORT.out.html
    versions = TFX_FINAL_REPORT.out.versions
}
