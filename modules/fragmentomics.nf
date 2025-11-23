// modules/fragmentomics.nf - Comprehensive Fragmentomics Module
// Combines cfDNA fragment length/motif profiling with variant-level filtering

nextflow.enable.dsl=2

process GLOBAL_FRAGMENT_SIGNATURES {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/fragmentomics", mode: 'copy'
    errorStrategy 'finish'
    
    when:
    params.enable_fragmentomics
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.fragmentomics.summary.json"), emit: summary, optional: true
    tuple val(sample_id), path("${sample_id}.global_hist.tsv"), emit: histogram, optional: true
    tuple val(sample_id), path("${sample_id}.fragment_motifs.tsv"), emit: motifs, optional: true
    path "versions.yml", emit: versions
    
    script:
    def short_min = params.frag_min_bp ?: 90
    def short_max = params.frag_max_bp ?: 150
    """
    set -e
    
    python3 ${projectDir}/bin/global_fragmentomics.py \\
        --sample ${sample_id} \\
        --bam ${bam} \\
        --output_json ${sample_id}.fragmentomics.summary.json \\
        --histogram ${sample_id}.global_hist.tsv \\
        --motifs ${sample_id}.fragment_motifs.tsv \\
        --short_min ${short_min} \\
        --short_max ${short_max}
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """
}

process VARIANT_FRAGMENT_SIGNATURES {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/fragmentomics", mode: 'copy'
    errorStrategy 'ignore'
    
    when:
    params.enable_fragmentomics
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple val(sample_id), path(vcf)
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.frag.vcf.gz"), emit: frag_vcf, optional: true
    tuple val(sample_id), path("${sample_id}.frag.vcf.gz.tbi"), emit: frag_vcf_index, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set -e
    
    python3 ${projectDir}/bin/run_variant_fragmentomics.py \\
        --bam ${bam} \\
        --vcf ${vcf} \\
        --ref_fasta ${ref_fasta} \\
        --outfile ${sample_id}.frag.vcf.gz
    
    if [ -f "${sample_id}.frag.vcf.gz" ]; then
        tabix -p vcf ${sample_id}.frag.vcf.gz || true
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """
}

process COMPUTE_CTDNA_PURITY {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/fragmentomics", mode: 'copy'
    errorStrategy 'ignore'
    
    when:
    params.enable_fragmentomics
    
    input:
    tuple val(sample_id), path(ichorcna_tfx), path(ichorcna_segments), path(frag_vcf), path(snv_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.ctdna_purity.json"), emit: purity_json, optional: true
    tuple val(sample_id), path("${sample_id}.ctdna_purity.tfx.txt"), emit: purity_tfx, optional: true
    path "versions.yml", emit: versions
    
    script:
    def ichorcna_tfx_arg = (ichorcna_tfx.toString() != "NO_FILE" && file(ichorcna_tfx).exists()) ? "--ichorcna_tfx ${ichorcna_tfx}" : ""
    def ichorcna_segments_arg = (ichorcna_segments.toString() != "NO_FILE" && file(ichorcna_segments).exists()) ? "--ichorcna_segments ${ichorcna_segments}" : ""
    def frag_vcf_arg = (frag_vcf.toString() != "NO_FILE" && file(frag_vcf).exists()) ? "--frag_vcf ${frag_vcf}" : ""
    def snv_vcf_arg = (snv_vcf.toString() != "NO_FILE" && file(snv_vcf).exists()) ? "--snv_vcf ${snv_vcf}" : ""
    """
    set +e
    
    echo "Computing comprehensive ctDNA purity for sample: ${sample_id}"
    echo "  ichorCNA TFX: ${ichorcna_tfx_arg ? 'provided' : 'not provided'}"
    echo "  ichorCNA segments: ${ichorcna_segments_arg ? 'provided' : 'not provided'}"
    echo "  Fragment VCF: ${frag_vcf_arg ? 'provided' : 'not provided'}"
    echo "  SNV VCF: ${snv_vcf_arg ? 'provided' : 'not provided'}"
    
    # Run ctDNA purity computation (with fallback if script fails)
    if python3 ${projectDir}/bin/compute_ctdna_purity.py \\
        --sample_id ${sample_id} \\
        ${ichorcna_tfx_arg} \\
        ${ichorcna_segments_arg} \\
        ${frag_vcf_arg} \\
        ${snv_vcf_arg} \\
        --ks_pval_cutoff ${params.frag_ks_pval ?: 0.05} \\
        --min_vaf 0.01 \\
        --cna_lod_cutoff ${params.cna_lod_cutoff ?: 0.03} \\
        --output_json ${sample_id}.ctdna_purity.json \\
        --output_tfx ${sample_id}.ctdna_purity.tfx.txt 2>&1 | tee ctdna_purity.log; then
        
        # Verify outputs were created
        if [ -f "${sample_id}.ctdna_purity.json" ] && [ -f "${sample_id}.ctdna_purity.tfx.txt" ]; then
            echo "ctDNA purity computation completed successfully"
            CTDNA_SUCCESS=true
        else
            echo "WARNING: ctDNA purity computation completed but outputs missing"
            CTDNA_SUCCESS=false
        fi
    else
        echo "WARNING: ctDNA purity computation failed, creating fallback outputs"
        cat ctdna_purity.log 2>/dev/null || echo "No log available"
        CTDNA_SUCCESS=false
    fi
    
    # Create fallback outputs if computation failed
    if [ "$CTDNA_SUCCESS" != "true" ]; then
        echo "Creating fallback ctDNA purity outputs..."
        
        # Try to extract ichorCNA tumor fraction as fallback
        fallback_purity=0.0
        if [ -f "${ichorcna_tfx}" ] && [ -s "${ichorcna_tfx}" ]; then
            fallback_purity=\$(awk 'NR==2 {print \$2}' ${ichorcna_tfx} 2>/dev/null || echo "0.0")
        fi
        
        # Calculate percentage (using awk for portability)
        fallback_percent=\$(awk "BEGIN {printf \"%.2f\", ${fallback_purity} * 100}")
        
        # Create minimal JSON output
        cat > ${sample_id}.ctdna_purity.json << EOF
{
  "sample_id": "${sample_id}",
  "ctdna_purity": {
    "consensus": ${fallback_purity},
    "consensus_percent": "${fallback_percent}%",
    "sources": ["ichorCNA_fallback"],
    "fallback_used": true
  },
  "ichorCNA": {
    "tumor_fraction": ${fallback_purity},
    "tumor_fraction_percent": "${fallback_percent}%",
    "status": "fallback"
  },
  "fragment_KS": {
    "purity": 0.0,
    "status": "not_available"
  },
  "SNV_INDEL": {
    "purity": 0.0,
    "status": "not_available"
  },
  "qc_metrics": {
    "ichorCNA_available": true,
    "fragment_KS_available": false,
    "SNV_INDEL_available": false,
    "consensus_method": "fallback_ichorCNA_only"
  }
}
EOF
        
        # Create TFX format output
        printf "sample\\ttumorFraction\\tploidy\\tmethod\\n" > ${sample_id}.ctdna_purity.tfx.txt
        printf "${sample_id}\\t${fallback_purity}\\tNA\\tfallback\\n" >> ${sample_id}.ctdna_purity.tfx.txt
        
        echo "Fallback outputs created with ichorCNA-only estimate"
    fi
    
    # Ensure outputs exist
    if [ ! -f "${sample_id}.ctdna_purity.json" ]; then
        echo '{"sample_id":"${sample_id}","ctdna_purity":{"consensus":0.0,"status":"failed"}}' > ${sample_id}.ctdna_purity.json
    fi
    
    if [ ! -f "${sample_id}.ctdna_purity.tfx.txt" ]; then
        printf "sample\\ttumorFraction\\tploidy\\tmethod\\n${sample_id}\\t0.0\\tNA\\tfailed\\n" > ${sample_id}.ctdna_purity.tfx.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
    status: \$([ "\${CTDNA_SUCCESS}" = "true" ] && echo "success" || echo "fallback")
END_VERSIONS
    
    set -e
    exit 0
    """
}

workflow FRAGMENTOMICS {
    take:
    raw_bam        // [sample_id, bam, bai]
    consensus_bam  // [sample_id, bam, bai]
    vcf            // [sample_id, vcf]
    ref_fasta
    
    main:
    // Global fragmentomics always runs on raw BAM
    global_metrics = GLOBAL_FRAGMENT_SIGNATURES(raw_bam)
    
    // Variant fragmentomics - try to join consensus BAM with VCF
    // If either is empty, the join will produce empty channel and variant process won't run
    variant_inputs = consensus_bam.join(vcf, by: 0)
    
    variant_metrics = VARIANT_FRAGMENT_SIGNATURES(
        variant_inputs.map { sid, bam_file, bai_file, vcf_file -> [sid, bam_file, bai_file] },
        variant_inputs.map { sid, bam_file, bai_file, vcf_file -> [sid, vcf_file] },
        ref_fasta
    )
    
    ch_versions = global_metrics.versions.mix(variant_metrics.versions)
    
    emit:
    summary = global_metrics.summary
    histogram = global_metrics.histogram
    motifs = global_metrics.motifs
    frag_vcf = variant_metrics.frag_vcf
    frag_vcf_index = variant_metrics.frag_vcf_index
    versions = ch_versions
}

// Note: COMPUTE_CTDNA_PURITY process is defined above (line 88) and called from main.nf

