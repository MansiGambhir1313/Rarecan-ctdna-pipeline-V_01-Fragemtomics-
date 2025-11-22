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

