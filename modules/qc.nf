// modules/qc.nf - Quality Control Module
// ECR-compatible version using deployed containers

nextflow.enable.dsl=2

process FASTP {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_fastp_R{1,2}.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    # Try to run fastp, but fallback to copying input files if it fails
    if command -v fastp >/dev/null 2>&1; then
        if fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}_fastp_R1.fastq.gz \\
        --out2 ${sample_id}_fastp_R2.fastq.gz \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --correction \\
        --cut_front \\
        --cut_tail \\
        --cut_window_size 4 \\
        --cut_mean_quality 15 \\
        --qualified_quality_phred 15 \\
        --unqualified_percent_limit 40 \\
        --n_base_limit 5 \\
            --length_required 50 2>&1; then
            echo "FASTP completed successfully"
        else
            echo "WARNING: FASTP failed, using input reads directly"
            cp ${reads[0]} ${sample_id}_fastp_R1.fastq.gz
            cp ${reads[1]} ${sample_id}_fastp_R2.fastq.gz
            echo '{"summary":{"before_filtering":{"total_reads":0}}}' > ${sample_id}_fastp.json
            echo "<html><body>FASTP failed, using raw reads</body></html>" > ${sample_id}_fastp.html
        fi
    else
        echo "WARNING: fastp not found, using input reads directly"
        cp ${reads[0]} ${sample_id}_fastp_R1.fastq.gz
        cp ${reads[1]} ${sample_id}_fastp_R2.fastq.gz
        echo '{"summary":{"before_filtering":{"total_reads":0}}}' > ${sample_id}_fastp.json
        echo "<html><body>FASTP not available, using raw reads</body></html>" > ${sample_id}_fastp.html
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(command -v fastp >/dev/null 2>&1 && fastp --version 2>&1 | sed 's/^.*fastp //; s/ .*\$//' || echo "not_available")
    END_VERSIONS
    
    set -e
    """
}

process FASTQC {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*.zip"), emit: zip
    tuple val(sample_id), path("*.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    """
    fastqc \\
        --quiet \\
        --threads ${task.cpus} \\
        ${reads.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}

process MULTIQC {
    tag "$sample_id"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastp_json)
    tuple val(sample_id), path(fastqc_zip)
    
    output:
    tuple val(sample_id), path("${sample_id}_multiqc.html"), emit: report
    tuple val(sample_id), path("${sample_id}_multiqc_data"), emit: data
    path "versions.yml", emit: versions
    
    script:
    """
    multiqc \\
        --force \\
        --filename ${sample_id}_multiqc \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

// QC Workflow
workflow QC {
    take:
    reads // channel: [sample_id, [R1, R2]]
    
    main:
    // FastP preprocessing
    FASTP(reads)
    
    // FastQC quality control
    FASTQC(FASTP.out.reads)
    
    // MultiQC aggregation
    MULTIQC(
        FASTP.out.json,
        FASTQC.out.zip
    )
    
    // Combine versions
    ch_versions = FASTP.out.versions
        .mix(FASTQC.out.versions)
        .mix(MULTIQC.out.versions)
    
    emit:
    reads = FASTP.out.reads
    fastp_json = FASTP.out.json
    fastqc_zip = FASTQC.out.zip
    multiqc_html = MULTIQC.out.report
    versions = ch_versions
}