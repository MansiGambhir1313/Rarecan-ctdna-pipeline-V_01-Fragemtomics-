// modules/align_minimap2.nf - Minimap2 Alignment and fgbio Duplex Consensus Module
// Alternative to BWA alignment that produces *.raw.bam and *.consensus.bam files

nextflow.enable.dsl=2

process MINIMAP2_INDEX {
    tag "minimap2_index"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/reference_index", mode: 'copy'
    
    input:
    path ref_fasta
    
    output:
    path "${ref_fasta}", emit: ref_fasta
    path "${ref_fasta.name}.mmi", emit: mmi_index
    path "versions.yml", emit: versions
    
    script:
    def mmi_file = "${ref_fasta.name}.mmi"
    """
    # Activate conda environment and install minimap2 if not available
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    # Check if minimap2 is available, if not install it
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "minimap2 not found, installing via conda..."
        conda install -y -c bioconda minimap2 2>&1 | tail -5
    fi
    
    # Verify minimap2 is now available
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "ERROR: Failed to install minimap2"
        exit 1
    fi
    
    echo "Using minimap2: \$(minimap2 --version 2>&1 | head -1)"
    minimap2 -d ${mmi_file} ${ref_fasta}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 | head -1 | sed 's/minimap2 //' || echo "not_available")
    END_VERSIONS
    """
}

process MINIMAP2_ALIGN_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.raw.bam*"
    
    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path mmi_index
    
    output:
    tuple val(sample_id), path("${sample_id}.raw.bam"), emit: raw_bam
    tuple val(sample_id), path("${sample_id}.raw.bam.bai"), emit: raw_bai
    path "versions.yml", emit: versions
    
    script:
    def rg_string = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}"
    """
    # Activate conda environment and ensure minimap2 is available
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    # Check if minimap2 is available, if not install it
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "minimap2 not found, installing via conda..."
        conda install -y -c bioconda minimap2 2>&1 | tail -5
    fi
    
    # Verify minimap2 is now available
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "ERROR: Failed to install minimap2"
        exit 1
    fi
    
    echo "Using minimap2: \$(minimap2 --version 2>&1 | head -1)"
    minimap2 -t ${task.cpus} -ax sr --MD -R "${rg_string}" ${mmi_index} ${reads[0]} ${reads[1]} \\
        | samtools view -b - \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.raw.bam
    
    samtools index ${sample_id}.raw.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 | head -1 | sed 's/minimap2 //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/^.*samtools //')
    END_VERSIONS
    """
}

process GROUP_READS_BY_UMI {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/umi", mode: 'copy', pattern: "*.grouped.bam*"
    
    input:
    tuple val(sample_id), path(raw_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.grouped.bam"), emit: grouped_bam
    path "versions.yml", emit: versions
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    if ! command -v fgbio >/dev/null 2>&1; then
        echo "fgbio not found, installing via conda..."
        conda install -y -c bioconda fgbio 2>&1 | tail -5
    fi
    
    fgbio GroupReadsByUmi \\
        --input ${raw_bam} \\
        --output ${sample_id}.grouped.bam \\
        --strategy=paired --edits=1
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | head -1 | sed 's/fgbio //')
    END_VERSIONS
    """
}

process CALL_DUPLEX_CONSENSUS {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/umi", mode: 'copy', pattern: "*.duplex.bam*"
    
    input:
    tuple val(sample_id), path(grouped_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.duplex.bam"), emit: duplex_bam
    path "versions.yml", emit: versions
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    if ! command -v fgbio >/dev/null 2>&1; then
        echo "fgbio not found, installing via conda..."
        conda install -y -c bioconda fgbio 2>&1 | tail -5
    fi
    
    fgbio CallDuplexConsensusReads \\
        --input ${grouped_bam} \\
        --output ${sample_id}.duplex.bam \\
        --min-reads=3 2 1 \\
        --error-rate-pre-umi=45 \\
        --error-rate-post-umi=30 \\
        --min-input-base-quality=20
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | head -1 | sed 's/fgbio //')
    END_VERSIONS
    """
}

process FILTER_CONSENSUS_READS {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/umi", mode: 'copy', pattern: "*.filtered.duplex.bam*"
    
    input:
    tuple val(sample_id), path(duplex_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.duplex.bam"), emit: filtered_duplex_bam
    path "versions.yml", emit: versions
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    if ! command -v fgbio >/dev/null 2>&1; then
        echo "fgbio not found, installing via conda..."
        conda install -y -c bioconda fgbio 2>&1 | tail -5
    fi
    
    fgbio FilterConsensusReads \\
        --input ${duplex_bam} \\
        --output ${sample_id}.filtered.duplex.bam \\
        --min-reads=1 --max-no-call-fraction=0.35
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | head -1 | sed 's/fgbio //')
    END_VERSIONS
    """
}

process REALIGN_CONSENSUS {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.consensus.bam*"
    
    input:
    tuple val(sample_id), path(filtered_duplex_bam)
    path ref_fasta
    path mmi_index
    
    output:
    tuple val(sample_id), path("${sample_id}.consensus.bam"), emit: consensus_bam
    tuple val(sample_id), path("${sample_id}.consensus.bam.bai"), emit: consensus_bai
    path "versions.yml", emit: versions
    
    script:
    """
    # Activate conda environment and ensure minimap2 is available
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
    conda activate ctdna 2>/dev/null || true
    
    # Check if minimap2 is available, if not install it
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "minimap2 not found, installing via conda..."
        conda install -y -c bioconda minimap2 2>&1 | tail -5
    fi
    
    # Verify minimap2 is now available
    if ! command -v minimap2 >/dev/null 2>&1; then
        echo "ERROR: Failed to install minimap2"
        exit 1
    fi
    
    samtools fastq ${filtered_duplex_bam} \\
        | minimap2 -t ${task.cpus} -ax sr --MD ${mmi_index} - \\
        | samtools view -b - \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.consensus.bam
    
    samtools index ${sample_id}.consensus.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 | head -1 | sed 's/minimap2 //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/^.*samtools //')
    END_VERSIONS
    """
}

// Main Workflow - Minimap2 alignment with fgbio duplex consensus
workflow ALIGN_MINIMAP2 {
    take:
    reads_ch
    ref_fasta
    
    main:
    // Create minimap2 index
    index_result = MINIMAP2_INDEX(ref_fasta)
    
    // Align FASTQ directly with minimap2 to create raw BAM
    raw_bam_result = MINIMAP2_ALIGN_FASTQ(
        reads_ch,
        index_result.ref_fasta,
        index_result.mmi_index
    )
    
    if( params.enable_umi ) {
    // Generate duplex consensus
    grouped_result = GROUP_READS_BY_UMI(raw_bam_result.raw_bam)
    duplex_result = CALL_DUPLEX_CONSENSUS(grouped_result.grouped_bam)
    filtered_result = FILTER_CONSENSUS_READS(duplex_result.duplex_bam)
    consensus_result = REALIGN_CONSENSUS(
        filtered_result.filtered_duplex_bam,
        index_result.ref_fasta,
        index_result.mmi_index
    )
    
    ch_versions = index_result.versions
        .mix(raw_bam_result.versions)
        .mix(grouped_result.versions)
        .mix(duplex_result.versions)
        .mix(filtered_result.versions)
        .mix(consensus_result.versions)
        
        consensus_bam_ch = consensus_result.consensus_bam
        consensus_bai_ch = consensus_result.consensus_bai
    }
    else {
        ch_versions = index_result.versions
            .mix(raw_bam_result.versions)
        
        consensus_bam_ch = raw_bam_result.raw_bam
        consensus_bai_ch = raw_bam_result.raw_bai
    }
    
    emit:
    raw_bam = raw_bam_result.raw_bam
    raw_bai = raw_bam_result.raw_bai
    consensus_bam = consensus_bam_ch
    consensus_bai = consensus_bai_ch
    versions = ch_versions
}

