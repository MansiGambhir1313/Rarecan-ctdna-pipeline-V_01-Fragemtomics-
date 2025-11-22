// modules/umi.nf - UMI Processing Module
// Fixed version with correct container and tool dependencies

nextflow.enable.dsl=2

process FASTQ_TO_BAM {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/umi", mode: 'copy'
    
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    
    input:
    tuple val(sample_id), path(reads)
    val umi_schema
    
    output:
    tuple val(sample_id), path("${sample_id}_unmapped.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    """
    # Convert FASTQ to unmapped BAM with UMI information
    # Using correct fgbio read structure format: 8M+T (8 bases molecular index + template)
    fgbio FastqToBam \\
        --input ${reads[0]} \\
        --input ${reads[1]} \\
        --output ${sample_id}_unmapped.bam \\
        --read-structures 8M+T +T \\
        --sample ${sample_id} \\
        --library ${sample_id} \\
        --platform illumina \\
        --read-group-id ${sample_id} \\
        --platform-unit ${sample_id} \\
        --sequencing-center "Clinical Lab" \\
        --predicted-insert-size 200 \\
        --sort false
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/^fgbio //' )
    END_VERSIONS
    """
}

process EXTRACT_UMIS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/umi", mode: 'copy'
    
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    
    input:
    tuple val(sample_id), path(unmapped_bam)
    val umi_schema
    
    output:
    tuple val(sample_id), path("${sample_id}_umi_extracted.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    """
    # Extract UMIs from unmapped BAM
    # CRITICAL FIX: Must provide --molecular-index-tags with exactly ONE tag 
    # for the single molecular index in read structure 8M+T +T
    fgbio ExtractUmisFromBam \\
        --input ${unmapped_bam} \\
        --output ${sample_id}_umi_extracted.bam \\
        --read-structure 8M+T +T \\
        --molecular-index-tags RX \\
        --annotate-read-names false
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/^fgbio //' )
    END_VERSIONS
    """
}

process ALIGN_UMIS {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    // CRITICAL FIX: Use fgbio container which has samtools, bwa-mem2, and fgbio
    // OR create a combined container with all three tools
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    
    input:
    tuple val(sample_id), path(umi_bam)
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned_umi.bam"), path("${sample_id}_aligned_umi.bam.bai"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    """
    # Check if required tools are available
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is not installed in container"; exit 1; }
    command -v bwa-mem2 >/dev/null 2>&1 || { echo "ERROR: bwa-mem2 is not installed in container"; exit 1; }
    command -v fgbio >/dev/null 2>&1 || { echo "ERROR: fgbio is not installed in container"; exit 1; }
    
    # Convert BAM to FASTQ for alignment
    samtools collate -u -O ${umi_bam} | \\
        samtools fastq -1 ${sample_id}_R1.fq -2 ${sample_id}_R2.fq -0 /dev/null -s /dev/null -n -
    
    # Align with BWA-MEM2
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}" \\
        ${ref_fasta} \\
        ${sample_id}_R1.fq ${sample_id}_R2.fq | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}_aligned_temp.bam -
    
    # Copy UMI tags from the original BAM to aligned BAM
    # Copy RX tag which contains the molecular index
    fgbio ZipperBam \\
        --unmapped ${umi_bam} \\
        --mapped ${sample_id}_aligned_temp.bam \\
        --output ${sample_id}_aligned_umi.bam \\
        --tags-to-copy RX
    
    # Index the final BAM
    samtools index ${sample_id}_aligned_umi.bam
    
    # Clean up temporary files
    rm -f ${sample_id}_R1.fq ${sample_id}_R2.fq ${sample_id}_aligned_temp.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/^fgbio //' )
    END_VERSIONS
    """
}

process GROUP_READS_BY_UMI {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/umi", mode: 'copy'
    
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_grouped.bam"), emit: bam
    path "${sample_id}_family_sizes.txt", emit: family_sizes
    path "versions.yml", emit: versions
    
    script:
    """
    fgbio GroupReadsByUmi \\
        --input ${bam} \\
        --output ${sample_id}_grouped.bam \\
        --strategy paired \\
        --edits 1 \\
        --family-size-histogram ${sample_id}_family_sizes.txt \\
        --min-map-q 20
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/^fgbio //' )
    END_VERSIONS
    """
}

process CALL_CONSENSUS_READS {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/umi", mode: 'copy'
    
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    
    input:
    tuple val(sample_id), path(grouped_bam)
    val enable_duplex
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_consensus.bam"), path("${sample_id}_consensus.bam.bai"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    def consensus_caller = enable_duplex ? "CallDuplexConsensusReads" : "CallMolecularConsensusReads"
    """
    # Call consensus reads
    if [ "${enable_duplex}" = "true" ]; then
        # Duplex consensus with proper min-reads parameter
        fgbio CallDuplexConsensusReads \\
            --input ${grouped_bam} \\
            --output ${sample_id}_consensus_raw.bam \\
            --error-rate-pre-umi 45 \\
            --error-rate-post-umi 40 \\
            --min-input-base-quality 20 \\
            --min-reads 3 2 1 \\
            --max-reads-per-strand 0
    else
        # Single-strand consensus
        fgbio CallMolecularConsensusReads \\
        --input ${grouped_bam} \\
        --output ${sample_id}_consensus_raw.bam \\
        --error-rate-pre-umi 45 \\
        --error-rate-post-umi 40 \\
        --min-input-base-quality 20 \\
        --min-reads 3 \\
        --max-reads-per-strand 0
    fi
    
    # Filter consensus reads
    fgbio FilterConsensusReads \\
        --input ${sample_id}_consensus_raw.bam \\
        --output ${sample_id}_consensus_filtered.bam \\
        --ref ${ref_fasta} \\
        --min-reads 3 \\
        --max-read-error-rate 0.05 \\
        --max-base-error-rate 0.1 \\
        --min-base-quality 20 \\
        --max-no-call-fraction 0.1
    
    # Rename to final output
    mv ${sample_id}_consensus_filtered.bam ${sample_id}_consensus.bam
    
    # Index the output BAM
    samtools index ${sample_id}_consensus.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(echo \$(fgbio --version 2>&1) | sed 's/^fgbio //' )
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

// UMI Processing Workflow
workflow UMI {
    take:
    reads            // channel: [sample_id, [R1, R2]]
    umi_schema       // string: UMI schema (e.g., "8M+T")
    enable_duplex    // boolean: enable duplex consensus calling
    ref_fasta        // path: reference FASTA
    
    main:
    // Step 1: Convert FASTQ to unmapped BAM with UMI structure
    FASTQ_TO_BAM(
        reads,
        umi_schema
    )
    
    // Step 2: Extract UMIs from unmapped BAM
    EXTRACT_UMIS(
        FASTQ_TO_BAM.out.bam,
        umi_schema
    )
    
    // Step 3: Align UMI-extracted reads
    ALIGN_UMIS(
        EXTRACT_UMIS.out.bam,
        ref_fasta
    )
    
    // Step 4: Group reads by UMI
    GROUP_READS_BY_UMI(ALIGN_UMIS.out.bam)
    
    // Step 5: Call consensus reads
    CALL_CONSENSUS_READS(
        GROUP_READS_BY_UMI.out.bam,
        enable_duplex,
        ref_fasta
    )
    
    // Combine versions
    ch_versions = FASTQ_TO_BAM.out.versions
        .mix(EXTRACT_UMIS.out.versions)
        .mix(ALIGN_UMIS.out.versions)
        .mix(GROUP_READS_BY_UMI.out.versions)
        .mix(CALL_CONSENSUS_READS.out.versions)
    
    emit:
    consensus_bam = CALL_CONSENSUS_READS.out.bam
    family_sizes = GROUP_READS_BY_UMI.out.family_sizes
    versions = ch_versions
}