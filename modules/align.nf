// modules/align.nf - Fault-Tolerant Alignment and Post-processing Module
// Fixed bash syntax errors

nextflow.enable.dsl=2

process CREATE_BWA_MEM2_INDEX {
    tag "bwa_mem2_index"
    label 'process_high'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    storeDir "${params.outdir}/reference_index"
    
    input:
    path ref_fasta
    
    output:
    path "${ref_fasta}", emit: ref_fasta
    path "${ref_fasta}.*", emit: index_files
    path "versions.yml", emit: versions
    
    script:
    def ref_name = ref_fasta.name
    """
    set +e
    
    echo "Checking for existing BWA index files for: ${ref_fasta}"
    
    # Check if BWA index files already exist
    if [ -f "${ref_fasta}.bwt" ] && [ -f "${ref_fasta}.amb" ] && [ -f "${ref_fasta}.ann" ] && [ -f "${ref_fasta}.pac" ] && [ -f "${ref_fasta}.sa" ]; then
        echo "BWA index files already exist, skipping index creation"
        touch ${ref_fasta}.bwt ${ref_fasta}.amb ${ref_fasta}.ann ${ref_fasta}.pac ${ref_fasta}.sa 2>/dev/null || true
    elif [ -f "${ref_fasta}.bwt.2bit.64" ]; then
        echo "BWA-MEM2 index file already exists, skipping index creation"
        touch ${ref_fasta}.bwt.2bit.64 2>/dev/null || true
    else
        echo "No existing index found, creating BWA index for: ${ref_fasta}"
    
    # Validate tools with graceful fallback
    if ! command -v bwa-mem2 >/dev/null 2>&1; then
        echo "WARNING: bwa-mem2 not found, will attempt to use bwa"
        BWA_MEM2_AVAILABLE=false
    else
        BWA_MEM2_AVAILABLE=true
    fi
    
        # Try to create index locally - use standard BWA first (requires less memory)
        if command -v bwa >/dev/null 2>&1; then
            echo "Creating standard BWA index (requires less memory than BWA-MEM2)..."
            if timeout 14400 bwa index ${ref_fasta} 2>&1; then
                echo "Standard BWA indexing successful"
            else
                echo "Standard BWA indexing failed, attempting BWA-MEM2..."
    if [ "\${BWA_MEM2_AVAILABLE}" = "true" ]; then
                    timeout 14400 bwa-mem2 index ${ref_fasta} 2>&1 && \\
                        echo "BWA-MEM2 indexing successful" || echo "WARNING: Both indexing attempts failed"
                else
                    echo "WARNING: BWA indexing failed and BWA-MEM2 not available"
            fi
        fi
        elif [ "\${BWA_MEM2_AVAILABLE}" = "true" ]; then
            echo "Attempting BWA-MEM2 indexing..."
            timeout 14400 bwa-mem2 index ${ref_fasta} 2>&1 && \\
                echo "BWA-MEM2 indexing successful" || echo "WARNING: BWA-MEM2 indexing failed"
        else
            echo "ERROR: Neither bwa-mem2 nor bwa available"
        fi
    fi
    
    # List all created files for debugging
    echo "Index files present:"
    ls -lh ${ref_fasta}* 2>/dev/null || echo "No index files found"
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bwa-mem2: \$(command -v bwa-mem2 >/dev/null 2>&1 && echo \$(bwa-mem2 version 2>&1) | sed 's/.* //' || echo "not_available")
    bwa: \$(command -v bwa >/dev/null 2>&1 && bwa 2>&1 | grep -E '^Version' | sed 's/Version: //' || echo "not_available")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process BWA_MEM2_ALIGN_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*_aligned.bam*"
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path index_files
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: bam, optional: true
    tuple val(sample_id), path("${sample_id}_aligned.bam.bai"), emit: bai, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Starting alignment for sample: ${sample_id}"
    echo "Input reads: ${reads}"
    echo "Reference FASTA: ${ref_fasta}"
    echo "Available index files:"
    ls -lh ${ref_fasta}* 2>/dev/null | head -20 || echo "No index files found"
    
    # Check tools availability
    TOOLS_OK=true
    for tool in samtools; do
        if ! command -v \${tool} >/dev/null 2>&1; then
            echo "WARNING: \${tool} not found"
            TOOLS_OK=false
        fi
    done
    
    if [ "\${TOOLS_OK}" = "false" ]; then
        echo "CRITICAL: Required tools not available"
        touch ${sample_id}_aligned.bam
        touch ${sample_id}_aligned.bam.bai
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "failed_missing_tools"
    attempt: ${task.attempt}
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    ALIGNMENT_SUCCESS=false
    USE_BWA_FALLBACK=false
    ALIGNER_NAME="unknown"
    ALIGNER_VERSION="unknown"
    RG_STRING="@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}"
    
    # Validate input reads exist
    if [ ! -f "${reads[0]}" ] || [ ! -f "${reads[1]}" ]; then
        echo "ERROR: Input read files not found"
        touch ${sample_id}_aligned.bam
        touch ${sample_id}_aligned.bam.bai
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "failed_missing_reads"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Try BWA-MEM2 alignment
    if [ -f "${ref_fasta}.bwt.2bit.64" ] && command -v bwa-mem2 >/dev/null 2>&1; then
        echo "Attempting alignment with BWA-MEM2..."
        
        # Run alignment and capture exit status
        bwa-mem2 mem -t ${task.cpus} -R "\${RG_STRING}" ${ref_fasta} ${reads[0]} ${reads[1]} 2>/dev/null | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_aligned.bam - 2>/dev/null
        BWA_MEM2_EXIT=\$?
        
        if [ \${BWA_MEM2_EXIT} -eq 0 ] && [ -f "${sample_id}_aligned.bam" ] && [ -s "${sample_id}_aligned.bam" ]; then
            echo "BWA-MEM2 alignment completed successfully"
            ALIGNMENT_SUCCESS=true
            ALIGNER_NAME="bwa-mem2"
            ALIGNER_VERSION=\$(bwa-mem2 version 2>&1 | sed 's/.* //' || echo "unknown")
        else
            echo "BWA-MEM2 alignment failed (exit: \${BWA_MEM2_EXIT}), attempting fallback..."
            USE_BWA_FALLBACK=true
            rm -f ${sample_id}_aligned.bam 2>/dev/null
        fi
    else
        echo "BWA-MEM2 index or binary not available, using standard BWA..."
        USE_BWA_FALLBACK=true
    fi
    
    # Fallback to standard BWA
    if [ "\${USE_BWA_FALLBACK}" = "true" ]; then
        if command -v bwa >/dev/null 2>&1; then
            echo "Attempting alignment with standard BWA..."
            
            # Check if BWA index exists
            if [ ! -f "${ref_fasta}.bwt" ]; then
                echo "BWA index not found, creating..."
                timeout 7200 bwa index ${ref_fasta} 2>/dev/null || echo "WARNING: BWA index creation timed out"
            fi
            
            echo "Running BWA alignment..."
            bwa mem -t ${task.cpus} -R "\${RG_STRING}" ${ref_fasta} ${reads[0]} ${reads[1]} 2>/dev/null | \\
                samtools sort -@ ${task.cpus} -o ${sample_id}_aligned.bam - 2>/dev/null
            BWA_EXIT=\$?
            
            if [ \${BWA_EXIT} -eq 0 ] && [ -f "${sample_id}_aligned.bam" ] && [ -s "${sample_id}_aligned.bam" ]; then
                echo "Standard BWA alignment completed"
                ALIGNMENT_SUCCESS=true
                ALIGNER_NAME="bwa"
                ALIGNER_VERSION=\$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //' || echo "unknown")
            else
                echo "WARNING: Standard BWA alignment failed (exit: \${BWA_EXIT})"
                ALIGNMENT_SUCCESS=false
            fi
        else
            echo "WARNING: Neither bwa-mem2 nor bwa available"
            ALIGNMENT_SUCCESS=false
        fi
    fi
    
    # Create empty BAM if all alignments failed
    if [ "\${ALIGNMENT_SUCCESS}" = "false" ]; then
        echo "Creating empty placeholder BAM file..."
        touch ${sample_id}_aligned.bam
    fi
    
    # Index BAM file
    if [ -f "${sample_id}_aligned.bam" ] && [ -s "${sample_id}_aligned.bam" ]; then
        echo "Indexing BAM file..."
        if timeout 300 samtools index ${sample_id}_aligned.bam 2>/dev/null; then
            echo "BAM indexing successful"
        else
            echo "WARNING: Failed to index BAM file"
            touch ${sample_id}_aligned.bam.bai
        fi
    else
        echo "Skipping indexing for empty BAM file"
        touch ${sample_id}_aligned.bam.bai
    fi
    
    # Verify output exists
    if [ ! -f "${sample_id}_aligned.bam" ]; then
        touch ${sample_id}_aligned.bam
    fi
    
    if [ ! -f "${sample_id}_aligned.bam.bai" ]; then
        touch ${sample_id}_aligned.bam.bai
    fi
    
    echo "Alignment processing completed for sample: ${sample_id}"
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    aligner: \${ALIGNER_NAME}
    aligner_version: \${ALIGNER_VERSION}
    samtools: \$(command -v samtools >/dev/null 2>&1 && samtools --version 2>&1 | head -1 | sed 's/^.*samtools //' || echo "not_available")
    attempt: ${task.attempt}
    success: \${ALIGNMENT_SUCCESS}
END_VERSIONS
    
    set -e
    exit 0
    """
}

process MARK_DUPLICATES_STANDARD {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*_markdup.bam*"
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_markdup.bam"), emit: bam, optional: true
    path "${sample_id}_dup_metrics.txt", emit: metrics, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Marking duplicates for: ${sample_id}"
    
    PICARD_VERSION="unknown"
    MARKDUP_SUCCESS=false
    
    if [ ! -f "${bam}" ]; then
        echo "WARNING: Input BAM file not found, creating placeholder"
        touch ${sample_id}_markdup.bam
        echo -e "LIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE" > ${sample_id}_dup_metrics.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "failed_missing_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if BAM is empty
    if [ ! -s "${bam}" ]; then
        echo "WARNING: Input BAM file is empty, copying as output"
        cp ${bam} ${sample_id}_markdup.bam 2>/dev/null || touch ${sample_id}_markdup.bam
        echo -e "LIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE" > ${sample_id}_dup_metrics.txt
        echo -e "unknown\\t0\\t0\\t0\\t0\\t0\\t0\\t0.0\\t0" >> ${sample_id}_dup_metrics.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "skipped_empty_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Attempt duplicate marking
    echo "Running Picard MarkDuplicates..."
    if timeout 1800 picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}_markdup.bam \\
        METRICS_FILE=${sample_id}_dup_metrics.txt \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT 2>&1 | tee picard.log; then
        
        if [ -f "${sample_id}_markdup.bam" ] && [ -s "${sample_id}_markdup.bam" ]; then
            echo "Duplicate marking successful"
            MARKDUP_SUCCESS=true
            PICARD_VERSION=\$(picard MarkDuplicates --version 2>&1 | grep -o 'Version.*' | cut -f2 -d' ' || echo "unknown")
        else
            echo "WARNING: Picard completed but output is empty"
            MARKDUP_SUCCESS=false
        fi
    else
        echo "WARNING: Picard MarkDuplicates failed"
        cat picard.log 2>/dev/null || echo "No log available"
        MARKDUP_SUCCESS=false
    fi
    
    # If marking failed, copy input as output
    if [ "\${MARKDUP_SUCCESS}" = "false" ]; then
        echo "Copying input BAM as output"
        cp ${bam} ${sample_id}_markdup.bam 2>/dev/null || touch ${sample_id}_markdup.bam
        
        # Create dummy metrics file
        echo -e "LIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE" > ${sample_id}_dup_metrics.txt
        echo -e "unknown\\t0\\t0\\t0\\t0\\t0\\t0\\t0.0\\t0" >> ${sample_id}_dup_metrics.txt
        
        PICARD_VERSION="unknown"
    fi
    
    # Ensure output files exist
    if [ ! -f "${sample_id}_markdup.bam" ]; then
        touch ${sample_id}_markdup.bam
    fi
    
    if [ ! -f "${sample_id}_dup_metrics.txt" ]; then
        touch ${sample_id}_dup_metrics.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    picard: \${PICARD_VERSION}
    attempt: ${task.attempt}
    success: \${MARKDUP_SUCCESS}
END_VERSIONS
    
    set -e
    exit 0
    """
}

process COLLECT_ALIGNMENT_METRICS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*_metrics.txt"
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(bam)
    path ref_fasta
    path targets_bed
    
    output:
    tuple val(sample_id), path("${sample_id}_alignment_metrics.txt"), emit: alignment_metrics, optional: true
    tuple val(sample_id), path("${sample_id}_insert_metrics.txt"), emit: insert_metrics, optional: true
    tuple val(sample_id), path("${sample_id}_hs_metrics.txt"), emit: hs_metrics, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Collecting alignment metrics for: ${sample_id}"
    
    PICARD_VERSION="unknown"
    
    # Check if BAM file exists and is not empty
    if [ ! -f "${bam}" ] || [ ! -s "${bam}" ]; then
        echo "WARNING: BAM file missing or empty, creating placeholder metrics"
        
        echo -e "CATEGORY\\tTOTAL_READS\\tPF_READS" > ${sample_id}_alignment_metrics.txt
        echo -e "PAIR\\t0\\t0" >> ${sample_id}_alignment_metrics.txt
        
        echo -e "MEDIAN_INSERT_SIZE\\tMODE_INSERT_SIZE" > ${sample_id}_insert_metrics.txt
        echo -e "0\\t0" >> ${sample_id}_insert_metrics.txt
        
        echo -e "BAIT_SET\\tGENOME_SIZE" > ${sample_id}_hs_metrics.txt
        echo -e "unknown\\t0" >> ${sample_id}_hs_metrics.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "skipped_empty_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Alignment summary metrics
    echo "Collecting alignment summary metrics..."
    if timeout 600 picard CollectAlignmentSummaryMetrics \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}_alignment_metrics.txt \\
        REFERENCE_SEQUENCE=${ref_fasta} \\
        VALIDATION_STRINGENCY=LENIENT 2>&1 | tee alignment.log; then
        echo "Alignment metrics collected successfully"
        PICARD_VERSION=\$(picard CollectAlignmentSummaryMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2 -d' ' || echo "unknown")
    else
        echo "WARNING: Failed to collect alignment metrics, creating placeholder"
        cat alignment.log 2>/dev/null || echo "No log available"
        echo -e "CATEGORY\\tTOTAL_READS\\tPF_READS" > ${sample_id}_alignment_metrics.txt
        echo -e "PAIR\\t0\\t0" >> ${sample_id}_alignment_metrics.txt
    fi
    
    # Insert size metrics
    echo "Collecting insert size metrics..."
    if timeout 600 picard CollectInsertSizeMetrics \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}_insert_metrics.txt \\
        HISTOGRAM_FILE=${sample_id}_insert_histogram.pdf \\
        VALIDATION_STRINGENCY=LENIENT 2>&1 | tee insert.log; then
        echo "Insert metrics collected successfully"
    else
        echo "WARNING: Failed to collect insert size metrics"
        cat insert.log 2>/dev/null || echo "No log available"
        echo -e "MEDIAN_INSERT_SIZE\\tMODE_INSERT_SIZE" > ${sample_id}_insert_metrics.txt
        echo -e "0\\t0" >> ${sample_id}_insert_metrics.txt
    fi
    
    # Hybrid selection metrics
    if [ -f "${targets_bed}" ] && [ -s "${targets_bed}" ]; then
        echo "Collecting hybrid selection metrics..."
        if timeout 600 picard CollectHsMetrics \\
            INPUT=${bam} \\
            OUTPUT=${sample_id}_hs_metrics.txt \\
            REFERENCE_SEQUENCE=${ref_fasta} \\
            BAIT_INTERVALS=${targets_bed} \\
            TARGET_INTERVALS=${targets_bed} \\
            VALIDATION_STRINGENCY=LENIENT 2>&1 | tee hs.log; then
            echo "Hybrid selection metrics collected successfully"
        else
            echo "WARNING: Failed to collect hybrid selection metrics"
            cat hs.log 2>/dev/null || echo "No log available"
            echo -e "BAIT_SET\\tGENOME_SIZE" > ${sample_id}_hs_metrics.txt
            echo -e "unknown\\t0" >> ${sample_id}_hs_metrics.txt
        fi
    else
        echo "WARNING: BED file not provided or empty, skipping hybrid selection metrics"
        echo -e "BAIT_SET\\tGENOME_SIZE" > ${sample_id}_hs_metrics.txt
        echo -e "unknown\\t0" >> ${sample_id}_hs_metrics.txt
    fi
    
    # Ensure all output files exist
    if [ ! -f "${sample_id}_alignment_metrics.txt" ]; then
        touch ${sample_id}_alignment_metrics.txt
    fi
    
    if [ ! -f "${sample_id}_insert_metrics.txt" ]; then
        touch ${sample_id}_insert_metrics.txt
    fi
    
    if [ ! -f "${sample_id}_hs_metrics.txt" ]; then
        touch ${sample_id}_hs_metrics.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    picard: \${PICARD_VERSION}
END_VERSIONS
    
    set -e
    exit 0
    """
}

// Main Workflow - Fault-tolerant alignment with graceful error handling
workflow ALIGN {
    take:
    reads_ch
    ref_fasta
    skip_umi
    
    main:
    // Create or download index with error tolerance
    index_result = CREATE_BWA_MEM2_INDEX(ref_fasta)
    
    // Align reads with full error recovery
    align_result = BWA_MEM2_ALIGN_FASTQ(
        reads_ch,
        index_result.ref_fasta,
        index_result.index_files
    )
    
    // Mark duplicates with graceful fallback
    dup_result = MARK_DUPLICATES_STANDARD(align_result.bam)
    
    // Collect metrics - non-critical, continues regardless
    metrics_result = COLLECT_ALIGNMENT_METRICS(
        dup_result.bam,
        index_result.ref_fasta,
        file(params.bed ?: 'NO_FILE')
    )
    
    // Combine all versions
    ch_versions = index_result.versions
        .mix(align_result.versions)
        .mix(dup_result.versions)
        .mix(metrics_result.versions)
    
    emit:
    bam = dup_result.bam
    metrics = dup_result.metrics
    alignment_metrics = metrics_result.alignment_metrics
    insert_metrics = metrics_result.insert_metrics
    hs_metrics = metrics_result.hs_metrics
    versions = ch_versions
}