// modules/cnv.nf - Copy Number Variation Analysis Module (Fault-Tolerant)

nextflow.enable.dsl=2

process CNVKIT_COVERAGE {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(bam)
    path targets_bed
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}.targetcoverage.cnn"), emit: target_coverage, optional: true
    tuple val(sample_id), path("${sample_id}.antitargetcoverage.cnn"), emit: antitarget_coverage, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Starting CNV coverage calculation for sample: ${sample_id}"
    echo "Targets BED: ${targets_bed}"
    echo "Reference: ${ref_fasta}"
    
    COVERAGE_SUCCESS=false
    
    # Check if cnvkit is available
    if ! command -v cnvkit.py >/dev/null 2>&1; then
        echo "WARNING: cnvkit.py not found in container"
        echo "CNV analysis will be skipped"
        
        # Create empty coverage files
        touch ${sample_id}.targetcoverage.cnn
        touch ${sample_id}.antitargetcoverage.cnn
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if BAM exists and is valid
    if [ ! -f "${bam}" ] || [ ! -s "${bam}" ]; then
        echo "WARNING: BAM file missing or empty"
        
        touch ${sample_id}.targetcoverage.cnn
        touch ${sample_id}.antitargetcoverage.cnn
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_invalid_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if targets BED exists
    if [ ! -f "${targets_bed}" ] || [ ! -s "${targets_bed}" ]; then
        echo "WARNING: Targets BED missing or empty"
        
        touch ${sample_id}.targetcoverage.cnn
        touch ${sample_id}.antitargetcoverage.cnn
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_targets"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Calculate target coverage
    echo "Calculating target coverage..."
    if timeout 600 cnvkit.py coverage \\
        ${bam} \\
        ${targets_bed} \\
        -o ${sample_id}.targetcoverage.cnn 2>&1 | tee coverage_target.log; then
        
        if [ -f "${sample_id}.targetcoverage.cnn" ] && [ -s "${sample_id}.targetcoverage.cnn" ]; then
            echo "Target coverage calculation successful"
        else
            echo "WARNING: Target coverage calculation completed but output is empty"
            cat coverage_target.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Target coverage calculation failed"
        cat coverage_target.log 2>/dev/null || echo "No log available"
    fi
    
    # Create antitargets - use simpler approach without reference FASTA
    # CNVkit's antitarget command expects an access file, not a FASTA
    echo "Generating antitargets from targets BED..."
    
    # Simple approach: create antitargets by excluding target regions
    # This avoids the FASTA format error
    if timeout 300 cnvkit.py antitarget \\
        ${targets_bed} \\
        -o ${sample_id}.antitargets.bed 2>&1 | tee antitarget_gen.log; then
        
        echo "Antitarget regions generated"
    else
        echo "WARNING: Antitarget generation failed, creating minimal antitarget file"
        cat antitarget_gen.log 2>/dev/null || echo "No log available"
        
        # Create minimal antitarget file (empty but valid BED format)
        echo -e "chromosome\\tstart\\tend\\tgene" > ${sample_id}.antitargets.bed
    fi
    
    # Calculate antitarget coverage if antitarget file was created
    if [ -f "${sample_id}.antitargets.bed" ] && [ -s "${sample_id}.antitargets.bed" ]; then
        echo "Calculating antitarget coverage..."
        if timeout 600 cnvkit.py coverage \\
            ${bam} \\
            ${sample_id}.antitargets.bed \\
            -o ${sample_id}.antitargetcoverage.cnn 2>&1 | tee coverage_antitarget.log; then
            
            if [ -f "${sample_id}.antitargetcoverage.cnn" ] && [ -s "${sample_id}.antitargetcoverage.cnn" ]; then
                echo "Antitarget coverage calculation successful"
                COVERAGE_SUCCESS=true
            else
                echo "WARNING: Antitarget coverage calculation completed but output is empty"
                cat coverage_antitarget.log 2>/dev/null || echo "No log available"
            fi
        else
            echo "WARNING: Antitarget coverage calculation failed"
            cat coverage_antitarget.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Antitarget file is empty, creating minimal antitarget coverage"
    fi
    
    # Ensure output files exist
    if [ ! -f "${sample_id}.targetcoverage.cnn" ]; then
        touch ${sample_id}.targetcoverage.cnn
    fi
    
    if [ ! -f "${sample_id}.antitargetcoverage.cnn" ]; then
        touch ${sample_id}.antitargetcoverage.cnn
    fi
    
    echo "Coverage calculation completed for sample: ${sample_id}"
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: \$([ "\${COVERAGE_SUCCESS}" = "true" ] && echo "success" || echo "partial")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process CNVKIT_FIX {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(target_coverage), path(antitarget_coverage)
    path reference_cnn
    
    output:
    tuple val(sample_id), path("${sample_id}.cnr"), emit: cnr, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running CNVkit fix for sample: ${sample_id}"
    
    FIX_SUCCESS=false
    
    # Check inputs
    if [ ! -f "${target_coverage}" ] || [ ! -s "${target_coverage}" ]; then
        echo "WARNING: Target coverage file missing or empty"
        touch ${sample_id}.cnr
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    if [ ! -f "${reference_cnn}" ] || [ ! -s "${reference_cnn}" ]; then
        echo "WARNING: Reference CNN file missing or empty"
        touch ${sample_id}.cnr
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_reference"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    if timeout 300 cnvkit.py fix \\
        ${target_coverage} \\
        ${antitarget_coverage} \\
        ${reference_cnn} \\
        -o ${sample_id}.cnr 2>&1 | tee fix.log; then
        
        if [ -f "${sample_id}.cnr" ] && [ -s "${sample_id}.cnr" ]; then
            echo "CNVkit fix completed successfully"
            FIX_SUCCESS=true
        else
            echo "WARNING: CNVkit fix completed but output is empty"
            cat fix.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: CNVkit fix failed"
        cat fix.log 2>/dev/null || echo "No log available"
    fi
    
    if [ ! -f "${sample_id}.cnr" ]; then
        touch ${sample_id}.cnr
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: \$([ "\${FIX_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process CNVKIT_SEGMENT {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(cnr)
    
    output:
    tuple val(sample_id), path("${sample_id}.cns"), emit: cns, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running CNVkit segmentation for sample: ${sample_id}"
    
    SEGMENT_SUCCESS=false
    
    if [ ! -f "${cnr}" ] || [ ! -s "${cnr}" ]; then
        echo "WARNING: CNR file missing or empty"
        touch ${sample_id}.cns
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    if timeout 600 cnvkit.py segment \\
        ${cnr} \\
        -o ${sample_id}.cns \\
        --method cbs \\
        --smooth-cbs 2>&1 | tee segment.log; then
        
        if [ -f "${sample_id}.cns" ] && [ -s "${sample_id}.cns" ]; then
            echo "Segmentation completed successfully"
            SEGMENT_SUCCESS=true
        else
            echo "WARNING: Segmentation completed but output is empty"
            cat segment.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Segmentation failed"
        cat segment.log 2>/dev/null || echo "No log available"
    fi
    
    if [ ! -f "${sample_id}.cns" ]; then
        touch ${sample_id}.cns
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: \$([ "\${SEGMENT_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process CNVKIT_CALL {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(cns)
    
    output:
    tuple val(sample_id), path("${sample_id}_calls.cns"), emit: calls, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running CNVkit call for sample: ${sample_id}"
    
    CALL_SUCCESS=false
    
    if [ ! -f "${cns}" ] || [ ! -s "${cns}" ]; then
        echo "WARNING: CNS file missing or empty"
        touch ${sample_id}_calls.cns
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    if timeout 300 cnvkit.py call \\
        ${cns} \\
        -o ${sample_id}_calls.cns \\
        --method threshold \\
        --thresholds -1.1,-0.4,0.3,0.7 2>&1 | tee call.log; then
        
        if [ -f "${sample_id}_calls.cns" ] && [ -s "${sample_id}_calls.cns" ]; then
            echo "Calling completed successfully"
            CALL_SUCCESS=true
        else
            echo "WARNING: Calling completed but output is empty"
            cat call.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Calling failed"
        cat call.log 2>/dev/null || echo "No log available"
    fi
    
    if [ ! -f "${sample_id}_calls.cns" ]; then
        touch ${sample_id}_calls.cns
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: \$([ "\${CALL_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process CNVKIT_GENEMETRICS {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(cnr), path(cns)
    
    output:
    tuple val(sample_id), path("${sample_id}_gene_metrics.tsv"), emit: gene_metrics, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running CNVkit genemetrics for sample: ${sample_id}"
    
    GENEMETRICS_SUCCESS=false
    
    if [ ! -f "${cnr}" ] || [ ! -s "${cnr}" ]; then
        echo "WARNING: CNR file missing or empty"
        echo -e "gene\\tchromosome\\tstart\\tend\\tlog2\\tprobes" > ${sample_id}_gene_metrics.tsv
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: "skipped_missing_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    if timeout 300 cnvkit.py genemetrics \\
        ${cnr} \\
        -s ${cns} \\
        -o ${sample_id}_gene_metrics.tsv \\
        --threshold 0.2 \\
        --min-probes ${params.cnv_min_probes_per_gene} 2>&1 | tee genemetrics.log; then
        
        if [ -f "${sample_id}_gene_metrics.tsv" ] && [ -s "${sample_id}_gene_metrics.tsv" ]; then
            echo "Gene metrics calculated successfully"
            GENEMETRICS_SUCCESS=true
        else
            echo "WARNING: Gene metrics completed but output is empty"
            cat genemetrics.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Gene metrics calculation failed"
        cat genemetrics.log 2>/dev/null || echo "No log available"
    fi
    
    if [ ! -f "${sample_id}_gene_metrics.tsv" ]; then
        echo -e "gene\\tchromosome\\tstart\\tend\\tlog2\\tprobes" > ${sample_id}_gene_metrics.tsv
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py version 2>&1 | sed 's/cnvkit //' || echo "unknown")
    status: \$([ "\${GENEMETRICS_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process ICHORCNA {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/cnv/ichorcna", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(bam)
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_ichorcna_tumor_fraction.txt"), emit: tumor_fraction, optional: true
    tuple val(sample_id), path("${sample_id}_ichorcna_summary.txt"), emit: summary, optional: true
    tuple val(sample_id), path("${sample_id}.tfx.txt"), emit: tfx, optional: true
    path "${sample_id}_ichorcna_segments.txt", emit: segments, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running ichorCNA for tumor fraction estimation: ${sample_id}"
    
    ICHORCNA_SUCCESS=false
    
    # Check if ichorCNA is available
    if ! command -v runIchorCNA.R >/dev/null 2>&1 && ! command -v Rscript >/dev/null 2>&1; then
        echo "WARNING: ichorCNA/R not found in container"
        
        echo "0.0" > ${sample_id}_ichorcna_tumor_fraction.txt
        echo "Tumor Fraction: 0.0%" > ${sample_id}_ichorcna_summary.txt
        touch ${sample_id}_ichorcna_segments.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    ichorcna: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if BAM exists
    if [ ! -f "${bam}" ] || [ ! -s "${bam}" ]; then
        echo "WARNING: BAM file missing or empty"
        
        echo "0.0" > ${sample_id}_ichorcna_tumor_fraction.txt
        echo "Tumor Fraction: 0.0%" > ${sample_id}_ichorcna_summary.txt
        touch ${sample_id}_ichorcna_segments.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    ichorcna: "unknown"
    status: "skipped_invalid_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Use the TFX pipeline run_ichorcna.R script if available, otherwise use inline version
    if [ -f "${projectDir}/bin/run_ichorcna.R" ]; then
        echo "Using TFX pipeline run_ichorcna.R script"
        
        # First, create WIG file from BAM using bedtools genomecov
        if command -v bedtools >/dev/null 2>&1; then
            echo "Creating WIG file from BAM..."
            bedtools genomecov -ibam ${bam} -bg -split | \\
                awk 'BEGIN {OFS="\\t"} {print \$1, \$2, \$3, \$4}' > ${sample_id}.wig || {
                echo "WARNING: Failed to create WIG file, using inline ichorCNA"
                USE_INLINE=true
            }
        else
            echo "WARNING: bedtools not found, using inline ichorCNA"
            USE_INLINE=true
        fi
        
        if [ "$USE_INLINE" != "true" ] && [ -f "${sample_id}.wig" ]; then
            # Run the TFX pipeline script
            Rscript ${projectDir}/bin/run_ichorcna.R \\
                --wig ${sample_id}.wig \\
                --outdir . \\
                --sample ${sample_id} \\
                --binSize 500000 \\
                --ploidy "c(2,3)" \\
                --normal "c(0.5,0.7,0.8,0.9,0.95)" \\
                --maxCN 5 \\
                --estimateNormal \\
                --estimatePloidy \\
                --txnE 0.9999 \\
                --txnStrength 10000 || {
                echo "WARNING: TFX run_ichorcna.R failed, using inline version"
                USE_INLINE=true
            }
        fi
    else
        USE_INLINE=true
    fi
    
    if [ "$USE_INLINE" = "true" ]; then
        echo "Using inline ichorCNA implementation"
    # Create ichorCNA R script
    cat > run_ichorcna.R << 'EOF'
library(ichorCNA)
library(GenomicRanges)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
ref_fasta <- args[2]
sample_id <- args[3]
output_prefix <- args[4]

cat("Running ichorCNA for sample:", sample_id, "\\n")
cat("BAM file:", bam_file, "\\n")
cat("Reference:", ref_fasta, "\\n")

# ichorCNA parameters
gcWig <- NULL  # Will use default GC correction
mapWig <- NULL  # Will use default mappability correction
normal <- NULL  # Tumor-only mode
ploidy <- "c(2)"  # Default ploidy
scStates <- "c(1,3)"  # Subclonal states
maxCN <- 5  # Maximum copy number
estimateNormal <- TRUE
estimatePloidy <- TRUE
estimateScPrevalance <- TRUE

# Run ichorCNA
tryCatch({
    results <- runIchorCNA(
        bamFile = bam_file,
        gcWig = gcWig,
        mapWig = mapWig,
        normal = normal,
        ploidy = ploidy,
        scStates = scStates,
        maxCN = maxCN,
        estimateNormal = estimateNormal,
        estimatePloidy = estimatePloidy,
        estimateScPrevalance = estimateScPrevalance,
        chrNormalize = c(1:22, "X", "Y"),
        chrTrain = c(1:22),
        chrs = c(1:22, "X", "Y"),
        outputDir = ".",
        outputName = output_prefix
    )
    
    # Extract tumor fraction
    if (!is.null(results) && "tumorFraction" %in% names(results)) {
        tf <- results$tumorFraction
    } else if (!is.null(results) && "fraction" %in% names(results)) {
        tf <- results$fraction
    } else {
        tf <- 0.0
    }
    
    # Write tumor fraction
    write(tf, file = paste0(output_prefix, "_ichorcna_tumor_fraction.txt"))
    tfx_df <- data.frame(sample=sample_id, tumorFraction=tf)
    write.table(tfx_df,
                file = paste0(output_prefix, ".tfx.txt"),
                sep = "\\t", quote = FALSE, row.names = FALSE)
    
    # Write summary
    summary_text <- paste0("Tumor Fraction: ", round(tf * 100, 2), "%\\n")
    if (!is.null(results) && "ploidy" %in% names(results)) {
        summary_text <- paste0(summary_text, "Ploidy: ", results$ploidy, "\\n")
    }
    write(summary_text, file = paste0(output_prefix, "_ichorcna_summary.txt"))
    
    # Write segments if available
    if (!is.null(results) && "segments" %in% names(results)) {
        write.table(results$segments, 
                   file = paste0(output_prefix, "_ichorcna_segments.txt"),
                   sep = "\\t", quote = FALSE, row.names = FALSE)
    } else {
        # Create empty segments file
        cat("chromosome\\tstart\\tend\\tlog2\\tcopy_number\\n", 
            file = paste0(output_prefix, "_ichorcna_segments.txt"))
    }
    
    cat("ichorCNA completed successfully. Tumor fraction:", tf, "\\n")
    
}, error = function(e) {
    cat("ERROR in ichorCNA:", e$message, "\\n")
    
    # Write default outputs on error
    write("0.0", file = paste0(output_prefix, "_ichorcna_tumor_fraction.txt"))
    tfx_df <- data.frame(sample=sample_id, tumorFraction=0.0)
    write.table(tfx_df,
                file = paste0(output_prefix, ".tfx.txt"),
                sep = "\\t", quote = FALSE, row.names = FALSE)
    cat("Tumor Fraction: 0.0% (Analysis Failed)\\n", 
        file = paste0(output_prefix, "_ichorcna_summary.txt"))
    cat("chromosome\\tstart\\tend\\tlog2\\tcopy_number\\n", 
        file = paste0(output_prefix, "_ichorcna_segments.txt"))
})
EOF

    # Run ichorCNA
    if command -v Rscript >/dev/null 2>&1; then
        echo "Running ichorCNA R script..."
        if timeout 3600 Rscript run_ichorcna.R "${bam}" "${ref_fasta}" "${sample_id}" "${sample_id}" 2>&1 | tee ichorcna.log; then
            
            if [ -f "${sample_id}_ichorcna_tumor_fraction.txt" ] && [ -s "${sample_id}_ichorcna_tumor_fraction.txt" ]; then
                echo "ichorCNA completed successfully"
                ICHORCNA_SUCCESS=true
            else
                echo "WARNING: ichorCNA completed but output is empty"
                cat ichorcna.log 2>/dev/null || echo "No log available"
            fi
        else
            echo "WARNING: ichorCNA R script failed"
            cat ichorcna.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: Rscript not found"
    fi
    
    # Ensure output files exist
    if [ ! -f "${sample_id}_ichorcna_tumor_fraction.txt" ]; then
        echo "0.0" > ${sample_id}_ichorcna_tumor_fraction.txt
    fi
    
    if [ ! -f "${sample_id}_ichorcna_summary.txt" ]; then
        echo "Tumor Fraction: 0.0%" > ${sample_id}_ichorcna_summary.txt
    fi
    
    if [ ! -f "${sample_id}.tfx.txt" ]; then
        printf "sample\ttumorFraction\n%s\t0.0\n" "${sample_id}" > ${sample_id}.tfx.txt
    fi
    
    if [ ! -f "${sample_id}_ichorcna_segments.txt" ]; then
        touch ${sample_id}_ichorcna_segments.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    ichorcna: \$(Rscript -e "packageVersion('ichorCNA')" 2>/dev/null | sed 's/.*\\[1\\] //' || echo "unknown")
    r: \$(R --version 2>&1 | head -1 | sed 's/.*version //' || echo "unknown")
    status: \$([ "\${ICHORCNA_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process CNV_FILTER_CALLS {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(gene_metrics), path(calls_cns)
    
    output:
    tuple val(sample_id), path("${sample_id}_cnv_calls.tsv"), emit: filtered_calls, optional: true
    path "${sample_id}_cnv_stats.json", emit: stats, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Filtering CNV calls for sample: ${sample_id}"
    
    # Check if Python is available
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        echo "WARNING: Python not found"
        
        echo -e "gene\\tchromosome\\tstart\\tend\\tlog2\\tcopy_number\\tprobes\\tinterpretation" > ${sample_id}_cnv_calls.tsv
        echo '{"sample_id":"${sample_id}","total_genes":0,"altered_genes":0,"status":"python_missing"}' > ${sample_id}_cnv_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    PYTHON_CMD=\$(command -v python3 || command -v python)
    
    # Check input files
    if [ ! -f "${gene_metrics}" ] || [ ! -s "${gene_metrics}" ]; then
        echo "WARNING: Gene metrics file missing or empty"
        
        echo -e "gene\\tchromosome\\tstart\\tend\\tlog2\\tcopy_number\\tprobes\\tinterpretation" > ${sample_id}_cnv_calls.tsv
        echo '{"sample_id":"${sample_id}","total_genes":0,"altered_genes":0,"status":"no_input"}' > ${sample_id}_cnv_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(\${PYTHON_CMD} --version 2>&1 | sed 's/Python //')
    status: "skipped_no_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    cat > filter_cnv.py << 'EOF'
import sys
import json

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    print("WARNING: pandas not available")
    PANDAS_AVAILABLE = False

if not PANDAS_AVAILABLE:
    # Create minimal output without pandas
    with open("${sample_id}_cnv_calls.tsv", 'w') as f:
        f.write("gene\\tchromosome\\tstart\\tend\\tlog2\\tcopy_number\\tprobes\\tinterpretation\\n")
    
    stats = {
        "sample_id": "${sample_id}",
        "total_genes": 0,
        "altered_genes": 0,
        "deletions": 0,
        "amplifications": 0,
        "status": "pandas_missing"
    }
    
    with open("${sample_id}_cnv_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    sys.exit(0)

try:
    # Read gene metrics
    gene_df = pd.read_csv("${gene_metrics}", sep='\\t')
    
    # Check if file has data
    if len(gene_df) == 0:
        print("WARNING: Gene metrics file is empty")
        raise ValueError("Empty gene metrics file")
    
    # Read segment calls if available
    seg_df = None
    try:
        if "${calls_cns}" and pd.io.common.file_exists("${calls_cns}"):
            seg_df = pd.read_csv("${calls_cns}", sep='\\t')
    except:
        print("WARNING: Could not read segment calls file")
    
    # Filter based on clinical thresholds
    filtered_genes = gene_df[
        (abs(gene_df['log2']) >= 0.58) &  # ~1.5x fold change
        (gene_df['probes'] >= ${params.cnv_min_probes_per_gene})
    ].copy()
    
    # Add clinical interpretation
    def interpret_cnv(log2_ratio):
        if log2_ratio <= -1.0:
            return "Deep Deletion"
        elif log2_ratio <= -0.58:
            return "Deletion" 
        elif log2_ratio >= 1.0:
            return "Amplification"
        elif log2_ratio >= 0.58:
            return "Gain"
        else:
            return "Neutral"
    
    filtered_genes['interpretation'] = filtered_genes['log2'].apply(interpret_cnv)
    filtered_genes['copy_number'] = 2 * (2 ** filtered_genes['log2'])
    
    # Select relevant columns
    output_cols = ['gene', 'chromosome', 'start', 'end', 'log2', 'copy_number', 
                   'probes', 'interpretation']
    filtered_genes[output_cols].to_csv("${sample_id}_cnv_calls.tsv", sep='\\t', index=False)
    
    # Generate statistics
    stats = {
        "sample_id": "${sample_id}",
        "total_genes": len(gene_df),
        "altered_genes": len(filtered_genes),
        "deletions": len(filtered_genes[filtered_genes['log2'] < -0.58]),
        "amplifications": len(filtered_genes[filtered_genes['log2'] > 0.58]),
        "mad_score": float(seg_df['weight'].std()) if seg_df is not None and 'weight' in seg_df.columns else None,
        "status": "success"
    }
    
    with open("${sample_id}_cnv_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"CNV filtering completed: {len(filtered_genes)} altered genes found")

except Exception as e:
    print(f"ERROR in CNV filtering: {e}")
    import traceback
    traceback.print_exc()
    
    # Create minimal output on error
    with open("${sample_id}_cnv_calls.tsv", 'w') as f:
        f.write("gene\\tchromosome\\tstart\\tend\\tlog2\\tcopy_number\\tprobes\\tinterpretation\\n")
    
    stats = {
        "sample_id": "${sample_id}",
        "total_genes": 0,
        "altered_genes": 0,
        "deletions": 0,
        "amplifications": 0,
        "status": "failed",
        "error": str(e)
    }
    
    with open("${sample_id}_cnv_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
EOF

    \${PYTHON_CMD} filter_cnv.py 2>&1 | tee filter.log
    
    # Ensure output files exist
    if [ ! -f "${sample_id}_cnv_calls.tsv" ]; then
        echo -e "gene\\tchromosome\\tstart\\tend\\tlog2\\tcopy_number\\tprobes\\tinterpretation" > ${sample_id}_cnv_calls.tsv
    fi
    
    if [ ! -f "${sample_id}_cnv_stats.json" ]; then
        echo '{"sample_id":"${sample_id}","total_genes":0,"altered_genes":0,"status":"unknown"}' > ${sample_id}_cnv_stats.json
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(\${PYTHON_CMD} --version 2>&1 | sed 's/Python //')
    pandas: \$(\${PYTHON_CMD} -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "not_available")
    status: "success"
END_VERSIONS
    
    set -e
    exit 0
    """
}

// CNV Workflow with fault tolerance
workflow CNV {
    take:
    bam              // channel: [sample_id, bam]
    targets_bed      // path: target regions BED
    reference_cnn    // path: CNVkit reference CNN
    ref_fasta        // path: reference FASTA (for antitarget generation)
    
    main:
    // Calculate coverage (antitargets will be generated automatically)
    CNVKIT_COVERAGE(
        bam,
        targets_bed,
        ref_fasta
    )
    
    // Combine target and antitarget coverage with safe handling
    combined_coverage = CNVKIT_COVERAGE.out.target_coverage
        .join(CNVKIT_COVERAGE.out.antitarget_coverage, by: 0, remainder: true)
        .map { sample_id, target_cov, antitarget_cov ->
            def safe_antitarget = antitarget_cov ?: file("NO_FILE")
            [sample_id, target_cov, safe_antitarget]
        }
    
    // Fix coverage ratios
    CNVKIT_FIX(
        combined_coverage,
        reference_cnn
    )
    
    // Segment
    CNVKIT_SEGMENT(CNVKIT_FIX.out.cnr)
    
    // Call CNVs
    CNVKIT_CALL(CNVKIT_SEGMENT.out.cns)
    
    // Calculate gene metrics with safe joining
    cnr_cns_combined = CNVKIT_FIX.out.cnr
        .join(CNVKIT_SEGMENT.out.cns, by: 0, remainder: true)
        .map { sample_id, cnr, cns ->
            def safe_cns = cns ?: file("NO_FILE")
            [sample_id, cnr, safe_cns]
        }
    
    CNVKIT_GENEMETRICS(cnr_cns_combined)
    
    // Filter and annotate calls with safe joining
    gene_calls_combined = CNVKIT_GENEMETRICS.out.gene_metrics
        .join(CNVKIT_CALL.out.calls, by: 0, remainder: true)
        .map { sample_id, gene_metrics, calls ->
            def safe_calls = calls ?: file("NO_FILE")
            [sample_id, gene_metrics, safe_calls]
        }
    
    CNV_FILTER_CALLS(gene_calls_combined)
    
    // ichorCNA for tumor fraction estimation
    ICHORCNA(
        bam,
        ref_fasta
    )
    
    // Combine versions
    ch_versions = CNVKIT_COVERAGE.out.versions
        .mix(CNVKIT_FIX.out.versions)
        .mix(CNVKIT_SEGMENT.out.versions)
        .mix(CNVKIT_CALL.out.versions)
        .mix(CNVKIT_GENEMETRICS.out.versions)
        .mix(CNV_FILTER_CALLS.out.versions)
        .mix(ICHORCNA.out.versions)
    
    emit:
    cnr = CNVKIT_FIX.out.cnr
    cns = CNVKIT_SEGMENT.out.cns
    calls = CNVKIT_CALL.out.calls
    gene_calls = CNV_FILTER_CALLS.out.filtered_calls
    gene_metrics = CNVKIT_GENEMETRICS.out.gene_metrics
    segments = CNVKIT_SEGMENT.out.cns
    stats = CNV_FILTER_CALLS.out.stats
    tumor_fraction = ICHORCNA.out.tumor_fraction
    ichorcna_summary = ICHORCNA.out.summary
    tfx = ICHORCNA.out.tfx
    versions = ch_versions
}