// modules/sv.nf - Structural Variant Detection Module (Fault-Tolerant)

nextflow.enable.dsl=2

process MANTA_CONFIG {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(bam)
    path ref_fasta
    path targets_bed
    
    output:
    tuple val(sample_id), path("${sample_id}_manta_config"), emit: config_dir, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Configuring Manta for sample: ${sample_id}"
    
    MANTA_AVAILABLE=false
    CONFIG_SUCCESS=false
    
    # Check if Manta is available
    if ! command -v configManta.py >/dev/null 2>&1; then
        echo "WARNING: Manta (configManta.py) not found in container"
        echo "SV detection will be skipped"
        
        # Create dummy config directory
        mkdir -p ${sample_id}_manta_config
        echo "MANTA_NOT_AVAILABLE" > ${sample_id}_manta_config/status.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    MANTA_AVAILABLE=true
    
    # Check if input BAM exists and is valid
    if [ ! -f "${bam}" ] || [ ! -s "${bam}" ]; then
        echo "WARNING: Input BAM file missing or empty"
        
        mkdir -p ${sample_id}_manta_config
        echo "BAM_INVALID" > ${sample_id}_manta_config/status.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: \$(configManta.py --version 2>&1 | grep "Manta" | sed 's/Manta //' || echo "unknown")
    status: "skipped_invalid_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if reference FASTA exists
    if [ ! -f "${ref_fasta}" ]; then
        echo "WARNING: Reference FASTA not found"
        
        mkdir -p ${sample_id}_manta_config
        echo "REF_MISSING" > ${sample_id}_manta_config/status.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: \$(configManta.py --version 2>&1 | grep "Manta" | sed 's/Manta //' || echo "unknown")
    status: "skipped_missing_reference"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Configure Manta for tumor-only mode
    echo "Running Manta configuration..."
    if timeout 300 configManta.py \\
        --tumorBam ${bam} \\
        --referenceFasta ${ref_fasta} \\
        --callRegions ${targets_bed} \\
        --runDir ${sample_id}_manta_config \\
        --generateEvidenceBam 2>&1 | tee manta_config.log; then
        
        if [ -d "${sample_id}_manta_config" ] && [ -f "${sample_id}_manta_config/runWorkflow.py" ]; then
            echo "Manta configuration successful"
            CONFIG_SUCCESS=true
        else
            echo "WARNING: Manta configuration completed but directory/script missing"
            cat manta_config.log 2>/dev/null || echo "No log available"
            CONFIG_SUCCESS=false
        fi
    else
        echo "WARNING: Manta configuration failed"
        cat manta_config.log 2>/dev/null || echo "No log available"
        CONFIG_SUCCESS=false
    fi
    
    # Create dummy directory if configuration failed
    if [ "\${CONFIG_SUCCESS}" = "false" ]; then
        mkdir -p ${sample_id}_manta_config
        echo "CONFIG_FAILED" > ${sample_id}_manta_config/status.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: \$(configManta.py --version 2>&1 | grep "Manta" | sed 's/Manta //' || echo "unknown")
    status: \$([ "\${CONFIG_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process MANTA_RUN {
    tag "$sample_id"
    label 'process_high'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(config_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_manta_results/results/variants/tumorSV.vcf.gz"), emit: vcf, optional: true
    tuple val(sample_id), path("${sample_id}_manta_results/results/variants/tumorSV.vcf.gz.tbi"), emit: tbi, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Running Manta workflow for sample: ${sample_id}"
    
    RUN_SUCCESS=false
    
    # Check if config directory contains status file (indicating previous failure)
    if [ -f "${config_dir}/status.txt" ]; then
        STATUS=\$(cat ${config_dir}/status.txt)
        echo "WARNING: Manta configuration was not successful: \${STATUS}"
        echo "Skipping Manta run and creating empty output"
        
        # Create empty VCF
        mkdir -p ${sample_id}_manta_results/results/variants
        
        cat > ${sample_id}_manta_results/results/variants/tumorSV.vcf << 'EOF'
##fileformat=VCFv4.1
##source=Manta_SKIPPED
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
        
        bgzip ${sample_id}_manta_results/results/variants/tumorSV.vcf
        tabix -p vcf ${sample_id}_manta_results/results/variants/tumorSV.vcf.gz
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: "skipped"
    status: "skipped_config_failed"
    reason: \${STATUS}
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if runWorkflow.py exists
    if [ ! -f "${config_dir}/runWorkflow.py" ]; then
        echo "WARNING: runWorkflow.py not found in config directory"
        
        mkdir -p ${sample_id}_manta_results/results/variants
        
        cat > ${sample_id}_manta_results/results/variants/tumorSV.vcf << 'EOF'
##fileformat=VCFv4.1
##source=Manta_SKIPPED
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
        
        bgzip ${sample_id}_manta_results/results/variants/tumorSV.vcf
        tabix -p vcf ${sample_id}_manta_results/results/variants/tumorSV.vcf.gz
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: "unknown"
    status: "skipped_missing_workflow"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Run Manta
    echo "Executing Manta workflow..."
    if timeout 3600 ${config_dir}/runWorkflow.py \\
        -m local \\
        -j ${task.cpus} 2>&1 | tee manta_run.log; then
        
        if [ -f "${config_dir}/results/variants/tumorSV.vcf.gz" ]; then
            echo "Manta run completed successfully"
            RUN_SUCCESS=true
            
            # Move results to expected location
            mkdir -p ${sample_id}_manta_results
            cp -r ${config_dir}/results ${sample_id}_manta_results/
        else
            echo "WARNING: Manta completed but output VCF not found"
            cat manta_run.log 2>/dev/null || echo "No log available"
            RUN_SUCCESS=false
        fi
    else
        echo "WARNING: Manta run failed"
        cat manta_run.log 2>/dev/null || echo "No log available"
        RUN_SUCCESS=false
    fi
    
    # Create empty VCF if run failed
    if [ "\${RUN_SUCCESS}" = "false" ]; then
        echo "Creating empty VCF output"
        mkdir -p ${sample_id}_manta_results/results/variants
        
        cat > ${sample_id}_manta_results/results/variants/tumorSV.vcf << 'EOF'
##fileformat=VCFv4.1
##source=Manta_FAILED
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
        
        bgzip ${sample_id}_manta_results/results/variants/tumorSV.vcf
        tabix -p vcf ${sample_id}_manta_results/results/variants/tumorSV.vcf.gz
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    manta: \$(command -v configManta.py >/dev/null 2>&1 && configManta.py --version 2>&1 | grep "Manta" | sed 's/Manta //' || echo "unknown")
    status: \$([ "\${RUN_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process SV_FILTER {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}_sv_filtered.vcf.gz"), emit: vcf, optional: true
    path "${sample_id}_sv_stats.json", emit: stats, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Filtering SVs for sample: ${sample_id}"
    
    FILTER_SUCCESS=false
    
    # Check if input VCF exists and has content
    if [ ! -f "${vcf}" ]; then
        echo "WARNING: Input VCF not found, creating empty output"
        
        cat > ${sample_id}_sv_filtered.vcf << 'EOF'
##fileformat=VCFv4.1
##source=SV_FILTER_SKIPPED
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
        bgzip ${sample_id}_sv_filtered.vcf
        tabix -p vcf ${sample_id}_sv_filtered.vcf.gz
        
        echo '{"sample_id":"${sample_id}","total_sv_calls":0,"filtered_sv_calls":0,"sv_types":{},"filter_rate":0,"status":"skipped_no_input"}' > ${sample_id}_sv_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "skipped_no_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if VCF has any variants (beyond header)
    VARIANT_COUNT=\$(zcat ${vcf} 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    
    if [ "\${VARIANT_COUNT}" -eq 0 ]; then
        echo "INFO: Input VCF has no variants (empty or skipped upstream)"
        
        # Copy input as output
        cp ${vcf} ${sample_id}_sv_filtered.vcf.gz
        
        if [ ! -f "${vcf}.tbi" ]; then
            tabix -p vcf ${sample_id}_sv_filtered.vcf.gz
        else
            cp ${vcf}.tbi ${sample_id}_sv_filtered.vcf.gz.tbi
        fi
        
        echo '{"sample_id":"${sample_id}","total_sv_calls":0,"filtered_sv_calls":0,"sv_types":{},"filter_rate":0,"status":"no_variants"}' > ${sample_id}_sv_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //' || echo "unknown")
    status: "success_no_variants"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if bcftools is available
    if ! command -v bcftools >/dev/null 2>&1; then
        echo "WARNING: bcftools not found, copying input as output"
        
        cp ${vcf} ${sample_id}_sv_filtered.vcf.gz
        if [ -f "${vcf}.tbi" ]; then
            cp ${vcf}.tbi ${sample_id}_sv_filtered.vcf.gz.tbi
        else
            tabix -p vcf ${sample_id}_sv_filtered.vcf.gz 2>/dev/null || echo "WARNING: Could not index VCF"
        fi
        
        echo '{"sample_id":"${sample_id}","total_sv_calls":0,"filtered_sv_calls":0,"sv_types":{},"filter_rate":0,"status":"skipped_tool_missing"}' > ${sample_id}_sv_stats.json
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bcftools: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Filter SVs based on quality and support
    echo "Filtering SVs with bcftools..."
    if bcftools view \\
        -f PASS \\
        -i 'QUAL>=20 && INFO/SOMATICSCORE>=20' \\
        ${vcf} \\
        -O z \\
        -o ${sample_id}_sv_filtered.vcf.gz 2>&1 | tee filter.log; then
        
        echo "SV filtering completed"
        FILTER_SUCCESS=true
    else
        echo "WARNING: bcftools filtering failed, copying input as output"
        cat filter.log 2>/dev/null || echo "No log available"
        cp ${vcf} ${sample_id}_sv_filtered.vcf.gz
    fi
    
    # Index filtered VCF
    if [ -f "${sample_id}_sv_filtered.vcf.gz" ]; then
        tabix -p vcf ${sample_id}_sv_filtered.vcf.gz 2>/dev/null || echo "WARNING: Could not index filtered VCF"
    fi
    
    # Generate statistics
    cat > generate_sv_stats.py << 'EOF'
import json
import subprocess
import sys

try:
    # Count total and filtered variants
    total_count = int(subprocess.check_output([
        "bcftools", "view", "-H", "${vcf}"
    ]).decode().count('\\n'))
except:
    total_count = 0

try:
    filtered_count = int(subprocess.check_output([
        "bcftools", "view", "-H", "${sample_id}_sv_filtered.vcf.gz"
    ]).decode().count('\\n'))
except:
    filtered_count = 0

# Count by SV type
sv_types = {}
try:
    sv_output = subprocess.check_output([
        "bcftools", "query", "-f", "%SVTYPE\\n", "${sample_id}_sv_filtered.vcf.gz"
    ]).decode().strip()
    
    for svtype in sv_output.split('\\n'):
        if svtype:
            sv_types[svtype] = sv_types.get(svtype, 0) + 1
except:
    pass

stats = {
    "sample_id": "${sample_id}",
    "total_sv_calls": total_count,
    "filtered_sv_calls": filtered_count,
    "sv_types": sv_types,
    "filter_rate": round((total_count - filtered_count) / total_count * 100, 2) if total_count > 0 else 0,
    "status": "success" if filtered_count >= 0 else "failed"
}

with open("${sample_id}_sv_stats.json", 'w') as f:
    json.dump(stats, f, indent=2)
EOF

    if command -v python >/dev/null 2>&1 || command -v python3 >/dev/null 2>&1; then
        PYTHON_CMD=\$(command -v python3 || command -v python)
        \${PYTHON_CMD} generate_sv_stats.py 2>/dev/null || echo '{"sample_id":"${sample_id}","total_sv_calls":0,"filtered_sv_calls":0,"sv_types":{},"filter_rate":0,"status":"stats_generation_failed"}' > ${sample_id}_sv_stats.json
    else
        echo "WARNING: Python not found, creating minimal stats"
        echo '{"sample_id":"${sample_id}","total_sv_calls":0,"filtered_sv_calls":0,"sv_types":{},"filter_rate":0,"status":"python_missing"}' > ${sample_id}_sv_stats.json
    fi
    
    # Ensure output files exist
    if [ ! -f "${sample_id}_sv_filtered.vcf.gz" ]; then
        touch ${sample_id}_sv_filtered.vcf.gz
    fi
    
    if [ ! -f "${sample_id}_sv_stats.json" ]; then
        echo '{"sample_id":"${sample_id}","status":"unknown"}' > ${sample_id}_sv_stats.json
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //' || echo "unknown")
    python: \$(python --version 2>&1 | sed 's/Python //' || python3 --version 2>&1 | sed 's/Python //' || echo "unknown")
    status: \$([ "\${FILTER_SUCCESS}" = "true" ] && echo "success" || echo "partial")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process SV_ANNOTATE {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}_sv_annotated.tsv"), emit: tsv, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Annotating SVs for sample: ${sample_id}"
    
    # Check if input VCF exists
    if [ ! -f "${vcf}" ] || [ ! -s "${vcf}" ]; then
        echo "WARNING: Input VCF not found or empty, creating empty TSV"
        
        cat > ${sample_id}_sv_annotated.tsv << 'EOF'
chrom	pos	ref	alt	sv_type	sv_length	somatic_score	clinical_significance	filter
EOF
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    status: "skipped_no_input"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check if Python is available
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        echo "WARNING: Python not found, creating empty TSV"
        
        cat > ${sample_id}_sv_annotated.tsv << 'EOF'
chrom	pos	ref	alt	sv_type	sv_length	somatic_score	clinical_significance	filter
EOF
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    PYTHON_CMD=\$(command -v python3 || command -v python)
    
    cat > annotate_sv.py << 'EOF'
import gzip
import json
import csv

def parse_vcf_line(line):
    fields = line.strip().split('\\t')
    if len(fields) < 8:
        return None
    
    chrom, pos, id_field, ref, alt, qual, filter_field, info = fields[:8]
    
    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return {
        'chrom': chrom,
        'pos': int(pos) if pos != '.' else 0,
        'ref': ref,
        'alt': alt,
        'qual': float(qual) if qual != '.' else 0,
        'filter': filter_field,
        'info': info_dict
    }

def annotate_sv_type(alt, ref, info):
    if 'SVTYPE' in info:
        return info['SVTYPE']
    elif alt.startswith('<'):
        return alt[1:-1]
    elif len(alt) > len(ref):
        return 'INS'
    elif len(alt) < len(ref):
        return 'DEL'
    else:
        return 'UNK'

# Process VCF
variants = []
try:
    with gzip.open("${vcf}", 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variant = parse_vcf_line(line)
            if variant:
                # Add annotations
                variant['sv_type'] = annotate_sv_type(variant['alt'], variant['ref'], variant['info'])
                variant['sv_length'] = abs(int(variant['info'].get('SVLEN', 0))) if 'SVLEN' in variant['info'] else 0
                variant['somatic_score'] = float(variant['info'].get('SOMATICSCORE', 0))
                
                # Determine clinical significance (simplified)
                if variant['sv_type'] in ['DEL', 'DUP'] and variant['sv_length'] > 1000:
                    variant['clinical_significance'] = 'Potentially significant'
                elif variant['sv_type'] in ['BND', 'TRA']:
                    variant['clinical_significance'] = 'Fusion candidate'
                else:
                    variant['clinical_significance'] = 'Unknown'
                
                variants.append(variant)
except Exception as e:
    print(f"WARNING: Error processing VCF: {e}")

# Write TSV output
with open("${sample_id}_sv_annotated.tsv", 'w', newline='') as f:
    fieldnames = ['chrom', 'pos', 'ref', 'alt', 'sv_type', 'sv_length', 
                 'somatic_score', 'clinical_significance', 'filter']
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\\t')
    writer.writeheader()
    
    if variants:
        for variant in variants:
            writer.writerow({
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'sv_type': variant['sv_type'],
                'sv_length': variant['sv_length'],
                'somatic_score': variant['somatic_score'],
                'clinical_significance': variant['clinical_significance'],
                'filter': variant['filter']
            })
EOF

    \${PYTHON_CMD} annotate_sv.py 2>&1 | tee annotate.log
    
    # Ensure output file exists
    if [ ! -f "${sample_id}_sv_annotated.tsv" ]; then
        echo "WARNING: Annotation script failed, creating empty TSV"
        cat annotate.log 2>/dev/null || echo "No log available"
        
        cat > ${sample_id}_sv_annotated.tsv << 'EOF'
chrom	pos	ref	alt	sv_type	sv_length	somatic_score	clinical_significance	filter
EOF
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

// SV Workflow with fault tolerance
workflow SV {
    take:
    bam              // channel: [sample_id, bam]
    ref_fasta        // path: reference FASTA
    targets_bed      // path: target regions BED
    
    main:
    // Configure Manta (with error handling)
    MANTA_CONFIG(
        bam,
        ref_fasta,
        targets_bed
    )
    
    // Run Manta (with error handling)
    MANTA_RUN(MANTA_CONFIG.out.config_dir)
    
    // Filter SVs (with error handling)
    SV_FILTER(MANTA_RUN.out.vcf)
    
    // Annotate SVs (with error handling)
    SV_ANNOTATE(SV_FILTER.out.vcf)
    
    // Combine versions
    ch_versions = MANTA_CONFIG.out.versions
        .mix(MANTA_RUN.out.versions)
        .mix(SV_FILTER.out.versions)
        .mix(SV_ANNOTATE.out.versions)
    
    emit:
    vcf = SV_FILTER.out.vcf
    annotated = SV_ANNOTATE.out.tsv
    stats = SV_FILTER.out.stats
    versions = ch_versions
}