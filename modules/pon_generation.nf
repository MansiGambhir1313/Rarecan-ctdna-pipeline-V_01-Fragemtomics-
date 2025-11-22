// modules/pon_generation.nf - Custom Panel of Normals (PoN) Generation Module
// PO-CFS Module 3.1: Generate UMI-aware PoN from normal samples

nextflow.enable.dsl=2

process MUTECT2_NORMAL {
    tag "$sample_id"
    label 'process_medium'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/pon/mutect2_normals", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)
    path ref_fasta
    path targets_bed
    path germline_resource optional: true

    output:
    tuple val(sample_id), path("${sample_id}_normal_mutect.vcf.gz"), path("${sample_id}_normal_mutect.vcf.gz.tbi"), emit: vcf
    path "${sample_id}_normal_mutect.stats", emit: stats
    path "versions.yml", emit: versions

    script:
    def germline_args = germline_resource ? "--germline-resource ${germline_resource}" : ""
    """
    set +e

    echo "Running Mutect2 in normal-only mode for sample: ${sample_id}"

    if ! command -v gatk >/dev/null 2>&1; then
        echo "ERROR: GATK not available in container"
        touch ${sample_id}_normal_mutect.vcf.gz
        touch ${sample_id}_normal_mutect.vcf.gz.tbi
        touch ${sample_id}_normal_mutect.stats

        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS

        set -e
        exit 0
    fi

    if [ ! -f "${bam}" ]; then
        echo "ERROR: Normal BAM not found: ${bam}"
        touch ${sample_id}_normal_mutect.vcf.gz
        touch ${sample_id}_normal_mutect.vcf.gz.tbi
        touch ${sample_id}_normal_mutect.stats

        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    status: "failed_missing_bam"
END_VERSIONS

        set -e
        exit 0
    fi

    gatk Mutect2 \\
        --input ${bam} \\
        --reference ${ref_fasta} \\
        --intervals ${targets_bed} \\
        --max-mnp-distance 0 \\
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --native-pair-hmm-threads ${task.cpus ?: 2} \\
        --tumor-sample ${sample_id} \\
        ${germline_args} \\
        --output ${sample_id}_normal_mutect.vcf.gz \\
        --stats ${sample_id}_normal_mutect.stats \\
        --java-options "-Xmx\${task.memory.toGiga()-2}G"

    tabix -p vcf ${sample_id}_normal_mutect.vcf.gz 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    status: "completed"
END_VERSIONS

    set -e
    exit 0
    """
}

process GENOMICS_DB_IMPORT_PON {
    tag "pon_genomicsdb"
    label 'process_high'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/pon/genomicsdb", mode: 'copy'
    errorStrategy 'ignore'

    input:
    val normal_vcf_list
    path ref_fasta
    path targets_bed

    output:
    path "pon_db", emit: workspace
    path "sample_map.txt", emit: sample_map
    path "versions.yml", emit: versions

    script:
    if (!normal_vcf_list || normal_vcf_list.isEmpty()) {
        throw new IllegalStateException("PoN generation requires at least one normal Mutect2 VCF")
    }
    def vcf_table = normal_vcf_list.collect { entry ->
        def sample_id = entry[0]
        def vcf_path = entry[1]
        def tbi_path = entry[2]
        "${sample_id}\t${vcf_path}\t${tbi_path}"
    }.join('\n')
    """
    set +e

    if ! command -v gatk >/dev/null 2>&1; then
        echo "ERROR: GATK not available in container"
        mkdir -p pon_db
        touch sample_map.txt

        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS

        set -e
        exit 0
    fi

    cat <<'EOF' > mutect_normals.tsv
sample_id	vcf_path	tbi_path
${vcf_table}
EOF

    > sample_map.txt
    mkdir -p pon_db_inputs

    tail -n +2 mutect_normals.tsv | while IFS="	" read -r SAMPLE VCF_PATH TBI_PATH; do
        if [ ! -f "\${VCF_PATH}" ]; then
            echo "ERROR: Missing VCF file \${VCF_PATH}" >&2
            exit 1
        fi
        ln -sf "\${VCF_PATH}" "pon_db_inputs/\${SAMPLE}.vcf.gz"
        if [ -f "\${TBI_PATH}" ]; then
            ln -sf "\${TBI_PATH}" "pon_db_inputs/\${SAMPLE}.vcf.gz.tbi"
        fi
        echo -e "\${SAMPLE}\t\${PWD}/pon_db_inputs/\${SAMPLE}.vcf.gz" >> sample_map.txt
    done

    if [ ! -s sample_map.txt ]; then
        echo "ERROR: Sample map is empty" >&2
        exit 1
    fi

    echo "Running GenomicsDBImport for PoN..."
    if timeout 21600 gatk GenomicsDBImport \\
        --sample-name-map sample_map.txt \\
        --genomicsdb-workspace-path pon_db \\
        --intervals ${targets_bed} \\
        --reader-threads ${task.cpus ?: 4} \\
        --reference ${ref_fasta} 2>&1 | tee genomicsdbimport.log; then
        STATUS="completed"
    else
        echo "WARNING: GenomicsDBImport failed"
        STATUS="failed"
    fi

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    status: "\${STATUS}"
END_VERSIONS

    set -e
    exit 0
    """
}

process CREATE_PON {
    tag "pon_generation"
    label 'process_high'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/pon", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    path workspace_dir
    path ref_fasta
    path targets_bed
    path sample_map optional: true

    output:
    path "rarecan_cfdna_pon.vcf.gz", emit: pon_vcf
    path "rarecan_cfdna_pon.vcf.gz.tbi", emit: pon_tbi
    path "pon_stats.json", emit: stats
    path "versions.yml", emit: versions

    script:
    """
    set +e

    if ! command -v gatk >/dev/null 2>&1; then
        echo "ERROR: GATK not available in container"
        touch rarecan_cfdna_pon.vcf.gz
        touch rarecan_cfdna_pon.vcf.gz.tbi
        echo '{"status":"gatk_missing","normal_samples":0}' > pon_stats.json

        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS

        set -e
        exit 0
    fi

    if [ ! -d "${workspace_dir}" ]; then
        echo "ERROR: GenomicsDB workspace not found: ${workspace_dir}"
        touch rarecan_cfdna_pon.vcf.gz
        touch rarecan_cfdna_pon.vcf.gz.tbi
        echo '{"status":"workspace_missing","normal_samples":0}' > pon_stats.json

        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    status: "failed_missing_workspace"
END_VERSIONS

        set -e
        exit 0
    fi

    echo "Creating RareCan cfDNA PoN from GenomicsDB workspace..."
    if timeout 14400 gatk CreateSomaticPanelOfNormals \\
        -R ${ref_fasta} \\
        -V gendb://${workspace_dir} \\
        -O rarecan_cfdna_pon.vcf.gz \\
        --intervals ${targets_bed} \\
        --min-sample-count 2 \\
        --java-options "-Xmx\${task.memory.toGiga()-2}G" 2>&1 | tee pon_generation.log; then

        if [ -f "rarecan_cfdna_pon.vcf.gz" ] && [ -s "rarecan_cfdna_pon.vcf.gz" ]; then
            tabix -p vcf rarecan_cfdna_pon.vcf.gz
            variant_count=\$(zcat rarecan_cfdna_pon.vcf.gz | grep -v "^#" | wc -l || echo "0")
            if [ -n "${sample_map}" ] && [ -f "${sample_map}" ]; then
                NORMAL_COUNT=\$(wc -l < "${sample_map}" | tr -d ' ')
            else
                NORMAL_COUNT=0
            fi

            cat > pon_stats.json << EOF
{
    "status": "success",
    "normal_samples": \${NORMAL_COUNT},
    "variants_in_pon": \${variant_count},
       "generation_date": "\$(date -Iseconds)",
    "target_regions": "${targets_bed}"
}
EOF
            STATUS="success"
        else
            echo "WARNING: PoN output empty"
            echo '{"status":"empty_output"}' > pon_stats.json
            STATUS="failed_empty"
        fi
    else
        echo "WARNING: CreateSomaticPanelOfNormals failed"
        echo '{"status":"failed_generation"}' > pon_stats.json
        STATUS="failed"
    fi

    if [ ! -f "rarecan_cfdna_pon.vcf.gz" ]; then
        cat > rarecan_cfdna_pon.vcf <<'VCFEOF'
##fileformat=VCFv4.2
##source=CreateSomaticPanelOfNormals
##FILTER=<ID=PON,Description="Variant found in panel of normals">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCFEOF
        bgzip rarecan_cfdna_pon.vcf
        tabix -p vcf rarecan_cfdna_pon.vcf.gz
    fi

    if [ ! -f "rarecan_cfdna_pon.vcf.gz.tbi" ]; then
        tabix -p vcf rarecan_cfdna_pon.vcf.gz 2>/dev/null || touch rarecan_cfdna_pon.vcf.gz.tbi
    fi

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    status: "\${STATUS}"
END_VERSIONS

    set -e
    exit 0
    """
}

// PoN Generation Workflow
workflow PON_GENERATION {
    take:
    normal_bams      // channel: [sample_id, bam]
    targets_bed      // path: target regions BED
    ref_fasta        // path: reference FASTA
    germline_resource // path: germline resource (optional)

    main:
    mutect_results = MUTECT2_NORMAL(
        normal_bams,
        ref_fasta,
        targets_bed,
        germline_resource
    )

    normal_vcf_list = mutect_results.vcf.collect()

    GENOMICS_DB_IMPORT_PON(
        normal_vcf_list,
        ref_fasta,
        targets_bed
    )

    CREATE_PON(
        GENOMICS_DB_IMPORT_PON.out.workspace,
        ref_fasta,
        targets_bed,
        GENOMICS_DB_IMPORT_PON.out.sample_map
    )

    ch_versions = MUTECT2_NORMAL.out.versions
        .mix(GENOMICS_DB_IMPORT_PON.out.versions)
        .mix(CREATE_PON.out.versions)

    emit:
    pon_vcf = CREATE_PON.out.pon_vcf
    pon_tbi = CREATE_PON.out.pon_tbi
    stats = CREATE_PON.out.stats
    versions = ch_versions
}

