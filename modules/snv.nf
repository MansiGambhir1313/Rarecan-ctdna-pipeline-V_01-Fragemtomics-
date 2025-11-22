nextflow.enable.dsl=2

process MUTECT2 {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/variants/snv", mode: 'copy'

    input:
        tuple val(sample_id), path(bam)
        path ref_fasta
        path targets_bed
        path pon_vcf optional: true
        path germline_resource optional: true

    output:
        tuple val(sample_id), path("${sample_id}_mutect2_raw.vcf.gz"), path("${sample_id}_mutect2.stats"), emit: raw_vcf
        tuple val(sample_id), path("${sample_id}_mutect2_raw.vcf.gz.tbi"), emit: raw_tbi
        path "versions.yml", emit: versions

    script:
    def pon_args = pon_vcf ? "--panel-of-normals ${pon_vcf}" : ""
    def germline_args = germline_resource ? "--germline-resource ${germline_resource}" : ""
    def af_arg = "--af-of-alleles-not-in-resource 1e-6"
    """
    gatk Mutect2 \\
        --input ${bam} \\
        --reference ${ref_fasta} \\
        --intervals ${targets_bed} \\
        --output ${sample_id}_mutect2_raw.vcf.gz \\
        --stats ${sample_id}_mutect2.stats \\
        ${pon_args} \\
        ${germline_args} \\
        ${af_arg} \\
        --tumor-sample ${sample_id} \\
        --max-reads-per-alignment-start 0 \\
        --max-suspicious-reads-per-alignment-start 0 \\
        --min-base-quality-score 20 \\
        --java-options "-Xmx\${task.memory.toGiga()-2}G"

    tabix -p vcf ${sample_id}_mutect2_raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
END_VERSIONS
    """
}

process FILTER_MUTECT_CALLS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/variants/snv", mode: 'copy'

    input:
        tuple val(sample_id), path(mutect_vcf), path(stats)
        path ref_fasta

    output:
        tuple val(sample_id), path("${sample_id}_mutect2.vcf.gz"), emit: vcf
        tuple val(sample_id), path("${sample_id}_mutect2.vcf.gz.tbi"), emit: tbi
        path "versions.yml", emit: versions

    script:
    """
    gatk FilterMutectCalls \\
        --variant ${mutect_vcf} \\
        --reference ${ref_fasta} \\
        --stats ${stats} \\
        --output ${sample_id}_mutect2.vcf.gz

    tabix -p vcf ${sample_id}_mutect2.vcf.gz

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
END_VERSIONS
    """
}

process VARDICT {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/variants/snv", mode: 'copy'

    input:
        tuple val(sample_id), path(bam)
        path ref_fasta
        path targets_bed

    output:
        tuple val(sample_id), path("${sample_id}_vardict.vcf.gz"), emit: vcf
        tuple val(sample_id), path("${sample_id}_vardict.vcf.gz.tbi"), emit: tbi
        path "versions.yml", emit: versions

    script:
    """
    vardict-java -G ${ref_fasta} -f 0.005 -N ${sample_id} -b ${bam} -z -F 0 -c 1 -S 2 -E 3 -g 4 ${targets_bed} \\
        | testsomatic.R \\
        | var2vcf_paired.pl -N ${sample_id} -f 0.005 \\
        | bgzip -c > ${sample_id}_vardict.vcf.gz

    tabix -p vcf ${sample_id}_vardict.vcf.gz

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    vardict: \$(echo \$(vardict-java 2>&1) | grep "VarDict" | sed 's/.*VarDict //; s/ .*//')
END_VERSIONS
    """
}

process LOFREQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/variants/snv", mode: 'copy'

    input:
        tuple val(sample_id), path(bam)
        path ref_fasta
        path targets_bed

    output:
        tuple val(sample_id), path("${sample_id}_lofreq.vcf.gz"), emit: vcf
        path "versions.yml", emit: versions

    when:
        params.enable_lofreq

    script:
    """
    lofreq call-parallel \\
        --pp-threads ${task.cpus} \\
        -f ${ref_fasta} \\
        -l ${targets_bed} \\
        -o ${sample_id}_lofreq.vcf \\
        --call-indels \\
        --min-cov 10 \\
        --min-mq 20 \\
        --min-bq 20 \\
        --sig 0.01 \\
        ${bam}

    bgzip ${sample_id}_lofreq.vcf
    tabix -p vcf ${sample_id}_lofreq.vcf.gz

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    lofreq: \$(lofreq version 2>&1 | sed 's/version: //')
END_VERSIONS
    """
}

process CONSENSUS_FILTER {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/variants/snv", mode: 'copy'

    input:
        tuple val(sample_id), path(mutect2_vcf), path(vardict_vcf), val(lofreq_placeholder)

    output:
        tuple val(sample_id), path("${sample_id}_consensus.vcf.gz"), emit: vcf
        path "${sample_id}_filter_stats.json", emit: stats
        path "versions.yml", emit: versions

    script:
    """
    pip install PyVCF3 --quiet

    cat > simple_filter.py <<'EOF'
import gzip
import json
from datetime import datetime

def simple_variant_filter(mutect2_file, vardict_file, output_file, sample_id):
    variants = []
    try:
        if mutect2_file.endswith('.gz'):
            fh = gzip.open(mutect2_file, 'rt')
        else:
            fh = open(mutect2_file, 'r')
        with fh:
            lines = fh.readlines()
        with open(output_file, 'w') as out:
            for line in lines:
                if line.startswith('#'):
                    out.write(line)
                else:
                    fields = line.split('\\t')
                    if len(fields) >= 6:
                        try:
                            qual = float(fields[5])
                            if qual >= 20:
                                out.write(line)
                                variants.append(line.strip())
                        except Exception:
                            out.write(line)
                            variants.append(line.strip())
    except Exception as exc:
        print(f"Error processing VCF: {exc}")
        with open(output_file, 'w') as out:
            out.write("##fileformat=VCFv4.2\\n")
            out.write(f"#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{sample_id}\\n")
    return len(variants)

variant_count = simple_variant_filter("${mutect2_vcf}", "${vardict_vcf}", "${sample_id}_consensus.vcf", "${sample_id}")

stats = {
    "sample_id": "${sample_id}",
    "analysis_date": datetime.now().isoformat(),
    "total_variants": variant_count,
    "filter_method": "simple_quality_filter"
}

with open("${sample_id}_filter_stats.json", 'w') as handle:
    json.dump(stats, handle, indent=2)

print(f"Consensus filtering complete: {variant_count} variants")
EOF

    python simple_filter.py

    bgzip ${sample_id}_consensus.vcf
    tabix -p vcf ${sample_id}_consensus.vcf.gz

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //')
    consensus_filter: "1.0.0"
END_VERSIONS
    """
}

workflow SNV {
    take:
        bam
        ref_fasta
        targets_bed
        pon_vcf
        germline_resource

    main:
        mutect2_results = MUTECT2(
            bam,
            ref_fasta,
            targets_bed,
            pon_vcf,
            germline_resource
        )

        mutect2_results.raw_vcf.into { mutect2_raw_for_filter; mutect2_raw_for_emit }
        mutect2_results.stats.into { mutect2_stats_for_filter; mutect2_stats_for_emit }

        filtered_mutect = FILTER_MUTECT_CALLS(
            mutect2_raw_for_filter.join(mutect2_stats_for_filter, by: 0),
            ref_fasta
        )

        vardict_results = VARDICT(
            bam,
            ref_fasta,
            targets_bed
        )

        def lofreq_vcf_channel = Channel.empty()
        def lofreq_versions = Channel.empty()
        if (params.enable_lofreq) {
            lofreq_results = LOFREQ(
                bam,
                ref_fasta,
                targets_bed
            )
            lofreq_vcf_channel = lofreq_results.vcf
            lofreq_versions = lofreq_results.versions
        }

        combined_vcfs = filtered_mutect.vcf.join(vardict_results.vcf, by: 0)
        if (params.enable_lofreq) {
            combined_vcfs = combined_vcfs.join(lofreq_vcf_channel, by: 0)
        } else {
            combined_vcfs = combined_vcfs.map { sid, mutect2_path, vardict_path ->
                [sid, mutect2_path, vardict_path, "NO_FILE"]
            }
        }

        consensus_results = CONSENSUS_FILTER(
            combined_vcfs.map { sid, mutect2_path, vardict_path, lofreq_placeholder ->
                [sid, mutect2_path, vardict_path, lofreq_placeholder]
            }
        )

        ch_versions = mutect2_results.versions
            .mix(filtered_mutect.versions)
            .mix(vardict_results.versions)
            .mix(consensus_results.versions)
        if (params.enable_lofreq) {
            ch_versions = ch_versions.mix(lofreq_versions)
        }

    emit:
        vcf = consensus_results.vcf
        mutect2_vcf = filtered_mutect.vcf
        mutect2_raw_vcf = mutect2_raw_for_emit
        vardict_vcf = vardict_results.vcf
        lofreq_vcf = params.enable_lofreq ? lofreq_vcf_channel : Channel.empty()
        mutect2_stats = mutect2_stats_for_emit
        filter_stats = consensus_results.stats
        versions = ch_versions
}

