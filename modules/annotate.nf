// modules/annotate.nf - Variant Annotation Module

nextflow.enable.dsl=2

process VEP_ANNOTATE {
    tag "$sample_id"
    label 'process_high'
    
    input:
    tuple val(sample_id), path(vcf)
    path vep_cache
    path ref_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_vep.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}_vep.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    def cache_args = vep_cache ? "--cache --dir_cache ${vep_cache}" : "--database"
    """
    # Run VEP annotation
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}_vep.vcf \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --stats_html ${sample_id}_vep.html \\
        --fork ${task.cpus} \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        ${cache_args} \\
        --fasta ${ref_fasta} \\
        --hgvs \\
        --symbol \\
        --numbers \\
        --domains \\
        --regulatory \\
        --canonical \\
        --protein \\
        --biotype \\
        --uniprot \\
        --tsl \\
        --appris \\
        --gene_phenotype \\
        --af \\
        --af_1kg \\
        --af_gnomad \\
        --max_af \\
        --pubmed \\
        --variant_class \\
        --sift b \\
        --polyphen b \\
        --humdiv \\
        --condel b \\
        --mastermind \\
        --clinvar \\
        --custom ${projectDir}/assets/kb/cosmic.vcf.gz,COSMIC,vcf,exact,0,CNT \\
        --plugin CADD,${projectDir}/assets/kb/whole_genome_SNVs.tsv.gz \\
        --plugin dbNSFP,${projectDir}/assets/kb/dbNSFP4.4a.txt.gz,REVEL_score,REVEL_rankscore
    
    # Index the output VCF
    tabix -p vcf ${sample_id}_vep.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}

process CLINICAL_ANNOTATION {
    tag "$sample_id"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(vep_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}_clinical_variants.tsv"), emit: clinical_tsv
    tuple val(sample_id), path("${sample_id}_annotation_summary.json"), emit: summary
    path "versions.yml", emit: versions
    
    script:
    """
    cat > clinical_annotation.py << 'EOF'
import gzip
import json
import csv
import re

def parse_vep_csq(csq_string, csq_header):
    \"\"\"Parse VEP CSQ annotation\"\"\"
    annotations = []
    
    if not csq_string or csq_string == '.':
        return annotations
    
    # Split CSQ header to get field names
    csq_fields = csq_header.split('|')
    
    # Parse each consequence
    for csq in csq_string.split(','):
        csq_values = csq.split('|')
        
        if len(csq_values) == len(csq_fields):
            annotation = dict(zip(csq_fields, csq_values))
            annotations.append(annotation)
    
    return annotations

def determine_clinical_significance(annotations, variant_info):
    \"\"\"Determine clinical significance based on annotations\"\"\"
    
    # Initialize significance
    max_significance = "Unknown"
    evidence_level = "5"  # Lowest evidence level
    
    for ann in annotations:
        consequence = ann.get('Consequence', '').lower()
        clinvar = ann.get('ClinVar_CLNSIG', '')
        cosmic = ann.get('COSMIC', '')
        sift = ann.get('SIFT', '')
        polyphen = ann.get('PolyPhen', '')
        
        # ClinVar significance
        if 'pathogenic' in clinvar.lower():
            if 'likely' in clinvar.lower():
                significance = "Likely Pathogenic"
                level = "2"
            else:
                significance = "Pathogenic"
                level = "1"
        elif 'benign' in clinvar.lower():
            significance = "Benign/Likely Benign"
            level = "4"
        elif 'uncertain' in clinvar.lower() or 'vus' in clinvar.lower():
            significance = "VUS"
            level = "3"
        
        # High impact consequences
        elif any(term in consequence for term in ['stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor']):
            significance = "High Impact"
            level = "2"
        
        # Moderate impact with damaging predictions
        elif 'missense' in consequence:
            if 'deleterious' in sift.lower() or 'damaging' in polyphen.lower():
                significance = "Moderate Impact - Damaging"
                level = "3"
            else:
                significance = "Moderate Impact"
                level = "4"
        
        # COSMIC presence
        elif cosmic and cosmic != '.':
            significance = "COSMIC Listed"
            level = "3"
        
        else:
            significance = "Unknown"
            level = "5"
        
        # Keep highest significance
        if int(level) < int(evidence_level):
            max_significance = significance
            evidence_level = level
    
    return max_significance, evidence_level

def extract_variant_info(fields):
    \"\"\"Extract basic variant information\"\"\"
    return {
        'chrom': fields[0],
        'pos': int(fields[1]),
        'ref': fields[3],
        'alt': fields[4],
        'qual': float(fields[5]) if fields[5] != '.' else 0,
        'filter': fields[6]
    }

def parse_info_field(info_string):
    \"\"\"Parse VCF INFO field\"\"\"
    info_dict = {}
    
    for item in info_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return info_dict

# Process VEP-annotated VCF
clinical_variants = []
csq_header = None

with gzip.open("${vep_vcf}", 'rt') as f:
    for line in f:
        if line.startswith('##INFO=<ID=CSQ'):
            # Extract CSQ header
            match = re.search(r'Format: (.+)">', line)
            if match:
                csq_header = match.group(1)
        
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\\t')
        if len(fields) < 8:
            continue
        
        variant_info = extract_variant_info(fields)
        info_dict = parse_info_field(fields[7])
        
        # Parse VEP annotations
        csq_string = info_dict.get('CSQ', '')
        annotations = parse_vep_csq(csq_string, csq_header) if csq_header else []
        
        if annotations:
            # Get the most severe annotation
            primary_annotation = annotations[0]
            
            # Determine clinical significance
            significance, evidence_level = determine_clinical_significance(annotations, variant_info)
            
            # Extract key information
            gene = primary_annotation.get('SYMBOL', '')
            consequence = primary_annotation.get('Consequence', '')
            hgvsp = primary_annotation.get('HGVSp', '')
            hgvsc = primary_annotation.get('HGVSc', '')
            
            # Calculate VAF from sample data if available
            vaf = 0.0
            if len(fields) > 9:
                sample_data = fields[9]
                # Simple VAF extraction (format-dependent)
                if ':' in sample_data:
                    sample_fields = sample_data.split(':')
                    # Look for AD field (allelic depth)
                    if len(sample_fields) > 1:
                        ad_str = sample_fields[1]
                        if ',' in ad_str:
                            ref_depth, alt_depth = map(int, ad_str.split(',')[:2])
                            total_depth = ref_depth + alt_depth
                            vaf = alt_depth / total_depth if total_depth > 0 else 0.0
            
            clinical_variant = {
                'gene': gene,
                'chromosome': variant_info['chrom'],
                'position': variant_info['pos'],
                'ref': variant_info['ref'],
                'alt': variant_info['alt'],
                'consequence': consequence,
                'hgvs_protein': hgvsp,
                'hgvs_cdna': hgvsc,
                'vaf': round(vaf, 4),
                'quality': variant_info['qual'],
                'filter': variant_info['filter'],
                'clinical_significance': significance,
                'evidence_level': evidence_level,
                'clinvar': primary_annotation.get('ClinVar_CLNSIG', ''),
                'cosmic': primary_annotation.get('COSMIC', ''),
                'sift': primary_annotation.get('SIFT', ''),
                'polyphen': primary_annotation.get('PolyPhen', ''),
                'gnomad_af': primary_annotation.get('gnomAD_AF', ''),
                'cadd_score': primary_annotation.get('CADD_PHRED', '')
            }
            
            clinical_variants.append(clinical_variant)

# Write clinical variants TSV
with open("${sample_id}_clinical_variants.tsv", 'w', newline='') as f:
    if clinical_variants:
        fieldnames = ['gene', 'chromosome', 'position', 'ref', 'alt', 'consequence',
                     'hgvs_protein', 'hgvs_cdna', 'vaf', 'quality', 'filter',
                     'clinical_significance', 'evidence_level', 'clinvar', 'cosmic',
                     'sift', 'polyphen', 'gnomad_af', 'cadd_score']
        
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\\t')
        writer.writeheader()
        writer.writerows(clinical_variants)
    else:
        # Empty file with header
        f.write("gene\\tchromosome\\tposition\\tref\\talt\\tconsequence\\thgvs_protein\\thgvs_cdna\\tvaf\\tquality\\tfilter\\tclinical_significance\\tevidence_level\\tclinvar\\tcosmic\\tsift\\tpolyphen\\tgnomad_af\\tcadd_score\\n")

# Create annotation summary
summary = {
    "sample_id": "${sample_id}",
    "total_variants": len(clinical_variants),
    "significance_counts": {},
    "consequence_counts": {},
    "gene_counts": {}
}

# Count by significance
for variant in clinical_variants:
    sig = variant['clinical_significance']
    summary['significance_counts'][sig] = summary['significance_counts'].get(sig, 0) + 1
    
    cons = variant['consequence']
    summary['consequence_counts'][cons] = summary['consequence_counts'].get(cons, 0) + 1
    
    gene = variant['gene']
    if gene:
        summary['gene_counts'][gene] = summary['gene_counts'].get(gene, 0) + 1

with open("${sample_id}_annotation_summary.json", 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Clinical Annotation Complete:")
print(f"  Sample: ${sample_id}")
print(f"  Total Variants: {len(clinical_variants)}")
print(f"  Pathogenic/Likely Pathogenic: {summary['significance_counts'].get('Pathogenic', 0) + summary['significance_counts'].get('Likely Pathogenic', 0)}")
EOF

    python clinical_annotation.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

// Annotation Workflow
workflow ANNOTATE {
    take:
    vcf              // channel: [sample_id, vcf]
    ref_fasta        // path: reference FASTA
    vep_cache        // path: VEP cache directory (optional)
    
    main:
    // Run VEP annotation
    VEP_ANNOTATE(
        vcf,
        vep_cache ?: [],
        ref_fasta
    )
    
    // Clinical annotation and interpretation
    CLINICAL_ANNOTATION(VEP_ANNOTATE.out.vcf)
    
    // Combine versions
    ch_versions = VEP_ANNOTATE.out.versions
        .mix(CLINICAL_ANNOTATION.out.versions)
    
    emit:
    vcf = VEP_ANNOTATE.out.vcf
    html = VEP_ANNOTATE.out.html
    clinical_variants = CLINICAL_ANNOTATION.out.clinical_tsv
    summary = CLINICAL_ANNOTATION.out.summary
    versions = ch_versions
}