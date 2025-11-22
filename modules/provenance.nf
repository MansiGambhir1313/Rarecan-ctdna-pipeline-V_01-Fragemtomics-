// modules/provenance.nf - Provenance Tracking Module
// Tracks all inputs, outputs, and processing steps for audit trail

nextflow.enable.dsl=2

process TRACK_PROVENANCE {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/provenance", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(input_files)
    path config_files
    val pipeline_version
    
    output:
    tuple val(sample_id), path("${sample_id}_provenance.json"), emit: provenance
    
    script:
    """
    python << 'EOF'
import json
import os
import hashlib
from datetime import datetime

def calculate_file_hash(filepath):
    \"\"\"Calculate SHA256 hash of file\"\"\"
    sha256_hash = hashlib.sha256()
    try:
        with open(filepath, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
    except Exception as e:
        return None

def get_file_info(filepath):
    \"\"\"Get file metadata\"\"\"
    try:
        if not os.path.exists(filepath):
            return {"path": filepath, "error": "File does not exist"}
        stat = os.stat(filepath)
        return {
            "path": filepath,
            "size": stat.st_size,
            "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
            "sha256": calculate_file_hash(filepath)
        }
    except Exception as e:
        return {"path": filepath, "error": f"Could not read file: {str(e)}"}

# Collect provenance information
provenance = {
    "sample_id": "${sample_id}",
    "pipeline_version": "${pipeline_version}",
    "analysis_date": datetime.now().isoformat(),
    "workflow_name": "RareCan-ctDNA-Pipeline",
    "inputs": [],
    "configurations": [],
    "environment": {
        "nextflow_version": os.environ.get("NXF_VER", "unknown"),
        "container_engine": os.environ.get("CONTAINER_ENGINE", "docker"),
        "executor": os.environ.get("NXF_EXECUTOR", "unknown")
    },
    "processing_steps": []
}

# Track input files - Nextflow stages files in work directory
# List all files in current directory and filter for input files
for item in os.listdir('.'):
    if os.path.isfile(item):
        # Check if it's an input file (FASTQ)
        if item.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
            provenance["inputs"].append(get_file_info(item))
        # Check if it's a config file (reference, BED)
        elif item.endswith(('.fa', '.fasta', '.bed', '.dict', '.fai')):
            provenance["configurations"].append(get_file_info(item))

# Write provenance JSON
with open("${sample_id}_provenance.json", 'w') as f:
    json.dump(provenance, f, indent=2)

print(f"Provenance tracking complete for ${sample_id}")
EOF
    """
}

// Provenance Tracking Workflow
workflow PROVENANCE {
    take:
    sample_id_ch
    input_files_ch
    config_files_ch
    pipeline_version_ch
    
    main:
    // For value channels, we need to extract values and create a single queue channel
    // Use map to extract values from value channels and create tuples
    sample_id_val = sample_id_ch.map { it }
    input_files_val = input_files_ch.map { it }
    config_files_val = config_files_ch.map { it }
    pipeline_version_val = pipeline_version_ch.map { it }
    
    // Combine the mapped channels (now they're queue channels)
    combined = sample_id_val
        .combine(input_files_val)
        .combine(config_files_val)
        .combine(pipeline_version_val)
    
    // Track provenance - create proper tuple structure
    TRACK_PROVENANCE(
        combined.map { sample_id, input_files, config_files, version -> [sample_id, input_files] },
        combined.map { sample_id, input_files, config_files, version -> config_files },
        combined.map { sample_id, input_files, config_files, version -> version }
    )
    
    emit:
    provenance = TRACK_PROVENANCE.out.provenance
}

