// modules/msi.nf - Microsatellite Instability Analysis Module (Fault-Tolerant)

nextflow.enable.dsl=2

process MSISENSOR_PRO_SCAN {
    tag "reference_scan"
    label 'process_medium'
    publishDir "${params.outdir}/msi/reference", mode: 'copy'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/msisensor-pro:latest'
    errorStrategy 'ignore'
    
    input:
    path ref_fasta
    
    output:
    path "microsatellites.list", emit: microsatellites, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Starting microsatellite scan of reference genome..."
    
    SCAN_SUCCESS=false
    
    # Check if msisensor-pro is available
    if ! command -v msisensor-pro >/dev/null 2>&1 && ! command -v msisensor_pro >/dev/null 2>&1; then
        echo "WARNING: msisensor-pro not found in container"
        echo "MSI analysis will be skipped"
        
        # Create empty microsatellites list
        touch microsatellites.list
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor-pro: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    MSISENSOR_CMD=\$(command -v msisensor-pro || command -v msisensor_pro || echo "msisensor-pro")
    
    # Check if reference exists
    if [ ! -f "${ref_fasta}" ]; then
        echo "WARNING: Reference FASTA not found"
        touch microsatellites.list
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor2: \$(msisensor2 --help 2>&1 | head -1 | sed 's/MSISensor2 wrapper //' || echo "unknown")
    status: "skipped_missing_reference"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Scan reference for microsatellites
    echo "Running msisensor-pro scan..."
    if timeout 1800 \${MSISENSOR_CMD} scan \\
        -d ${ref_fasta} \\
        -o microsatellites.list 2>&1 | tee scan.log; then
        
        if [ -f "microsatellites.list" ] && [ -s "microsatellites.list" ]; then
            line_count=\$(wc -l < microsatellites.list)
            echo "Found \$line_count microsatellite loci"
            
            if [ "\$line_count" -lt 10 ]; then
                echo "WARNING: Very few microsatellites detected (\$line_count)"
            else
                SCAN_SUCCESS=true
            fi
        else
            echo "WARNING: Scan completed but output file is empty or missing"
            cat scan.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: msisensor2 scan failed"
        cat scan.log 2>/dev/null || echo "No log available"
    fi
    
    # Create empty file if scan failed
    if [ "\${SCAN_SUCCESS}" = "false" ] || [ ! -f "microsatellites.list" ]; then
        echo "Creating empty microsatellites list"
        touch microsatellites.list
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor-pro: \$(\${MSISENSOR_CMD} --version 2>&1 | head -1 | sed 's/MSISensor-pro //' || echo "unknown")
    status: \$([ "\${SCAN_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    echo "Microsatellite scan completed"
    
    set -e
    exit 0
    """
}

process MSISENSOR_PRO_MSI {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/msi/raw", mode: 'copy'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/msisensor-pro:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(bam)
    path microsatellites
    path msi_bed
    
    output:
    tuple val(sample_id), path("${sample_id}_msi"), emit: msi_output, optional: true
    tuple val(sample_id), path("${sample_id}_msi_dis"), emit: msi_dis, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    echo "Starting MSI analysis for sample: ${sample_id}"
    echo "BAM file: ${bam}"
    echo "Microsatellites file: ${microsatellites}"
    echo "MSI regions file: ${msi_bed}"
    
    MSI_SUCCESS=false
    
    # Check if msisensor-pro is available
    MSISENSOR_CMD=\$(command -v msisensor-pro || command -v msisensor_pro || echo "msisensor-pro")
    
    if ! command -v \${MSISENSOR_CMD} >/dev/null 2>&1; then
        echo "WARNING: msisensor-pro not found, creating empty output"
        
        echo -e "Total_Number_of_Sites\\tNumber_of_Somatic_Sites\\tPercentage\\tMSI_Score" > ${sample_id}_msi
        echo -e "0\\t0\\t0.0\\t0.0000" >> ${sample_id}_msi
        touch ${sample_id}_msi_dis
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor-pro: "not_available"
    samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Validate input BAM
    if [ ! -f "${bam}" ] || [ ! -s "${bam}" ]; then
        echo "WARNING: BAM file not found or empty"
        
        echo -e "Total_Number_of_Sites\\tNumber_of_Somatic_Sites\\tPercentage\\tMSI_Score" > ${sample_id}_msi
        echo -e "0\\t0\\t0.0\\t0.0000" >> ${sample_id}_msi
        touch ${sample_id}_msi_dis
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor2: \$(msisensor2 --help 2>&1 | head -1 | sed 's/MSISensor2 wrapper //' || echo "unknown")
    samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    status: "skipped_invalid_bam"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Check microsatellites file
    if [ ! -f "${microsatellites}" ] || [ ! -s "${microsatellites}" ]; then
        echo "WARNING: Microsatellites file not found or empty (likely from failed scan)"
        
        echo -e "Total_Number_of_Sites\\tNumber_of_Somatic_Sites\\tPercentage\\tMSI_Score" > ${sample_id}_msi
        echo -e "0\\t0\\t0.0\\t0.0000" >> ${sample_id}_msi
        touch ${sample_id}_msi_dis
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor2: \$(msisensor2 --help 2>&1 | head -1 | sed 's/MSISensor2 wrapper //' || echo "unknown")
    samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    status: "skipped_no_microsatellites"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    # Create BAM index if needed
    if [ ! -f "${bam}.bai" ]; then
        echo "Creating BAM index..."
        if command -v samtools >/dev/null 2>&1; then
            samtools index ${bam} 2>/dev/null || echo "WARNING: Could not create BAM index"
        fi
    fi
    
    # Run MSI analysis
    echo "Running MSISensor-pro MSI analysis..."
    if timeout 1800 \${MSISENSOR_CMD} msi \\
        -M ${microsatellites} \\
        -t ${bam} \\
        -e ${msi_bed} \\
        -o ${sample_id} \\
        -b ${task.cpus} 2>&1 | tee msi.log; then
        
        if [ -f "${sample_id}_msi" ] && [ -s "${sample_id}_msi" ]; then
            echo "MSI analysis completed successfully"
            MSI_SUCCESS=true
        else
            echo "WARNING: MSI analysis completed but output file is empty or missing"
            cat msi.log 2>/dev/null || echo "No log available"
        fi
    else
        echo "WARNING: MSI analysis failed"
        cat msi.log 2>/dev/null || echo "No log available"
    fi
    
    # Create minimal output if analysis failed
    if [ "\${MSI_SUCCESS}" = "false" ] || [ ! -f "${sample_id}_msi" ]; then
        echo "Creating minimal output for downstream processing..."
        echo -e "Total_Number_of_Sites\\tNumber_of_Somatic_Sites\\tPercentage\\tMSI_Score" > ${sample_id}_msi
        echo -e "0\\t0\\t0.0\\t0.0000" >> ${sample_id}_msi
    fi
    
    if [ ! -f "${sample_id}_msi_dis" ]; then
        touch ${sample_id}_msi_dis
    fi
    
    # Display results summary
    echo "MSI analysis completed for ${sample_id}"
    if [ -f "${sample_id}_msi" ]; then
        echo "MSI results:"
        cat ${sample_id}_msi
    fi
    
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msisensor-pro: \$(\${MSISENSOR_CMD} --version 2>&1 | head -1 | sed 's/MSISensor-pro //' || echo "unknown")
    samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    status: \$([ "\${MSI_SUCCESS}" = "true" ] && echo "success" || echo "failed")
END_VERSIONS
    
    set -e
    exit 0
    """
}

process MSI_INTERPRET {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/msi/results", mode: 'copy'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(msi_output), path(msi_dis)
    
    output:
    tuple val(sample_id), path("${sample_id}_msi_result.json"), emit: result, optional: true
    tuple val(sample_id), path("${sample_id}_msi_summary.txt"), emit: summary, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    set +e
    
    # Check if Python is available
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        echo "WARNING: Python not found, creating minimal output"
        
        echo '{"sample_id":"${sample_id}","msi_status":"Failed","interpretation":"Python not available","analysis_success":false}' > ${sample_id}_msi_result.json
        echo "MSI Analysis Failed: Python not available" > ${sample_id}_msi_summary.txt
        
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: "not_available"
    status: "skipped_tool_missing"
END_VERSIONS
        
        set -e
        exit 0
    fi
    
    PYTHON_CMD=\$(command -v python3 || command -v python)
    
    cat > interpret_msi.py << 'EOF'
import json
import os
import sys
from datetime import datetime

def parse_msi_output(msi_file):
    \"\"\"Parse MSISensor2 output file\"\"\"
    print(f"Parsing MSI output file: {msi_file}")
    
    if not os.path.exists(msi_file):
        print(f"ERROR: MSI output file does not exist: {msi_file}")
        return None
    
    try:
        with open(msi_file, 'r') as f:
            lines = f.readlines()
        
        print(f"Read {len(lines)} lines from MSI output")
        
        # Find the data line (skip headers and comments)
        data_lines = []
        for i, line in enumerate(lines):
            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
                
            # Skip header line
            if any(word in line.upper() for word in ['TOTAL_NUMBER', 'SITES', 'PERCENTAGE', 'SCORE']):
                print(f"Skipping header line: {line}")
                continue
                
            data_lines.append(line)
        
        if not data_lines:
            print("No data lines found in MSI output")
            return None
        
        # Use the last data line
        result_line = data_lines[-1]
        print(f"Using result line: {result_line}")
        
        fields = result_line.split('\\t')
        print(f"Parsed {len(fields)} fields: {fields}")
        
        if len(fields) >= 3:
            # Parse fields with error handling
            total_sites = 0
            somatic_sites = 0
            msi_score = 0.0
            
            try:
                total_sites = int(float(fields[0])) if fields[0].replace('.', '').replace('-', '').isdigit() else 0
            except (ValueError, IndexError):
                pass
                
            try:
                somatic_sites = int(float(fields[1])) if fields[1].replace('.', '').replace('-', '').isdigit() else 0
            except (ValueError, IndexError):
                pass
            
            # MSI score is typically in the last field or 4th field
            for field in reversed(fields):
                try:
                    if '.' in field or field.replace('.', '').replace('-', '').isdigit():
                        score_candidate = float(field)
                        if 0 <= score_candidate <= 1:  # MSI scores are typically 0-1
                            msi_score = score_candidate
                            break
                except ValueError:
                    continue
            
            print(f"Parsed values - Total sites: {total_sites}, Somatic sites: {somatic_sites}, MSI score: {msi_score}")
            
            return {
                'total_sites': total_sites,
                'somatic_sites': somatic_sites,
                'msi_score': msi_score
            }
        
        print(f"Insufficient fields in result line: {len(fields)}")
        return None
        
    except Exception as e:
        print(f"Error parsing MSI output: {e}")
        import traceback
        traceback.print_exc()
        return None

def interpret_msi_status(msi_score, total_sites, somatic_sites):
    \"\"\"Interpret MSI status based on score and clinical thresholds\"\"\"
    
    # Get thresholds from environment or use defaults
    min_loci = int(os.environ.get('MSI_MIN_LOCI', '${params.msi_min_informative_loci ?: 20}'))
    threshold_high = float(os.environ.get('MSI_THRESHOLD_HIGH', '0.2'))
    threshold_low = float(os.environ.get('MSI_THRESHOLD_LOW', '0.1'))
    
    print(f"Using thresholds - Min loci: {min_loci}, High: {threshold_high}, Low: {threshold_low}")
    print(f"Sample values - Score: {msi_score}, Total sites: {total_sites}, Somatic sites: {somatic_sites}")
    
    # Check if we have sufficient data
    if total_sites < min_loci:
        return {
            "status": "Indeterminate",
            "interpretation": f"Insufficient informative loci ({total_sites} < {min_loci})",
            "confidence": "Low",
            "clinical_significance": "Cannot determine MSI status - insufficient data"
        }
    
    # Determine MSI status based on score
    if msi_score >= threshold_high:
        return {
            "status": "MSI-H",
            "interpretation": "High microsatellite instability",
            "confidence": "High" if total_sites >= 50 else "Medium",
            "clinical_significance": "May benefit from immunotherapy (e.g., pembrolizumab, nivolumab)"
        }
    elif msi_score >= threshold_low:
        return {
            "status": "MSI-L", 
            "interpretation": "Low microsatellite instability",
            "confidence": "Medium",
            "clinical_significance": "Intermediate MSI status - clinical significance uncertain"
        }
    else:
        return {
            "status": "MSS",
            "interpretation": "Microsatellite stable",
            "confidence": "High" if total_sites >= 30 else "Medium",
            "clinical_significance": "Unlikely to benefit from MSI-targeted immunotherapy"
        }

def generate_quality_metrics(msi_data, total_sites):
    \"\"\"Generate quality control metrics\"\"\"
    
    if not msi_data:
        return {
            "analysis_success": False,
            "data_quality": "Failed",
            "coverage_adequacy": "Unknown",
            "reliability": "Low"
        }
    
    # Assess data quality
    if total_sites >= 50:
        coverage_adequacy = "Excellent"
        reliability = "High"
    elif total_sites >= 30:
        coverage_adequacy = "Good"
        reliability = "High"
    elif total_sites >= 20:
        coverage_adequacy = "Adequate"
        reliability = "Medium"
    else:
        coverage_adequacy = "Poor"
        reliability = "Low"
    
    return {
        "analysis_success": True,
        "data_quality": "Pass",
        "coverage_adequacy": coverage_adequacy,
        "reliability": reliability,
        "informative_loci_count": total_sites
    }

# Main analysis
print("="*60)
print("Starting MSI interpretation...")
print("="*60)

# Parse MSI output
msi_data = parse_msi_output("${msi_output}")

if msi_data and (msi_data['total_sites'] > 0 or msi_data['msi_score'] > 0):
    # Interpret MSI status
    interpretation = interpret_msi_status(
        msi_data['msi_score'], 
        msi_data['total_sites'], 
        msi_data['somatic_sites']
    )
    
    # Generate quality metrics
    quality_metrics = generate_quality_metrics(msi_data, msi_data['total_sites'])
    
    # Create comprehensive result
    result = {
        "sample_id": "${sample_id}",
        "analysis_date": datetime.now().isoformat(),
        "msi_score": round(msi_data['msi_score'], 4),
        "total_sites": msi_data['total_sites'],
        "somatic_sites": msi_data['somatic_sites'],
        "msi_status": interpretation['status'],
        "interpretation": interpretation['interpretation'],
        "confidence": interpretation['confidence'],
        "clinical_significance": interpretation['clinical_significance'],
        "method": "MSISensor2",
        "version": "0.1.0_wrapper",
        "thresholds": {
            "high": threshold_high,
            "low": threshold_low,
            "min_informative_loci": min_loci
        },
        "quality_metrics": quality_metrics,
        "analysis_success": True
    }
    
else:
    # Handle failed analysis
    print("MSI analysis failed or returned no results - creating failure result")
    
    result = {
        "sample_id": "${sample_id}",
        "analysis_date": datetime.now().isoformat(),
        "msi_score": 0.0,
        "total_sites": 0,
        "somatic_sites": 0,
        "msi_status": "Failed",
        "interpretation": "MSI analysis failed - insufficient data or processing error",
        "confidence": "N/A",
        "clinical_significance": "Cannot determine - analysis failed",
        "method": "MSISensor2",
        "version": "0.1.0_wrapper",
        "thresholds": {
            "high": 0.2,
            "low": 0.1,
            "min_informative_loci": int(os.environ.get('MSI_MIN_LOCI', '${params.msi_min_informative_loci ?: 20}'))
        },
        "quality_metrics": {
            "analysis_success": False,
            "data_quality": "Failed",
            "coverage_adequacy": "Unknown",
            "reliability": "Low"
        },
        "analysis_success": False
    }

# Write JSON result
with open("${sample_id}_msi_result.json", 'w') as f:
    json.dump(result, f, indent=2, sort_keys=True)

# Write human-readable summary
with open("${sample_id}_msi_summary.txt", 'w') as f:
    f.write("MSI Analysis Summary\\n")
    f.write("=" * 50 + "\\n")
    f.write(f"Sample ID: {result['sample_id']}\\n")
    f.write(f"Analysis Date: {result['analysis_date']}\\n")
    f.write(f"Method: {result['method']}\\n")
    f.write("\\n")
    f.write("Results:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"MSI Status: {result['msi_status']}\\n")
    f.write(f"MSI Score: {result['msi_score']:.4f}\\n")
    f.write(f"Confidence: {result['confidence']}\\n")
    f.write(f"Total Informative Sites: {result['total_sites']}\\n")
    f.write(f"Somatic Sites: {result['somatic_sites']}\\n")
    f.write("\\n")
    f.write("Interpretation:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"{result['interpretation']}\\n")
    f.write("\\n")
    f.write("Clinical Significance:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"{result['clinical_significance']}\\n")
    f.write("\\n")
    f.write("Quality Metrics:\\n")
    f.write("-" * 20 + "\\n")
    f.write(f"Analysis Success: {result['quality_metrics']['analysis_success']}\\n")
    f.write(f"Data Quality: {result['quality_metrics']['data_quality']}\\n")
    f.write(f"Coverage Adequacy: {result['quality_metrics']['coverage_adequacy']}\\n")
    f.write(f"Reliability: {result['quality_metrics']['reliability']}\\n")

# Console output
print("="*60)
print("MSI Analysis Complete")
print("="*60)
print(f"Sample: {result['sample_id']}")
print(f"MSI Score: {result['msi_score']:.4f}")
print(f"Status: {result['msi_status']}")
print(f"Confidence: {result['confidence']}")
print(f"Informative Loci: {result['total_sites']}")
print(f"Analysis Success: {result['analysis_success']}")
print("="*60)
EOF

    # Set environment variables
    export MSI_MIN_LOCI=${params.msi_min_informative_loci ?: 20}
    export MSI_THRESHOLD_HIGH=0.2
    export MSI_THRESHOLD_LOW=0.1
    
    # Run interpretation
    \${PYTHON_CMD} interpret_msi.py 2>&1 | tee interpret.log
    
    # Validate outputs - create minimal files if missing
    if [ ! -f "${sample_id}_msi_result.json" ]; then
        echo "WARNING: MSI result JSON was not created, creating minimal output"
        cat interpret.log 2>/dev/null || echo "No log available"
        echo '{"sample_id":"${sample_id}","msi_status":"Failed","interpretation":"Interpretation script failed","analysis_success":false}' > ${sample_id}_msi_result.json
    fi
    
    if [ ! -f "${sample_id}_msi_summary.txt" ]; then
        echo "WARNING: MSI summary was not created, creating minimal output"
        echo "MSI Analysis Failed" > ${sample_id}_msi_summary.txt
        echo "Sample: ${sample_id}" >> ${sample_id}_msi_summary.txt
        echo "Status: Failed" >> ${sample_id}_msi_summary.txt
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

// MSI Workflow with fault tolerance
workflow MSI {
    take:
    bam              // channel: [sample_id, bam]
    ref_fasta        // path: reference FASTA
    msi_bed          // path: MSI loci BED file
    
    main:
    // Initialize channels
    ch_versions = Channel.empty()
    
    // Scan reference for microsatellites (only needs to be done once)
    MSISENSOR_PRO_SCAN(ref_fasta)
    ch_versions = ch_versions.mix(MSISENSOR_PRO_SCAN.out.versions)
    
    // Run MSI analysis on each sample
    MSISENSOR_PRO_MSI(
        bam,
        MSISENSOR_PRO_SCAN.out.microsatellites,
        msi_bed
    )
    ch_versions = ch_versions.mix(MSISENSOR_PRO_MSI.out.versions)
    
    // Combine MSI outputs for interpretation
    msi_combined = MSISENSOR_PRO_MSI.out.msi_output
        .join(MSISENSOR_PRO_MSI.out.msi_dis, by: 0, remainder: true)
        .map { sample_id, msi_output, msi_dis ->
            // Handle cases where msi_dis might be null
            def safe_msi_dis = msi_dis ?: file("NO_FILE")
            [sample_id, msi_output, safe_msi_dis]
        }
    
    // Interpret MSI results
    MSI_INTERPRET(msi_combined)
    ch_versions = ch_versions.mix(MSI_INTERPRET.out.versions)
    
    emit:
    result = MSI_INTERPRET.out.result           // [sample_id, json_result]
    summary = MSI_INTERPRET.out.summary         // [sample_id, summary_txt]
    raw_output = MSISENSOR_PRO_MSI.out.msi_output  // [sample_id, raw_msi_file]
    microsatellites = MSISENSOR_PRO_SCAN.out.microsatellites  // microsatellites.list
    versions = ch_versions
}