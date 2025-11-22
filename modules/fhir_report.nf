// modules/fhir_report.nf - FHIR-Compliant Reporting Module
// Generates FHIR R4 Observation resources for clinical reporting

nextflow.enable.dsl=2

process GENERATE_FHIR_REPORT {
    tag "$sample_id"
    label 'process_low'
    container '965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest'
    publishDir "${params.outdir}/fhir", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(clinical_report)
    tuple val(sample_id), path(tumor_fraction)
    tuple val(sample_id), path(msi_result)
    tuple val(sample_id), path(tmb_result)
    
    output:
    tuple val(sample_id), path("${sample_id}_fhir_bundle.json"), emit: fhir_bundle
    
    script:
    """
    python << 'EOF'
import json
import sys
from datetime import datetime
from uuid import uuid4

def create_fhir_observation(code, value, unit=None, interpretation=None):
    \"\"\"Create FHIR Observation resource\"\"\"
    obs = {
        "resourceType": "Observation",
        "id": str(uuid4()),
        "status": "final",
        "code": {
            "coding": [code]
        },
        "valueQuantity": {
            "value": value
        }
    }
    
    if unit:
        obs["valueQuantity"]["unit"] = unit
        obs["valueQuantity"]["system"] = "http://unitsofmeasure.org"
        obs["valueQuantity"]["code"] = unit
    
    if interpretation:
        obs["interpretation"] = [{"coding": [interpretation]}]
    
    return obs

# Load clinical data
clinical_data = {}
try:
    with open("${clinical_report}", 'r') as f:
        clinical_data = json.load(f)
except:
    pass

tumor_fraction_value = 0.0
try:
    with open("${tumor_fraction}", 'r') as f:
        tumor_fraction_value = float(f.read().strip())
except:
    pass

msi_data = {}
try:
    with open("${msi_result}", 'r') as f:
        msi_data = json.load(f)
except:
    pass

tmb_data = {}
try:
    with open("${tmb_result}", 'r') as f:
        tmb_data = json.load(f)
except:
    pass

# Create FHIR Bundle
bundle = {
    "resourceType": "Bundle",
    "id": str(uuid4()),
    "type": "collection",
    "timestamp": datetime.now().isoformat(),
    "entry": []
}

# Tumor Fraction Observation
tf_obs = create_fhir_observation(
    {
        "system": "http://loinc.org",
        "code": "33747-0",
        "display": "Circulating tumor DNA fraction"
    },
    tumor_fraction_value * 100,  # Convert to percentage
    unit="%"
)
bundle["entry"].append({"resource": tf_obs})

# MSI Status Observation
if msi_data.get("msi_status"):
    msi_obs = create_fhir_observation(
        {
            "system": "http://loinc.org",
            "code": "94035-7",
            "display": "Microsatellite instability status"
        },
        1 if msi_data["msi_status"] == "MSI-H" else 0,
        interpretation={
            "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
            "code": "POS" if msi_data["msi_status"] == "MSI-H" else "NEG",
            "display": msi_data["msi_status"]
        }
    )
    bundle["entry"].append({"resource": msi_obs})

# TMB Observation
if tmb_data.get("tmb_per_mb"):
    tmb_obs = create_fhir_observation(
        {
            "system": "http://loinc.org",
            "code": "94034-0",
            "display": "Tumor mutational burden"
        },
        tmb_data["tmb_per_mb"],
        unit="mutations/Mb"
    )
    bundle["entry"].append({"resource": tmb_obs})

# Variant Observations
if clinical_data.get("small_variants", {}).get("variants"):
    for variant in clinical_data["small_variants"]["variants"][:10]:  # Limit to top 10
        variant_obs = {
            "resourceType": "Observation",
            "id": str(uuid4()),
            "status": "final",
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": "48019-4",
                    "display": "Genetic variant assessment"
                }]
            },
            "valueString": f"{variant.get('gene', 'Unknown')}: {variant.get('hgvs_protein', 'Unknown')}",
            "interpretation": [{
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                    "code": "POS" if "pathogenic" in variant.get("clinical_significance", "").lower() else "IND",
                    "display": variant.get("clinical_significance", "Unknown")
                }]
            }]
        }
        bundle["entry"].append({"resource": variant_obs})

# Write FHIR Bundle
with open("${sample_id}_fhir_bundle.json", 'w') as f:
    json.dump(bundle, f, indent=2)

print(f"FHIR Bundle generated for ${sample_id}")
print(f"  Observations: {len(bundle['entry'])}")
EOF
    """
}

// FHIR Reporting Workflow
workflow FHIR_REPORT {
    take:
    clinical_report    // channel: [sample_id, clinical_report_json]
    tumor_fraction     // channel: [sample_id, tumor_fraction_file]
    msi_result         // channel: [sample_id, msi_result_json]
    tmb_result         // channel: [sample_id, tmb_result_json]
    
    main:
    // Combine all inputs
    combined = clinical_report
        .join(tumor_fraction, by: 0)
        .join(msi_result, by: 0)
        .join(tmb_result, by: 0)
    
    // Generate FHIR report
    GENERATE_FHIR_REPORT(
        combined.map { sample_id, report, tf, msi, tmb -> [sample_id, report] },
        combined.map { sample_id, report, tf, msi, tmb -> [sample_id, tf] },
        combined.map { sample_id, report, tf, msi, tmb -> [sample_id, msi] },
        combined.map { sample_id, report, tf, msi, tmb -> [sample_id, tmb] }
    )
    
    emit:
    fhir_bundle = GENERATE_FHIR_REPORT.out.fhir_bundle
}

