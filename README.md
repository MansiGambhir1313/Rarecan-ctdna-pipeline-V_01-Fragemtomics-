# RareCan v1.0 Clinical-Grade ctDNA Pipeline

A comprehensive Nextflow-based pipeline for clinical-grade circulating tumor DNA (ctDNA) analysis with UMI-based error correction, designed for AWS HealthOmics deployment.

## 🎯 Project Status

**Current Status:** ~75% Complete - Core Pipeline Functional, Fragmentomics Implemented

**See:** [Project Analysis](docs/PROJECT_ANALYSIS.md) for detailed status

## 📋 Features

### ✅ Implemented
- Quality Control (FastP, FastQC, MultiQC)
- Alignment (Minimap2 with BWA fallback)
- Variant Calling (Mutect2, VarDict, LoFreq)
- CNV Analysis (CNVkit)
- Structural Variant Detection (Manta)
- MSI Analysis (MSIsensor2)
- TMB Calculation
- Variant Annotation (VEP)
- Clinical Reporting (JSON/HTML)
- **Fragmentomics Analysis** (Global + Variant signatures)
- **TFX Clinical Report Generation**
- **FHIR-Compliant Reporting**

### ⚠️ In Progress / Missing
- Complete UMI Consensus Pipeline (partial implementation)
- ichorCNA Integration (tumor fraction estimation)
- CHIP Filtering (module exists, needs integration)
- Validation Framework

## 📁 Project Structure

```
Rarecan_final/
├── main.nf                 # Main Nextflow workflow
├── nextflow.config         # Pipeline configuration
├── modules/                # Nextflow modules
│   ├── qc.nf              # Quality control
│   ├── umi.nf              # UMI processing
│   ├── align.nf            # Alignment (BWA)
│   ├── align_minimap2.nf   # Alignment (Minimap2)
│   ├── snv.nf              # SNV/INDEL calling
│   ├── cnv.nf              # Copy number analysis
│   ├── sv.nf               # Structural variants
│   ├── msi.nf              # MSI analysis
│   ├── tmb.nf              # TMB calculation
│   ├── annotate.nf         # Variant annotation
│   ├── fragmentomics.nf   # Fragmentomics analysis
│   ├── tfx_report.nf       # TFX clinical reports
│   ├── fhir_report.nf      # FHIR reporting
│   └── report.nf           # Clinical reporting
├── bin/                    # Analysis scripts
│   ├── global_fragmentomics.py      # Fragment length/motif analysis
│   ├── run_variant_fragmentomics.py # Variant-level fragmentomics
│   ├── run_final_report.py          # TFX report generator
│   ├── run_ichorcna.R               # ichorCNA script
│   └── [other analysis scripts]
├── assets/                 # Reference files, databases
├── config/                 # Configuration files
├── docs/                   # Documentation
│   ├── PROJECT_ANALYSIS.md
│   └── IMPROVEMENTS_NEEDED.md
└── validation/             # Validation datasets and workflows
```

## 🚀 Quick Start

### Prerequisites
- Nextflow >= 22.10.0
- AWS HealthOmics access (for cloud execution)
- Docker or Singularity (for container execution)

### Basic Usage

```bash
nextflow run main.nf \
    --read1 sample_R1.fastq.gz \
    --read2 sample_R2.fastq.gz \
    --sample_id SAMPLE001 \
    --ref s3://references/hg38.fasta \
    --bed targets.bed \
    --pon_vcf pon.vcf.gz \
    --cnv_pon reference.cnn \
    --outdir /output
```

### Parameters

**Required:**
- `--read1`, `--read2`: Input FASTQ files
- `--sample_id`: Sample identifier
- `--ref`: Reference genome FASTA
- `--bed`: Target regions BED file
- `--pon_vcf`: Panel of Normals VCF
- `--cnv_pon`: CNV reference CNN file

**Optional:**
- `--enable_lofreq`: Enable LoFreq caller (default: false)
- `--enable_duplex`: Enable duplex UMI consensus (default: true)
- `--enable_fragmentomics`: Enable fragmentomics analysis (default: false)
- `--enable_umi`: Enable UMI processing (default: true)
- `--enable_cnv`: Enable CNV analysis (default: true)
- `--enable_sv`: Enable structural variant detection (default: true)
- `--msi_bed`: MSI loci BED file
- `--vep_cache`: VEP cache directory
- `--frag_min_bp`: Minimum fragment length for short fragments (default: 90)
- `--frag_max_bp`: Maximum fragment length for short fragments (default: 150)

## 🔧 Configuration

### AWS HealthOmics

The pipeline is configured for AWS HealthOmics deployment. Update `nextflow.config` with your:
- ECR registry: `965747689553.dkr.ecr.eu-west-2.amazonaws.com/ctdna-universal:latest`
- AWS Batch queue
- S3 bucket paths

### Container

Currently uses a single universal container. See [IMPROVEMENTS_NEEDED.md](docs/IMPROVEMENTS_NEEDED.md) for planned one-process-one-container strategy.

## 📊 Output Structure

```
output/
├── qc/                    # Quality control reports
├── umi/                   # UMI processing outputs
├── alignment/             # Alignment BAMs and metrics
├── fragmentomics/         # Fragmentomics analysis
│   ├── *.fragmentomics.summary.json  # Summary metrics
│   ├── *.global_hist.tsv            # Fragment length histogram
│   └── *.fragment_motifs.tsv        # End motif frequencies
├── variants/
│   ├── snv/              # SNV/INDEL VCFs
│   ├── cnv/              # CNV calls
│   └── sv/               # Structural variants
├── msi/                  # MSI analysis results
├── tmb/                  # TMB calculations
├── annotation/            # Annotated variants
└── reports/              # Clinical reports (JSON/HTML)
```

## ⚠️ Known Issues & Limitations

1. **UMI Processing**: Module exists but not fully integrated into main workflow
2. **ichorCNA**: Module exists but requires proper container configuration
3. **CHIP Filtering**: Module exists, needs full integration
4. **Reference Genome**: Currently using placeholder reference - requires full hg38 for production
5. **Licensing**: Manta, MSIsensor2 require commercial licenses for production use

See [IMPROVEMENTS_NEEDED.md](docs/IMPROVEMENTS_NEEDED.md) for detailed list.

## 📚 Documentation

- [Project Analysis](docs/PROJECT_ANALYSIS.md) - Comprehensive status and gap analysis
- [Improvements Needed](docs/IMPROVEMENTS_NEEDED.md) - Critical fixes and enhancements

## 🔬 Clinical Validation

**Status:** Validation framework not yet implemented

**Planned:**
- Tier 1: Analytical validation (SafeMut, HG002)
- Tier 2: Scientific validation (MSK-IMPACT concordance)
- Tier 3: Clinical validation (EGA-4847 dataset)

## 📝 License & Compliance

**Current Status:** ⚠️ License compliance issues exist

- Manta: Polyform Strict License (requires Illumina negotiation)
- MSIsensor2: GPL-3.0 (should use MSIsensor-pro for commercial)
- ichorCNA: GPL-3.0 (requires isolation)

**Action Required:** Resolve licensing before commercial deployment

## 🤝 Contributing

This is a clinical-grade pipeline under active development. See project documentation for contribution guidelines.

## 📧 Contact

For questions or issues, refer to project documentation or contact the development team.

## 🔄 Version History

- **v1.0.0** (Current): Foundation implementation with core modules
  - Basic variant calling pipeline
  - Partial UMI implementation
  - Clinical reporting framework

---

**Note:** This pipeline is designed for clinical use and requires proper validation and regulatory compliance before deployment in clinical settings.

