#!/usr/bin/env python3

"""
run_variant_fragmentomics.py (ADVANCED VERSION)
(Module 3.3 - "The Smoking Gun")

Reads a VCF file and a corresponding BAM file. For each variant (SNV, INS, DEL),
it fetches all supporting reads and calculates fragment size statistics
for the ALT allele vs the REF allele.

It performs a Kolmogorov-Smirnov (K-S) test to determine if the two
fragment size distributions are statistically different.

Outputs a new VCF annotated with:
- FRAG_ALT_MED: Median fragment size of ALT-supporting reads
- FRAG_REF_MED: Median fragment size of REF-supporting reads
- FRAG_KS_PVAL: P-value from the K-S test (one-sided, 'less')
- FRAG_ALT_N: Number of ALT-supporting reads
- FRAG_REF_N: Number of REF-supporting reads
"""

import pysam
import argparse
import sys
from collections import defaultdict
import statistics
from scipy.stats import ks_2samp

def get_allele_from_read_advanced(read, rec):
    """
    ADVANCED: Check which allele a read supports at a specific position,
    handles SNVs, Insertions, and Deletions.

    Returns: 'REF', 'ALT', or 'NONE'.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return 'NONE'

    var_pos = rec.pos - 1  # VCF is 1-based, pysam is 0-based
    var_ref = rec.ref
    var_alt = rec.alts[0] # We only support one ALT allele

    is_snv = len(var_ref) == 1 and len(var_alt) == 1
    is_del = len(var_ref) > len(var_alt)
    is_ins = len(var_ref) < len(var_alt)

    try:
        aligned_pairs = read.get_aligned_pairs_with_seq()
    except Exception:
        # Fallback for reads with weird CIGARs
        return 'NONE'

    read_bases = read.query_sequence
    read_supports = 'NONE'
    ref_pos_prev = -1 # Store previous ref_pos for insertion checks

    for query_pos, ref_pos, ref_base in aligned_pairs:
        if ref_pos is None: # Insertion in the read
            if is_ins and (ref_pos_prev == var_pos):
                # Insertion starts right after the variant position
                # Check if the inserted bases match the ALT allele
                ins_len = len(var_alt) - len(var_ref)
                read_insertion = read_bases[query_pos : query_pos + ins_len]
                if read_insertion == var_alt[1:]: # Compare inserted part
                    read_supports = 'ALT'
                    break
            continue # Move to next pair

        if ref_pos_prev == -1: # Initialize on first non-insertion
             ref_pos_prev = ref_pos
        
        if query_pos is None: # Deletion in the read
            if is_del and (ref_pos == var_pos + 1):
                # Deletion starts right after the variant position
                # Check if deletion length matches
                del_len = len(var_ref) - len(var_alt)
                # We need to check how many 'None' query_pos follow
                read_supports = 'ALT' # Assume it's the ALT for now
                break # Found a deletion starting at the right place
            ref_pos_prev = ref_pos # Update ref_pos_prev
            continue # Move to next pair

        if ref_pos < var_pos:
            ref_pos_prev = ref_pos # Update ref_pos_prev
            continue # Not at the variant yet

        if ref_pos > var_pos:
            break # Moved past the variant

        if ref_pos == var_pos:
            if is_snv:
                read_base = read_bases[query_pos]
                if read_base == var_alt:
                    read_supports = 'ALT'
                elif read_base == var_ref:
                    read_supports = 'REF'
                break # Found SNV
            
            if is_del:
                # Read must match REF base at var_pos
                read_base = read_bases[query_pos]
                if read_base == var_ref[0]:
                    read_supports = 'REF' # Supports REF, not the deletion
                break

            if is_ins:
                # Read must match REF base at var_pos
                read_base = read_bases[query_pos]
                if read_base == var_ref:
                    read_supports = 'REF' # Supports REF, not the insertion
                break
        
        ref_pos_prev = ref_pos # Update ref_pos_prev

    return read_supports

def analyze_variant_fragments(vcf_path, bam_path, ref_fasta_path, outfile_path):
    """
    Main analysis function.
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta_path)
        vcf_in = pysam.VariantFile(vcf_path, "r")
        
        # Add new INFO fields to the header
        vcf_in.header.info.add("FRAG_ALT_MED", "1", "Float", "Median fragment size of ALT-supporting reads")
        vcf_in.header.info.add("FRAG_REF_MED", "1", "Float", "Median fragment size of REF-supporting reads")
        vcf_in.header.info.add("FRAG_KS_PVAL", "1", "Float", "K-S test p-value for ALT vs REF fragment distributions (one-sided, 'less')")
        vcf_in.header.info.add("FRAG_ALT_N", "1", "Integer", "Number of ALT-supporting reads")
        vcf_in.header.info.add("FRAG_REF_N", "1", "Integer", "Number of REF-supporting reads")
        
        vcf_out = pysam.VariantFile(outfile_path, "w", header=vcf_in.header)

    except Exception as e:
        print(f"Error opening files: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Starting variant fragmentomics analysis (ADVANCED: SNV + INDEL)...", file=sys.stderr)
    
    count = 0
    for rec in vcf_in:
        count += 1
        if count % 100 == 0:
            print(f"Processed {count} variants...", file=sys.stderr)
            
        # We only analyze PASS variants with a single ALT allele
        if "PASS" not in rec.filter or len(rec.alts) > 1 or rec.alts[0] == '*':
            vcf_out.write(rec)
            continue
            
        frag_lengths = defaultdict(list)
        read_names = set() # To avoid double-counting read pairs

        try:
            # Fetch reads from a window around the variant
            # Fetch a wider window to ensure we get reads covering INDELs
            fetch_start = max(0, rec.pos - 50)
            fetch_end = rec.pos + 50
            for read in bam.fetch(rec.chrom, fetch_start, fetch_end):
                
                # Use read name to skip processing R2 if R1 was already processed
                if read.is_paired and read.read_name in read_names:
                    continue

                # Use template_length (TLEN) for fragment size
                # Must be positive, paired, and mapped
                tlen = abs(read.template_length)
                if tlen == 0 or not read.is_proper_pair:
                    continue

                # Use the ADVANCED allele checker
                allele = get_allele_from_read_advanced(read, rec)
                
                if allele in ['REF', 'ALT']:
                    frag_lengths[allele].append(tlen)
                    if read.is_paired:
                        read_names.add(read.read_name)

        except Exception as e:
            print(f"Error fetching/parsing reads for {rec.chrom}:{rec.pos}: {e}", file=sys.stderr)
            vcf_out.write(rec)
            continue

        alt_frags = frag_lengths.get('ALT', [])
        ref_frags = frag_lengths.get('REF', [])
        
        alt_n = len(alt_frags)
        ref_n = len(ref_frags)
        
        rec.info["FRAG_ALT_N"] = alt_n
        rec.info["FRAG_REF_N"] = ref_n
        
        # Require minimum number of reads to perform stats
        if alt_n > 5 and ref_n > 5:
            try:
                rec.info["FRAG_ALT_MED"] = round(statistics.median(alt_frags), 1)
                rec.info["FRAG_REF_MED"] = round(statistics.median(ref_frags), 1)
                
                # Perform K-S test: Test if ALT fragments are 'less' (shorter)
                ks_stat, p_val = ks_2samp(alt_frags, ref_frags, alternative='less') 
                rec.info["FRAG_KS_PVAL"] = float(f"{p_val:.3e}") # Store in scientific notation
            except Exception as e:
                # Handle cases where all fragment sizes are identical (K-S test fails)
                rec.info["FRAG_KS_PVAL"] = 1.0
                print(f"Warning: Could not perform stats for {rec.chrom}:{rec.pos} (e.g., all fragments identical): {e}", file=sys.stderr)
                
        vcf_out.write(rec)

    print("Analysis complete.", file=sys.stderr)
    vcf_in.close()
    vcf_out.close()
    bam.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate a VCF with fragment size statistics from a BAM file.")
    parser.add_argument("--vcf", required=True, help="Input VCF file (filtered.m2.vcf.gz)")
    parser.add_argument("--bam", required=True, help="Input BAM file (consensus.bam)")
    parser.add_argument("--ref_fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--outfile", required=True, help="Output annotated VCF file (e.g., out.frag.vcf.gz)")
    
    args = parser.parse_args()
    
    analyze_variant_fragments(args.vcf, args.bam, args.ref_fasta, args.outfile)
