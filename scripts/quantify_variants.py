#!/usr/bin/env python3
"""
Quantify variants from VCF files and MSI results
"""

import pandas as pd
import pysam
import os
from pathlib import Path

def parse_vcf_variants(vcf_file):
    """Parse VCF file and extract variant information"""
    variants = []
    
    if not os.path.exists(vcf_file):
        return variants
    
    try:
        vcf = pysam.VariantFile(vcf_file)
        
        for record in vcf.fetch():
            # Extract basic variant information
            variant_info = {
                'chromosome': record.chrom,
                'position': record.pos,
                'reference': record.ref,
                'alternate': ','.join([alt for alt in record.alts]) if record.alts else '',
                'quality': record.qual,
                'filter': ','.join(record.filter) if record.filter else 'PASS'
            }
            
            # Extract sample-specific information
            if record.samples:
                sample = list(record.samples.keys())[0]
                sample_data = record.samples[sample]
                
                # Get depth and allele frequency if available
                if 'DP' in sample_data:
                    variant_info['depth'] = sample_data['DP']
                else:
                    variant_info['depth'] = 0
                
                if 'AD' in sample_data and sample_data['AD']:
                    ref_depth = sample_data['AD'][0]
                    alt_depth = sum(sample_data['AD'][1:])
                    total_depth = ref_depth + alt_depth
                    if total_depth > 0:
                        variant_info['allele_frequency'] = alt_depth / total_depth
                    else:
                        variant_info['allele_frequency'] = 0
                else:
                    variant_info['allele_frequency'] = 0
                
                # Get genotype
                if 'GT' in sample_data:
                    variant_info['genotype'] = sample_data['GT']
                else:
                    variant_info['genotype'] = './.'
            
            # Determine variant type
            if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts):
                variant_info['type'] = 'SNV'
            else:
                variant_info['type'] = 'INDEL'
            
            variants.append(variant_info)
            
    except Exception as e:
        print(f"Error parsing VCF file {vcf_file}: {e}")
    
    return variants

def parse_sv_variants(sv_vcf_file):
    """Parse structural variant VCF file"""
    sv_variants = []
    
    if not os.path.exists(sv_vcf_file):
        return sv_variants
    
    try:
        vcf = pysam.VariantFile(sv_vcf_file)
        
        for record in vcf.fetch():
            sv_info = {
                'chromosome': record.chrom,
                'position': record.pos,
                'sv_type': record.info.get('SVTYPE', 'UNKNOWN'),
                'sv_length': record.info.get('SVLEN', 0),
                'quality': record.qual,
                'filter': ','.join(record.filter) if record.filter else 'PASS'
            }
            
            # Get end position if available
            if 'END' in record.info:
                sv_info['end_position'] = record.info['END']
            else:
                sv_info['end_position'] = record.pos
            
            sv_variants.append(sv_info)
            
    except Exception as e:
        print(f"Error parsing SV VCF file {sv_vcf_file}: {e}")
    
    return sv_variants

def parse_msi_results(msi_file):
    """Parse MSI results file"""
    msi_info = {
        'msi_score': 0,
        'msi_status': 'Unknown',
        'total_sites': 0,
        'unstable_sites': 0
    }
    
    if not os.path.exists(msi_file):
        return msi_info
    
    try:
        with open(msi_file, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            line = line.strip()
            if line.startswith('Total_Number_of_Sites'):
                msi_info['total_sites'] = int(line.split('\t')[1])
            elif line.startswith('Number_of_Somatic_Sites'):
                msi_info['unstable_sites'] = int(line.split('\t')[1])
            elif line.startswith('MSI_Score'):
                msi_info['msi_score'] = float(line.split('\t')[1])
            elif line.startswith('MSI_Status'):
                msi_info['msi_status'] = line.split('\t')[1]
                
    except Exception as e:
        print(f"Error parsing MSI file {msi_file}: {e}")
    
    return msi_info

def main():
    # Get input files from snakemake
    snv_vcf = snakemake.input.snv_vcf
    sv_vcf = snakemake.input.sv_vcf
    msi_file = snakemake.input.msi_txt
    output_file = snakemake.output.report
    
    # Extract sample name from file path
    sample_name = Path(snv_vcf).stem.replace('_variants', '')
    
    # Parse variant files
    snv_variants = parse_vcf_variants(snv_vcf)
    sv_variants = parse_sv_variants(sv_vcf)
    msi_info = parse_msi_results(msi_file)
    
    # Create summary report
    report_data = []
    
    # Add SNV/INDEL summary
    snv_count = len([v for v in snv_variants if v['type'] == 'SNV'])
    indel_count = len([v for v in snv_variants if v['type'] == 'INDEL'])
    
    report_data.append({
        'sample': sample_name,
        'marker_type': 'SNV',
        'count': snv_count,
        'details': f"Total SNVs detected: {snv_count}"
    })
    
    report_data.append({
        'sample': sample_name,
        'marker_type': 'INDEL',
        'count': indel_count,
        'details': f"Total INDELs detected: {indel_count}"
    })
    
    # Add SV summary
    sv_count = len(sv_variants)
    report_data.append({
        'sample': sample_name,
        'marker_type': 'SV',
        'count': sv_count,
        'details': f"Total SVs detected: {sv_count}"
    })
    
    # Add MSI summary
    report_data.append({
        'sample': sample_name,
        'marker_type': 'MSI',
        'count': msi_info['msi_score'],
        'details': f"MSI Status: {msi_info['msi_status']}, Score: {msi_info['msi_score']}"
    })
    
    # Save report
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df = pd.DataFrame(report_data)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Quantification report saved to: {output_file}")

if __name__ == "__main__":
    main()
