#!/usr/bin/env python3
"""
Generate HTML summary report for all samples
"""

import pandas as pd
import os
from pathlib import Path
from jinja2 import Template

def load_all_reports(report_files):
    """Load and combine all quantification reports"""
    all_data = []
    
    for report_file in report_files:
        if os.path.exists(report_file):
            df = pd.read_csv(report_file, sep='\t')
            all_data.append(df)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df
    else:
        return pd.DataFrame()

def generate_summary_stats(df):
    """Generate summary statistics"""
    if df.empty:
        return {}
    
    stats = {
        'total_samples': df['sample'].nunique(),
        'total_snvs': df[df['marker_type'] == 'SNV']['count'].sum(),
        'total_indels': df[df['marker_type'] == 'INDEL']['count'].sum(),
        'total_svs': df[df['marker_type'] == 'SV']['count'].sum(),
        'msi_positive_samples': len(df[(df['marker_type'] == 'MSI') & (df['details'].str.contains('MSI-H|MSI-L'))]),
        'avg_snvs_per_sample': df[df['marker_type'] == 'SNV']['count'].mean(),
        'avg_indels_per_sample': df[df['marker_type'] == 'INDEL']['count'].mean(),
        'avg_svs_per_sample': df[df['marker_type'] == 'SV']['count'].mean()
    }
    
    return stats

def create_html_report(df, stats, output_file):
    """Create HTML summary report"""
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>UMI-Aware Targeted DNA-Sequencing Pipeline Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
            .stats { margin: 20px 0; }
            .stat-box { display: inline-block; margin: 10px; padding: 15px; background-color: #e8f4f8; border-radius: 5px; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #4CAF50; color: white; }
            .sample-group { margin: 20px 0; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>UMI-Aware Targeted DNA-Sequencing Pipeline Summary</h1>
            <p>Generated on: {{ timestamp }}</p>
        </div>
        
        <div class="stats">
            <h2>Overall Statistics</h2>
            <div class="stat-box">
                <strong>Total Samples:</strong> {{ stats.total_samples }}
            </div>
            <div class="stat-box">
                <strong>Total SNVs:</strong> {{ stats.total_snvs }}
            </div>
            <div class="stat-box">
                <strong>Total INDELs:</strong> {{ stats.total_indels }}
            </div>
            <div class="stat-box">
                <strong>Total SVs:</strong> {{ stats.total_svs }}
            </div>
            <div class="stat-box">
                <strong>MSI Positive Samples:</strong> {{ stats.msi_positive_samples }}
            </div>
        </div>
        
        <div class="stats">
            <h2>Average per Sample</h2>
            <div class="stat-box">
                <strong>Avg SNVs:</strong> {{ "%.1f"|format(stats.avg_snvs_per_sample) }}
            </div>
            <div class="stat-box">
                <strong>Avg INDELs:</strong> {{ "%.1f"|format(stats.avg_indels_per_sample) }}
            </div>
            <div class="stat-box">
                <strong>Avg SVs:</strong> {{ "%.1f"|format(stats.avg_svs_per_sample) }}
            </div>
        </div>
        
        <h2>Sample Details</h2>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Marker Type</th>
                    <th>Count</th>
                    <th>Details</th>
                </tr>
            </thead>
            <tbody>
                {% for index, row in df.iterrows() %}
