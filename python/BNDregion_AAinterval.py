#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def extract_intervals_to_bed(input_path, aa_barcode=None, band_width=1000, refbuild='GRCh37'):
    """
    Extract intervals from cycle.txt files and safely expand them within chromosome bounds.
    
    Args:
        input_path (str): Path to directory containing cycle.txt files
        aa_barcode (str, optional): Specific barcode to filter files
        band_width (int): Width to expand intervals on each side
        refbuild (str): Reference genome build
    
    Returns:
        pandas.DataFrame: DataFrame of extracted and expanded intervals
    """
    # Chromosome lengths for GRCh37 and GRCh38
    GRCh37_lengths = {
        'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
        'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
        'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
        'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
        'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566
    }
    
    GRCh38_lengths = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    # Find all cycle.txt files
    cycle_files = []
    for root, dirs, files in os.walk(input_path):
        cycle_files.extend([
            os.path.join(root, file) 
            for file in files 
            if file.endswith('_cycles.txt')
        ])
    
    # Filter files by barcode if specified
    if aa_barcode:
        cycle_files = [f for f in cycle_files if aa_barcode in f]
    
    # Initialize list to store intervals
    all_intervals = []
    
    # Process each cycle file
    for file in cycle_files:
        try:
            # Read file lines
            with open(file, 'r') as f:
                lines = f.readlines()
            
            # Extract interval lines
            interval_lines = [
                line.strip() 
                for line in lines 
                if line.startswith('Interval')
            ]
            
            # Parse interval lines
            for line in interval_lines:
                parts = line.split('\t')
                chr = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                
                # Safely expand intervals within chromosome bounds
                chr_lengths = GRCh38_lengths if refbuild == 'GRCh38' else GRCh37_lengths
                chr_length = chr_lengths.get(chr, float('inf'))
                
                # Ensure expanded start is not negative and does not exceed chromosome length
                expanded_start = max(1, start - band_width)
                
                # Ensure expanded end does not exceed chromosome length
                expanded_end = min(chr_length, end + band_width)
                
                all_intervals.append({
                    'barcode': os.path.basename(os.path.dirname(file)),
                    'chr': chr,
                    'original_start': start,
                    'original_end': end,
                    'start': expanded_start,
                    'end': expanded_end,
                    'interval_index': int(parts[1])
                })
        
        except Exception as e:
            print(f"Error processing file {file}: {e}")
    
    # Convert to DataFrame
    return pd.DataFrame(all_intervals)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Extract and safely expand intervals from cycle.txt files')
    parser.add_argument('--aa_out_path', 
                        type=str, 
                        required=True, 
                        help='Path to directory containing cycle.txt files')
    parser.add_argument('--aa_barcode', 
                        type=str, 
                        default=None, 
                        help='Specific AA barcode to process')
    parser.add_argument('-o', '--out_dir', 
                        type=str, 
                        required=True, 
                        help='Directory where output BED file will be written')
    parser.add_argument('--refbuild', 
                        type=str, 
                        default='GRCh37', 
                        choices=['GRCh37', 'GRCh38'], 
                        help='Reference genome build')
    parser.add_argument('--band_width', 
                        type=int, 
                        default=1000, 
                        help='Width to expand intervals on each side (default: 1000)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Extract intervals
    intervals_df = extract_intervals_to_bed(
        input_path=args.aa_out_path, 
        aa_barcode=args.aa_barcode,
        band_width=args.band_width,
        refbuild=args.refbuild
    )
    
    # Ensure output directory exists
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Prepare output filename
    if args.aa_barcode:
        output_filename = os.path.join(args.out_dir, f"{args.aa_barcode}_AA_interval_{args.band_width}_BND.bed")
    else:
        output_filename = os.path.join(args.out_dir, "all_expanded_intervals.bed")
    
    # Write BED file (0-based start)
    bed_data = pd.DataFrame({
        'chr': intervals_df['chr'],
        'start': intervals_df['start'] - 1,  # Convert to 0-based start
        'end': intervals_df['end']
    })
    
    # Save to BED file
    bed_data.to_csv(
        output_filename, 
        sep='\t', 
        header=False, 
        index=False
    )
    
    # Print summary
    print(f"Reference build: {args.refbuild}")
    print(f"Band width: ±{args.band_width} bp")
    print(f"Extracted {len(intervals_df)} intervals.")
    print(f"Output written to: {output_filename}")
    
    return intervals_df

if __name__ == '__main__':
    main()
