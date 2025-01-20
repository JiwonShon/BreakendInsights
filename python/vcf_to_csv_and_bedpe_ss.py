"""
=========================================================
    Script Name   : vcf_to_csv_and_bedpe_ss.py
    Author        : Jiwon Shon
    Created Date  : 2024-11-20
    Description   : Transform svaba.vcf to bedpe
=========================================================
"""
import pandas as pd
import sys
import os
import argparse

def process_vcf_to_csv_and_bedpe(svaba_file, object_id):
    #=== STEP 1: Read VCF and convert to DataFrame
    print(f">> Reading VCF file: {svaba_file}")
    if not os.path.exists(svaba_file) or os.stat(svaba_file).st_size == 0:
        svaba_bedpe_filename = f"{svaba_file}.bedpe"
        svaba_bedpe_csv_filename = f"{svaba_file}.bedpe.csv"
        with open(svaba_bedpe_filename, 'w') as f:
            pass  # Do nothing, just create the file
        with open(svaba_bedpe_csv_filename, 'w') as f:
            pass  # Do nothing, just create the file
        return

    with open(svaba_file, 'r') as f:
        lines = f.readlines()
    # Extract header and body
    header_line = ''
    for i, line in enumerate(lines):
        if line.startswith('#CHROM'):
            header_line = line.strip()
            body_start_idx = i + 1
            break

    columns = header_line.split('\t')
    columns[0] = 'CHROM'
    body_lines = [line.strip().split("\t") for line in lines[body_start_idx:]]
    
    vcf_df = pd.DataFrame(body_lines, columns=columns)

    vcf_df = vcf_df.rename(
        columns={
            columns[-2]: "NM_sample",
            columns[-1]: "TM_sample"
        }
    )
    print(f"VCF body shape: {vcf_df.shape}")


    #=== STEP2: Split INFO field into individual columns
    print(">> Splitting INFO fields...")
    for idx, row in vcf_df.iterrows():
        info_fields = row['INFO'].split(';')
        
        for field in info_fields:
            if '=' in field:
                key, value = field.split('=', 1)  # '=' 기준으로 key, value로 나눔
                column_name = f'INFO_{key}'
                
                if column_name not in vcf_df.columns:
                    # 새로운 컬럼을 생성하고 NaN으로 채움
                    vcf_df[column_name] = pd.NA
                    
                # 해당 컬럼에 값 할당
                vcf_df.at[idx, column_name] = value

    # NEW STEP 스킵 처리 추가
    if vcf_df.shape[0] == 0:
        print(f"!! No data in {svaba_file}. Skipping this file.")
        svaba_bedpe_filename = f"{svaba_file}.bedpe"
        svaba_bedpe_csv_filename = f"{svaba_file}.bedpe.csv"
        with open(svaba_bedpe_filename, 'w') as f:
            pass  # Do nothing, just create the file
        with open(svaba_bedpe_csv_filename, 'w') as f:
            pass  # Do nothing, just create the file
        return  # 이 파일을 건너뛰고 종료


    vcf_df = vcf_df.drop(columns=['INFO'])  # INFO 컬럼 삭제

    # Check if INFO_HOMSEQ exists, if not create it with NaN values
    if 'INFO_HOMSEQ' not in vcf_df.columns:
        print("INFO_HOMSEQ column not found. Creating column with NaN values.")
        vcf_df['INFO_HOMSEQ'] = pd.NA

    # Calculate INFO_HOMLEN based on INFO_HOMSEQ
    vcf_df['INFO_HOMLEN'] = vcf_df['INFO_HOMSEQ'].apply(lambda x: len(x) if pd.notna(x) else None).astype('Int64')

    # Ensure INFO_EVDNC is string for comparison
    vcf_df['INFO_EVDNC'] = vcf_df['INFO_EVDNC'].astype(str)  # Convert to string for comparison

    # Apply the logic: if INFO_EVDNC is 'ASDIS' and INFO_HOMSEQ is NA, set INFO_HOMLEN to 0
    vcf_df.loc[(vcf_df['INFO_EVDNC'] == 'ASDIS') & (vcf_df['INFO_HOMSEQ'].isna()), 'INFO_HOMLEN'] = 0  

    # Remove INFO_HOMLEN and insert it next to INFO_HOMSEQ
    columns = vcf_df.columns.tolist()
    info_homseq_index = columns.index('INFO_HOMSEQ')
    columns.remove('INFO_HOMLEN')
    columns.insert(info_homseq_index + 1, 'INFO_HOMLEN')
    vcf_df = vcf_df[columns]


    #=== STEP3: Split FORMAT field into individual columns by TM_sample and NM_sample
    format_cols = vcf_df['FORMAT'].iloc[0].split(':')  # FORMAT 컬럼의 첫 번째 값에서 분리
    tm_sample = 'TM_sample'
    nm_sample = 'NM_sample'

    # TM_sample
    tm_values = vcf_df[tm_sample].str.split(':', expand=True)
    tm_values.columns = [f"FORMAT_{name}_tm" for name in format_cols]  # 컬럼 이름 변경

    # NM_sample
    nm_values = vcf_df[nm_sample].str.split(':', expand=True)
    nm_values.columns = [f"FORMAT_{name}_nm" for name in format_cols]  # 컬럼 이름 변경

    # 기존 데이터프레임과 새로운 컬럼 결합
    vcf_df = pd.concat([vcf_df, tm_values, nm_values], axis=1)
    # 기존 FORMAT, TM_sample, NM_sample 컬럼 삭제
    vcf_df = vcf_df.drop(columns=['FORMAT', tm_sample, nm_sample])

    # Save CSV
    svaba_csv_filename = f"{svaba_file}.csv"
    vcf_df.to_csv(svaba_csv_filename, index=False)
    print(f"CSV saved: {svaba_csv_filename}")

    #=== STEP4: Convert to BEDPE-like format
    print(">> Converting to BEDPE-like format...")
    bedpe_rows = []

    # Process each row in the VCF DataFrame
    for _, row in vcf_df.iterrows():
        # Check if the ID ends with ':1' or ':2'
        if ":1" in row['ID']:
            mate_id = row['INFO_MATEID']
            pair_rows = vcf_df[vcf_df['ID'] == mate_id]  # Find the mate pair
            if not pair_rows.empty:
                mate_row = pair_rows.iloc[0]  # Assuming there is only one matching row

                # Combine the two rows into BEDPE format
                bedpe_row = {
                    "CHROM_1": row['CHROM'],
                    "POS_1": row['POS'],
                    "ID_1": row['ID'],
                    "REF_1": row['REF'],
                    "ALT_1": row['ALT'],
                    "CHROM_2": mate_row['CHROM'],
                    "POS_2": mate_row['POS'],
                    "ID_2": mate_row['ID'],
                    "REF_2": mate_row['REF'],
                    "ALT_2": mate_row['ALT']
                }

                # Add common columns from the original VCF
                for col in vcf_df.columns:
                    if col not in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'MATEID']:
                        bedpe_row[col] = row[col]  # Common columns

                bedpe_rows.append(bedpe_row)

    bedpe_df = pd.DataFrame(bedpe_rows)

    # Save BEDPE file
    svaba_bedpe_filename = f"{svaba_file}.bedpe"
    bedpe_df.to_csv(svaba_bedpe_filename, index=False)
    print(f"BEDPE-like file saved: {svaba_bedpe_filename}")

    svaba_bedpe_csv_filename = f"{svaba_file}.bedpe.csv"
    bedpe_df.to_csv(svaba_bedpe_csv_filename, index=False)
    print(f"BEDPE-like file as csv saved: {svaba_bedpe_csv_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VCF to CSV and BEDPE format")
    parser.add_argument("svaba_file", help="Input VCF file")
    parser.add_argument("object_id", help="Object ID for output filenames")
    args = parser.parse_args()

    process_vcf_to_csv_and_bedpe(args.svaba_file, args.object_id)
