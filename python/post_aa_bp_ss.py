"""
=========================================================
    Script Name   : post_aa_bp_ss.py
    Author        : Jiwon Shon
    Created Date  : 2025-01-13
    last edit     : 2025-01-15
    Description   : Integrates breakend data from Amplicon Architect, SvABA, and BP_Seq_Seek.
=========================================================
"""
import os
import pandas as pd
import numpy as np
import argparse
import re

def create_empty_files(aliquot_barcode, filenames):
    for filename in filenames:
        filepath = os.path.join(f"{aliquot_barcode}{filename}.csv")
        with open(filepath, 'w') as f:
            pass  # Just create an empty file
        print(f"Created empty file: {filepath}")

def split_and_align(df, cols_to_split):
    for col in cols_to_split:
        df[col] = df[col].apply(lambda x: x.split(',') if isinstance(x, str) else [])
    max_length = df[cols_to_split].applymap(len).max(axis=1)
    # max_length = df[cols_to_split].map(lambda col: col.apply(lambda x: len(x) if isinstance(x, str) else 0)).max(axis=1)
    for col in cols_to_split:
        df[col] = df.apply(lambda row: row[col] + [None] * (max_length[row.name] - len(row[col])), axis=1)
    df = df.explode(cols_to_split, ignore_index=True)
    return df

def process_info_sctg(df):
    result_rows = [] 
    for _, row in df.iterrows():
        sctg_1bp = row['SCTG_1bp']
        sctg_10bp = row['SCTG_10bp']
        # Case 1: 1bp와 10bp 값이 동일할 경우
        if sctg_1bp == sctg_10bp and pd.notna(sctg_1bp):
            result_rows.append({
                "INFO_SCTG": sctg_1bp,
                "INFO_SCTG_bp": "1bp",
                "aa_barcode": row['aa_barcode'],
                "amplicon_idx": row['amplicon_idx'],
                "IsCyclicPath": row['IsCyclicPath'],
                "CycleClass":row['CycleClass']
            })
        # Case 2: 10bp에만 값이 있고 1bp는 없는 경우
        elif pd.isna(sctg_1bp) and pd.notna(sctg_10bp):
            result_rows.append({
                "INFO_SCTG": sctg_10bp,
                "INFO_SCTG_bp": "10bp",
                "aa_barcode": row['aa_barcode'],
                "amplicon_idx": row['amplicon_idx'],
                "IsCyclicPath": row['IsCyclicPath'],
                "CycleClass":row['CycleClass']
            })
        # Case 3: 1bp에만 값이 있고 10bp는 없는 경우
        elif pd.notna(sctg_1bp) and pd.isna(sctg_10bp):
            result_rows.append({
                "INFO_SCTG": sctg_1bp,
                "INFO_SCTG_bp": "1bp",
                "aa_barcode": row['aa_barcode'],
                "amplicon_idx": row['amplicon_idx'],
                "IsCyclicPath": row['IsCyclicPath'],
                "CycleClass":row['CycleClass']
            })
        # Case 4: 1bp와 10bp가 모두 다를 경우
        elif pd.notna(sctg_1bp) and pd.notna(sctg_10bp):
            result_rows.append({
                "INFO_SCTG": sctg_1bp,
                "INFO_SCTG_bp": "1bp",
                "aa_barcode": row['aa_barcode'],
                "amplicon_idx": row['amplicon_idx'],
                "IsCyclicPath": row['IsCyclicPath'],
                "CycleClass":row['CycleClass']
            })
            result_rows.append({
                "INFO_SCTG": sctg_10bp,
                "INFO_SCTG_bp": "10bp",
                "aa_barcode": row['aa_barcode'],
                "amplicon_idx": row['amplicon_idx'],
                "IsCyclicPath": row['IsCyclicPath'],
                "CycleClass":row['CycleClass']
            })
    return pd.DataFrame(result_rows)

def match_segedge_index(svaba, segedges, span, pos_col, chrom_col):
    # print(svaba)
    matching_rows_start = segedges[
        (segedges["seqnames"] == svaba[chrom_col]) &
        ((segedges["start"] - span <= svaba[pos_col]) & (svaba[pos_col] <= segedges["start"] + span))
    ]
    matching_rows_end = segedges[
        (segedges["seqnames"] == svaba[chrom_col]) &
        ((segedges["end"] - span <= svaba[pos_col]) & (svaba[pos_col] <= segedges["end"] + span))
    ]
    # print(f"matching_rows_start, {matching_rows_start}")
    # print(f"matching_rows_end, {matching_rows_end}")
    
    matching_rows_start = matching_rows_start.copy()
    matching_rows_end = matching_rows_end.copy()
    matching_rows_start["SegEdges_idx"] += "_start"
    matching_rows_end["SegEdges_idx"] += "_end"
    matching_rows = pd.concat([matching_rows_start, matching_rows_end])

    if not matching_rows.empty:
        return matching_rows["SegEdges_idx"].tolist()
    return []

def extract_unique_amplicons(row, columns):
    amplicons = set()
    for col in columns:
        col_values = row[col] if isinstance(row[col], list) else []  
        for value in col_values:  
            amp = value.split('_')[0]  
            amplicons.add(amp)
    return ','.join(sorted(amplicons)) if amplicons else None  

def update_amplicon_idx(row):
    if pd.isna(row['amplicon_idx']):  
        return row['SegEdge_amplicon'] if row['SegEdge_amplicon'] is not None else None  
    return row['amplicon_idx']  

def update_nan_values(df):
    df.loc[df['amplicon_idx'].notna() & df['amplicon_decomposition_class'].isna(), 'amplicon_decomposition_class'] = 'No amp/Invalid'
    # Update 'ecDNA' and 'BFB' columns
    for col in ['ecDNA', 'BFB']:
        df.loc[df['amplicon_idx'].notna() & df[col].isna(), col] = 'None detected'
    # Update patient-related columns with unique values
    columns_to_fill = ['patient_barcode', 'cohort_barcode', 'project', 'primaryTumorLocation', 'source', 'tm_cov',
                        'tumor_purity', 'all_cutoff_passed']
    for col in columns_to_fill:
        unique_value = df[col].dropna().unique()
        if len(unique_value) == 1:
            df[col].fillna(unique_value[0], inplace=True)
    return df

def find_overlaps(row,tsi_l):
    if tsi_l is None or tsi_l.empty:
        return pd.Series({
            "overlaped_TSI_L": None,
            "num_overlaped_TSI_L": 0
        })

    overlaps = tsi_l[
        ((tsi_l["CHROM_1"] == row["CHROM_1"]) & (tsi_l["POS_1"] == row["POS_1"])) |
        ((tsi_l["CHROM_2"] == row["CHROM_2"]) & (tsi_l["POS_2"] == row["POS_2"])) |
        ((tsi_l["CHROM_1"] == row["CHROM_2"]) & (tsi_l["POS_1"] == row["POS_2"])) |
        ((tsi_l["CHROM_2"] == row["CHROM_1"]) & (tsi_l["POS_2"] == row["POS_1"])) 
    ]
    if overlaps.empty:
        return pd.Series({
            "overlaped_TSI_L": None,
            "num_overlaped_TSI_L": 0
        })

    ids = []
    for _, overlap_row in overlaps.iterrows():
        if (overlap_row["CHROM_1"] == row["CHROM_1"] and overlap_row["POS_1"] == row["POS_1"]) or \
           (overlap_row["CHROM_1"] == row["CHROM_2"] and overlap_row["POS_1"] == row["POS_2"]):
            ids.append(overlap_row["ID_1"])
        elif (overlap_row["CHROM_2"] == row["CHROM_1"] and overlap_row["POS_2"] == row["POS_1"]) or \
             (overlap_row["CHROM_2"] == row["CHROM_2"] and overlap_row["POS_2"] == row["POS_2"]):
            ids.append(overlap_row["ID_2"])

    return pd.Series({
        "overlaped_TSI_L": ",".join(ids) if ids else None,
        "num_overlaped_TSI_L": len(ids)
    })


def process_files(aliquot_barcode, bedpe_path, bps_path, bpss_output, gr4_path, amp_path):
    # load data
    output_filenames = [
        ".svaba.somatic.sv.bedpe.aa.unfiltered.breakend",
        ".svaba.somatic.sv.bedpe.aa.TSI_G_filtered.breakend",
        ".svaba.somatic.sv.bedpe.aa.breakend",
    ]
    input_files = {
        "bedpe_path": bedpe_path,
        "bps_path": bps_path,
        # "bpss_output": bpss_output,
        "gr4_path": gr4_path,
        "amp_path": amp_path,
    }
    for name, path in input_files.items():
        if not os.path.exists(path) or os.stat(path).st_size == 0:
            print(f"Skipping {name}: {path} (File is empty or does not exist)")
            create_empty_files(aliquot_barcode, output_filenames)
            return

    bedpe = pd.read_csv(bedpe_path)
    bps_txt=pd.read_csv(bps_path,sep='\t')
    # bpss_output=pd.read_csv(bpss_output,sep='\t')
    gr4=pd.read_csv(gr4_path,sep='\t')
    amp=pd.read_csv(amp_path,sep='\t')
    

    # ==== step 1: bspp_output reproducing
    if not os.path.exists(bpss_output) or os.stat(bpss_output).st_size == 0:
        print(f"Warning: {bpss_output} is empty or does not exist.")
        bpss_output = pd.DataFrame()
        bpss_output_filtered = pd.DataFrame({
            "INFO_SCTG": [None],
            "INFO_SCTG_bp": [None],
            "aa_barcode": [None],
            "amplicon_idx": [None],
            "IsCyclicPath": [None],   
            "CycleClass": [None],  
            "AA_discordants": [None],  
        })
        
    else: 
        bpss_output=pd.read_csv(bpss_output,sep='\t')
        
        bpss_output['aa_barcode'] = bpss_output['file_name'].apply(lambda x: x.split('/')[-1].split('_amplicon')[0])
        bpss_output['amplicon_idx'] = bpss_output['file_name'].apply(lambda x: x.split('/')[-1].split('_')[-3])
        bpss_output_extract = bpss_output[['SCTG_1bp', 'SCTG_10bp', 'aa_barcode', 'amplicon_idx','IsCyclicPath','CycleClass']].drop_duplicates()
        bpss_output_extract['IsCyclicPath'] = bpss_output_extract['IsCyclicPath'].str.split('=').str[1]
        bpss_output_extract['CycleClass'] = bpss_output_extract['CycleClass'].str.split('=').str[1]
        print(bpss_output.shape, bpss_output_extract.shape)

        sctg_bps= split_and_align(bpss_output_extract, ["SCTG_1bp", "SCTG_10bp"])
        
        processed_sctg_bps = process_info_sctg(sctg_bps)

        if processed_sctg_bps.empty:
            print("processed_sctg_bps_drop_dup is empty. Creating a default DataFrame.")
            bpss_output_filtered = pd.DataFrame({
                "INFO_SCTG": [None],
                "INFO_SCTG_bp": [None],
                "aa_barcode": [None],
                "amplicon_idx": [None],
                "IsCyclicPath": [None],   
                "CycleClass": [None],  
                "AA_discordants": [None],  
            })
        else:
            processed_sctg_bps_drop_dup=processed_sctg_bps.drop_duplicates()

            bpss_output_filtered = processed_sctg_bps_drop_dup.sort_values(by="IsCyclicPath", ascending=False).drop_duplicates(subset=["INFO_SCTG"], keep="first")
            bpss_output_filtered["AA_discordants"] = 'True'



    # ==== Step 2: bps_txt strand annotation to bedpe_path
    print(f"bedpe: {bedpe.shape} \n bps_txt: {bps_txt.shape}")
    df1=bps_txt[['chr1', 'pos1', 'strand1', 'chr2', 'pos2', 'strand2','contig']].drop_duplicates()
    df1['chr1'] = df1['chr1'].astype('str')
    df1['chr2'] = df1['chr2'].astype('str')
    df1['pos1'] = df1['pos1'].astype('Int64')
    df1['pos2'] = df1['pos2'].astype('Int64')
    df2=bedpe.copy()
    df2['CHROM_2'] = df2['CHROM_2'].astype('str')
    df2['CHROM_1'] = df2['CHROM_1'].astype('str')

    merged_df=pd.merge(df1,df2,left_on=['chr1', 'pos1','chr2', 'pos2','contig'],right_on=['CHROM_1', 'POS_1', 'CHROM_2', 'POS_2', 'INFO_SCTG'], how='right')
    merged_df = merged_df.rename(columns={'strand1': 'STRAND_1', 'strand2': 'STRAND_2'})
    columns = merged_df.columns.tolist()
    pos_1_index = columns.index('POS_1')
    print('pos_1_index', pos_1_index)
    columns.insert(pos_1_index, columns.pop(columns.index('STRAND_1')))
    pos_2_index = columns.index('POS_2')
    columns.insert(pos_2_index, columns.pop(columns.index('STRAND_2')))

    merged_df = merged_df[columns]
    merged_df = merged_df.drop(columns=['chr1', 'pos1', 'chr2', 'pos2','contig'])

    bedpe_update=merged_df.copy()
    print(merged_df.shape)
    print(bedpe_update.shape)


    # ==== Step3:  bedpe + bp annotation
    bedpe_bpss = pd.merge(bedpe_update, bpss_output_filtered, on='INFO_SCTG', how='left')
    bedpe_bpss["AA_discordants"] = bedpe_bpss["AA_discordants"].apply(
        lambda x: "True" if x == "True" else "False"
    )
    print(bedpe_update.shape)
    print(bedpe_bpss.shape)


    # ==== Step4: gr4 -> segment annotation (True/False)
    gr4_pair = gr4[gr4['aa_barcode'] == aliquot_barcode].copy()
    if gr4_pair.shape[0] == 0:
        print(f"Skipping sample with barcode {aliquot_barcode} as no matching entries were found in GR4.")
        create_empty_files(aliquot_barcode, output_filenames)
        return
    gr4_pair['SegEdges_idx'] = gr4_pair['amplicon_idx'] + "_" + gr4_pair['idx'].astype(str)
    gr4_sub = gr4_pair[['seqnames', 'start', 'end', 'SegEdges_idx']].copy()
    gr4_sub['seqnames'] = gr4_sub['seqnames'].astype('str')
    gr4_sub['start'] = gr4_sub['start'].astype('Int64')
    gr4_sub['end'] = gr4_sub['end'].astype('Int64')

    for span in [1, 10]:
        bedpe_bpss[f"POS_1_SegEdges_idx_{span}bp"] = bedpe_bpss.apply(lambda row: match_segedge_index(row, gr4_sub, span, "POS_1", "CHROM_1"), axis=1)
        bedpe_bpss[f"POS_2_SegEdges_idx_{span}bp"] = bedpe_bpss.apply(lambda row: match_segedge_index(row, gr4_sub, span, "POS_2", "CHROM_2"), axis=1)

    columns_to_process = [
        "POS_1_SegEdges_idx_1bp", 
        "POS_1_SegEdges_idx_10bp", 
        "POS_2_SegEdges_idx_1bp", 
        "POS_2_SegEdges_idx_10bp"
    ]

    bedpe_bpss["SegEdge_amplicon"] = bedpe_bpss.apply(lambda row: extract_unique_amplicons(row, columns_to_process), axis=1)
    bedpe_bpss['AA_SegEdge'] = bedpe_bpss['SegEdge_amplicon'].apply(lambda x: 'True' if pd.notna(x) else 'False')
    bedpe_bpss['amplicon_idx'] = bedpe_bpss.apply(update_amplicon_idx, axis=1)
    bedpe_bpss['amplicon_idx'].replace('', pd.NA, inplace=True)
    bedpe_bpss['aa_barcode'] = bedpe_bpss['aa_barcode'].fillna(aliquot_barcode) 

    # ==== Step5: Process Amplicon Profiles
    amp_pair = amp[amp['aa_barcode'] == aliquot_barcode].copy()
    if amp_pair.shape[0] == 0:
        print(f"Skipping sample with barcode {aliquot_barcode} as no matching entries were found in amplicon wize.")
        create_empty_files(aliquot_barcode, output_filenames)
        return
    amp_pair=amp_pair[['aa_barcode', 'amplicon_index','amplicon_decomposition_class', 'ecDNA', 'BFB', 'patient_barcode',
            'cohort_barcode', 'project', 'primaryTumorLocation', 'source', 'tm_cov', 'tumor_purity', 'all_cutoff_passed']]        
    amp_pair = amp_pair.rename(columns={"amplicon_index": "amplicon_idx"})
    bedpe_bpss_aa = pd.merge(bedpe_bpss, amp_pair, on=['aa_barcode', 'amplicon_idx'], how='left')

    bedpe_bpss_aa_update = update_nan_values(bedpe_bpss_aa)
    
    bedpe_bpss_aa_update["amplicon_number"] = bedpe_bpss_aa_update["amplicon_idx"].str.extract(r'(\d+)$')
    bedpe_bpss_aa_update["amplicon_number"] = pd.to_numeric(bedpe_bpss_aa_update["amplicon_number"], errors="coerce")

    bedpe_bpss_aa_update = bedpe_bpss_aa_update.sort_values(
        by=["amplicon_number", "AA_SegEdge", "AA_discordants"],
        ascending=[True, False, False],  # True → 위, False → 아래, NaN → 맨 아래
        key=lambda col: col.map({True: 0, False: 1, np.nan: 2}) if col.name in ["AA_SegEdge", "AA_discordants"] else col
    )
    bedpe_bpss_aa_update = bedpe_bpss_aa_update.drop(columns=["amplicon_number"])
    
    bedpe_bpss_aa_update['INFO_MATEID'] = bedpe_bpss_aa_update['INFO_MATEID'].str.extract(r'^([^:]+)')

    # ==== Final updates and export
    print("Final updates and exporting")
    print("Output_file_1: raw file ==================")
    output_file = os.path.join(f"{aliquot_barcode}.svaba.somatic.sv.bedpe.aa.unfiltered.breakend.csv")
    bedpe_bpss_aa_update.to_csv(output_file, index=False)
    print(f"Process completed. Output saved to {output_file}")

    print("Output_file_2: raw file =================== ~ ['INFO_EVDNC'] != 'TSI_G' & ['QUAL'] != 0") 
    columns_to_drop = [
        "POS_1_SegEdges_idx_1bp", "POS_1_SegEdges_idx_10bp",
        "POS_2_SegEdges_idx_1bp", "POS_2_SegEdges_idx_10bp",
        "SegEdge_amplicon"
    ]    
    filtered_df = bedpe_bpss_aa_update[bedpe_bpss_aa_update['QUAL'] != 0]
    filtered_df = filtered_df.drop(columns=columns_to_drop)
    
    filtered_df2 = filtered_df[filtered_df['INFO_EVDNC'] != 'TSI_G']
    output_file2 = os.path.join(f"{aliquot_barcode}.svaba.somatic.sv.bedpe.aa.breakend.csv")    
    filtered_df2.to_csv(output_file2, index=False)
    print(f"Process completed. Output saved to {output_file2}")

    print("Output_file_3: ==============================only TSI_G")
    tsi_g = filtered_df[filtered_df["INFO_EVDNC"] == "TSI_G"].copy()
    tsi_l = filtered_df[filtered_df["INFO_EVDNC"] == "TSI_L"].copy()
    

    
    if tsi_g.shape[0] == 0:
        print(f"Skipping sample with barcode {aliquot_barcode} has no TGI_G evidence")
        output_file3 = os.path.join(f"{aliquot_barcode}.svaba.somatic.sv.bedpe.aa.TSI_G.breakend.csv")
        with open(output_file3, 'w') as f:
            pass  
        print(f"Created empty file: {output_file3}")
        return

    tsi_g[["overlaped_TSI_L", "num_overlaped_TSI_L"]] = tsi_g.apply(find_overlaps, axis=1, args=(tsi_l,))


    output_file3 = os.path.join(f"{aliquot_barcode}.svaba.somatic.sv.bedpe.aa.TSI_G.breakend.csv")
    tsi_g.to_csv(output_file3, index=False)
    print(f"Process completed. Output saved to {output_file3}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files for VCF to CSV and BEDPE conversions")
    parser.add_argument("aliquot_barcode", help="Path to the bedpe file")
    parser.add_argument("bedpe_path", help="Path to the svaba.sv.bedpe file")
    parser.add_argument("bps_path", help="Path to the bps.txt")
    parser.add_argument("bpss_output", help="Path to the matched_annotated_aa file")
    parser.add_argument("gr4_path", help="Path to the GR4 file")
    parser.add_argument("amp_path", help="Path to the Amplicon profiles file")
    args = parser.parse_args()

    process_files(
        aliquot_barcode=args.aliquot_barcode,
        bedpe_path=args.bedpe_path,
        bps_path=args.bps_path,
        bpss_output=args.bpss_output,
        gr4_path=args.gr4_path,
        amp_path=args.amp_path,
    )
