import pandas as pd
import os
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# === CONFIGURATION ===
base_dir = os.getcwd()
os.chdir(base_dir)

file_daf = "Master_Genomic_Scan_FST_DAF.csv"
file_pbs = "Master_PBS_Results.csv"

def main():
    print("[INFO] Loading master CSV files...")
    df_daf = pd.read_csv(file_daf)
    df_pbs = pd.read_csv(file_pbs)

    print("[INFO] Merging datasets...")
    df_pbs.rename(columns={'CHROM': 'CHR', 'POS': 'BP'}, inplace=True)
    df_merged = pd.merge(df_daf, df_pbs, on=['CHR', 'BP'], how='inner')

    # --- WHOLE ASSAY LEVEL CALCULATION ---
    print(f"[INFO] Calculating P-values and FDR for all {len(df_merged)} SNPs...")
    
    # 1. Get P-values for every single merged SNP
    df_merged['p_BPS'] = norm.sf(df_merged['Z_PBS_MLY'])
    
    # 2. Apply Benjamini-Hochberg to the ENTIRE genome-wide set
    # This accounts for the true number of tests performed (N)
    _, q_values, _, _ = multipletests(df_merged['p_BPS'], alpha=0.05, method='fdr_bh')
    df_merged['q_FDR'] = q_values

    print("[INFO] Filtering for target regions (hg19) AFTER FDR correction...")
    chr_col = df_merged['CHR'].astype(str)

    cond_3q29 = (chr_col == '3') & (df_merged['BP'] >= 195000000) & (df_merged['BP'] <= 198022430)
    cond_7p11_2 = (chr_col == '7') & (df_merged['BP'] >= 54800000) & (df_merged['BP'] <= 58000000)
    cond_17p13_1 = (chr_col == '17') & (df_merged['BP'] >= 6600000) & (df_merged['BP'] <= 8500000)

    # Apply the regional filter to the already-corrected data
    df_regions = df_merged[cond_3q29 | cond_7p11_2 | cond_17p13_1].copy()

    # --- TASK 1: Save all variants in target regions (with global FDR) ---
    output1 = "Target_Regions_All_Variants_Global_FDR.csv"
    df_regions.to_csv(output1, index=False)
    print(f"[SUCCESS] Saved {len(df_regions)} variants to: {output1}")

    # --- TASK 2: Save SIGNIFICANT variants (Global FDR q < 0.05) ---
    print("[INFO] Filtering for Global FDR q < 0.05...")
    df_significant = df_regions[df_regions['q_FDR'] < 0.05].copy()

    # Select and rename columns
    cols_to_keep = ['CHR', 'SNP', 'BP', 'PBS_MLY', 'Z_PBS_MLY', 'p_BPS', 'q_FDR', 'DeltaDAF_EAS', 'Z_DeltaDAF_EAS']
    df_significant = df_significant[cols_to_keep]
    df_significant.rename(columns={
        'BP': 'POS',
        'PBS_MLY': 'BPS',
        'Z_PBS_MLY': 'ZBPS',
        'DeltaDAF_EAS': 'DAF_EAS',
        'Z_DeltaDAF_EAS': 'Z_DAF_EAS'
    }, inplace=True)

    output2 = "Target_Regions_Global_FDR_Significant.csv"
    df_significant.to_csv(output2, index=False)
    print(f"[SUCCESS] Saved {len(df_significant)} variants with Global FDR < 0.05 to: {output2}")

if __name__ == "__main__":
    main()
