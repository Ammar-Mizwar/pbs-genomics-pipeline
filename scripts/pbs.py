import pandas as pd
import numpy as np
import os
import glob

# === CONFIGURATION ===
base_dir = os.getcwd()
os.chdir(base_dir)

target = "MLY"
ref = "EAS"
outgroup = "AFR"


# Function to calculate Wright's FST directly from allele frequencies
def calc_wright_fst(p1, p2):
    h1 = 2 * p1 * (1 - p1)
    h2 = 2 * p2 * (1 - p2)
    hs = (h1 + h2) / 2
    p_mean = (p1 + p2) / 2
    ht = 2 * p_mean * (1 - p_mean)

    # Calculate FST, avoiding division by zero if Ht is 0
    fst = np.where(ht > 0, (ht - hs) / ht, 0)
    return np.clip(fst, 0, 0.9999)  # Clip any statistical noise


# Detect processed FST files
chr_files = glob.glob(f"results_fst_chr*_{outgroup}.weir.fst")
chromosomes = sorted([int(f.split('_chr')[1].split('_')[0]) for f in chr_files])

print(f"Detected processed FST files for {len(chromosomes)} chromosomes.")

all_pbs_data = []

for chrom in chromosomes:
    print(f"\n{'=' * 40}\nPROCESSING CHROMOSOME {chrom}\n{'=' * 40}")

    fst_mly_eas_file = f"results_fst_chr{chrom}_{ref}.weir.fst"
    fst_mly_afr_file = f"results_fst_chr{chrom}_{outgroup}.weir.fst"
    freq_file = f"freq_data/freq_chr{chrom}.csv"

    if not os.path.exists(freq_file):
        print(f"[ERROR] Missing frequency file {freq_file}. Skipping.")
        continue

    print(f"[INFO] Calculating mathematically derived FST for {ref} vs {outgroup}...")

    # 1. Load Frequencies (Directly from your compiled CSVs)
    df_freq = pd.read_csv(freq_file)

    # Rename columns to match the standard VCFtools headers (CHROM, POS)
    df_freq.rename(columns={'CHR': 'CHROM', 'BP': 'POS', 'Freq_EAS': 'MAF_EAS', 'Freq_AFR': 'MAF_AFR'}, inplace=True)

    # Calculate FST directly from the frequencies
    df_freq['FST_EAS_AFR'] = calc_wright_fst(df_freq['MAF_EAS'], df_freq['MAF_AFR'])

    # 2. Load VCFtools FSTs
    try:
        df_mly_eas = pd.read_csv(fst_mly_eas_file, sep='\t')
        df_mly_afr = pd.read_csv(fst_mly_afr_file, sep='\t')
    except Exception as e:
        print(f"[ERROR] Could not read FST files for Chr {chrom}. Skipping.")
        continue

    df_mly_eas.rename(columns={'WEIR_AND_COCKERHAM_FST': 'FST_MLY_EAS'}, inplace=True)
    df_mly_afr.rename(columns={'WEIR_AND_COCKERHAM_FST': 'FST_MLY_AFR'}, inplace=True)

    # 3. Merge all three FST values on exactly shared SNPs using CHROM and POS
    df_merged = df_mly_eas.merge(df_mly_afr[['CHROM', 'POS', 'FST_MLY_AFR']], on=['CHROM', 'POS'], how='inner')
    df_merged = df_merged.merge(df_freq[['CHROM', 'POS', 'FST_EAS_AFR']], on=['CHROM', 'POS'], how='inner')
    df_merged.dropna(subset=['FST_MLY_EAS', 'FST_MLY_AFR', 'FST_EAS_AFR'], inplace=True)

    # 4. Calculate T-values and PBS
    print("[INFO] Calculating PBS scores...")
    for col in ['FST_MLY_EAS', 'FST_MLY_AFR', 'FST_EAS_AFR']:
        df_merged[col] = np.clip(df_merged[col], 0, 0.9999)

    t_mly_eas = -np.log(1 - df_merged['FST_MLY_EAS'])
    t_mly_afr = -np.log(1 - df_merged['FST_MLY_AFR'])
    t_eas_afr = -np.log(1 - df_merged['FST_EAS_AFR'])

    df_merged['PBS_MLY'] = (t_mly_eas + t_mly_afr - t_eas_afr) / 2
    df_merged['PBS_MLY'] = np.maximum(df_merged['PBS_MLY'], 0)  # Floor negative PBS to 0

    all_pbs_data.append(df_merged)

# === FINAL COMPILATION ===
if all_pbs_data:
    print("\n" + "=" * 40)
    print("COMPILING MASTER DATASET")
    print("=" * 40)

    master_df = pd.concat(all_pbs_data, ignore_index=True)

    # Standardize to Z-PBS
    print("[INFO] Calculating Standardized Z-PBS scores...")
    pbs_mean = master_df['PBS_MLY'].mean()
    pbs_std = master_df['PBS_MLY'].std()
    master_df['Z_PBS_MLY'] = (master_df['PBS_MLY'] - pbs_mean) / pbs_std

    # Clean up and export
    cols = ['CHROM', 'POS', 'FST_MLY_EAS', 'FST_MLY_AFR', 'FST_EAS_AFR', 'PBS_MLY', 'Z_PBS_MLY']
    master_df = master_df[cols]

    output_file = "Master_PBS_Results.csv"
    master_df.to_csv(output_file, index=False)

    print(f"[SUCCESS] Analyzed {len(master_df)} SNPs.")
    print(f"[SUCCESS] Master file saved to: {output_file}")
    print("=" * 40)
else:
    print("\n[ERROR] No data was processed. Double-check your freq_data folder.")