import os
import subprocess
import pandas as pd
import numpy as np
import sys
import shutil
import glob
import re

# ==========================================
# CONFIGURATION
# ==========================================
WD = os.getcwd()
TEMP_DIR = os.path.join(WD, "temp_work")
FREQ_DIR = os.path.join(WD, "freq_data")
POP_DIR = os.path.join(WD, "pop_structure")
OS_ENV = os.environ.copy()

# File Paths
LOCAL_PLINK_PREFIX = "Plink_Output.cleaned.top.snpqc"

# URLs
KGP_PANEL_URL = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
KGP_VCF_BASE = "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"

# Superpopulations
SUPERPOPS = ['AFR', 'EAS']

# ==========================================
# HELPER FUNCTIONS
# ==========================================
def run_cmd(cmd):
    """Runs a shell command with error checking."""
    print(f"[CMD] {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, env=OS_ENV, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}")
        raise e

def count_snps(prefix):
    if os.path.exists(f"{prefix}.bim"):
        with open(f"{prefix}.bim") as f:
            return sum(1 for _ in f)
    return 0

def count_individuals(prefix):
    if os.path.exists(f"{prefix}.fam"):
        with open(f"{prefix}.fam") as f:
            return sum(1 for _ in f)
    return 0

def get_chrom_num(filename):
    match = re.search(r'chr(\d+)', filename)
    if match: return int(match.group(1))
    return 999

def download_file(url, output_path, min_size_bytes=1024):
    if os.path.exists(output_path):
        if os.path.getsize(output_path) < min_size_bytes:
            print(f"[WARN] File {output_path} is too small. Re-downloading...")
            os.remove(output_path)
        else:
            print(f"[INFO] File {output_path} exists. Skipping.")
            return
    print(f"[DOWNLOAD] {url} -> {output_path}")
    try:
        run_cmd(f"curl -f -L -o \"{output_path}\" \"{url}\"")
    except Exception:
        print(f"[ERROR] Failed to download {url}")
        sys.exit(1)

def check_tools():
    tools = ['plink', 'bcftools', 'vcftools', 'curl']
    for t in tools:
        if shutil.which(t) is None:
            print(f"[ERROR] Tool '{t}' not found. Ensure it is in PATH.")
            sys.exit(1)

def get_superpop(s_id, mapping_dict):
    if s_id in mapping_dict: return mapping_dict[s_id]
    parts = s_id.split('_')
    if len(parts) > 1 and parts[-1] in mapping_dict: return mapping_dict[parts[-1]]
    if len(parts) > 1 and parts[0] in mapping_dict: return mapping_dict[parts[0]]
    return None

def create_polarization_rules(aa_file, ref_alt_file, output_rule_file):
    aa_dict = {}
    with open(aa_file, 'r') as f:
        for line in f:
            p = line.strip().split()
            if len(p) >= 2:
                aa = p[1].split('|')[0].upper()
                if aa in ['A', 'C', 'G', 'T']: aa_dict[p[0]] = aa

    with open(ref_alt_file, 'r') as f_in, open(output_rule_file, 'w') as f_out:
        for line in f_in:
            p = line.strip().split()
            if len(p) < 3: continue
            snp, ref, alt = p[0], p[1], p[2]
            if snp not in aa_dict: continue
            anc = aa_dict[snp]
            if anc == ref:
                f_out.write(f"{snp}\tKEEP\n")
            elif anc == alt:
                f_out.write(f"{snp}\tFLIP\n")

# ==========================================
# POPULATION STRUCTURE (PCA ONLY)
# ==========================================
def run_pca_only():
    print("\n" + "=" * 50)
    print("RUNNING GLOBAL PCA (AUTOSOMES 1-22)")
    print("=" * 50)
    
    os.chdir(POP_DIR)
    
    if not os.path.exists("qc_chr1.bed"):
        print("[ERROR] qc_chr1.bed not found in pop_structure. PCA cannot proceed.")
        os.chdir(WD)
        return

    merge_list = "merge_list.txt"
    valid_chroms = []
    with open(merge_list, 'w') as f:
        for i in range(2, 23):
            if os.path.exists(f"qc_chr{i}.bed"):
                f.write(f"qc_chr{i}.bed qc_chr{i}.bim qc_chr{i}.fam\n")
                valid_chroms.append(i)
    
    global_prefix = "autosomes_merged"
    
    if not os.path.exists(f"{global_prefix}.bed"):
        print(f"--- Merging {len(valid_chroms)+1} Autosomes for Whole-Genome PCA ---")
        try:
            run_cmd(f"plink --bfile qc_chr1 --merge-list {merge_list} --make-bed --out {global_prefix}")
        except subprocess.CalledProcessError:
            print("[WARN] Merge failed. Applying automated fix (.missnp exclusion)...")
            missnp_file = f"{global_prefix}-merge.missnp"
            
            if os.path.exists(missnp_file):
                run_cmd(f"plink --bfile qc_chr1 --exclude {missnp_file} --make-bed --out qc_chr1_clean")
                with open(merge_list, 'w') as f:
                    for i in valid_chroms:
                        run_cmd(f"plink --bfile qc_chr{i} --exclude {missnp_file} --make-bed --out qc_chr{i}_clean")
                        f.write(f"qc_chr{i}_clean.bed qc_chr{i}_clean.bim qc_chr{i}_clean.fam\n")
                run_cmd(f"plink --bfile qc_chr1_clean --merge-list {merge_list} --make-bed --out {global_prefix}")
            else:
                print("[ERROR] Merge failed but no .missnp file was generated. PCA skipped.")
                os.chdir(WD)
                return

    pruned_prefix = f"{global_prefix}_pruned"
    if os.path.exists(f"{global_prefix}.bed") and not os.path.exists(f"{pruned_prefix}.bed"):
        print("--- LD Pruning for PCA ---")
        run_cmd(f"plink --bfile {global_prefix} --indep-pairwise 50 10 0.1 --out pop_prune")
        run_cmd(f"plink --bfile {global_prefix} --extract pop_prune.prune.in --make-bed --out {pruned_prefix}")

    if os.path.exists(f"{pruned_prefix}.bed") and not os.path.exists("plink_pca.eigenvec"):
        print("--- Running PCA (10 Components) ---")
        run_cmd(f"plink --bfile {pruned_prefix} --pca 10 --out plink_pca")

    os.chdir(WD)

# ==========================================
# NORMALIZATION & CSV GENERATION
# ==========================================
def normalize_and_generate_csv():
    print("\n" + "=" * 50)
    print("CALCULATING Z-SCORES & GENERATING MASTER CSV")
    print("=" * 50)

    freq_files = sorted(glob.glob(os.path.join(FREQ_DIR, "freq_chr*.csv")), key=get_chrom_num)
    if not freq_files:
        print("[ERROR] No frequency files found. Run pipeline first.")
        return

    meta_df = pd.concat((pd.read_csv(f) for f in freq_files), ignore_index=True)
    meta_df['CHR'] = meta_df['CHR'].astype(str)
    meta_df['SNP'] = meta_df['SNP'].astype(str)

    pops = ['MLY'] + SUPERPOPS
    for pop in pops:
        col_raw = f"Freq_{pop}"
        col_daf = f"Derived_Freq_{pop}"
        if col_raw in meta_df.columns:
            meta_df[col_daf] = np.where(meta_df['Rule'] == 'KEEP', meta_df[col_raw], 1.0 - meta_df[col_raw])

    for sp in SUPERPOPS:
        fst_files = sorted(glob.glob(f"results_fst_chr*_{sp}.weir.fst"), key=get_chrom_num)
        fst_dfs = []
        for f in fst_files:
            d = pd.read_csv(f, sep='\t')
            d.rename(columns={'WEIR_AND_COCKERHAM_FST': f'FST_{sp}', 'POS': 'BP', 'CHROM': 'CHR'}, inplace=True)
            fst_dfs.append(d[['CHR', 'BP', f'FST_{sp}']])

        if fst_dfs:
            fst_all = pd.concat(fst_dfs)
            fst_all['CHR'] = fst_all['CHR'].astype(str)
            meta_df = pd.merge(meta_df, fst_all, on=['CHR', 'BP'], how='left')
            meta_df[f'Z_FST_{sp}'] = (meta_df[f'FST_{sp}'] - meta_df[f'FST_{sp}'].mean()) / meta_df[f'FST_{sp}'].std()

        col_daf_mly = "Derived_Freq_MLY"
        col_daf_sp = f"Derived_Freq_{sp}"
        col_delta = f"DeltaDAF_{sp}"

        if col_daf_sp in meta_df.columns:
            meta_df[col_delta] = meta_df[col_daf_mly] - meta_df[col_daf_sp]
            meta_df[f'Z_DeltaDAF_{sp}'] = (meta_df[col_delta] - meta_df[col_delta].mean()) / meta_df[col_delta].std()

    final_cols = ['CHR', 'BP', 'SNP', 'AA', 'Derived_Freq_MLY']
    for sp in SUPERPOPS:
        final_cols.extend([f'FST_{sp}', f'Z_FST_{sp}', f'DeltaDAF_{sp}', f'Z_DeltaDAF_{sp}'])

    exist_cols = [c for c in final_cols if c in meta_df.columns]
    out_name = "Master_Genomic_Scan_FST_DAF.csv"
    print(f"--- Writing {out_name} ---")
    meta_df[exist_cols].to_csv(out_name, index=False, float_format='%.5f')

# ==========================================
# MAIN PIPELINE
# ==========================================
def main():
    if not os.path.exists(WD):
        print(f"[ERROR] Working directory {WD} does not exist.")
        sys.exit(1)

    os.chdir(WD)
    os.makedirs(TEMP_DIR, exist_ok=True)
    os.makedirs(FREQ_DIR, exist_ok=True)
    os.makedirs(POP_DIR, exist_ok=True)
    check_tools()

    print("\n=== STEP 0: PRE-PROCESSING ===")
    download_file(KGP_PANEL_URL, "1kgp_panel.txt")
    panel = pd.read_csv("1kgp_panel.txt", sep='\t')
    kgp_map = dict(zip(panel['sample'], panel['super_pop']))

    qc1_out = f"{LOCAL_PLINK_PREFIX}.strictly_cleaned.nonpruned"
    prune = f"{LOCAL_PLINK_PREFIX}.pruning_files"
    qc_out = f"{LOCAL_PLINK_PREFIX}.strictly_cleaned.pruned.nonrelated"
    
    if not os.path.exists(f"{qc_out}.bed"):
        run_cmd(f"plink --bfile \"{LOCAL_PLINK_PREFIX}\" --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 --biallelic-only strict --snps-only --keep-allele-order --make-bed --out \"{qc1_out}\"")
        #run_cmd(f"plink --bfile \"{qc1_out}\" --indep-pairwise 50 5 0.2 --out \"{prune}\"")
        run_cmd(f"plink --bfile \"{qc1_out}\" --rel-cutoff 0.125 --make-bed --out \"{qc_out}\"")

    # --- CHROMOSOME LOOP (AUTOSOMES ONLY) ---
    for chrom in range(1, 23):
        print(f"\n\n{'=' * 40}\nPROCESSING CHROMOSOME {chrom}\n{'=' * 40}")
        vcf_filename = f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        merged_vcf = os.path.join(TEMP_DIR, f"merged_chr{chrom}.vcf")

        try:
            # 1. Download & Prep Subsets
            download_file(f"{KGP_VCF_BASE}{vcf_filename}", vcf_filename, 1000000)
            
            run_cmd(f"plink --bfile \"{qc_out}\" --chr {chrom} --keep-allele-order --make-bed --out \"{TEMP_DIR}/mly_chr{chrom}\"")
            run_cmd(f"plink --vcf \"{vcf_filename}\" --biallelic-only strict --snps-only --double-id --keep-allele-order --make-bed --out \"{TEMP_DIR}/kgp_chr{chrom}\"")

            # 2. Extract Common SNPs to avoid mismatch errors early on
            with open(os.path.join(TEMP_DIR, "common_snps.txt"), "w") as f:
                mly_bim = pd.read_csv(f"{TEMP_DIR}/mly_chr{chrom}.bim", sep=r'\s+', header=None)[1]
                kgp_bim = pd.read_csv(f"{TEMP_DIR}/kgp_chr{chrom}.bim", sep=r'\s+', header=None)[1]
                common = set(mly_bim).intersection(set(kgp_bim))
                for snp in common: f.write(f"{snp}\n")

            run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_chr{chrom}\" --extract \"{TEMP_DIR}/common_snps.txt\" --make-bed --out \"{TEMP_DIR}/mly_subset\"")
            run_cmd(f"plink --bfile \"{TEMP_DIR}/kgp_chr{chrom}\" --extract \"{TEMP_DIR}/common_snps.txt\" --make-bed --out \"{TEMP_DIR}/kgp_subset\"")

            # 3. Robust 3-Strike Merging Strategy
            try:
                run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_subset\" --bmerge \"{TEMP_DIR}/kgp_subset\" --keep-allele-order --make-bed --out \"{TEMP_DIR}/merged_chr{chrom}_raw\"")
            except subprocess.CalledProcessError:
                print("[WARN] First merge failed (Likely allele flip). Attempting to flip strands...")
                missnp = f"{TEMP_DIR}/merged_chr{chrom}_raw-merge.missnp"
                
                if os.path.exists(missnp):
                    run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_subset\" --flip \"{missnp}\" --make-bed --out \"{TEMP_DIR}/mly_flipped\"")
                    try:
                        run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_flipped\" --bmerge \"{TEMP_DIR}/kgp_subset\" --keep-allele-order --make-bed --out \"{TEMP_DIR}/merged_chr{chrom}_raw\"")
                    except subprocess.CalledProcessError:
                        print("[WARN] Flip failed. Excluding stubborn multi-allelics...")
                        missnp2 = f"{TEMP_DIR}/merged_chr{chrom}_raw-merge.missnp"
                        run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_flipped\" --exclude \"{missnp2}\" --make-bed --out \"{TEMP_DIR}/mly_cleaned\"")
                        run_cmd(f"plink --bfile \"{TEMP_DIR}/kgp_subset\" --exclude \"{missnp2}\" --make-bed --out \"{TEMP_DIR}/kgp_cleaned\"")
                        run_cmd(f"plink --bfile \"{TEMP_DIR}/mly_cleaned\" --bmerge \"{TEMP_DIR}/kgp_cleaned\" --keep-allele-order --make-bed --out \"{TEMP_DIR}/merged_chr{chrom}_raw\"")

            # 4. Filter Merged Result & Convert to VCF
            run_cmd(f"plink --bfile \"{TEMP_DIR}/merged_chr{chrom}_raw\" --maf 0.01 --hwe 1e-6 --geno 0.1 --snps-only --keep-allele-order --make-bed --out \"{TEMP_DIR}/qc_chr{chrom}\"")
            run_cmd(f"plink --bfile \"{TEMP_DIR}/qc_chr{chrom}\" --recode vcf --out \"{TEMP_DIR}/merged_chr{chrom}\"")

            # 5. Populations Setup
            fam_file = os.path.join(TEMP_DIR, f"qc_chr{chrom}.fam")
            fam_df = pd.read_csv(fam_file, sep=r'\s+', header=None, usecols=[0, 1], names=['FID', 'IID'])
            fam_df['FID'] = fam_df['FID'].astype(str)
            fam_df['IID'] = fam_df['IID'].astype(str)

            mly_fam = fam_df[fam_df['IID'].apply(lambda x: get_superpop(x, kgp_map) is None)]
            if len(mly_fam) == 0: raise ValueError("No target/local samples found after merging.")

            mly_fam.to_csv(os.path.join(TEMP_DIR, "pop_MLY.txt"), sep='\t', header=False, index=False)

            sp_counts = {}
            for sp in SUPERPOPS:
                sp_fam = fam_df[fam_df['IID'].apply(lambda x: get_superpop(x, kgp_map) == sp)]
                sp_fam.to_csv(os.path.join(TEMP_DIR, f"pop_{sp}.txt"), sep='\t', header=False, index=False)
                sp_counts[sp] = len(sp_fam)

            sample_file = os.path.join(TEMP_DIR, "samples.txt")
            with open(sample_file, "w") as f_out:
                subprocess.run(["bcftools", "query", "-l", merged_vcf], stdout=f_out, check=True)

            vcf_samples = [l.strip() for l in open(sample_file)]
            with open(os.path.join(TEMP_DIR, "pop_MLY_bcf.txt"), "w") as f:
                for s in [s for s in vcf_samples if get_superpop(s, kgp_map) is None]: f.write(f"{s}\n")

            for sp in SUPERPOPS:
                with open(os.path.join(TEMP_DIR, f"pop_{sp}_bcf.txt"), "w") as f:
                    for s in [s for s in vcf_samples if get_superpop(s, kgp_map) == sp]: f.write(f"{s}\n")

            # 6. Polarization & Frequencies
            f_aa = os.path.join(TEMP_DIR, f"aa_chr{chrom}.txt")
            f_ra = os.path.join(TEMP_DIR, f"ref_alt_chr{chrom}.txt")
            f_rules = os.path.join(TEMP_DIR, f"rules_chr{chrom}.txt")
            f_pos = os.path.join(TEMP_DIR, f"pos_chr{chrom}.txt")

            run_cmd(f"bcftools query -f '%ID\t%INFO/AA\n' \"{vcf_filename}\" > \"{f_aa}\"")
            run_cmd(f"bcftools query -f '%ID\t%REF\t%ALT\n' \"{merged_vcf}\" > \"{f_ra}\"")
            create_polarization_rules(f_aa, f_ra, f_rules)
            run_cmd(f"bcftools query -f '%ID\t%CHROM\t%POS\n' \"{merged_vcf}\" > \"{f_pos}\"")

            df_rules = pd.read_csv(f_rules, sep='\t', header=None, names=['SNP', 'Rule'])
            df_pos = pd.read_csv(f_ra, sep='\t', header=None, names=['SNP', 'Ref', 'Alt'])
            df_bp = pd.read_csv(f_pos, sep='\t', header=None, names=['SNP', 'CHR', 'BP'])
            df_bp['CHR'], df_bp['SNP'] = df_bp['CHR'].astype(str), df_bp['SNP'].astype(str)
            df_meta = pd.merge(df_bp, df_rules, on='SNP')

            def process_plink_freq(freq_file, df_ra, col):
                if not os.path.exists(freq_file): return None
                d = pd.read_csv(freq_file, sep=r'\s+')
                t = pd.merge(d, df_ra, on='SNP')
                t[col] = np.where(t['A1'] == t['Alt'], t['MAF'], 1.0 - t['MAF'])
                return t[['SNP', col]]

            run_cmd(f"plink --bfile \"{TEMP_DIR}/qc_chr{chrom}\" --allow-no-sex --keep \"{TEMP_DIR}/pop_MLY.txt\" --freq --out \"{TEMP_DIR}/freq_mly_chr{chrom}\"")
            df_mly = process_plink_freq(f"{TEMP_DIR}/freq_mly_chr{chrom}.frq", df_pos, "Freq_MLY")
            if df_mly is not None: df_meta = pd.merge(df_meta, df_mly, on='SNP')

            for sp in SUPERPOPS:
                if sp_counts[sp] > 0:
                    run_cmd(f"plink --bfile \"{TEMP_DIR}/qc_chr{chrom}\" --allow-no-sex --keep \"{TEMP_DIR}/pop_{sp}.txt\" --freq --out \"{TEMP_DIR}/freq_{sp}_chr{chrom}\"")
                    df_sp = process_plink_freq(f"{TEMP_DIR}/freq_{sp}_chr{chrom}.frq", df_pos, f"Freq_{sp}")
                    if df_sp is not None: df_meta = pd.merge(df_meta, df_sp, on='SNP', how='left')
                else:
                    df_meta[f"Freq_{sp}"] = np.nan

            df_meta.to_csv(os.path.join(FREQ_DIR, f"freq_chr{chrom}.csv"), index=False)

            # 7. FST Calculation
            for sp in SUPERPOPS:
                if sp_counts[sp] > 0 and not os.path.exists(f"results_fst_chr{chrom}_{sp}.weir.fst"):
                    run_cmd(f"vcftools --vcf \"{merged_vcf}\" --weir-fst-pop \"{TEMP_DIR}/pop_MLY_bcf.txt\" --weir-fst-pop \"{TEMP_DIR}/pop_{sp}_bcf.txt\" --out results_fst_chr{chrom}_{sp}")

            # Stash autosomes for Global PCA before wiping TEMP_DIR
            for ext in ['.bed', '.bim', '.fam']:
                src = os.path.join(TEMP_DIR, f"qc_chr{chrom}{ext}")
                dst = os.path.join(POP_DIR, f"qc_chr{chrom}{ext}")
                if os.path.exists(src): shutil.copy(src, dst)

            print(f"--- Cleaning Temp Files for Chr {chrom} ---")
            shutil.rmtree(TEMP_DIR)
            os.makedirs(TEMP_DIR, exist_ok=True)

        except Exception as e:
            print(f"\n[CRITICAL ERROR] Chromosome {chrom} failed: {e}")
            import traceback
            traceback.print_exc()
            print("Skipping to next, but check logs.\n")
            continue

    run_pca_only()
    normalize_and_generate_csv()
    print("\n=== PIPELINE FINISHED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
