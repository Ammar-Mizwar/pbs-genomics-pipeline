import pandas as pd
import numpy as np
import os
from plotnine import *
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

class GenomicAnalysisSuite:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        os.chdir(base_dir)
        # Professional Serif Theme (matches your R setup)
        self.theme_style = theme_minimal(base_family='serif', base_size=12) + \
                           theme(panel_grid_major=element_blank(),
                                 panel_grid_minor=element_blank(),
                                 panel_border=element_rect(color='black', fill=None, size=0.5),
                                 plot_title=element_text(weight='bold', hjust=0.5),
                                 strip_background=element_rect(fill='grey90'))
        
        # Color palette for 22 autosomes
        self.chr_colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                           "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", 
                           "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", 
                           "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
                           "#FDBF6F", "#FF7F00"]

    # ==========================================
    # DATA PROCESSING & FDR CORRECTION
    # ==========================================
    # ==========================================
    # DATA PROCESSING & FDR CORRECTION
    # ==========================================
    def process_fdr_and_regions(self, file_daf, file_pbs):
        """Calculates Global FDR and filters for target regions."""
        print("[INFO] Processing Global FDR (Whole Assay Level)...")
        df_daf = pd.read_csv(file_daf)
        df_pbs = pd.read_csv(file_pbs).rename(columns={'CHROM': 'CHR', 'POS': 'BP'})
        
        df_merged = pd.merge(df_daf, df_pbs, on=['CHR', 'BP'], how='inner')
        
        # Whole Assay Level P-value and FDR calculation
        df_merged['p_BPS'] = norm.sf(df_merged['Z_PBS_MLY'])
        _, q_values, _, _ = multipletests(df_merged['p_BPS'], alpha=0.05, method='fdr_bh')
        df_merged['q_FDR'] = q_values

        # Filter for hg19 Target Regions
        chr_str = df_merged['CHR'].astype(str)
        cond_3q29 = (chr_str == '3') & (df_merged['BP'].between(195000000, 198022430))
        cond_7p11_2 = (chr_str == '7') & (df_merged['BP'].between(54800000, 58000000))
        cond_17p13_1 = (chr_str == '17') & (df_merged['BP'].between(6600000, 8500000))
        
        # FIX: Assign the Locus to the full dataframe first to avoid length mismatches
        df_merged['Locus'] = np.select(
            [cond_3q29, cond_7p11_2, cond_17p13_1], 
            ['3q29', '7p11.2', '17p13.1'], default='Other'
        )
        
        # Now slice out only our targeted regions
        df_regions = df_merged[df_merged['Locus'] != 'Other'].copy()
        
        df_regions.to_csv("Target_Regions_All_Variants_Global_FDR.csv", index=False)
        return df_merged, df_regions

    # ==========================================
    # POPULATION STRUCTURE PLOTS
    # ==========================================
    def plot_pca(self, pca_file="pop_structure/plink_pca.eigenvec"):
        """Plots PCA Results (PC1 vs PC2)."""
        print("[INFO] Plotting PCA...")
        if not os.path.exists(pca_file):
            print(f"[WARN] {pca_file} not found. Skipping PCA plot.")
            return

        pca_df = pd.read_csv(pca_file, sep=r'\s+', header=None, 
                             names=['FID', 'IID'] + [f'PC{i}' for i in range(1, 11)])
        
        p = (ggplot(pca_df, aes(x='PC1', y='PC2'))
             + geom_point(alpha=0.6, color="#377EB8")
             + self.theme_style
             + labs(title="Principal Component Analysis (PCA)", x="PC1", y="PC2"))
        
        p.save("PCA_PC1_PC2.png", width=8, height=6, dpi=300)

    def plot_admixture(self, pruned_prefix="genome_wide_merged_pruned"):
        """Plots Admixture results for the best K found during the pipeline."""
        print("[INFO] Plotting Admixture results...")
        summary_file = "pop_structure/best_k_summary.txt"
        
        # Safety check: skip if Admixture wasn't run
        if not os.path.exists(summary_file):
            print(f"[WARN] {summary_file} not found. Skipping Admixture plot.")
            return

        # Automatically find the best K from summary file generated in pipeline
        with open(summary_file, "r") as f:
            best_k = [line for line in f if "Optimal K:" in line][0].split()[-1]
        
        q_file = f"pop_structure/{pruned_prefix}.{best_k}.Q"
        if not os.path.exists(q_file):
            print(f"[WARN] {q_file} not found. Skipping Admixture plot.")
            return

        df_q = pd.read_csv(q_file, sep=r'\s+', header=None)
        df_q['ID'] = range(len(df_q))
        df_long = df_q.melt(id_vars='ID', var_name='Ancestry', value_name='Fraction')
        
        p = (ggplot(df_long, aes(x='factor(ID)', y='Fraction', fill='factor(Ancestry)'))
             + geom_bar(stat='identity', width=1)
             + scale_fill_brewer(type='qual', palette='Set1')
             + self.theme_style
             + theme(axis_text_x=element_blank(), legend_title=element_blank())
             + labs(title=f"Admixture Ancestry Proportions (K={best_k})", x="Individuals", y="Ancestry Fraction"))
        
        p.save(f"Admixture_K{best_k}.png", width=12, height=4, dpi=300)

    # ==========================================
    # SELECTION SCAN PLOTS
    # ==========================================
    def plot_manhattan_fdr(self, df_merged):
        """Genome-wide Manhattan plot using -log10(q_FDR)."""
        print("[INFO] Plotting Manhattan (FDR Threshold)...")
        # Prepare cumulative positions for plotting
        df = df_merged[df_merged['CHR'] != '23'].copy()
        df['CHR'] = pd.to_numeric(df['CHR'])
        df = df.sort_values(['CHR', 'BP'])
        
        df['-logQ'] = -np.log10(df['q_FDR'].replace(0, 1e-300))
        
        # Calculate cumulative BP for the x-axis
        df['pos_cum'] = 0
        s = 0
        centers = []
        for ch, group in df.groupby('CHR'):
            df.loc[group.index, 'pos_cum'] = group['BP'] + s
            centers.append(s + group['BP'].max()/2)
            s += group['BP'].max()

        p = (ggplot(df, aes(x='pos_cum', y='-logQ', color='factor(CHR)'))
             + geom_point(alpha=0.8, size=1)
             + scale_color_manual(values=self.chr_colors)
             + geom_hline(yintercept=-np.log10(0.05), color="red", linetype="dashed")
             + self.theme_style
             + theme(legend_position="none", 
                     axis_text_x=element_blank(), 
                     axis_ticks_major_x=element_blank())
             + labs(title="Genome-Wide Selection Scan (FDR corrected)", 
                    x="Chromosome", y="-log10(q-value)"))
        
        p.save("Manhattan_FDR.png", width=18, height=6, dpi=600)

# === EXECUTION ===
if __name__ == "__main__":
    suite = GenomicAnalysisSuite(os.getcwd())
    
    # 1. Process Data & FDR
    merged_data, regional_data = suite.process_fdr_and_regions(
        "Master_Genomic_Scan_FST_DAF.csv", "Master_PBS_Results.csv"
    )
    
    # 2. Generate Plots
    suite.plot_pca()
    suite.plot_admixture() 
    suite.plot_manhattan_fdr(merged_data)
