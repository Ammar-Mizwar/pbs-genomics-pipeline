import pandas as pd
import requests
import time
import sys

# === CONFIGURATION ===
bim_file = "Plink_Output.cleaned.top.snpqc.bim"
variants_file = "VEP_Annotations_Results.csv"
output_file = "Missense_PolyPhen_SIFT_Results.csv"

# CRITICAL: Using the GRCh37 (hg19) specific server
server = "https://grch37.rest.ensembl.org"

print("[INFO] Loading PLINK .bim file...")
# 1. Load the .bim file into a dictionary for blazing fast lookups
try:
    bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                         names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])

    # Create a dictionary mapping (CHR, BP) -> A1 (Alternate Allele)
    bim_dict = {(str(row['CHR']), str(row['BP'])): row['A1'] for _, row in bim_df.iterrows()}
    print(f"[INFO] Loaded {len(bim_dict)} markers from the .bim file.")
except FileNotFoundError:
    print(f"[ERROR] Could not find {bim_file}. Check your directory.")
    sys.exit(1)

print("[INFO] Loading targeted variants...")
# 2. Load your significant variants
try:
    df_variants = pd.read_csv(variants_file)
except FileNotFoundError:
    print(f"[ERROR] Could not find {variants_file}.")
    sys.exit(1)

# Filter for only missense variants
df_missense = df_variants[df_variants['Consequence'].str.contains('missense', case=False, na=False)].copy()
print(f"[INFO] Found {len(df_missense)} missense variants to analyze.")

# ---> THE FIX: Create an empty file for Nextflow before exiting <---
if len(df_missense) == 0:
    print(f"[INFO] No missense variants found. Creating empty {output_file} to satisfy Nextflow and exiting.")
    df_missense.to_csv(output_file, index=False)
    sys.exit(0)

# 3. Prepare to query the API
headers = {"Content-Type": "application/json"}
results_list = []

print("[INFO] Querying Ensembl for PolyPhen & SIFT scores...")

for index, row in df_missense.iterrows():
    # Handle column names flexibly
    chrom = str(row.get('CHR', row.get('CHROM', '')))
    pos = str(row.get('BP', row.get('POS', '')))
    rsid = row.get('rsID', f"{chrom}:{pos}")
    gene = row.get('Gene_Symbol', 'Unknown')

    # Lookup the Alternate Allele (A1) from the .bim file using CHR and POS
    alt_allele = bim_dict.get((chrom, pos))

    if not alt_allele:
        print(f"[WARNING] Could not find {chrom}:{pos} in the .bim file. Skipping.")
        continue

    # Endpoint for a specific region and allele, requesting SIFT and PolyPhen data
    endpoint = f"/vep/human/region/{chrom}:{pos}-{pos}/{alt_allele}?sift=b&polyphen=b"

    try:
        response = requests.get(server + endpoint, headers=headers)

        if response.ok:
            data = response.json()

            sift_score, sift_pred = "N/A", "N/A"
            poly_score, poly_pred = "N/A", "N/A"

            # Extract the scores from the transcript consequences
            if len(data) > 0 and 'transcript_consequences' in data[0]:
                for transcript in data[0]['transcript_consequences']:
                    # We want to match the gene symbol if possible
                    if transcript.get('gene_symbol') == gene or gene == 'Unknown':
                        sift_score = transcript.get('sift_score', sift_score)
                        sift_pred = transcript.get('sift_prediction', sift_pred)
                        poly_score = transcript.get('polyphen_score', poly_score)
                        poly_pred = transcript.get('polyphen_prediction', poly_pred)
                        break  # Found the primary transcript, stop looking

            results_list.append({
                "rsID": rsid,
                "CHR": chrom,
                "BP": pos,
                "Gene": gene,
                "Alt_Allele_from_BIM": alt_allele,
                "SIFT_Prediction": sift_pred,
                "SIFT_Score": sift_score,
                "PolyPhen_Prediction": poly_pred,
                "PolyPhen_Score": poly_score
            })
            print(f" -> Processed {rsid} ({gene}): SIFT={sift_pred}, PolyPhen={poly_pred}")

        else:
            print(f"[WARNING] Failed API request for {rsid}: HTTP {response.status_code}")

    except Exception as e:
        print(f"[ERROR] Connection issue for {rsid}: {e}")

    # Be polite to the REST API
    time.sleep(0.3)

# 4. Save and merge results
df_scores = pd.DataFrame(results_list)

# Standardize the column name if it was changed to 'POS'
df_missense.rename(columns={'POS': 'BP'}, inplace=True, errors='ignore')

# Force the columns to be strings so they match perfectly
df_missense['CHR'] = df_missense['CHR'].astype(str)
df_missense['BP'] = df_missense['BP'].astype(str)
df_scores['CHR'] = df_scores['CHR'].astype(str)
df_scores['BP'] = df_scores['BP'].astype(str)

# Now they will merge without crashing!
df_final = pd.merge(df_missense, df_scores, on=["rsID", "CHR", "BP"], how="left", suffixes=('', '_drop'))

# Clean up duplicate gene column if it exists
if 'Gene_drop' in df_final.columns:
    df_final.drop(columns=['Gene_drop'], inplace=True)

# Rename it back to POS right before saving
df_final.rename(columns={'BP': 'POS'}, inplace=True, errors='ignore')

df_final.to_csv(output_file, index=False)
print(f"\n[SUCCESS] PolyPhen and SIFT analysis complete! Saved to {output_file}")
