import pandas as pd
import requests
import time
import sys

# === CONFIGURATION ===
# Using the merged file from your VEP step
input_file = "VEP_Annotations_Results.csv"
output_file = "Intergenic_ENCODE_Annotations.csv"

# CRITICAL: Using the GRCh37 (hg19) specific server to match your LocusZoom data
server = "https://grch37.rest.ensembl.org"

print("[INFO] Loading VEP annotated variants...")

try:
    df = pd.read_csv(input_file)
except FileNotFoundError:
    print(f"[ERROR] Could not find {input_file}. Check your working directory.")
    sys.exit(1)

# 1. Filter for intergenic variants
# (Adjust 'Consequence' if your column is named differently)
if 'Consequence' not in df.columns:
    print("[ERROR] Could not find the 'Consequence' column from the VEP output.")
    sys.exit(1)

# Keep only rows where the consequence includes the word 'intergenic'
df_intergenic = df[df['Consequence'].str.contains('intergenic', case=False, na=False)].copy()

print(f"[INFO] Found {len(df_intergenic)} intergenic variants to check against ENCODE data.")

if len(df_intergenic) == 0:
    print("[INFO] No intergenic variants found. Exiting.")
    sys.exit(0)

# 2. Prepare the API Request
headers = {
    "Content-Type": "application/json"
}

results_list = []

print("[INFO] Querying ENCODE/Regulatory database for overlapping features...")

# 3. Iterate through each intergenic variant
# Assuming your columns for chromosome and position are named 'CHR' and 'BP' (or 'CHROM' and 'POS')
chrom_col = 'CHR' if 'CHR' in df.columns else 'CHROM'
pos_col = 'BP' if 'BP' in df.columns else 'POS'

for index, row in df_intergenic.iterrows():
    chrom = row[chrom_col]
    pos = row[pos_col]
    rsid = row.get('SNP', f"{chrom}:{pos}")

    # Endpoint for overlapping regulatory features
    endpoint = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=regulatory"

    try:
        response = requests.get(server + endpoint, headers=headers)

        if response.ok:
            data = response.json()

            # If the API returns a list with items, we found an ENCODE feature!
            if len(data) > 0:
                for feature in data:
                    results_list.append({
                        "rsID": rsid,
                        "CHR": chrom,
                        "BP": pos,
                        "Regulatory_Type": feature.get('description', 'Unknown Regulatory Feature'),
                        "Epigenome_Activity": feature.get('activity', 'N/A')
                    })
            else:
                # No regulatory feature found at this exact base pair
                results_list.append({
                    "rsID": rsid,
                    "CHR": chrom,
                    "BP": pos,
                    "Regulatory_Type": "None Found",
                    "Epigenome_Activity": "N/A"
                })

        else:
            print(f"[WARNING] Failed to fetch data for {chrom}:{pos} - Status: {response.status_code}")

    except Exception as e:
        print(f"[ERROR] Connection issue for {chrom}:{pos}: {e}")

    # Be polite to the server to avoid getting blocked
    time.sleep(0.3)

# 4. Save the results
print("[INFO] Formatting and saving results...")
df_regulatory = pd.DataFrame(results_list)

# ---> THE FIX: Standardize the column name if it was changed to 'POS' <---
df_intergenic.rename(columns={'POS': 'BP'}, inplace=True, errors='ignore')
df_regulatory.rename(columns={'POS': 'BP'}, inplace=True, errors='ignore')

# (This is your original line 99)
df_final = pd.merge(df_intergenic, df_regulatory, on=["rsID", "CHR", "BP"], how="left")
df_final.to_csv(output_file, index=False)

print(f"[SUCCESS] ENCODE analysis complete! Results saved to {output_file}")