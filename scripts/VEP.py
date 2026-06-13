import pandas as pd
import requests
import json
import time
import sys

# === CONFIGURATION ===
input_file = "Target_Regions_Global_FDR_Significant.csv"
output_file = "VEP_Annotations_Results.csv"
server = "https://rest.ensembl.org"
endpoint = "/vep/human/id"

print("[INFO] Loading significant variants...")

# 1. Load the CSV
try:
    df = pd.read_csv(input_file)
except FileNotFoundError:
    print(f"[ERROR] Could not find {input_file}. Check your working directory.")
    sys.exit(1)

# 2. Find the column containing the variant IDs (rsIDs)
# We will guess common column names, but you can hardcode this if yours is different!
id_col = None
for col in ['rsID', 'SNP', 'ID', 'Variant_ID']:
    if col in df.columns:
        id_col = col
        break

if not id_col:
    print("[ERROR] Could not find a column named 'rsID', 'SNP', or 'ID'. Please specify the column name.")
    sys.exit(1)

# Get unique list of IDs, dropping any NaNs
variant_ids = df[id_col].dropna().unique().tolist()
print(f"[INFO] Found {len(variant_ids)} unique variants to query.")

# If you have zero variants, stop here.
if len(variant_ids) == 0:
    print("[INFO] No variants to process. Exiting.")
    sys.exit(0)

# 3. Prepare the API Request
headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

# Ensembl allows up to 200 variants per request.
# We will batch them just in case you have more than that.
batch_size = 200
all_vep_results = []

print("[INFO] Querying Ensembl VEP API...")

for i in range(0, len(variant_ids), batch_size):
    batch = variant_ids[i:i + batch_size]
    data = json.dumps({"ids": batch})

    # Send the POST request to Ensembl
    response = requests.post(server + endpoint, headers=headers, data=data)

    if not response.ok:
        print(f"[ERROR] API request failed with status {response.status_code}")
        print(response.text)
        sys.exit(1)

    # Parse the JSON response
    results = response.json()

    # Extract the juicy biological info we care about
    for variant in results:
        v_id = variant.get('id', 'Unknown')
        most_severe_consequence = variant.get('most_severe_consequence', 'Unknown')

        # A variant can map to multiple transcripts (genes).
        # We will pull the first overlapping gene name if it exists.
        gene_symbol = "Intergenic"
        biotype = "N/A"

        if 'transcript_consequences' in variant:
            # Just grabbing the first transcript's gene symbol for simplicity
            for transcript in variant['transcript_consequences']:
                if 'gene_symbol' in transcript:
                    gene_symbol = transcript['gene_symbol']
                    biotype = transcript.get('biotype', 'N/A')
                    break  # Stop after finding the primary gene

        all_vep_results.append({
            "rsID": v_id,
            "Gene_Symbol": gene_symbol,
            "Consequence": most_severe_consequence,
            "Gene_Biotype": biotype
        })

    # Be polite to the server, pause briefly between batches
    time.sleep(1)

# 4. Save the results
print("[INFO] Formatting and saving results...")
df_vep = pd.DataFrame(all_vep_results)

# Merge the new annotations back into your original significant variants dataframe
df_final = pd.merge(df, df_vep, left_on=id_col, right_on="rsID", how="left")
df_final.to_csv(output_file, index=False)

print(f"[SUCCESS] VEP analysis complete! Results saved to {output_file}")
