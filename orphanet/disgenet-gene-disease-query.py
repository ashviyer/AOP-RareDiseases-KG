import warnings
warnings.filterwarnings("ignore")

import requests
import pandas as pd
import time

# ── Config ────────────────────────────────────────────────────────────────────
API_KEY = "68f26dab-c675-4d4e-97e8-dc3f982f7e1f"

HEADERS = {
    "Authorization": API_KEY,
    "accept":        "application/json",
}

BATCH_SIZE = 100

# ── Load UMLS CUIs ────────────────────────────────────────────────────────────
df_umls   = pd.read_csv("orphacode_umls_only.tsv", sep="\t", dtype=str)
df_umls   = df_umls.dropna(subset=["UMLS_CUI"])
cui_list  = df_umls["UMLS_CUI"].tolist()
print(f"Total UMLS CUIs : {len(cui_list):,}")
print(f"Batch size      : {BATCH_SIZE}")
print(f"Total batches   : {(len(cui_list) + BATCH_SIZE - 1) // BATCH_SIZE}")

# Build lookup: CUI → (OrphaCode, Disorder_Name)
cui_to_orpha = dict(zip(df_umls["UMLS_CUI"], df_umls["OrphaCode"]))
cui_to_name  = dict(zip(df_umls["UMLS_CUI"], df_umls["Disorder_Name"]))

# ── Batch query ───────────────────────────────────────────────────────────────
print("\nQuerying DisGeNET in batches...")
all_results = []

batches = [cui_list[i:i + BATCH_SIZE]
           for i in range(0, len(cui_list), BATCH_SIZE)]

for b_idx, batch in enumerate(batches):
    # Format: UMLS_C0020179 (prefix required)
    disease_param = ",".join(f"UMLS_{cui}" for cui in batch)

    params = {
        "disease":     disease_param,
        "page_number": "0",
    }

    try:
        response = requests.get(
            "https://api.disgenet.com/api/v1/gda/summary",
            params=params,
            headers=HEADERS,
            verify=False,
            timeout=60,
        )

        if response.status_code == 429:
            print(f"  [429] Rate limited — waiting 60s...")
            time.sleep(60)
            response = requests.get(
                "https://api.disgenet.com/api/v1/gda/summary",
                params=params,
                headers=HEADERS,
                verify=False,
                timeout=60,
            )

        data = response.json()

        if data.get("status") != "OK":
            print(f"  [Batch {b_idx+1}] Failed: "
                  f"{data.get('payload', {}).get('details', data.get('status'))}")
            time.sleep(1)
            continue

        payload = data.get("payload", [])

        for row in payload:
            cui = row.get("disease_id", "").replace("UMLS_", "")
            row["OrphaCode"]     = cui_to_orpha.get(cui, "")
            row["Disorder_Name"] = cui_to_name.get(cui, "")
        all_results.extend(payload)

    except Exception as e:
        print(f"  [Batch {b_idx+1}] Error: {e}")

    if (b_idx + 1) % 10 == 0:
        print(f"  Batch {b_idx+1:,}/{len(batches):,} — "
              f"{len(all_results):,} GDAs collected")
        # Checkpoint save
        if all_results:
            pd.DataFrame(all_results).to_csv(
                "disgenet_orphanet_gda_checkpoint.tsv", sep="\t", index=False)

    time.sleep(1)

# ── Save final ────────────────────────────────────────────────────────────────
df = pd.DataFrame(all_results)
print(f"\nTotal GDAs      : {len(df):,}")
print(f"Columns         : {df.columns.tolist()}")
print(f"Unique diseases : {df['disease_id'].nunique():,}" if 'disease_id' in df.columns else "")
print(f"Unique genes    : {df['gene_symbol'].nunique():,}" if 'gene_symbol' in df.columns else "")
print(df.head(10).to_string())

df.to_csv("disgenet_orphanet_gda.tsv", sep="\t", index=False)
print("\nSaved → disgenet_orphanet_gda.tsv")