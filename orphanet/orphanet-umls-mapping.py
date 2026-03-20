import warnings
warnings.filterwarnings("ignore")

import requests
import xml.etree.ElementTree as ET
import pandas as pd

# ── Download en_product1.xml (disorder cross-references) ─────────────────────
# This is the Orphanet file with UMLS, ICD-10, MeSH etc at disorder level
URL = "https://www.orphadata.com/data/xml/en_product1.xml"

print(f"Downloading {URL} ...")
r = requests.get(URL, timeout=300, verify=False)
r.raise_for_status()

with open("en_product1.xml", "wb") as f:
    f.write(r.content)
print(f"Downloaded {len(r.content) / 1e6:.1f} MB → en_product1.xml")

# ── Parse ─────────────────────────────────────────────────────────────────────
print("Parsing...")
with open("en_product1.xml", "rb") as f:
    tree = ET.parse(f)
root = tree.getroot()

rows = []

for disorder in root.findall(".//Disorder"):
    orpha_code    = disorder.findtext("OrphaCode")
    disorder_name = disorder.findtext("Name")
    disorder_id   = disorder.get("id")

    # Collect all cross-references at disorder level
    umls_cuis = []
    icd10     = []
    omim      = []
    mesh      = []

    for xref in disorder.findall(".//ExternalReference"):
        source = xref.findtext("Source", "").strip()
        ref    = xref.findtext("Reference", "").strip()
        if source == "UMLS":
            umls_cuis.append(ref)
        elif source in ("ICD-10", "ICD10"):
            icd10.append(ref)
        elif source == "OMIM":
            omim.append(ref)
        elif source in ("MeSH", "MESH"):
            mesh.append(ref)

    # One row per UMLS CUI (some disorders have multiple)
    if umls_cuis:
        for cui in umls_cuis:
            rows.append({
                "OrphaCode":      orpha_code,
                "Disorder_ID":    disorder_id,
                "Disorder_Name":  disorder_name,
                "UMLS_CUI":       cui,
                "ICD10":          "|".join(icd10) if icd10 else None,
                "OMIM":           "|".join(omim)  if omim  else None,
                "MeSH":           "|".join(mesh)  if mesh  else None,
            })
    else:
        rows.append({
            "OrphaCode":     orpha_code,
            "Disorder_ID":   disorder_id,
            "Disorder_Name": disorder_name,
            "UMLS_CUI":      None,
            "ICD10":         "|".join(icd10) if icd10 else None,
            "OMIM":          "|".join(omim)  if omim  else None,
            "MeSH":          "|".join(mesh)  if mesh  else None,
        })

df = pd.DataFrame(rows)

print(f"\nTotal disorders         : {df['OrphaCode'].nunique():,}")
print(f"With UMLS CUI           : {df['UMLS_CUI'].notna().sum():,}")
print(f"Without UMLS CUI        : {df['UMLS_CUI'].isna().sum():,}")
print()
print(df[df["UMLS_CUI"].notna()].head(10).to_string())

# ── Save ──────────────────────────────────────────────────────────────────────
df.to_csv("orphacode_umls_mapping.tsv", sep="\t", index=False)
print("\nSaved → orphacode_umls_mapping.tsv")

# Also save UMLS-only rows for DisGeNET queries
df_umls = df[df["UMLS_CUI"].notna()][["OrphaCode", "Disorder_Name", "UMLS_CUI"]]
df_umls.to_csv("orphacode_umls_only.tsv", sep="\t", index=False)
print(f"Saved → orphacode_umls_only.tsv ({len(df_umls):,} rows)")