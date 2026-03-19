import xml.etree.ElementTree as ET
import pandas as pd

# ── Parse Orphanet XML ────────────────────────────────────────────────────────
with open("orphanet-list.xml", "rb") as f:
    tree = ET.parse(f)
root = tree.getroot()

# ── Extract disease-gene associations ─────────────────────────────────────────
rows = []
for disorder in root.findall(".//Disorder"):
    disorder_id   = disorder.get("id", "N/A")
    disorder_name = disorder.findtext("Name", "N/A").strip()
    orpha_code    = disorder.findtext("OrphaCode", "N/A").strip()

    genes = disorder.findall(".//Gene")

    if not genes:
        rows.append({
            "Disorder Name":    disorder_name,
            "OrphaCode":        orpha_code,
            "Disorder ID":      disorder_id,
            "Gene Name":        None,
            "Gene Symbol":      None,
            "Ensembl Reference": None,
        })
        continue

    for gene in genes:
        gene_name   = gene.findtext("Name",   "N/A").strip()
        gene_symbol = gene.findtext("Symbol", "N/A").strip()

        # Extract Ensembl ID from ExternalReference list
        ensembl_id = None
        for ext_ref in gene.findall(".//ExternalReference"):
            source = ext_ref.findtext("Source", "").strip()
            ref    = ext_ref.findtext("Reference", "").strip()
            if source == "Ensembl":
                ensembl_id = ref
                break

        rows.append({
            "Disorder Name":     disorder_name,
            "OrphaCode":         orpha_code,
            "Disorder ID":       disorder_id,
            "Gene Name":         gene_name,
            "Gene Symbol":       gene_symbol,
            "Ensembl Reference": ensembl_id,
        })

# ── Build DataFrame ───────────────────────────────────────────────────────────
df = pd.DataFrame(rows)

print(f"Shape: {df.shape}")
print(f"Unique diseases : {df['OrphaCode'].nunique():,}")
print(f"Unique genes    : {df['Gene Symbol'].nunique():,}")
print()
print(df.head(10).to_string())

# ── Save ───────────────────────────────────────────────────────────────────────
df.to_csv("orphanet_gene_disease.csv", index=False)
print("\nSaved → orphanet_gene_disease.csv")