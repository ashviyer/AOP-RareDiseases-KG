import pandas as pd

# ── Load the three dataframes ─────────────────────────────────────────────────
df_disgenet = pd.read_csv("disgenet_orphanet_gda_subset.tsv", sep="\t", dtype=str)
df_orphanet = pd.read_csv("orphanet_gene_disease.csv", dtype=str)
df_mapping  = pd.read_csv("orphacode_umls_mapping.tsv", sep="\t", dtype=str)

print(f"DisGeNET GDA     : {df_disgenet.shape}")
print(f"Orphanet GDA     : {df_orphanet.shape}")
print(f"UMLS mapping     : {df_mapping.shape}")

# ── Standardise column names ──────────────────────────────────────────────────
# DisGeNET: symbolOfGene, diseaseUMLSCUI, geneEnsemblIDs, diseaseName
df_disgenet = df_disgenet.rename(columns={
    "symbolOfGene":  "Gene_Symbol",
    "diseaseUMLSCUI":"UMLS_CUI",
    "geneEnsemblIDs":"Ensembl_ID",
    "diseaseName":   "Disease_Name",
})

# Orphanet: Disorder Name, OrphaCode, Disorder ID, Gene Name, Gene Symbol, Ensembl Reference
df_orphanet = df_orphanet.rename(columns={
    "Disorder Name":    "Disease_Name",
    "OrphaCode":        "OrphaCode",
    "Disorder ID":      "Disorder_ID",
    "Gene Name":        "Gene_Name",
    "Gene Symbol":      "Gene_Symbol",
    "Ensembl Reference":"Ensembl_ID",
})

# Mapping: OrphaCode, Disorder_ID, Disorder_Name, UMLS_CUI, ICD10, OMIM, MeSH
df_mapping = df_mapping.rename(columns={"Disorder_Name": "Disease_Name"})

# ── Merge DisGeNET with UMLS mapping to get OrphaCode ────────────────────────
df_disgenet_mapped = df_disgenet.merge(
    df_mapping[["OrphaCode", "Disorder_ID", "UMLS_CUI", "ICD10", "OMIM"]],
    on="UMLS_CUI",
    how="left",
)
df_disgenet_mapped["source"] = "DisGeNET"

# ── Prepare Orphanet GDA ──────────────────────────────────────────────────────
df_orphanet_prep = df_orphanet.merge(
    df_mapping[["OrphaCode", "UMLS_CUI", "ICD10", "OMIM"]],
    on="OrphaCode",
    how="left",
)
df_orphanet_prep["source"] = "Orphanet"

# ── Select common columns and combine ─────────────────────────────────────────
COLS = ["Disease_Name", "OrphaCode", "Disorder_ID", "UMLS_CUI",
        "ICD10", "OMIM", "Gene_Symbol", "Ensembl_ID", "source"]

df_combined = pd.concat([
    df_disgenet_mapped.reindex(columns=COLS),
    df_orphanet_prep.reindex(columns=COLS),
], ignore_index=True)

# ── Deduplicate ───────────────────────────────────────────────────────────────
df_combined = df_combined.drop_duplicates(
    subset=["OrphaCode", "UMLS_CUI", "Gene_Symbol"]
)

df_combined = df_combined.sort_values(["Disease_Name", "Gene_Symbol"])

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\nCombined shape   : {df_combined.shape}")
print(f"Unique diseases  : {df_combined['Disease_Name'].nunique():,}")
print(f"Unique OrphaCode : {df_combined['OrphaCode'].nunique():,}")
print(f"Unique UMLS CUIs : {df_combined['UMLS_CUI'].nunique():,}")
print(f"Unique genes     : {df_combined['Gene_Symbol'].nunique():,}")
print(f"\nSources:")
print(df_combined["source"].value_counts().to_string())
print()
print(df_combined.head(10).to_string())

# ── Save ──────────────────────────────────────────────────────────────────────
df_combined.to_csv("disease_gene_combined.tsv", sep="\t", index=False)
print("\nSaved → disease_gene_combined.tsv")

print(f"Shape: {df_combined.shape}")
print(f"Unique diseases : {df_combined['OrphaCode'].nunique():,}")
print(f"Unique genes    : {df_combined['Gene_Symbol'].nunique():,}")
print()