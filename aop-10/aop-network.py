import sys
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import io
import requests
import json
import re
from SPARQLWrapper import SPARQLWrapper, CSV

# ── SPARQL query: Key Events → AOPs ───────────────────────────────────────────
endpoint_url = "https://aopwiki.rdf.bigcat-bioinformatics.org/sparql"

query = """
PREFIX aopo: <http://aopkb.org/aop_ontology#>
PREFIX dc:   <http://purl.org/dc/elements/1.1/>

SELECT DISTINCT ?ke ?keLabel ?aop ?aopLabel
WHERE {
  ?aop a aopo:AdverseOutcomePathway ;
       aopo:has_key_event ?ke ;
       dc:title ?aopLabel .
  ?ke dc:title ?keLabel .
}
ORDER BY ?keLabel
"""

print("Querying AOP-Wiki SPARQL endpoint...")
sparql = SPARQLWrapper(endpoint_url)
sparql.setQuery(query)
sparql.setReturnFormat(CSV)

raw = sparql.query().convert()
df  = pd.read_csv(io.BytesIO(raw))
print(f"Shape: {df.shape}")
df.to_csv("aop_ke_list.tsv", sep="\t", index=False)
print("Saved → aop_ke_list.tsv")

# ── Filter to AOP 10 ──────────────────────────────────────────────────────────
df_aop10 = df[df["aop"] == "https://identifiers.org/aop/10"].copy()

# Extract numeric KE ID from URL e.g. https://identifiers.org/aop.events/667 → 667
df_aop10["ke_id"] = df_aop10["ke"].str.extract(r"/(\d+)$")
print(f"\nAOP 10 KEs: {df_aop10['ke_id'].tolist()}")

df_aop10.to_csv("aop10_ke_list.tsv", sep="\t", index=False)
print("Saved → aop10_ke_list.tsv")

# ── Fetch KE details from AOP-Wiki API ────────────────────────────────────────
def fetch_ke(ke_id):
    url  = f"https://aopwiki.org/events/{ke_id}.json"
    resp = requests.get(url, timeout=10)
    if resp.status_code != 200 or not resp.text.strip():
        print(f"  [WARN] No data for KE {ke_id} (status {resp.status_code})")
        return None
    data = resp.json()

    components = []
    for ec in data.get("event_components", []):
        components.append({
            "ec_id":        ec.get("id"),
            "process_term": ec.get("process", {}).get("term"),
            "process_id":   ec.get("process", {}).get("source_id"),
            "object_term":  ec.get("object",  {}).get("term"),
            "object_id":    ec.get("object",  {}).get("source_id"),
            "action_term":  ec.get("action",  {}).get("term"),
            "action_id":    ec.get("action",  {}).get("source_id"),
        })

    return {
        "ke_id":        data.get("id"),
        "title":        data.get("title"),
        "short_name":   data.get("short_name"),
        "bio_org":      data.get("biological_organization"),
        "cell_term":    data.get("cell_term",  {}).get("official_name"),
        "cell_term_id": data.get("cell_term",  {}).get("source_id"),
        "organ_term":   data.get("organ_term", {}).get("official_name"),
        "organ_term_id":data.get("organ_term", {}).get("source_id"),
        "event_components": components,
    }

# ── Loop over individual KE IDs (fix: ke_id not ke_ids) ──────────────────────
ke_ids = df_aop10["ke_id"].tolist()   # numeric IDs: ['667', '682', ...]

rows = []
for ke_id in ke_ids:                  # FIX 1: iterate ke_id
    print(f"Fetching KE {ke_id}...")
    result = fetch_ke(ke_id)          # FIX 2: pass ke_id not ke_ids
    if result is None:
        continue
    if result["event_components"]:
        for ec in result["event_components"]:
            rows.append({**{k: v for k, v in result.items() if k != "event_components"}, **ec})
    else:
        rows.append({k: v for k, v in result.items() if k != "event_components"})

df_ke = pd.DataFrame(rows)
print(f"\nKE details shape: {df_ke.shape}")
print(df_ke.to_string())
df_ke.to_csv("aop10_ke_details.tsv", sep="\t", index=False)
print("Saved → aop10_ke_details.tsv")

##
df_ke.subset = df_ke[["ke_id", "process_id"]]
df_ke.subset = df_ke.subset.drop_duplicates()
print(df_ke.subset)
  
from gseapy import Biomart
import pandas as pd

bm = Biomart()

# Get unique GO IDs, drop NBO (not a GO term)
go_ids = df_ke["process_id"].dropna()
go_ids = go_ids[go_ids.str.startswith("GO:")].unique().tolist()
print(f"GO IDs to query: {go_ids}")

# filters must be a dict: {filter_name: value_or_list}
results = bm.query(
    dataset="hsapiens_gene_ensembl",
    attributes=["ensembl_gene_id", "external_gene_name", "entrezgene_id", "go_id"],
    filters={"go": go_ids},   # "go" is the Biomart filter name for GO IDs
)

print(results.shape)
print(results.tail(10))

results.to_csv("aop10_ke_genes_biomart.tsv", sep="\t", index=False)
print("Saved → aop10_ke_genes_biomart.tsv")

filters=bm.get_filters(dataset="hsapiens_gene_ensembl")
filters.to_csv("filters.tsv", sep="\t", index=False)
print("Saved → filters.tsv")

# Keep rows WITH these GO IDs
df_subset = results[results["go_id"].isin(go_ids)]
df_subset.to_csv("aop10_ke_genes_biomart_subset.tsv", sep="\t", index=False)
print("Saved → aop10_ke_genes_biomart_subset.tsv")

df_subset["go_id"].value_counts()

# df_ke.subset has: ke_id, process_id
# df_go has: ensembl_gene_id, external_gene_name, entrezgene_id, go_id

df_merged = df_subset.merge(
    df_ke.subset[["ke_id", "process_id"]],
    left_on="go_id",
    right_on="process_id",
    how="left"
)

# Drop redundant process_id column
df_merged = df_merged.drop(columns=["process_id"])
new_row = pd.DataFrame([{
    "ensembl_gene_id":    "",
    "external_gene_name": "",
    "entrezgene_id":      "",
    "go_id":              "",
    "ke_id":              613
}])

df_merged = pd.concat([df_merged, new_row], ignore_index=True)
print(df_merged.tail(3))

print(df_merged.shape)
print(df_merged.head(10).to_string())

df_merged.to_csv("aop10_go_genes_with_ke.tsv", sep="\t", index=False)
print("Saved → aop10_go_genes_with_ke.tsv")


