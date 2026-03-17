import requests
import zipfile
import os

URL      = "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip"
ZIP_PATH = "human_intact.zip"
OUT_DIR  = "."

# Download
print("Downloading IntAct human PPI...")
r = requests.get(URL, stream=True, timeout=300)
r.raise_for_status()

with open(ZIP_PATH, "wb") as f:
    downloaded = 0
    for chunk in r.iter_content(chunk_size=1 << 16):
        f.write(chunk)
        downloaded += len(chunk)
        print(f"\r  {downloaded / 1e6:.1f} MB downloaded", end="", flush=True)

print(f"\nSaved: {ZIP_PATH} ({os.path.getsize(ZIP_PATH) / 1e6:.1f} MB)")

# Validate
if not zipfile.is_zipfile(ZIP_PATH):
    raise RuntimeError("Downloaded file is not a valid ZIP.")

# Unzip
print("Unzipping...")
with zipfile.ZipFile(ZIP_PATH, "r") as zf:
    all_files = zf.namelist()
    print(f"  Archive contains: {all_files}")
    zf.extractall(OUT_DIR)

# Remove everything except human.txt
print("Cleaning up — keeping only human.txt...")
for fname in all_files:
    fpath = os.path.join(OUT_DIR, fname)
    if fname != "human.txt" and os.path.isfile(fpath):
        os.remove(fpath)
        print(f"  Removed: {fname}")

print(f"Done. Kept: human.txt ({os.path.getsize('human.txt') / 1e6:.1f} MB)")

import pandas as pd

# Load
print("Loading human.txt...")
df = pd.read_csv("human.txt", sep="\t", header=0, low_memory=False, dtype=str)
print(f"Raw shape: {df.shape}")

df.columns = [c.strip() for c in df.columns]

ID_A  = "#ID(s) interactor A"
ID_B  = "ID(s) interactor B"
TAX_A = "Taxid interactor A"
TAX_B = "Taxid interactor B"
SCORE = "Confidence value(s)"

def get_uniprot(field):
    if not isinstance(field, str):
        return None
    for part in field.split("|"):
        if part.startswith("uniprotkb:"):
            return part[10:]
    return None

def get_taxid(field):
    if not isinstance(field, str):
        return None
    for part in field.split("|"):
        if part.startswith("taxid:"):
            return part[6:].split("(")[0]
    return None

def get_score(field):
    if not isinstance(field, str):
        return None
    for part in field.split("|"):
        if part.startswith("intact-miscore:"):
            try:
                return float(part[15:])
            except ValueError:
                pass
    return None

# Build slim dataframe
print("Extracting PPI + score...")
result = pd.DataFrame({
    "uniprotA": df[ID_A].apply(get_uniprot),
    "uniprotB": df[ID_B].apply(get_uniprot),
    "taxA":     df[TAX_A].apply(get_taxid),
    "taxB":     df[TAX_B].apply(get_taxid),
    "mi_score": df[SCORE].apply(get_score),
})

# Keep only human-human (taxid 9606) with valid UniProt IDs
mask = (
    (result["taxA"] == "9606") &
    (result["taxB"] == "9606") &
    result["uniprotA"].notna() &
    result["uniprotB"].notna()
)
result = result[mask][["uniprotA", "uniprotB", "mi_score"]].reset_index(drop=True)

print(f"\nHuman-human PPIs : {len(result):,}")
print(f"Unique proteins  : {pd.concat([result['uniprotA'], result['uniprotB']]).nunique():,}")
print(f"Score available  : {result['mi_score'].notna().sum():,} / {len(result):,}")
print()
print(result.head(10).to_string())

# Save
result.to_csv("intact_human_ppi_scored.tsv", sep="\t", index=False)
print("\nSaved → intact_human_ppi_scored.tsv")

import pandas as pd
import mygene

# ── Load PPI table ─────────────────────────────────────────────────────────────
print("Loading intact_human_ppi_scored.tsv...")
ppi = pd.read_csv("intact_human_ppi_scored.tsv", sep="\t", dtype=str)
ppi["mi_score"] = pd.to_numeric(ppi["mi_score"], errors="coerce")
print(f"Loaded: {len(ppi):,} interactions")

# ── Collect unique UniProt IDs ─────────────────────────────────────────────────
uniprot_ids = list(
    set(ppi["uniprotA"].dropna().tolist() + ppi["uniprotB"].dropna().tolist())
)
print(f"Unique UniProt IDs to map: {len(uniprot_ids):,}")

# ── Query MyGene.info ──────────────────────────────────────────────────────────
mg = mygene.MyGeneInfo()
print("Querying MyGene.info (batch, may take ~1 min)...")

hits = mg.querymany(
    uniprot_ids,
    scopes="uniprot",
    fields="entrezgene,ensembl.gene,symbol",
    species="human",
    as_dataframe=False,
    returnall=False,
)

# ── Build mapping dict ─────────────────────────────────────────────────────────
def parse_ensembl(val):
    if isinstance(val, dict):
        return val.get("gene")
    if isinstance(val, list) and val:
        return val[0].get("gene")
    return None

mapping = {}
for hit in hits:
    query = hit.get("query")
    if not query or hit.get("notfound"):
        continue
    ensembl_raw = hit.get("ensembl", hit.get("ensembl.gene"))
    mapping[query] = {
        "symbol":  hit.get("symbol"),
        "entrez":  str(int(hit["entrezgene"])) if "entrezgene" in hit else None,
        "ensembl": parse_ensembl(ensembl_raw),
    }

mapped_n = sum(1 for v in mapping.values() if v["ensembl"] is not None)
print(f"Ensembl mapped: {mapped_n:,} / {len(uniprot_ids):,} UniProt IDs")

# ── Annotate PPI table — explicit column names ─────────────────────────────────
ppi["geneA_symbol"]  = ppi["uniprotA"].map(lambda u: mapping.get(u, {}).get("symbol"))
ppi["geneA_entrez"]  = ppi["uniprotA"].map(lambda u: mapping.get(u, {}).get("entrez"))
ppi["geneA_ensembl"] = ppi["uniprotA"].map(lambda u: mapping.get(u, {}).get("ensembl"))

ppi["geneB_symbol"]  = ppi["uniprotB"].map(lambda u: mapping.get(u, {}).get("symbol"))
ppi["geneB_entrez"]  = ppi["uniprotB"].map(lambda u: mapping.get(u, {}).get("entrez"))
ppi["geneB_ensembl"] = ppi["uniprotB"].map(lambda u: mapping.get(u, {}).get("ensembl"))

# Reorder
ppi = ppi[[
    "uniprotA", "geneA_symbol", "geneA_entrez", "geneA_ensembl",
    "uniprotB", "geneB_symbol", "geneB_entrez", "geneB_ensembl",
    "mi_score",
]]

# ── Summary ────────────────────────────────────────────────────────────────────
both_ensembl = ppi["geneA_ensembl"].notna() & ppi["geneB_ensembl"].notna()
print(f"\nInteractions — both Ensembl IDs : {both_ensembl.sum():,} / {len(ppi):,}")
print(f"Interactions — at least one NaN : {(~both_ensembl).sum():,}")
print()
print(ppi.head(10).to_string())

# ── Save ───────────────────────────────────────────────────────────────────────
ppi.to_csv("intact_human_ppi_annotated.tsv", sep="\t", index=False)
print("\nSaved → intact_human_ppi_annotated.tsv")