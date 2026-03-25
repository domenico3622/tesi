"""
create_unified_dataset.py
-------------------------
Builds a unified dataset starting from uniprot_yeast.tsv:
  - Replaces "Gene Ontology (cellular component)" with the macro-zone columns
    from protein_zones.csv (nucleus, cytoplasm, …, zones)
  - Adds a "sequence" column with amino-acid sequences from the FASTA file,
    matched by UniProt Entry ID.

Output: unified_protein_dataset.csv
"""

import re
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
UNIPROT_TSV   = "uniprot_yeast.tsv"
ZONES_CSV     = "protein_zones.csv"
FASTA_FILE    = "../fasta_lievito/UP000002311_559292.fasta"
OUTPUT_CSV    = "unified_protein_filter_ds.csv"

# ── 1. Parse FASTA → dict { entry_id : sequence } ────────────────────────────
print("Parsing FASTA file …")
fasta_dict = {}          # UniProt accession → sequence
fasta_by_name = {}       # Entry name (e.g. THDH_YEAST) → sequence

header_re = re.compile(r"^>(\w+)\|([A-Za-z0-9_]+)\|(\w+)")

current_id = None
current_name = None
seq_parts = []

with open(FASTA_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Save previous sequence
            if current_id and seq_parts:
                seq = "".join(seq_parts)
                fasta_dict[current_id] = seq
                if current_name:
                    fasta_by_name[current_name] = seq
            # Parse new header:  >sp|P00927|THDH_YEAST ...
            m = header_re.match(line)
            if m:
                current_id   = m.group(2)   # e.g. P00927
                current_name = m.group(3)    # e.g. THDH_YEAST
            else:
                # Fallback: try to grab anything after >
                parts = line[1:].split("|")
                current_id   = parts[1] if len(parts) > 1 else parts[0].split()[0]
                current_name = parts[2].split()[0] if len(parts) > 2 else None
            seq_parts = []
        else:
            seq_parts.append(line)
    # Last entry
    if current_id and seq_parts:
        seq = "".join(seq_parts)
        fasta_dict[current_id] = seq
        if current_name:
            fasta_by_name[current_name] = seq

print(f"  {len(fasta_dict):,} sequences loaded (by accession).")
print(f"  {len(fasta_by_name):,} sequences loaded (by entry name).\n")

# ── 2. Load uniprot_yeast.tsv ────────────────────────────────────────────────
print("Loading uniprot_yeast.tsv …")
df = pd.read_csv(UNIPROT_TSV, sep="\t", dtype=str, keep_default_na=False)
df.columns = [c.strip() for c in df.columns]
print(f"  {len(df):,} proteins.\n")

# ── 3. Load protein_zones.csv ────────────────────────────────────────────────
print("Loading protein_zones.csv …")
zones_df = pd.read_csv(ZONES_CSV, dtype=str, keep_default_na=False)
zone_cols = ["nucleus", "cytoplasm", "mitochondrion", "plasma_membrane",
             "endoplasmic_reticulum", "golgi_apparatus", "vacuole", "zones"]
print(f"  {len(zones_df):,} rows.\n")

# ── 4. Merge zones onto uniprot data (by Entry) ─────────────────────────────
# protein_zones has "entry" which matches uniprot "Entry"
zones_to_merge = zones_df[["entry"] + zone_cols].copy()
zones_to_merge = zones_to_merge.rename(columns={"entry": "Entry"})

merged = df.merge(zones_to_merge, on="Entry", how="left")

# Drop the original GO column
if "Gene Ontology (cellular component)" in merged.columns:
    merged.drop(columns=["Gene Ontology (cellular component)"], inplace=True)

print(f"  After merge: {len(merged):,} rows.")

# ── 5. Add sequence column ──────────────────────────────────────────────────
def lookup_sequence(row):
    """Try matching by Entry (accession), then by Entry Name."""
    entry = row.get("Entry", "")
    name  = row.get("Entry Name", "")
    if entry in fasta_dict:
        return fasta_dict[entry]
    if name in fasta_by_name:
        return fasta_by_name[name]
    return ""

print("Mapping sequences …")
merged["sequence"] = merged.apply(lookup_sequence, axis=1)

found     = (merged["sequence"] != "").sum()
not_found = (merged["sequence"] == "").sum()
print(f"  Sequences found : {found:,}")
print(f"  Sequences missing: {not_found:,}\n")

# ── 6. Filter by sequence length (50 – 800 aa) ──────────────────────────────
print("Filtering proteins by length (50 – 800 aa) …")
before = len(merged)
merged["Length"] = pd.to_numeric(merged["Length"], errors="coerce")
merged = merged[(merged["Length"] >= 50) & (merged["Length"] <= 800)].copy()
print(f"  Kept {len(merged):,} / {before:,} proteins after length filter.\n")

# ── 7. Save ─────────────────────────────────────────────────────────────────
merged.to_csv(OUTPUT_CSV, index=False)
print(f"Saved → {OUTPUT_CSV}  ({len(merged):,} rows, {len(merged.columns)} columns)")
print(f"Columns: {list(merged.columns)}")
print("\nDone.")
