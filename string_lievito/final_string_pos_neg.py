"""
final_string_pos_neg.py
-----------------------
1. Builds a gene-name → UniProt Entry mapping from uniprot_yeast.tsv
2. Translates STRING protein IDs (systematic gene names) to UniProt accessions
3. Keeps only pairs whose proteins are present in unified_protein_filter_ds.csv
   (i.e. sequence available and length 50-800 aa)
4. Generates negative pairs at a 10:1 (neg:pos) ratio by random pairing,
   ensuring no pair appears in the positive set
5. Saves  string_final_dataset.csv  with columns: protein1, protein2, label
"""

import random
import pandas as pd

random.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
STRING_CSV   = "string_high_confidence.csv"
UNIPROT_TSV  = "../uniprot_lievito/uniprot_yeast.tsv"
FILTERED_DS  = "../uniprot_lievito/unified_protein_filter_ds.csv"
OUTPUT_CSV   = "string_final_dataset.csv"

# ── 1. Build gene-name token → UniProt Entry mapping ─────────────────────────
print("Building gene-name → UniProt Entry map …")
uni = pd.read_csv(UNIPROT_TSV, sep="\t", dtype=str, keep_default_na=False)
uni.columns = [c.strip() for c in uni.columns]

gene_to_entry = {}   # any token in Gene Names → UniProt accession
for _, row in uni.iterrows():
    entry = row["Entry"]
    for token in row["Gene Names"].split():
        t = token.strip().upper()
        if t:
            gene_to_entry[t] = entry

print(f"  {len(gene_to_entry):,} gene-name tokens mapped to {uni['Entry'].nunique():,} UniProt entries.\n")

# ── 2. Load STRING positives ──────────────────────────────────────────────────
print("Loading STRING high-confidence positives …")
df_str = pd.read_csv(STRING_CSV, dtype=str)
print(f"  {len(df_str):,} positive pairs before translation.\n")

# ── 3. Translate STRING gene names → UniProt accessions ──────────────────────
print("Translating STRING IDs …")
df_str["entry1"] = df_str["protein1"].str.upper().map(gene_to_entry)
df_str["entry2"] = df_str["protein2"].str.upper().map(gene_to_entry)

# Drop pairs where either protein couldn't be translated
before = len(df_str)
df_str = df_str.dropna(subset=["entry1", "entry2"])
print(f"  Translated: {len(df_str):,} / {before:,} pairs ({before - len(df_str):,} dropped – no UniProt match).\n")

# ── 4. Keep only proteins present in the filtered unified dataset ─────────────
print("Filtering to proteins in unified_protein_filter_ds.csv …")
ds = pd.read_csv(FILTERED_DS, dtype=str, keep_default_na=False)
valid_proteins = set(ds["Entry"].dropna())
print(f"  Valid proteins in filtered dataset: {len(valid_proteins):,}")

before = len(df_str)
df_str = df_str[
    df_str["entry1"].isin(valid_proteins) &
    df_str["entry2"].isin(valid_proteins)
].copy()
print(f"  Positive pairs after filter: {len(df_str):,} / {before:,}\n")

# ── 5. Build positive set ─────────────────────────────────────────────────────
# Canonical form: always (min, max) to catch duplicates regardless of order
pos_pairs_set = set(
    tuple(sorted([r.entry1, r.entry2]))
    for r in df_str.itertuples()
)
df_positives = pd.DataFrame(
    [(a, b, 1) for a, b in pos_pairs_set],
    columns=["protein1", "protein2", "label"]
)
n_pos = len(df_positives)
print(f"Unique positive pairs: {n_pos:,}")

# ── 6. Generate negatives at 10:1 ratio ──────────────────────────────────────
n_neg_target = n_pos * 10
print(f"Generating {n_neg_target:,} negative pairs (ratio 10:1) …")

protein_pool = list(valid_proteins)
neg_pairs = set()

max_attempts = n_neg_target * 20
attempts = 0
while len(neg_pairs) < n_neg_target and attempts < max_attempts:
    a, b = random.sample(protein_pool, 2)
    pair = tuple(sorted([a, b]))
    if pair not in pos_pairs_set and pair not in neg_pairs:
        neg_pairs.add(pair)
    attempts += 1

print(f"  Generated: {len(neg_pairs):,} negative pairs.\n")

df_negatives = pd.DataFrame(
    [(a, b, 0) for a, b in neg_pairs],
    columns=["protein1", "protein2", "label"]
)

# ── 7. Combine, shuffle, save ─────────────────────────────────────────────────
df_final = pd.concat([df_positives, df_negatives], ignore_index=True)
df_final = df_final.sample(frac=1, random_state=42).reset_index(drop=True)

df_final.to_csv(OUTPUT_CSV, index=False)

print(f"── Final dataset ────────────────────────────────────────────────────")
print(f"  Positive pairs (label=1): {(df_final['label']==1).sum():,}")
print(f"  Negative pairs (label=0): {(df_final['label']==0).sum():,}")
print(f"  Total pairs             : {len(df_final):,}")
print(f"\nSaved → {OUTPUT_CSV}")
