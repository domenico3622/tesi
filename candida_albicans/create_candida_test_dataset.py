"""
create_candida_test_dataset.py
------------------------------
Script unico per Candida albicans che:
1) Legge il file STRING physical links detailed
2) Pulisce gli ID rimuovendo il prefisso tassonomico "237561."
3) Mostra la distribuzione della colonna experimental
4) Filtra le interazioni positive ad alta confidenza (experimental >= 700)
5) Mantiene solo le coppie con proteine presenti nel FASTA filtrato
6) Costruisce il dataset finale con positivi (label=1) e negativi (label=0)
   con rapporto 10:1 (neg:pos), evitando overlap con i positivi
7) Salva il dataset finale in CSV
"""

import random
from pathlib import Path

import pandas as pd


random.seed(42)

# ── Configurazione ────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent

STRING_TXT = BASE_DIR / "237561.protein.physical.links.detailed.v12.0.txt"
FILTERED_FASTA = BASE_DIR / "candida_albicans_filter.fasta"
HIGH_CONF_CSV = BASE_DIR / "string_high_confidence.csv"
OUTPUT_CSV = BASE_DIR / "string_final_test_dataset.csv"

TAXON_PREFIX = "237561."
CONFIDENCE_THRESHOLD = 700
NEG_POS_RATIO = 10


def parse_fasta_ids(fasta_path: Path) -> set[str]:
	valid_ids = set()
	with fasta_path.open("r", encoding="utf-8") as handle:
		for line in handle:
			if line.startswith(">"):
				header = line[1:].strip()
				if header:
					valid_ids.add(header.split()[0].strip().upper())
	return valid_ids


print("Caricamento del database STRING originale...")
df_string = pd.read_csv(STRING_TXT, sep=r"\s+")

print("Pulizia degli ID (rimozione prefisso tassonomico)...")
df_string["protein1"] = df_string["protein1"].astype(str).str.replace(TAXON_PREFIX, "", regex=False)
df_string["protein2"] = df_string["protein2"].astype(str).str.replace(TAXON_PREFIX, "", regex=False)

print("\n── Distribuzione colonna 'experimental' ──────────────────────────────")
total = len(df_string)
bins = {
	"= 0": df_string["experimental"] == 0,
	"1 – 300": (df_string["experimental"] >= 1) & (df_string["experimental"] <= 300),
	"301 – 699": (df_string["experimental"] >= 301) & (df_string["experimental"] <= 699),
	">= 700": df_string["experimental"] >= 700,
}
for label, mask in bins.items():
	count = int(mask.sum())
	print(f"  {label:<12s}  {count:>8,}  ({100 * count / total:.1f} %)")
print(f"  {'TOTALE':<12s}  {total:>8,}\n")

print(f"Filtraggio positivi ad alta confidenza (experimental >= {CONFIDENCE_THRESHOLD})...")
df_positives_raw = df_string[df_string["experimental"] >= CONFIDENCE_THRESHOLD].copy()

cols_to_drop = [c for c in ["database", "textmining", "combined_score"] if c in df_positives_raw.columns]
if cols_to_drop:
	df_positives_raw.drop(columns=cols_to_drop, inplace=True)
	print(f"Colonne rimosse: {cols_to_drop}")

df_positives_raw.to_csv(HIGH_CONF_CSV, index=False)
print(f"Interazioni positive ad alta confidenza salvate in: {HIGH_CONF_CSV.name}")

print("\nCaricamento proteine valide dal FASTA filtrato...")
valid_proteins = parse_fasta_ids(FILTERED_FASTA)
print(f"  Proteine valide nel FASTA filtrato: {len(valid_proteins):,}")

df_positives_raw["protein1"] = df_positives_raw["protein1"].str.upper()
df_positives_raw["protein2"] = df_positives_raw["protein2"].str.upper()

before = len(df_positives_raw)
df_positives_raw = df_positives_raw[
	df_positives_raw["protein1"].isin(valid_proteins)
	& df_positives_raw["protein2"].isin(valid_proteins)
	& (df_positives_raw["protein1"] != df_positives_raw["protein2"])
].copy()
print(f"  Positivi dopo filtro FASTA: {len(df_positives_raw):,} / {before:,}")

pos_pairs_set = {
	tuple(sorted((row.protein1, row.protein2)))
	for row in df_positives_raw.itertuples()
}

df_positives = pd.DataFrame(
	[(a, b, 1) for a, b in pos_pairs_set],
	columns=["protein1", "protein2", "label"],
)
n_pos = len(df_positives)
print(f"  Coppie positive uniche: {n_pos:,}")

protein_pool = sorted(valid_proteins)
max_neg_possible = (len(protein_pool) * (len(protein_pool) - 1) // 2) - n_pos
n_neg_target = n_pos * NEG_POS_RATIO

if n_neg_target > max_neg_possible:
	print(
		f"[ATTENZIONE] Richieste {n_neg_target:,} negative ma massimo possibile è {max_neg_possible:,}. "
		f"Uso il massimo disponibile."
	)
	n_neg_target = max_neg_possible

print(f"Generazione negative: target={n_neg_target:,} (rapporto richiesto {NEG_POS_RATIO}:1)")
neg_pairs = set()
attempts = 0
max_attempts = max(n_neg_target * 20, 1000)

while len(neg_pairs) < n_neg_target and attempts < max_attempts:
	a, b = random.sample(protein_pool, 2)
	pair = tuple(sorted((a, b)))
	if pair not in pos_pairs_set and pair not in neg_pairs:
		neg_pairs.add(pair)
	attempts += 1

if len(neg_pairs) < n_neg_target:
	print(
		f"[ATTENZIONE] Generate {len(neg_pairs):,} negative su {n_neg_target:,} target "
		f"(tentativi massimi raggiunti: {max_attempts:,})."
	)

df_negatives = pd.DataFrame(
	[(a, b, 0) for a, b in neg_pairs],
	columns=["protein1", "protein2", "label"],
)

df_final = pd.concat([df_positives, df_negatives], ignore_index=True)
df_final = df_final.sample(frac=1, random_state=42).reset_index(drop=True)
df_final.to_csv(OUTPUT_CSV, index=False)

print("\n── Final dataset Candida albicans ───────────────────────────────────")
print(f"  Positive pairs (label=1): {(df_final['label'] == 1).sum():,}")
print(f"  Negative pairs (label=0): {(df_final['label'] == 0).sum():,}")
print(f"  Totale coppie           : {len(df_final):,}")
print(f"Salvato in: {OUTPUT_CSV.name}")
