import pandas as pd

TSV_INPUT = "fasta_drosophila/drosophila_filter.tsv"
CSV_OUTPUT = "fasta_drosophila/drosophila_filter.csv"

# Carica il TSV (senza header)
df = pd.read_csv(TSV_INPUT, sep='\t', header=None, names=['protein1', 'protein2', 'label'])

# Converti label da float a int
df['label'] = df['label'].astype(int)

# Salva in CSV con header
df.to_csv(CSV_OUTPUT, index=False)

print(f"File salvato: {CSV_OUTPUT}")
print(f"Righe totali: {len(df)}")
print(f"  - Positive (1): {(df['label'] == 1).sum()}")
print(f"  - Negative (0): {(df['label'] == 0).sum()}")
print(f"\nPrime 5 righe:")
print(df.head())
