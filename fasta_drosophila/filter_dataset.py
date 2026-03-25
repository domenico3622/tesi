import random
import pandas as pd
from Bio import SeqIO

FASTA_INPUT = "fasta_drosophila\\drosophila.fasta"
TSV_INPUT = "fasta_drosophila\\drosophila.tsv"
TSV_OUTPUT = "fasta_drosophila\\drosophila_filter.tsv"
FASTA_OUTPUT = "fasta_drosophila\\drosophila_filter.fasta"
N_PROTEINE = 10000

# 1. Caricare tutte le proteine dal FASTA e campionarne 10000
records = list(SeqIO.parse(FASTA_INPUT, "fasta"))
print(f"Proteine nel FASTA originale: {len(records)}")

random.seed(42)
sampled = random.sample(records, N_PROTEINE)
sampled_ids = {rec.id for rec in sampled}
print(f"Proteine campionate: {len(sampled_ids)}")

# 2. Salvare il nuovo FASTA con 10000 proteine
with open(FASTA_OUTPUT, "w") as f:
    for rec in sampled:
        f.write(f">{rec.id}\n{str(rec.seq)}\n")

# 3. Caricare il TSV e filtrare solo interazioni dove entrambe le proteine sono nel FASTA
df = pd.read_csv(TSV_INPUT, sep='\t', header=None, names=['p1', 'p2', 'label'])
filtered_df = df[df['p1'].isin(sampled_ids) & df['p2'].isin(sampled_ids)]
filtered_df.to_csv(TSV_OUTPUT, sep='\t', index=False, header=False)

positivi = filtered_df[filtered_df['label'] == 1.0]
negativi = filtered_df[filtered_df['label'] == 0.0]

print(f"\nInterazioni nel TSV originale: {len(df)}")
print(f"Interazioni nel TSV filtrato:  {len(filtered_df)}")
print(f"  - Positive (1): {len(positivi)}")
print(f"  - Negative (0): {len(negativi)}")
