import pandas as pd

# 1. Carica il file STRING (usa sep=' ' perché STRING di solito usa lo spazio)
print("Caricamento del database STRING originale...")
df_string = pd.read_csv('4932.protein.physical.links.detailed.v12.0.txt', sep=' ')

# 2. Pulisci gli ID di STRING (rimuovi il prefisso "4932.")
print("Pulizia degli ID...")
df_string['protein1'] = df_string['protein1'].str.replace('4932.', '', regex=False)
df_string['protein2'] = df_string['protein2'].str.replace('4932.', '', regex=False)

# 3. Metriche sulla colonna experimental prima del filtraggio
print("\n── Distribuzione colonna 'experimental' ──────────────────────────────")
total = len(df_string)
bins = {
    "= 0":       df_string['experimental'] == 0,
    "1 – 300":   (df_string['experimental'] >= 1)   & (df_string['experimental'] <= 300),
    "301 – 699": (df_string['experimental'] >= 301)  & (df_string['experimental'] <= 699),
    ">= 700":    df_string['experimental'] >= 700,
}
for label, mask in bins.items():
    count = mask.sum()
    print(f"  {label:<12s}  {count:>8,}  ({100 * count / total:.1f} %)")
print(f"  {'TOTALE':<12s}  {total:>8,}")
print()

# 4. Filtra solo le interazioni fisiche sperimentali ad ALTA CONFIDENZA (score >= 700)
print("Filtraggio per confidenza sperimentale (score >= 700)...")
df_positives = df_string[df_string['experimental'] >= 700].copy()

# 5. Rimuovi le colonne non necessarie
cols_to_drop = [c for c in ['database', 'textmining', 'combined_score'] if c in df_positives.columns]
if cols_to_drop:
    df_positives.drop(columns=cols_to_drop, inplace=True)
    print(f"Colonne rimosse: {cols_to_drop}")

# 6. SALVATAGGIO DEL NUOVO DB PULITO
output_filename = 'string_high_confidence.csv'
df_positives.to_csv(output_filename, index=False)

print(f"\nFinito! Interazioni positive finali: {len(df_positives)}")
print(f"Colonne nel file di output: {list(df_positives.columns)}")
print(f"Il tuo nuovo DB pulito è stato salvato nel file: '{output_filename}'")