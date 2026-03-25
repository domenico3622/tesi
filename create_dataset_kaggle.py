import pandas as pd
import re
from sklearn.model_selection import train_test_split

print("1. Lettura mappa dei Cluster da CD-HIT...")
cluster_map = {}
with open("fasta_lievito\proteins_cdhit.clstr", "r") as f:
    current_cluster = ""
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster = line.strip().replace(">", "")
        else:
            # Estrae l'ID (es. P00927) e gli assegna il cluster
            match = re.search(r'>(.*?)\.\.\.', line)
            if match:
                cluster_map[match.group(1)] = current_cluster

print(f"   Cluster letti               : {len(set(cluster_map.values())):,}")
print(f"   Proteine mappate da .clstr  : {len(cluster_map):,}")

print("2. Filtraggio del dataset per rimuovere la ridondanza...")
df = pd.read_csv("string_lievito\string_final_dataset.csv")

# --- Metriche iniziali dataset ---
n_interazioni_in = len(df)
pos_in = int((df["label"] == 1).sum())
neg_in = int((df["label"] == 0).sum())
proteine_in = set(df["protein1"]).union(set(df["protein2"]))

print("\n--- Metriche input (pre-cleaning cluster) ---")
print(f"Interazioni totali             : {n_interazioni_in:,}")
print(f"Interazioni positive           : {pos_in:,}")
print(f"Interazioni negative           : {neg_in:,}")
print(f"Proteine uniche coinvolte      : {len(proteine_in):,}")

# Copertura proteine nel file .clstr
proteine_mancanti_clstr = [p for p in proteine_in if p not in cluster_map]
print(f"Proteine NON presenti in .clstr: {len(proteine_mancanti_clstr):,}")
if proteine_mancanti_clstr:
    preview = ", ".join(sorted(proteine_mancanti_clstr)[:20])
    print(f"  Esempi (max 20): {preview}")

# Costruzione metrica su coppie cluster-cluster prima del filtro effettivo
df_metrics = df.copy()
df_metrics["cluster1"] = df_metrics["protein1"].map(lambda p: cluster_map.get(p, f"Missing_{p}"))
df_metrics["cluster2"] = df_metrics["protein2"].map(lambda p: cluster_map.get(p, f"Missing_{p}"))
df_metrics["cluster_pair"] = df_metrics.apply(
    lambda r: tuple(sorted([r["cluster1"], r["cluster2"]])), axis=1
)

cluster_pair_sizes = df_metrics.groupby("cluster_pair").size()
n_cluster_pairs_uniche = int(cluster_pair_sizes.shape[0])
n_interazioni_duplicate_clusterpair = int(n_interazioni_in - n_cluster_pairs_uniche)

label_conflicts = (
    df_metrics.groupby("cluster_pair")["label"].nunique() > 1
).sum()

print("\n--- Metriche ridondanza cluster-pair ---")
print(f"Cluster-pair uniche            : {n_cluster_pairs_uniche:,}")
print(f"Interazioni duplicate scartabili: {n_interazioni_duplicate_clusterpair:,}")
print(f"Cluster-pair con label conflittuali (0/1): {int(label_conflicts):,}")

visti = set()
righe_valide = []

for _, row in df.iterrows():
    p1 = row['protein1']
    p2 = row['protein2']
    
    # Troviamo i cluster. Se una proteina non c'è, le diamo un cluster fittizio unico
    c1 = cluster_map.get(p1, f"Missing_{p1}")
    c2 = cluster_map.get(p2, f"Missing_{p2}")
    
    # Ordiniamo i cluster per non avere duplicati speculari
    coppia_cluster = tuple(sorted([c1, c2]))
    
    # Se non abbiamo mai visto un'interazione tra queste due famiglie di proteine, la teniamo
    if coppia_cluster not in visti:
        visti.add(coppia_cluster)
        righe_valide.append(row)

df_filtrato = pd.DataFrame(righe_valide)

# --- Metriche dopo cleaning cluster ---
n_interazioni_out = len(df_filtrato)
pos_out = int((df_filtrato["label"] == 1).sum())
neg_out = int((df_filtrato["label"] == 0).sum())

proteine_out = set(df_filtrato["protein1"]).union(set(df_filtrato["protein2"]))
proteine_scartate = proteine_in - proteine_out

print("\n--- Metriche output (post-cleaning cluster) ---")
print(f"Interazioni mantenute          : {n_interazioni_out:,}")
print(f"Interazioni scartate           : {n_interazioni_in - n_interazioni_out:,}")
print(f"Retention interazioni          : {100 * n_interazioni_out / n_interazioni_in:.2f}%")
print(f"Positive mantenute             : {pos_out:,}  (da {pos_in:,})")
print(f"Negative mantenute             : {neg_out:,}  (da {neg_in:,})")
print(f"Proteine uniche mantenute      : {len(proteine_out):,}  (da {len(proteine_in):,})")
print(f"Proteine uniche scartate       : {len(proteine_scartate):,}")
if proteine_scartate:
    preview = ", ".join(sorted(proteine_scartate)[:20])
    print(f"  Esempi proteine scartate (max 20): {preview}")

print("\n--- Controllo coerenza filtro ---")
print(f"Cluster-pair viste nel filtro  : {len(visti):,}")
print(f"Interazioni finali             : {n_interazioni_out:,}")
print("(Devono coincidere: una riga per cluster-pair)")

print("3. Split 80/20 per il Training...")
df_train, df_test = train_test_split(
    df_filtrato, 
    test_size=0.20, 
    random_state=42, 
    stratify=df_filtrato['label']
)

df_train.to_csv("dscript_train.csv", index=False)
df_test.to_csv("dscript_test.csv", index=False)

print("\n--- Metriche split train/test ---")
print(f"Train size                     : {len(df_train):,}")
print(f"Test size                      : {len(df_test):,}")
print(f"Train positive                 : {(df_train['label'] == 1).sum():,}")
print(f"Train negative                 : {(df_train['label'] == 0).sum():,}")
print(f"Test positive                  : {(df_test['label'] == 1).sum():,}")
print(f"Test negative                  : {(df_test['label'] == 0).sum():,}")

print("\nFile pronti per Kaggle: 'dscript_train.csv' e 'dscript_test.csv' creati!")
