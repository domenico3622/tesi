# Questo script è pensato per essere trasformato in un Notebook su Kaggle.
# SU KAGGLE: pip install dscript torch pandas numpy tqdm scikit-learn

import os
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
from tqdm import tqdm
import dscript

# --- PARAMETRI KAGGLE ---
# Assicurati di cambiare questi path su Kaggle (es: /kaggle/input/...)
DIR_DATA = "./" 
TRAIN_CSV = os.path.join(DIR_DATA, "dscript_train.csv")
TEST_CSV = os.path.join(DIR_DATA, "dscript_test.csv")
FASTA_FILE = os.path.join(DIR_DATA, "fasta_lievito", "cdhit_ready_yeast.fasta")
UNIPROT_FILTER_CSV = os.path.join(DIR_DATA, "uniprot_lievito", "unified_protein_filter_ds.csv")

# Parametri Iper-addestramento
BATCH_SIZE = 16
EPOCHS = 5
LEARNING_RATE = 1e-4
PENALTY_LAMBDA = 0.5  # Peso della penalità zone
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print(f"Device in uso: {DEVICE}")

# --- 1. PREPARAZIONE DATI E "ZONE" ---

print("Caricamento mappa Zone...")
df_zones = pd.read_csv(UNIPROT_FILTER_CSV)
# Creiamo un dizionario Entry -> Zone
# Sostituiamo i NaN con stringa vuota per evitare errori
df_zones['zones'] = df_zones['zones'].fillna("")
zone_map = pd.Series(df_zones.zones.values, index=df_zones.Entry).to_dict()

# Caricamento Fasta
print("Caricamento sequenze FASTA...")
# dscript.fasta.parse_fasta restituisce una tupla (name, sequence) o un dizionario
# se preferisci, usiamo biopython per semplicità, ma dscript ha le sue utilities.
with open(FASTA_FILE, "r") as f:
    fasta_dict = {}
    curr_name = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # estrae solo l'id (es: P00927) rimuovendo la descrizione
            curr_name = line[1:].split()[0]
            fasta_dict[curr_name] = ""
        else:
            fasta_dict[curr_name] += line

# ---- 2. CARICAMENTO MODELLO D-SCRIPT ---
print("Caricamento modello pre-addestrato D-SCRIPT Human...")
# Instanzia D-SCRIPT: Il modo corretto è caricare i pre-trained weights da un file se scaricati,
# Kaggle potrebbe averlo. In dscript puoi usare i modelli built-in se esisti, 
# altrimenti devi caricarli dal percorso locale se scaricati.
# Assumiamo che possiedi "human_v1.pt" caricato nel dataset kaggle.
# path_modello_human = "/kaggle/input/dscript-models/human_v1.pt"
# Per ora, come placeholder simuliamo l'uso. La libreria D-SCRIPT lo inizializza dal file .pt
try:
    # Se il file pt esiste nel PATH
    # model = torch.load(path_modello_human)
    # model.eval()
    pass
except Exception as e:
    print("Modello non caricato in questo script mockup, su Kaggle assicurati di caricare il .pt corretto.")

# Nota: per importare un modello D-SCRIPT vero, la reference è:
# model = dscript.models.interaction.Model.load("path_al_modello.pt")

# --- 3. CUSTOM DATASET PYTORCH ---
class YeastPPIDataset(Dataset):
    def __init__(self, tsv_file, fasta_dict, vocab):
        self.data = pd.read_csv(tsv_file)
        self.fasta_dict = fasta_dict
        self.vocab = vocab  # dizionario aminoacidi fornito da dscript
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        row = self.data.iloc[idx]
        p1_name = row['protein1']
        p2_name = row['protein2']
        label = float(row['label'])
        
        # Recupera sequenze
        seq1 = self.fasta_dict.get(p1_name, "")
        seq2 = self.fasta_dict.get(p2_name, "")
        
        # In un'applicazione reale, dovresti filtrare le sequenze mancanti prima
        # e tokenizzare la sequenza in numeri basandoti sul Vocabolario di D-SCRIPT.
        return p1_name, p2_name, seq1, seq2, label

# Per un vero dataloader:
# vocab = dscript.alphabets.Uniprot21() # O simile
# train_dataset = YeastPPIDataset(TRAIN_CSV, fasta_dict, vocab=None)
# train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)


# --- 4. TRAINING LOOP CON LOSS CUSTOM ---
def fine_tune_dscript(model, train_loader, zone_map, epochs=EPOCHS, lambda_penalty=PENALTY_LAMBDA):
    # Freezing the Language Model layers and Projection modules
    # In D-SCRIPT, l'architettura tipicamente è composta da un 'embedding' (LM) 
    # e un modulo 'contact' e 'interaction'.
    # Dipende se la versione corrente di dscript su pt consente iterazioni sui layer
    for name, param in model.named_parameters():
        if "embedding" in name or "projection" in name:
            param.requires_grad = False
        else:
            # Scongela solo Contact Map e Interaction Map
            param.requires_grad = True

    model.to(DEVICE)
    model.train()
    
    optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=LEARNING_RATE)
    criterion_bce = nn.BCEWithLogitsLoss()
    
    for epoch in range(epochs):
        epoch_loss = 0.0
        
        # tqdm per la barra di progresso su Kaggle
        for batch in tqdm(train_loader, desc=f"Epoch {epoch+1}/{epochs}"):
            p1_names, p2_names, seqs1, seqs2, labels = batch
            
            # Qui andrebbero convertite seqs in tensori in base a come gradisce il modello,
            # lo ommettiamo nel boilerplate
            # tensor_seq1, tensor_seq2 = tokenize_and_pad(seqs1), tokenize_and_pad(seqs2)
            # labels = labels.to(DEVICE)
            
            optimizer.zero_grad()
            
            # 1. Forward pass del Modello
            # prob_interazione, ... = model(tensor_seq1, tensor_seq2)
            # LOGITS! D-Script potrebbe uscire dal livello finale (senza sigmoid) per BCEWithLogits
            # preds_logits = prob_interazione.squeeze()
            
            # --- SIMULAZIONE CALCOLO LOSS ---
            # bce_loss = criterion_bce(preds_logits, labels)
            
            # 2. Calcolo Penalty "Zone"
            penalty = 0.0
            # Convertiamo logs in probabilità per ponderare la loss
            # pred_probs = torch.sigmoid(preds_logits)
            
            for i in range(len(p1_names)):
                z1 = zone_map.get(p1_names[i], "")
                z2 = zone_map.get(p2_names[i], "")
                
                # Applica penalità SE le zone sono note E le proteine sono in zone diverse (e non sono multiple)
                # (nota: 'zones' nel dataset a volte contiene stringhe combinate, occorre splittare o verificare la differenza netta)
                if z1 != "" and z2 != "" and z1 != z2:
                    # Penalizziamo proporzionalmente a quanto il modello crede che interagiscano
                    # penalty += pred_probs[i] 
                    pass
            
            # penalty_tensor = (penalty / len(p1_names)) * lambda_penalty
            # total_loss = bce_loss + penalty_tensor
            
            # Backward pass
            # total_loss.backward()
            optimizer.step()
            
            # epoch_loss += total_loss.item()
            
        print(f"Epoch {epoch+1} completed. Loss media (Simulata): ???")
        
    print("Fine-tuning terminato.")
    return model

# ESEMPIO CHIAMATA:
# model_finetuned = fine_tune_dscript(model, train_loader, zone_map)
# torch.save(model_finetuned.state_dict(), "dscript_yeast_finetuned.pt")
print("Script template generato! Trasformami in un blocco di codice per Jupyter Kaggle.")
