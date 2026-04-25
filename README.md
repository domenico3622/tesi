# Tesi: Predizione di Interazioni Proteina-Proteina con ESM-2 e Modello Custom

## Panoramica
Questo repository contiene il lavoro di tesi sulla predizione delle interazioni proteina-proteina (PPI), con confronto tra:

- baseline D-SCRIPT pre-addestrata
- modello custom basato su embedding ESM-2 + modulo di cross-attention bidirezionale + classificazione logistica

Il flusso sperimentale principale e' documentato nel notebook:

- [appunti/esm2-v2 (4).ipynb](appunti/esm2-v2%20(4).ipynb)

L'idea centrale e' migliorare la capacita' di generalizzazione cross-specie sfruttando embedding proteici contestuali (ESM-2) e una contact map predetta con attention multi-head.

## Obiettivo della tesi
Obiettivo del progetto e' costruire e valutare una pipeline completa che:

1. prepara dataset PPI bilanciati/sbilanciati in modo controllato
2. genera embedding proteici da sequenza con ESM-2
3. addestra un modello custom per inferire la probabilita' di interazione
4. confronta le prestazioni con D-SCRIPT su piu' organismi

## Architettura del modello
La figura seguente riassume l'architettura implementata nel notebook.

![Diagramma generale architettura](appunti/diagramma_generale.png)

### Lettura del diagramma
- Ogni proteina viene codificata da ESM-2 in una sequenza di vettori (embedding per residuo).
- I due tensori Z1 e Z2 passano in blocchi Conv1D dilatati paralleli (dilatazioni 1, 2, 4) per catturare pattern locali su scale diverse.
- Le feature passano in due MLP (q_net e k_net) e vengono splittate in multi-head.
- Viene applicata una cross-attention bidirezionale tra le due proteine.
- Le mappe vengono fuse e attivate per ottenere una Predicted Contact Map.
- Un modulo di pooling adattivo + attivazione logistica produce la probabilita' finale di interazione.

Nel notebook sono presenti due varianti principali:

- CustomAttentionPPI (matrice di affinita' spaziale statica)
- CustomAttentionPPI_Auto (matrice di affinita' 7x7 apprendibile)

## Pipeline ETL (Extract, Transform, Load)
La fase ETL e' distribuita in piu' script e prepara i dati per training e valutazione.

### 1) Extract
Sorgenti principali:

- STRING physical links (lievito, pombe, candida)
- UniProt TSV e FASTA
- dataset TSV/CSV organism-specific (es. drosophila)

File/sorgenti coinvolte nella codebase:

- [string_lievito/4932.protein.physical.links.detailed.v12.0.txt](string_lievito/4932.protein.physical.links.detailed.v12.0.txt)
- [fasta_pombe/284812.protein.physical.links.detailed.v12.0.txt](fasta_pombe/284812.protein.physical.links.detailed.v12.0.txt)
- [candida_albicans/237561.protein.physical.links.detailed.v12.0.txt](candida_albicans/237561.protein.physical.links.detailed.v12.0.txt)
- [uniprot_lievito/uniprot_yeast.tsv](uniprot_lievito/uniprot_yeast.tsv)
- [fasta_lievito/UP000002311_559292.fasta](fasta_lievito/UP000002311_559292.fasta)

### 2) Transform
Trasformazioni principali:

- pulizia ID STRING (rimozione prefisso tassonomico)
- selezione positivi ad alta confidenza (experimental >= 700)
- costruzione negativi con campionamento casuale controllato (tipicamente rapporto 10:1)
- filtraggio sequenze per lunghezza
- deduplicazione/riduzione ridondanza tramite CD-HIT e cluster map
- mapping gene/protein ID verso accession UniProt
- arricchimento con localizzazione subcellulare (7 macro-zone GO)

Script principali:

- [string_lievito/cleaning_string.py](string_lievito/cleaning_string.py)
- [string_lievito/final_string_pos_neg.py](string_lievito/final_string_pos_neg.py)
- [create_dataset_kaggle.py](create_dataset_kaggle.py)
- [uniprot_lievito/create_protein_zones.py](uniprot_lievito/create_protein_zones.py)
- [uniprot_lievito/create_unified_dataset.py](uniprot_lievito/create_unified_dataset.py)
- [fasta_lievito/cd-hit.py](fasta_lievito/cd-hit.py)
- [fasta_pombe/create_pombe_training_dataset.py](fasta_pombe/create_pombe_training_dataset.py)
- [candida_albicans/create_candida_test_dataset.py](candida_albicans/create_candida_test_dataset.py)
- [fasta_drosophila/filter_dataset.py](fasta_drosophila/filter_dataset.py)

Output ETL tipici:

- dataset train/test per il fine-tuning: [dscript_train.csv](dscript_train.csv), [dscript_test.csv](dscript_test.csv)
- dataset organism-specific di test
- FASTA filtrati
- file unificato con zone GO e sequenze: [uniprot_lievito/unified_protein_filter_ds.csv](uniprot_lievito/unified_protein_filter_ds.csv)

### 3) Load
Il caricamento avviene nel notebook e negli script di training:

- conversione sequenze in embedding ESM-2 e salvataggio H5
- caricamento embeddings in RAM con dataset PyTorch custom
- costruzione DataLoader con padding dinamico e maschere di lunghezza

Nel notebook il dataset class e':

- YeastPPIDataset

La funzione di collate crea:

- tensori padded per entrambe le proteine
- maschere booleane per ignorare il padding
- tensori di localizzazione subcellulare

## Flusso sperimentale nel notebook
Nel notebook [appunti/esm2-v2 (4).ipynb](appunti/esm2-v2%20(4).ipynb) il workflow e' organizzato in blocchi:

1. installazione dipendenze
2. generazione embedding ESM-2 in H5
3. definizione dataset e dataloader
4. definizione architetture CustomAttentionPPI e CustomAttentionPPI_Auto
5. training multi-configurazione (diversi pesi della loss)
6. valutazione su baseline D-SCRIPT e modelli custom
7. confronto su organismi multipli (S. cerevisiae, Drosophila, S. pombe, Candida albicans)

## Loss e informazione biologica
La loss combina contributi multipli:

- BCE per classificazione interazione/non-interazione
- termine sulla magnitudine/struttura della contact map
- penalizzazione spaziale basata su compatibilita' di localizzazione subcellulare

Nella variante Auto, la matrice di affinita' tra compartimenti cellulari e' apprendibile.
Nella variante Spatial, la matrice e' statica (hard-coded).

## Dataset e organismi
Il repository contiene dati e script per:

- Saccharomyces cerevisiae (pipeline principale)
- Drosophila
- Schizosaccharomyces pombe
- Candida albicans

Cartelle principali:

- [fasta_lievito](fasta_lievito)
- [string_lievito](string_lievito)
- [uniprot_lievito](uniprot_lievito)
- [fasta_drosophila](fasta_drosophila)
- [fasta_pombe](fasta_pombe)
- [candida_albicans](candida_albicans)
- [grafici](grafici)
- [appunti](appunti)

## Valutazione e grafici
La valutazione confronta baseline e modello custom con metriche:

- AUPR
- AUROC
- Precision
- Recall
- tempo di inferenza/testing

Script di supporto ai grafici:

- [grafici/compare_custom_vs_dscript.py](grafici/compare_custom_vs_dscript.py)

## Requisiti software
Dipendenze usate nel notebook (Kaggle/Colab-like):

- torch
- transformers
- h5py
- pandas
- numpy
- biopython
- tqdm
- scikit-learn
- dscript

## Come leggere il repository
Percorso consigliato:

1. apri il notebook principale [appunti/esm2-v2 (4).ipynb](appunti/esm2-v2%20(4).ipynb)
2. leggi la pipeline ETL nei moduli di preprocessing in [string_lievito](string_lievito), [uniprot_lievito](uniprot_lievito), [fasta_lievito](fasta_lievito)
3. verifica i dataset finali [dscript_train.csv](dscript_train.csv) e [dscript_test.csv](dscript_test.csv)
4. consulta i grafici finali in [grafici](grafici)

## Nota finale
Questo progetto implementa una pipeline end-to-end per PPI prediction che integra:

- conoscenza biologica (localizzazione cellulare)
- rappresentazioni proteiche da language model (ESM-2)
- modulo neurale custom con cross-attention bidirezionale

con validazione comparativa multi-organismo contro baseline D-SCRIPT.