"""
Questo codice fa il preprocessing (pulizia) del file FASTA grezzo per renderlo 
compatibile con i passaggi successivi:
- Pulisce i nomi: Prende le intestazioni lunghe e complesse di UniProt 
 (es. >sp|P00927|THDH...) e le taglia per tenere solo la Chiave Primaria (es. >P00927).
- Applica i filtri: Conta la lunghezza di ogni sequenza. Salva solo 
  quelle comprese tra 80 e 400 amminoacidi (il vincolo imposto dagli autori).
- Output: Genera un file FASTA snello e pulito (cdhit_ready_yeast.fasta), perfetto 
  per essere letto dall'algoritmo successivo.

Cos'è CD-HIT e perché gli autori lo hanno usato:
Algoritmo bioinformatico che confronta migliaia di sequenze FASTA e raggruppa 
(clusterizza) quelle che si assomigliano molto.
L'Obiettivo degli autori (Rimozione della Ridondanza): Gli autori lo usano 
impostando una soglia di somiglianza del 40%.
Perché lo fanno: Se non usassero CD-HIT, il modello potrebbe imparare a memoria 
coppie di proteine facilissime e quasi identiche tra loro, imbrogliando durante 
il test. CD-HIT assicura che il dataset contenga solo famiglie di proteine molto
diverse tra loro, costringendo la rete neurale a imparare le vere regole fisiche 
delle interazioni e non a memorizzare i dati."""
from Bio import SeqIO

input_fasta = "candida_albicans\\candida_albicans.fasta"
output_fasta = "candida_albicans\\candida_albicans_filter.fasta"

proteine_valide = 0
proteine_scartate = 0

print(f"Lettura del file grezzo {input_fasta} in corso...")

# Apriamo il nuovo file in modalità scrittura
with open(output_fasta, "w") as out_file:
    # SeqIO.parse legge il FASTA blocco per blocco senza saturare la RAM
    for record in SeqIO.parse(input_fasta, "fasta"):
        
        # 1. Estraiamo l'ID puro
        # record.id di solito è "sp|P00927|THDH_YEAST". 
        # Lo splittiamo usando il carattere "|" e prendiamo il secondo elemento (indice 1)
        parti_id = record.id.split('|')
        
        if len(parti_id) >= 2:
            id_puro = parti_id[1] # Questo diventa "P00927"
        else:
            id_puro = record.id # Fallback nel caso il formato sia strano
            
        # 2. Misuriamo la lunghezza della sequenza (Filtro D-SCRIPT)
        lunghezza_seq = len(record.seq)
        
        if 80 <= lunghezza_seq <= 700:
            # 3. Scriviamo il FASTA pulito se passa il filtro
            out_file.write(f">{id_puro}\n{str(record.seq)}\n")
            proteine_valide += 1
        else:
            proteine_scartate += 1

print("\n--- RISULTATI ---")
print(f"Proteine valide salvate per CD-HIT (80-700 aa): {proteine_valide}")
print(f"Proteine scartate (troppo corte o troppo lunghe): {proteine_scartate}")
print(f"File pronto generato: '{output_fasta}'")


# IMPORTANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# una volta eseguito questo script, bisogna eseguire CD-HIT da terminale:
# !cd-hit -i cdhit_ready_yeast.fasta -o proteins_cdhit -c 0.4 -n 2 -M 0 -d 0

"""
1. proteins_cdhit.clstr
è la mappa dei cluster. È l'elenco di chi è finito in quale gruppo.
>Cluster 0: È il nome del primo gruppo creato.
0 500aa, >P32606... *: Ti dice che la proteina P32606 fa parte del Cluster 0.

Cosa significa l'asterisco (*)? Significa che questa proteina è il "Leader" 
di questo gruppo. Siccome nel tuo estratto c'è solo lei nel Cluster 0, vuol dire 
che CD-HIT non ha trovato nessun'altra proteina nel lievito che le assomigli per 
più del 40%. È unica nel suo genere.

Se ci fossero state proteine simili, avrei visto una seconda riga sotto il Cluster 0, 
senza asterisco, con scritto ad esempio at 45%, a indicare che un'altra proteina 
è stata "assorbita" da questo leader.

"""

"""
PERCHE COME PARAMENTRO DI CD-HIT HANNO SCELTO 0.4 (40%)?

Perché proprio il 40%? (Da dire in sede di tesi)

In bioinformatica strutturale, la soglia circa 40% è nota come "Twilight Zone"
dell'identità di sequenza.
Se due proteine hanno un'identità > 40%, è molto probabile che si ripieghino 
nella stessa forma 3D. È facile per un'IA "indovinare" le loro interazioni per 
analogia.
Sotto il 40%, le sequenze sembrano totalmente diverse a occhio nudo, anche se 
magari mantengono una struttura spaziale simile.
Impostando il filtro al 40%, gli autori di D-SCRIPT non possono più affidarsi a 
semplici somiglianze di sequenza per prevedere le interazioni, ma a imparare le 
vere regole chimico-fisiche del legame invece di imparare a memoria le sequenze 
simili.
"""