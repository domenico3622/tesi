"""
run_etl_pipeline.py
===================
Pipeline ETL unificata che esegue in sequenza tutti i 6 step:

  Step 1: create_protein_zones    → protein_zones.csv
  Step 2: create_unified_dataset  → unified_protein_filter_ds.csv
  Step 3: cleaning_string         → string_high_confidence.csv
  Step 4: final_string_pos_neg    → string_final_dataset.csv
  Step 5: cd-hit prep             → cdhit_ready_yeast.fasta  (+ comando CD-HIT manuale)
  Step 6: create_dataset_kaggle   → dscript_train.csv, dscript_test.csv

Uso:
  python run_etl_pipeline.py --output-dir ./output
  python run_etl_pipeline.py --go-obo path/go-basic.obo --uniprot-tsv path/uniprot_yeast.tsv ...
"""

import argparse
import os
import re
import random
import subprocess
import sys

# Forza l'encoding UTF-8 per evitare errori con le emoji su Windows
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding='utf-8')

import pandas as pd
from Bio import SeqIO
from goatools.obo_parser import GODag
from sklearn.model_selection import train_test_split

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich import print as rprint

console = Console()

random.seed(42)

# ═══════════════════════════════════════════════════════════════════════════════
#  ARGPARSE
# ═══════════════════════════════════════════════════════════════════════════════
def parse_args():
    description = "Pipeline ETL unificata per il progetto D-SCRIPT sul lievito."
    epilog = """
[bold cyan]Esempi di utilizzo:[/bold cyan]

  [yellow]# Esecuzione standard con tutti i default[/yellow]
  python run_etl_pipeline.py --output-dir ./output_test

  [yellow]# Esecuzione saltando gli step iniziali e definendo soglie personalizzate[/yellow]
  python run_etl_pipeline.py --skip-steps 1,2 --exp-threshold 900 --neg-ratio 5

  [yellow]# Esecuzione con binario CD-HIT specifico e avvio automatico del clustering[/yellow]
  python run_etl_pipeline.py --run-cdhit --cdhit-bin C:\\path\\to\\cd-hit.exe

  [yellow]# Modifica dei vincoli di lunghezza sequenza[/yellow]
  python run_etl_pipeline.py --min-len 50 --max-len 800
    """

    p = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False  # Disabilitiamo l'help automatico per gestirlo con Rich senza errori di encoding
    )
    # Aggiungiamo manualmente -h/--help
    p.add_argument("-h", "--help", action="store_true", help="Mostra questo messaggio di aiuto ed esce")
    # ── file di input ──
    p.add_argument("--go-obo",       default="uniprot_lievito/go-basic.obo",
                   help="File GO ontology (go-basic.obo)")
    p.add_argument("--uniprot-tsv",  default="uniprot_lievito/uniprot_yeast.tsv",
                   help="File TSV UniProt del lievito")
    p.add_argument("--fasta",        default="fasta_lievito/UP000002311_559292.fasta",
                   help="File FASTA proteoma del lievito")
    p.add_argument("--string-txt",   default="string_lievito/4932.protein.physical.links.detailed.v12.0.txt",
                   help="File STRING physical links")
    p.add_argument("--cdhit-clstr",  default="fasta_lievito/proteins_cdhit.clstr",
                   help="File .clstr generato da CD-HIT")

    # ── directory di output ──
    p.add_argument("--output-dir",   default="output_etl",
                   help="Directory dove salvare tutti gli output")

    # ── parametri configurabili ──
    p.add_argument("--exp-threshold",   type=int, default=700,
                   help="Soglia minima score sperimentale STRING (default: 700)")
    p.add_argument("--neg-ratio",       type=int, default=10,
                   help="Rapporto negativi:positivi (default: 10)")
    p.add_argument("--min-len",         type=int, default=80,
                   help="Lunghezza minima amminoacidi per filtraggio proteine (default: 80)")
    p.add_argument("--max-len",         type=int, default=400,
                   help="Lunghezza massima amminoacidi per filtraggio proteine (default: 400)")
    p.add_argument("--test-size",       type=float, default=0.20,
                   help="Frazione test split (default: 0.20)")
    p.add_argument("--run-cdhit",       action="store_true",
                   help="Esegui automaticamente il comando cd-hit dopo lo step 5")
    p.add_argument("--cdhit-threshold", type=float, default=0.4,
                   help="Soglia di similarita' per CD-HIT (default: 0.4)")
    p.add_argument("--cdhit-bin",       default="cd-hit",
                   help="Percorso dell'eseguibile cd-hit (default: cd-hit)")
    p.add_argument("--skip-steps",      type=str, default="",
                   help="Step da saltare, separati da virgola (es. '1,5')")
    p.add_argument("-i", "--interactive", action="store_true",
                   help="Avvia la pipeline in modalita' interattiva (chiede i parametri)")

    return p, p.parse_args()


def banner(step_num, title):
    console.print(f"\n")
    console.print(Panel(
        f"[bold cyan]STEP {step_num}: {title}[/bold cyan]",
        border_style="cyan",
        expand=False,
        padding=(0, 2)
    ))


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 1 — create_protein_zones
# ═══════════════════════════════════════════════════════════════════════════════
MACRO_ZONES = {
    "nucleus":               "GO:0005634",
    "cytoplasm":             "GO:0005737",
    "mitochondrion":         "GO:0005739",
    "plasma_membrane":       "GO:0005886",
    "endoplasmic_reticulum": "GO:0005783",
    "golgi_apparatus":       "GO:0005794",
    "vacuole":               "GO:0005773",
}

def get_all_ancestors(godag, go_id):
    if go_id not in godag:
        return set()
    term = godag[go_id]
    ancestors = set()
    queue = list(term.parents)
    while queue:
        parent = queue.pop()
        if parent.item_id not in ancestors:
            ancestors.add(parent.item_id)
            queue.extend(parent.parents)
    return ancestors

GO_PATTERN = re.compile(r"GO:\d{7}")

def step1_create_protein_zones(go_obo_path, uniprot_tsv_path, output_dir):
    banner(1, "CREATE PROTEIN ZONES")
    output_path = os.path.join(output_dir, "protein_zones.csv")

    print(f"Loading GO ontology da {go_obo_path} ...")
    godag = GODag(go_obo_path, optional_attrs={"relationship"})
    print(f"  Loaded {len(godag):,} GO terms.\n")

    print(f"Reading {uniprot_tsv_path} ...")
    df = pd.read_csv(uniprot_tsv_path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    print(f"  {len(df):,} proteins found.\n")

    zone_names  = list(MACRO_ZONES.keys())
    zone_go_ids = list(MACRO_ZONES.values())
    results = []

    for _, row in df.iterrows():
        cc_field = row.get("Gene Ontology (cellular component)", "")
        protein_go_ids = GO_PATTERN.findall(cc_field)
        hits = {z: 0 for z in zone_names}

        for go_id in protein_go_ids:
            if go_id not in godag:
                continue
            ancestors = get_all_ancestors(godag, go_id)
            ancestors.add(go_id)
            for zn, zgo in zip(zone_names, zone_go_ids):
                if zgo in ancestors:
                    hits[zn] = 1

        matched = [z for z in zone_names if hits[z] == 1]
        record = {
            "entry":      row.get("Entry", ""),
            "entry_name": row.get("Entry Name", ""),
            "gene_names": row.get("Gene Names", ""),
        }
        record.update(hits)
        record["zones"] = ",".join(matched)
        results.append(record)

    out_df = pd.DataFrame(results, columns=["entry","entry_name","gene_names"]+zone_names+["zones"])
    out_df.to_csv(output_path, index=False)
    print(f"Saved {output_path}  ({len(out_df):,} rows)\n")

    for z in zone_names:
        c = out_df[z].sum()
        print(f"  {z:<30s} {c:>5,} proteins  ({100*c/len(out_df):.1f} %)")
    no_zone = (out_df["zones"] == "").sum()
    multi   = (out_df["zones"].str.contains(",")).sum()
    print(f"\n  NO zone: {no_zone:,}  |  MULTIPLE zones: {multi:,}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 2 — create_unified_dataset
# ═══════════════════════════════════════════════════════════════════════════════
def step2_create_unified_dataset(uniprot_tsv, zones_csv, fasta_file, output_dir,
                                  min_len, max_len):
    banner(2, "CREATE UNIFIED DATASET")
    output_path = os.path.join(output_dir, "unified_protein_filter_ds.csv")

    # Parse FASTA
    print(f"Parsing FASTA {fasta_file} ...")
    fasta_dict, fasta_by_name = {}, {}
    header_re = re.compile(r"^>(\w+)\|([A-Za-z0-9_]+)\|(\w+)")
    current_id = current_name = None
    seq_parts = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and seq_parts:
                    seq = "".join(seq_parts)
                    fasta_dict[current_id] = seq
                    if current_name:
                        fasta_by_name[current_name] = seq
                m = header_re.match(line)
                if m:
                    current_id, current_name = m.group(2), m.group(3)
                else:
                    parts = line[1:].split("|")
                    current_id = parts[1] if len(parts) > 1 else parts[0].split()[0]
                    current_name = parts[2].split()[0] if len(parts) > 2 else None
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_id and seq_parts:
            seq = "".join(seq_parts)
            fasta_dict[current_id] = seq
            if current_name:
                fasta_by_name[current_name] = seq

    print(f"  {len(fasta_dict):,} sequences (by accession).\n")

    # Load TSV + zones
    df = pd.read_csv(uniprot_tsv, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    zones_df = pd.read_csv(zones_csv, dtype=str, keep_default_na=False)
    zone_cols = ["nucleus","cytoplasm","mitochondrion","plasma_membrane",
                 "endoplasmic_reticulum","golgi_apparatus","vacuole","zones"]
    zones_to_merge = zones_df[["entry"]+zone_cols].rename(columns={"entry":"Entry"})
    merged = df.merge(zones_to_merge, on="Entry", how="left")
    if "Gene Ontology (cellular component)" in merged.columns:
        merged.drop(columns=["Gene Ontology (cellular component)"], inplace=True)

    def lookup_seq(row):
        e, n = row.get("Entry",""), row.get("Entry Name","")
        return fasta_dict.get(e, fasta_by_name.get(n, ""))

    merged["sequence"] = merged.apply(lookup_seq, axis=1)
    found = (merged["sequence"] != "").sum()
    print(f"  Sequences found: {found:,} / {len(merged):,}")

    before = len(merged)
    merged["Length"] = pd.to_numeric(merged["Length"], errors="coerce")
    merged = merged[(merged["Length"] >= min_len) & (merged["Length"] <= max_len)].copy()
    print(f"  After length filter ({min_len}-{max_len}): {len(merged):,} / {before:,}\n")

    merged.to_csv(output_path, index=False)
    print(f"Saved → {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 3 — cleaning_string
# ═══════════════════════════════════════════════════════════════════════════════
def step3_cleaning_string(string_txt, output_dir, exp_threshold):
    banner(3, "CLEANING STRING")
    output_path = os.path.join(output_dir, "string_high_confidence.csv")

    print(f"Loading STRING da {string_txt} ...")
    df = pd.read_csv(string_txt, sep=' ')
    df['protein1'] = df['protein1'].str.replace('4932.', '', regex=False)
    df['protein2'] = df['protein2'].str.replace('4932.', '', regex=False)

    total = len(df)
    print(f"  Totale interazioni: {total:,}")
    print(f"  Filtraggio experimental >= {exp_threshold} ...")

    df_pos = df[df['experimental'] >= exp_threshold].copy()
    cols_drop = [c for c in ['database','textmining','combined_score'] if c in df_pos.columns]
    if cols_drop:
        df_pos.drop(columns=cols_drop, inplace=True)

    df_pos.to_csv(output_path, index=False)
    print(f"  Interazioni positive finali: {len(df_pos):,}")
    print(f"Saved → {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 4 — final_string_pos_neg
# ═══════════════════════════════════════════════════════════════════════════════
def step4_final_string_pos_neg(string_hc_csv, uniprot_tsv, filtered_ds, output_dir,
                                neg_ratio):
    banner(4, "FINAL STRING POS/NEG")
    output_path = os.path.join(output_dir, "string_final_dataset.csv")

    # Gene-name map
    uni = pd.read_csv(uniprot_tsv, sep="\t", dtype=str, keep_default_na=False)
    uni.columns = [c.strip() for c in uni.columns]
    gene_to_entry = {}
    for _, row in uni.iterrows():
        entry = row["Entry"]
        for token in row["Gene Names"].split():
            t = token.strip().upper()
            if t:
                gene_to_entry[t] = entry
    print(f"  Gene-name tokens: {len(gene_to_entry):,}")

    # Translate STRING
    df_str = pd.read_csv(string_hc_csv, dtype=str)
    df_str["entry1"] = df_str["protein1"].str.upper().map(gene_to_entry)
    df_str["entry2"] = df_str["protein2"].str.upper().map(gene_to_entry)
    before = len(df_str)
    df_str = df_str.dropna(subset=["entry1","entry2"])
    print(f"  Translated: {len(df_str):,} / {before:,}")

    # Filter valid proteins
    ds = pd.read_csv(filtered_ds, dtype=str, keep_default_na=False)
    valid_proteins = set(ds["Entry"].dropna())
    df_str = df_str[df_str["entry1"].isin(valid_proteins) & df_str["entry2"].isin(valid_proteins)].copy()

    pos_pairs_set = set(tuple(sorted([r.entry1, r.entry2])) for r in df_str.itertuples())
    df_pos = pd.DataFrame([(a,b,1) for a,b in pos_pairs_set], columns=["protein1","protein2","label"])
    n_pos = len(df_pos)
    print(f"  Unique positive pairs: {n_pos:,}")

    # Negatives
    n_neg = n_pos * neg_ratio
    print(f"  Generating {n_neg:,} negatives (ratio {neg_ratio}:1) ...")
    pool = list(valid_proteins)
    neg_pairs = set()
    attempts, max_att = 0, n_neg * 20
    while len(neg_pairs) < n_neg and attempts < max_att:
        a, b = random.sample(pool, 2)
        pair = tuple(sorted([a, b]))
        if pair not in pos_pairs_set and pair not in neg_pairs:
            neg_pairs.add(pair)
        attempts += 1

    df_neg = pd.DataFrame([(a,b,0) for a,b in neg_pairs], columns=["protein1","protein2","label"])
    df_final = pd.concat([df_pos, df_neg], ignore_index=True).sample(frac=1, random_state=42).reset_index(drop=True)
    df_final.to_csv(output_path, index=False)
    print(f"  Pos: {(df_final['label']==1).sum():,}  Neg: {(df_final['label']==0).sum():,}  Tot: {len(df_final):,}")
    print(f"Saved → {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 5 — cd-hit prep
# ═══════════════════════════════════════════════════════════════════════════════
def step5_cdhit_prep(fasta_file, output_dir, fasta_min, fasta_max,
                     run_cdhit, cdhit_threshold, cdhit_bin="cd-hit"):
    banner(5, "CD-HIT PREP")
    output_fasta = os.path.join(output_dir, "cdhit_ready_yeast.fasta")

    valid, skipped = 0, 0
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            parts = record.id.split('|')
            pure_id = parts[1] if len(parts) >= 2 else record.id
            length = len(record.seq)
            if fasta_min <= length <= fasta_max:
                out.write(f">{pure_id}\n{str(record.seq)}\n")
                valid += 1
            else:
                skipped += 1

    print(f"  Valide ({fasta_min}-{fasta_max} aa): {valid}  |  Scartate: {skipped}")
    print(f"Saved → {output_fasta}")

    cdhit_output = os.path.join(output_dir, "proteins_cdhit")
    if run_cdhit:
        # Use quotes around paths to handle spaces
        cmd = f'"{cdhit_bin}" -i "{output_fasta}" -o "{cdhit_output}" -c {cdhit_threshold} -n 2 -M 0 -d 0'
        print(f"\n  Esecuzione CD-HIT: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        if result.returncode != 0:
            raise RuntimeError(
                f"CD-HIT ha restituito codice {result.returncode}.\n"
                f"Verifica che il percorso del binario sia corretto: '{cdhit_bin}'\n"
                f"Stderr: {result.stderr.strip()}"
            )
        print("  CD-HIT completato!")
    else:
        print(f"\n  ⚠ Per completare lo step, esegui manualmente:")
        print(f"    \"{cdhit_bin}\" -i \"{output_fasta}\" -o \"{cdhit_output}\" -c {cdhit_threshold} -n 2 -M 0 -d 0")

    return output_fasta


# ═══════════════════════════════════════════════════════════════════════════════
#  STEP 6 — create_dataset_kaggle
# ═══════════════════════════════════════════════════════════════════════════════
def step6_create_dataset_kaggle(cdhit_clstr, string_final_csv, output_dir, test_size):
    banner(6, "CREATE DATASET KAGGLE")

    # Parse cluster map
    cluster_map = {}
    with open(cdhit_clstr, "r") as f:
        current_cluster = ""
        for line in f:
            if line.startswith(">Cluster"):
                current_cluster = line.strip().replace(">", "")
            else:
                match = re.search(r'>(.*?)\.\.\.', line)
                if match:
                    cluster_map[match.group(1)] = current_cluster

    print(f"  Cluster: {len(set(cluster_map.values())):,}  |  Proteine mappate: {len(cluster_map):,}")

    df = pd.read_csv(string_final_csv)
    n_in = len(df)

    # Cluster-pair deduplication
    seen = set()
    valid_rows = []
    for _, row in df.iterrows():
        c1 = cluster_map.get(row['protein1'], f"Missing_{row['protein1']}")
        c2 = cluster_map.get(row['protein2'], f"Missing_{row['protein2']}")
        cp = tuple(sorted([c1, c2]))
        if cp not in seen:
            seen.add(cp)
            valid_rows.append(row)

    df_filt = pd.DataFrame(valid_rows)
    print(f"  Post-cleaning: {len(df_filt):,} / {n_in:,}")

    # Split
    df_train, df_test = train_test_split(
        df_filt, test_size=test_size, random_state=42, stratify=df_filt['label']
    )
    train_path = os.path.join(output_dir, "dscript_train.csv")
    test_path  = os.path.join(output_dir, "dscript_test.csv")
    df_train.to_csv(train_path, index=False)
    df_test.to_csv(test_path, index=False)

    print(f"  Train: {len(df_train):,} (pos={int((df_train['label']==1).sum()):,})")
    print(f"  Test:  {len(df_test):,} (pos={int((df_test['label']==1).sum()):,})")
    print(f"Saved → {train_path}, {test_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    parser, args = parse_args()
    
    # Se l'utente chiede help, usiamo Rich per renderizzare i tag nell'epilog
    if args.help:
        console.print(Panel(parser.description, style="bold blue", expand=False))
        import io
        f = io.StringIO()
        parser.print_help(f)
        console.print(f.getvalue())
        sys.exit(0)

    if args.interactive:
        from rich.prompt import Prompt, IntPrompt, FloatPrompt
        rprint(Panel("[bold]MODALITÀ INTERATTIVA[/bold]\nPremi Invio per mantenere i valori di default.", border_style="yellow"))
        
        args.output_dir = Prompt.ask("Directory di output", default=args.output_dir)
        args.exp_threshold = IntPrompt.ask("Soglia sperimentale STRING", default=args.exp_threshold)
        args.neg_ratio = IntPrompt.ask("Rapporto negativi:positivi", default=args.neg_ratio)
        args.min_len = IntPrompt.ask("Lunghezza minima proteina", default=args.min_len)
        args.max_len = IntPrompt.ask("Lunghezza massima proteina", default=args.max_len)
        
        if Prompt.ask("Vuoi saltare qualche step?", choices=["y", "n"], default="n") == "y":
            args.skip_steps = Prompt.ask("Inserisci step da saltare (es. 1,2,5)", default="")
            
        if Prompt.ask("Avviare CD-HIT automaticamente?", choices=["y", "n"], default="y" if args.run_cdhit else "n") == "y":
            args.run_cdhit = True
            args.cdhit_bin = Prompt.ask("Percorso binario CD-HIT", default=args.cdhit_bin)

    skip = set(args.skip_steps.split(",")) if args.skip_steps else set()

    os.makedirs(args.output_dir, exist_ok=True)
    
    rprint(Panel(
        f"[bold blue]D-SCRIPT Yeast ETL Pipeline[/bold blue]\n"
        f"Output Directory: [green]{args.output_dir}[/green]\n"
        f"Config: STRING Thresh=[yellow]{args.exp_threshold}[/yellow], "
        f"Length=[yellow]{args.min_len}-{args.max_len}[/yellow]",
        title="[bold]🚀 Avvio Pipeline[/bold]",
        border_style="blue"
    ))

    # Paths per gli output intermedi
    zones_csv   = os.path.join(args.output_dir, "protein_zones.csv")
    unified_csv = os.path.join(args.output_dir, "unified_protein_filter_ds.csv")
    string_hc   = os.path.join(args.output_dir, "string_high_confidence.csv")
    string_fin  = os.path.join(args.output_dir, "string_final_dataset.csv")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        
        # STEP 1
        if "1" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 1: Protein Zones...", total=None)
            step1_create_protein_zones(args.go_obo, args.uniprot_tsv, args.output_dir)
        else:
            rprint("\n[yellow]⏭ Step 1 saltato.[/yellow]")

        # STEP 2
        if "2" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 2: Unified Dataset...", total=None)
            step2_create_unified_dataset(
                args.uniprot_tsv, zones_csv, args.fasta, args.output_dir,
                args.min_len, args.max_len
            )
        else:
            rprint("\n[yellow]⏭ Step 2 saltato.[/yellow]")

        # STEP 3
        if "3" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 3: Cleaning STRING...", total=None)
            step3_cleaning_string(args.string_txt, args.output_dir, args.exp_threshold)
        else:
            rprint("\n[yellow]⏭ Step 3 saltato.[/yellow]")

        # STEP 4
        if "4" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 4: Final STRING Pos/Neg...", total=None)
            step4_final_string_pos_neg(
                string_hc, args.uniprot_tsv, unified_csv, args.output_dir, args.neg_ratio
            )
        else:
            rprint("\n[yellow]⏭ Step 4 saltato.[/yellow]")

        # STEP 5
        if "5" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 5: CD-HIT Prep...", total=None)
            step5_cdhit_prep(
                args.fasta, args.output_dir,
                args.min_len, args.max_len,
                args.run_cdhit, args.cdhit_threshold, args.cdhit_bin
            )
        else:
            rprint("\n[yellow]⏭ Step 5 saltato.[/yellow]")

        # STEP 6
        if "6" not in skip:
            progress.add_task(description="[cyan]Esecuzione Step 6: Create Dataset Kaggle...", total=None)
            step6_create_dataset_kaggle(
                args.cdhit_clstr, string_fin, args.output_dir, args.test_size
            )
        else:
            rprint("\n[yellow]⏭ Step 6 saltato.[/yellow]")

    rprint("\n")
    rprint(Panel(
        "[bold green]✅ PIPELINE ETL COMPLETATA CON SUCCESSO![/bold green]\n"
        "Tutti i dataset finali sono pronti per il training.",
        border_style="green"
    ))


if __name__ == "__main__":
    main()
