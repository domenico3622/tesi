"""
create_pombe_training_dataset.py
--------------------------------
Pipeline completa per Schizosaccharomyces pombe:
1) Filtra il FASTA e tiene solo proteine con lunghezza 50-700 aa
2) Pulisce il file STRING physical links
3) Crea un CSV etichettato (label=1 se experimental >= 700, altrimenti 0)
4) Tiene solo coppie con entrambe le proteine nel FASTA filtrato
5) Costruisce il dataset finale con rapporto negativi:positivi = 10:1

Output:
- pombe_filtered_50_700.fasta
- pombe_filtered_proteins_50_700.csv
- pombe_pairs_filtered_labeled.csv
- pombe_training_10to1.csv
"""

from pathlib import Path
import random

import pandas as pd
from Bio import SeqIO


MIN_LEN = 50
MAX_LEN = 700
SCORE_THRESHOLD = 700
NEG_POS_RATIO = 10
RANDOM_SEED = 42


BASE_DIR = Path(__file__).resolve().parent
FASTA_IN = BASE_DIR / "pombe.fasta"
STRING_TXT = BASE_DIR / "284812.protein.physical.links.detailed.v12.0.txt"

FASTA_OUT = BASE_DIR / "pombe_filtered_50_700.fasta"
PROTEIN_CSV_OUT = BASE_DIR / "pombe_filtered_proteins_50_700.csv"
LABELED_PAIRS_OUT = BASE_DIR / "pombe_pairs_filtered_labeled.csv"
TRAIN_OUT = BASE_DIR / "pombe_training_10to1.csv"


random.seed(RANDOM_SEED)


def extract_accession(record_id: str) -> str:
    parts = record_id.split("|")
    if len(parts) >= 2:
        return parts[1].strip()
    return record_id.strip()


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def filter_fasta_by_length() -> dict[str, int]:
    print(f"Lettura FASTA: {FASTA_IN}")

    kept_sequences: dict[str, str] = {}
    kept_lengths: dict[str, int] = {}

    total = 0
    discarded = 0

    for rec in SeqIO.parse(str(FASTA_IN), "fasta"):
        total += 1
        acc = extract_accession(rec.id)
        seq = str(rec.seq)
        length = len(seq)

        if MIN_LEN <= length <= MAX_LEN:
            if acc not in kept_sequences:
                kept_sequences[acc] = seq
                kept_lengths[acc] = length
        else:
            discarded += 1

    with FASTA_OUT.open("w", encoding="utf-8") as f:
        for acc, seq in kept_sequences.items():
            f.write(f">{acc}\n{seq}\n")

    protein_df = pd.DataFrame(
        [(acc, length) for acc, length in kept_lengths.items()],
        columns=["protein", "length"],
    ).sort_values("protein")
    protein_df.to_csv(PROTEIN_CSV_OUT, index=False)

    print("\n--- FASTA filtrato ---")
    print(f"Proteine lette           : {total:,}")
    print(f"Proteine mantenute       : {len(kept_sequences):,}")
    print(f"Proteine scartate        : {discarded:,}")
    print(f"Salvato FASTA filtrato   : {FASTA_OUT.name}")
    print(f"Salvato CSV lunghezze    : {PROTEIN_CSV_OUT.name}\n")

    return kept_lengths


def load_and_label_string(valid_proteins: set[str]) -> pd.DataFrame:
    print(f"Lettura STRING: {STRING_TXT}")

    df = pd.read_csv(STRING_TXT, sep=r"\s+", dtype={"experimental": int})

    df["protein1"] = df["protein1"].astype(str).str.split(".", n=1).str[-1]
    df["protein2"] = df["protein2"].astype(str).str.split(".", n=1).str[-1]

    df = df[df["protein1"] != df["protein2"]].copy()

    mask_valid = df["protein1"].isin(valid_proteins) & df["protein2"].isin(valid_proteins)
    df = df[mask_valid].copy()

    p1 = df[["protein1", "protein2"]].min(axis=1)
    p2 = df[["protein1", "protein2"]].max(axis=1)
    df["protein1"] = p1
    df["protein2"] = p2

    pair_max = (
        df.groupby(["protein1", "protein2"], as_index=False)["experimental"]
        .max()
        .rename(columns={"experimental": "experimental_max"})
    )

    pair_max["label"] = (pair_max["experimental_max"] >= SCORE_THRESHOLD).astype(int)

    pair_max.to_csv(LABELED_PAIRS_OUT, index=False)

    n_pos = int((pair_max["label"] == 1).sum())
    n_neg = int((pair_max["label"] == 0).sum())

    print("--- Coppie filtrate ed etichettate ---")
    print(f"Coppie uniche totali      : {len(pair_max):,}")
    print(f"Positivi (>= {SCORE_THRESHOLD})     : {n_pos:,}")
    print(f"Negativi (< {SCORE_THRESHOLD})      : {n_neg:,}")
    print(f"Salvato CSV etichettato   : {LABELED_PAIRS_OUT.name}\n")

    return pair_max


def build_train_10_to_1(labeled_pairs: pd.DataFrame, valid_proteins: list[str]) -> pd.DataFrame:
    pos_df = labeled_pairs[labeled_pairs["label"] == 1][["protein1", "protein2"]].copy()
    neg_df = labeled_pairs[labeled_pairs["label"] == 0][["protein1", "protein2"]].copy()

    n_pos = len(pos_df)
    if n_pos == 0:
        raise RuntimeError("Nessuna coppia positiva trovata con la soglia impostata.")

    target_neg = n_pos * NEG_POS_RATIO

    pos_set = {canonical_pair(a, b) for a, b in pos_df.itertuples(index=False)}
    neg_set = {canonical_pair(a, b) for a, b in neg_df.itertuples(index=False)}

    if len(neg_set) < target_neg:
        print(
            f"Negativi da STRING insufficienti ({len(neg_set):,} < {target_neg:,}). "
            "Genero negativi aggiuntivi casuali."
        )

        attempts = 0
        max_attempts = target_neg * 30

        while len(neg_set) < target_neg and attempts < max_attempts:
            a, b = random.sample(valid_proteins, 2)
            pair = canonical_pair(a, b)
            if pair not in pos_set and pair not in neg_set:
                neg_set.add(pair)
            attempts += 1

    neg_list = list(neg_set)
    if len(neg_list) >= target_neg:
        neg_sampled = random.sample(neg_list, target_neg)
    else:
        neg_sampled = neg_list

    pos_final = pd.DataFrame(sorted(pos_set), columns=["protein1", "protein2"])
    pos_final["label"] = 1

    neg_final = pd.DataFrame(neg_sampled, columns=["protein1", "protein2"])
    neg_final["label"] = 0

    final_df = pd.concat([pos_final, neg_final], ignore_index=True)
    final_df = final_df.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)

    final_df.to_csv(TRAIN_OUT, index=False)

    print("--- Dataset training finale ---")
    print(f"Positivi finali           : {(final_df['label'] == 1).sum():,}")
    print(f"Negativi finali           : {(final_df['label'] == 0).sum():,}")
    print(f"Rapporto neg:pos          : {NEG_POS_RATIO}:1 (target)")
    print(f"Totale righe              : {len(final_df):,}")
    print(f"Salvato training CSV      : {TRAIN_OUT.name}\n")

    return final_df


def main() -> None:
    if not FASTA_IN.exists():
        raise FileNotFoundError(f"File FASTA non trovato: {FASTA_IN}")
    if not STRING_TXT.exists():
        raise FileNotFoundError(f"File STRING non trovato: {STRING_TXT}")

    kept_lengths = filter_fasta_by_length()
    valid_proteins = set(kept_lengths.keys())

    labeled_pairs = load_and_label_string(valid_proteins)
    build_train_10_to_1(labeled_pairs, sorted(valid_proteins))

    print("Pipeline completata con successo.")


if __name__ == "__main__":
    main()
