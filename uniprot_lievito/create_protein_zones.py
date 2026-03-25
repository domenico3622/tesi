"""
create_protein_zones.py
-----------------------
Reads uniprot_yeast.tsv, maps each protein's GO cellular-component terms to
the 7 macro-zones defined below (using ancestor traversal via goatools), and
writes protein_zones.csv.

Macro-zones
-----------
nucleus            GO:0005634
cytoplasm          GO:0005737
mitochondrion      GO:0005739
plasma_membrane    GO:0005886
endoplasmic_reticulum GO:0005783
golgi_apparatus    GO:0005794
vacuole            GO:0005773

Output columns
--------------
entry              UniProt accession
entry_name
gene_names
nucleus            1/0
cytoplasm          1/0
mitochondrion      1/0
plasma_membrane    1/0
endoplasmic_reticulum 1/0
golgi_apparatus    1/0
vacuole            1/0
zones              comma-separated list of matched zone names (empty = unknown)
"""

import re
import pandas as pd
from goatools.obo_parser import GODag

# ── 1. Macro-zone definitions ────────────────────────────────────────────────
MACRO_ZONES = {
    "nucleus":                 "GO:0005634",
    "cytoplasm":               "GO:0005737",
    "mitochondrion":           "GO:0005739",
    "plasma_membrane":         "GO:0005886",
    "endoplasmic_reticulum":   "GO:0005783",
    "golgi_apparatus":         "GO:0005794",
    "vacuole":                 "GO:0005773",
}

# ── 2. Load GO DAG ───────────────────────────────────────────────────────────
print("Loading GO ontology from go-basic.obo …")
godag = GODag("go-basic.obo", optional_attrs={"relationship"})
print(f"  Loaded {len(godag):,} GO terms.\n")

# Pre-build ancestor sets for each macro-zone so we only compute once.
# A GO term T belongs to zone Z if:
#   T == Z  OR  Z is in ancestors(T)
def get_all_ancestors(go_id: str) -> set:
    """Return the set of all ancestor GO IDs for a given term (not including itself)."""
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

# ── 3. Parse the UniProt TSV ─────────────────────────────────────────────────
print("Reading uniprot_yeast.tsv …")
df = pd.read_csv(
    "uniprot_yeast.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False,
)
df.columns = [c.strip() for c in df.columns]
print(f"  {len(df):,} proteins found.\n")

# ── 4. Extract GO term IDs from the cellular-component column ────────────────
GO_PATTERN = re.compile(r"GO:\d{7}")

def extract_go_ids(cell_value: str) -> list:
    """Return list of GO IDs found in a free-text annotation cell."""
    return GO_PATTERN.findall(cell_value)

# ── 5. Map each protein to macro-zones ──────────────────────────────────────
zone_names   = list(MACRO_ZONES.keys())
zone_go_ids  = list(MACRO_ZONES.values())

# Cache: for each protein GO term, store which zones it hits
# We'll compute: for each protein GO term T,
#   ancestors_of_T = get_all_ancestors(T) | {T}
#   then check if zone_id in ancestors_of_T   (T is under or equal to zone)

results = []

for _, row in df.iterrows():
    cc_field = row.get("Gene Ontology (cellular component)", "")
    protein_go_ids = extract_go_ids(cc_field)

    # For each protein GO term, pre-compute its ancestor set (include self)
    hits = {z: 0 for z in zone_names}

    for go_id in protein_go_ids:
        if go_id not in godag:
            continue
        # ancestors of this term + the term itself
        ancestors = get_all_ancestors(go_id)
        ancestors.add(go_id)

        for zone_name, zone_go in zip(zone_names, zone_go_ids):
            if zone_go in ancestors:
                hits[zone_name] = 1

    matched_zones = [z for z in zone_names if hits[z] == 1]
    zones_str = ",".join(matched_zones)

    record = {
        "entry":      row.get("Entry",      ""),
        "entry_name": row.get("Entry Name", ""),
        "gene_names": row.get("Gene Names", ""),
    }
    record.update(hits)
    record["zones"] = zones_str
    results.append(record)

# ── 6. Build and save the output DataFrame ───────────────────────────────────
out_df = pd.DataFrame(results, columns=["entry", "entry_name", "gene_names"] + zone_names + ["zones"])

output_path = "protein_zones.csv"
out_df.to_csv(output_path, index=False)
print(f"Saved {output_path}  ({len(out_df):,} rows)\n")

# ── 7. Quick summary ─────────────────────────────────────────────────────────
print("── Zone distribution ──────────────────────────────────────────────────")
for z in zone_names:
    count = out_df[z].sum()
    print(f"  {z:<30s} {count:>5,} proteins  ({100*count/len(out_df):.1f} %)")

no_zone = (out_df["zones"] == "").sum()
multi   = (out_df["zones"].str.contains(",")).sum()
print(f"\n  Proteins with NO matched zone : {no_zone:,}")
print(f"  Proteins in MULTIPLE zones    : {multi:,}")
print("\nDone.")


"""   
Output stats (6,743 proteins):
Zone	Count	%
nucleus	1,985	29.4%
cytoplasm	2,113	31.3%
mitochondrion	1,150	17.1%
plasma_membrane	554	8.2%
ER	595	8.8%
Golgi	148	2.2%
vacuole	314	4.7%
no zone matched	2,080	30.8%
multi-zone	1,860	27.6%
"""