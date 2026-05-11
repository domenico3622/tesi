"""
etl_ui.py - UI grafica per la pipeline ETL D-SCRIPT (lievito)
Lancia con: python etl_ui.py
"""
import os, sys, threading, tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

# Fix per UI sfocata su Windows (High DPI awareness)
if sys.platform == "win32":
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

# Import step functions from the CLI pipeline
from run_etl_pipeline import (
    step1_create_protein_zones, step2_create_unified_dataset,
    step3_cleaning_string, step4_final_string_pos_neg,
    step5_cdhit_prep, step6_create_dataset_kaggle,
)

# ── Costanti UI ──────────────────────────────────────────────────────────────
BG       = "#1e1e2e"
BG2      = "#2a2a3d"
FG       = "#cdd6f4"
ACCENT   = "#89b4fa"
ACCENT2  = "#a6e3a1"
WARN_CLR = "#fab387"
ERR_CLR  = "#f38ba8"
FONT     = ("Segoe UI", 10)
FONT_B   = ("Segoe UI", 10, "bold")
FONT_H   = ("Segoe UI", 14, "bold")
FONT_S   = ("Segoe UI", 9)

FILE_EXAMPLES = {
    "go_obo": (
        "GO Ontology (go-basic.obo)",
        "File OBO dell'ontologia Gene Ontology.\n"
        "Struttura:\n"
        "  [Term]\n"
        "  id: GO:0005634\n"
        "  name: nucleus\n"
        "  namespace: cellular_component\n"
        "  is_a: GO:0043231 ! intracellular ...",
        False,
    ),
    "uniprot_tsv": (
        "UniProt TSV (uniprot_yeast.tsv)",
        "File TSV scaricato da UniProt con le proteine del lievito.\n"
        "Colonne: Entry | Entry Name | Gene Names | Length | Gene Ontology (...) | ...\n"
        "Esempio riga:\n"
        "  P00927  THDH_YEAST  PDC1 ...  350  nucleus [GO:0005634]; ...",
        False,
    ),
    "fasta": (
        "Proteoma FASTA (.fasta)",
        "File FASTA del proteoma del lievito (es. UP000002311_559292.fasta).\n"
        "Struttura:\n"
        "  >sp|P00927|THDH_YEAST Pyruvate decarboxylase ...\n"
        "  MSEIERLNANDTLVGMPAGAKQLHQALTSGVSAIPNAAENSYQYVAT...\n"
        "  >sp|P00445|SODC_YEAST Superoxide dismutase ...\n"
        "  MVQAVAVLK...",
        False,
    ),
    "string_txt": (
        "STRING Physical Links (.txt)",
        "File STRING con le interazioni fisiche.\n"
        "Colonne separate da spazio:\n"
        "  protein1 protein2 experimental database textmining combined_score\n"
        "  4932.Q0045 4932.YGR222W 292 0 0 292\n"
        "  4932.Q0045 4932.YML030W 734 0 0 734",
        False,
    ),
    "cdhit_clstr": (
        "CD-HIT Clusters (.clstr)  ⚡ OPZIONALE",
        "File .clstr generato dal comando cd-hit.\n"
        "Se NON viene fornito, lo script eseguira' cd-hit\n"
        "automaticamente (cd-hit deve essere installato).\n\n"
        "Struttura:\n"
        "  >Cluster 0\n"
        "  0  500aa, >P32606... *\n"
        "  >Cluster 1\n"
        "  0  398aa, >P00927... *\n"
        "  1  312aa, >Q12345... at 45%",
        True,
    ),
    "cdhit_bin": (
        "Eseguibile CD-HIT  ⚡ OPZIONALE",
        "Percorso completo dell'eseguibile cd-hit.\n"
        "Lascia vuoto se cd-hit e' gia' nel PATH di sistema.\n\n"
        "Esempi:\n"
        "  Windows: C:\\cd-hit\\cd-hit.exe\n"
        "  Linux/Mac: /usr/local/bin/cd-hit\n\n"
        "Se non hai cd-hit installato, scaricalo da:\n"
        "  https://github.com/weizhongli/cdhit/releases",
        True,
    ),
}

PARAM_INFO = {
    "exp_thresh": (
        "Soglia Sperimentale STRING",
        "Score minimo della colonna 'experimental' nel file STRING.\n"
        "Serve a filtrare solo le interazioni con prove sperimentali solide.\n"
        "Valore consigliato: >= 700 (Alta confidenza)."
    ),
    "neg_ratio": (
        "Rapporto Negativi:Positivi",
        "Il numero di coppie di proteine non interagenti (label=0) da generare\n"
        "per ogni coppia interagente (label=1).\n"
        "Esempio: 10 significa che per ogni positivo ci sono 10 negativi."
    ),
    "min_prot_len": (
        "Lunghezza Minima Proteina",
        "Numero minimo di amminoacidi che una proteina deve avere.\n"
        "Proteine piu' corte vengono scartate dal dataset.\n"
        "Default: 80 (vincolo imposto dagli autori D-SCRIPT)."
    ),
    "max_prot_len": (
        "Lunghezza Massima Proteina",
        "Numero massimo di amminoacidi che una proteina puo' avere.\n"
        "Proteine piu' lunghe vengono scartate dal dataset.\n"
        "Default: 400 (vincolo imposto dagli autori D-SCRIPT)."
    ),
    "test_size": (
        "Frazione Dataset di Test",
        "Frazione del dataset finale da riservare per il test (valutazione).\n"
        "Il valore deve essere compreso tra 0 e 1.\n"
        "Esempio: 0.20 significa 20% test e 80% training."
    ),
    "cdhit_thresh": (
        "Soglia Similarità CD-HIT",
        "Soglia di identita' di sequenza usata da CD-HIT per il clustering.\n"
        "0.4 significa 40% di identita' (Twilight Zone).\n"
        "Sopra il 40% le proteine sono troppo simili per un test equo."
    )
}

DEFAULTS = {
    "go_obo":      "uniprot_lievito/go-basic.obo",
    "uniprot_tsv": "uniprot_lievito/uniprot_yeast.tsv",
    "fasta":       "fasta_lievito/UP000002311_559292.fasta",
    "string_txt":  "string_lievito/4932.protein.physical.links.detailed.v12.0.txt",
    "cdhit_clstr": "",
    "cdhit_bin":   "",
}


class ETLApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("D-SCRIPT Yeast ETL Pipeline")
        self.configure(bg=BG)
        self.geometry("950x1100")
        self.minsize(900, 900)
        self.running = False

        # ── Variabili ──
        self.file_vars = {k: tk.StringVar(value=v) for k, v in DEFAULTS.items()}
        self.output_dir = tk.StringVar(value="output_etl")
        self.exp_thresh = tk.IntVar(value=700)
        self.neg_ratio  = tk.IntVar(value=10)
        self.min_prot_len = tk.IntVar(value=80)
        self.max_prot_len = tk.IntVar(value=400)
        self.test_size  = tk.DoubleVar(value=0.20)
        self.cdhit_thresh = tk.DoubleVar(value=0.4)

        self._build_ui()

    # ── Costruzione UI ───────────────────────────────────────────────────────
    def _build_ui(self):
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TFrame", background=BG)
        style.configure("Card.TFrame", background=BG2)
        style.configure("TLabel", background=BG, foreground=FG, font=FONT)
        style.configure("Card.TLabel", background=BG2, foreground=FG, font=FONT)
        style.configure("H.TLabel", background=BG, foreground=ACCENT, font=FONT_H)
        style.configure("Opt.TLabel", background=BG2, foreground=WARN_CLR, font=FONT_S)
        style.configure("TButton", font=FONT_B)
        style.configure("Run.TButton", font=("Segoe UI", 12, "bold"))
        style.configure("Small.TButton", font=("Segoe UI", 8), padding=1)
        style.configure("Info.TButton", font=("Segoe UI", 8))

        # Scrollable canvas
        canvas = tk.Canvas(self, bg=BG, highlightthickness=0)
        vsb = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        self.main = ttk.Frame(canvas, style="TFrame")
        canvas.create_window((0, 0), window=self.main, anchor="nw")
        self.main.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(-1*(e.delta//120), "units"))

        pad = {"padx": 14, "pady": 4}

        # ── Header ──
        ttk.Label(self.main, text="🧬  D-SCRIPT Yeast — ETL Pipeline", style="H.TLabel").pack(padx=14, pady=(14,2))
        ttk.Label(self.main, text="Seleziona i file di input, configura i parametri e avvia la pipeline.",
                  style="TLabel").pack(padx=14, pady=4)

        # ── Input Files ──
        self._section("📂  File di Input", [
            ("go_obo",      "GO Ontology (.obo)"),
            ("uniprot_tsv", "UniProt TSV"),
            ("fasta",       "Proteoma FASTA"),
            ("string_txt",  "STRING Physical Links"),
        ])

        # ── Optional Files ──
        self._section("📂  File Opzionali", [
            ("cdhit_clstr", "CD-HIT Clusters (.clstr)"),
            ("cdhit_bin",   "Eseguibile CD-HIT"),
        ], optional=True)

        # ── Output Dir ──
        fr_out = ttk.LabelFrame(self.main, text="📁  Directory di Output", style="Card.TFrame")
        fr_out.pack(fill="x", padx=14, pady=(8,4))
        fr_out.configure(labelwidget=self._colored_label(fr_out, "📁  Directory di Output"))
        row = ttk.Frame(fr_out, style="Card.TFrame")
        row.pack(fill="x", padx=10, pady=6)
        ttk.Entry(row, textvariable=self.output_dir, font=FONT, width=30).pack(side="left")
        ttk.Button(row, text="📁", width=3, style="Small.TButton",
                   command=lambda: self._browse_dir(self.output_dir)).pack(side="left", padx=(6,0))

        # ── Info Panel ──
        fr_info = ttk.LabelFrame(self.main, text="📖  Dettagli File", style="Card.TFrame")
        fr_info.pack(fill="x", padx=14, pady=(8,4))
        fr_info.configure(labelwidget=self._colored_label(fr_info, "📖  Dettagli File"))

        self.info_box = scrolledtext.ScrolledText(fr_info, height=7, bg=BG, fg=FG,
                                                 font=("Consolas", 9), relief="flat",
                                                 highlightthickness=0, state="disabled")
        self.info_box.pack(fill="x", padx=10, pady=8)
        self.info_box.tag_configure("header", foreground=ACCENT, font=("Consolas", 9, "bold"))
        self.info_box.tag_configure("opt", foreground=WARN_CLR, font=("Consolas", 9, "bold"))

        self._update_info_box("Seleziona '?' per visualizzare i dettagli del file e la sua struttura.")

        # ── Parameters ──
        self._build_params()

        # ── Run Button ──
        self.run_btn = ttk.Button(self.main, text="▶  AVVIA PIPELINE", style="Run.TButton",
                                   command=self._on_run)
        self.run_btn.pack(pady=12, ipady=6, ipadx=20)

        # ── Progress ──
        self.progress = ttk.Progressbar(self.main, mode="determinate", maximum=6)
        self.progress.pack(fill="x", padx=14, pady=(0,6))

        # ── Log ──
        self.log = scrolledtext.ScrolledText(self.main, height=14, bg="#11111b", fg=FG,
                                              font=("Consolas", 9), insertbackground=FG,
                                              state="disabled", wrap="word")
        self.log.pack(fill="both", expand=True, padx=14, pady=(0,14))
        self.log.tag_configure("ok", foreground=ACCENT2)
        self.log.tag_configure("warn", foreground=WARN_CLR)
        self.log.tag_configure("err", foreground=ERR_CLR)

    def _colored_label(self, parent, text):
        lbl = tk.Label(parent, text=text, bg=BG2, fg=ACCENT, font=FONT_B)
        return lbl

    def _section(self, title, items, optional=False):
        lf = ttk.LabelFrame(self.main, style="Card.TFrame")
        lf.pack(fill="x", padx=14, pady=(8,4))
        lf.configure(labelwidget=self._colored_label(lf, title))

        if optional:
            note = ttk.Label(lf, style="Opt.TLabel",
                text="⚡ Se non inseriti, verranno generati automaticamente (cd-hit deve essere installato).")
            note.pack(anchor="w", padx=10, pady=(4,0))

        for key, label in items:
            self._file_row(lf, key, label)

    def _file_row(self, parent, key, label):
        row = ttk.Frame(parent, style="Card.TFrame")
        row.pack(fill="x", padx=10, pady=4)

        ttk.Label(row, text=label, style="Card.TLabel", width=24, anchor="w").pack(side="left")
        ttk.Entry(row, textvariable=self.file_vars[key], font=FONT, width=30).pack(side="left", padx=(4,0))
        ttk.Button(row, text="📁", width=2, style="Small.TButton",
                   command=lambda k=key: self._browse_file(k)).pack(side="left", padx=(4,0))
        ttk.Button(row, text="?", width=1, style="Small.TButton",
                   command=lambda k=key: self._show_info(k)).pack(side="left", padx=(2,0))

    def _build_params(self):
        lf = ttk.LabelFrame(self.main, style="Card.TFrame")
        lf.pack(fill="x", padx=14, pady=(8,4))
        lf.configure(labelwidget=self._colored_label(lf, "⚙  Parametri Configurabili"))

        ttk.Label(lf, style="Opt.TLabel",
            text="Valori di default gia' impostati. Modifica solo se necessario."
        ).pack(anchor="w", padx=10, pady=(4,0))

        params = [
            ("Soglia sperimentale STRING",  self.exp_thresh,   50,  999, "exp_thresh"),
            ("Rapporto neg:pos",            self.neg_ratio,     1,   50, "neg_ratio"),
            ("Min lunghezza proteina (aa)", self.min_prot_len,  10,  500, "min_prot_len"),
            ("Max lunghezza proteina (aa)", self.max_prot_len, 100, 5000, "max_prot_len"),
        ]
        for label, var, lo, hi, key in params:
            self._param_row(lf, label, var, lo, hi, key)

        # Float params
        fr = ttk.Frame(lf, style="Card.TFrame")
        fr.pack(fill="x", padx=10, pady=4)
        ttk.Label(fr, text="Test split fraction", style="Card.TLabel", width=28, anchor="w").pack(side="left")
        ttk.Entry(fr, textvariable=self.test_size, font=FONT, width=8).pack(side="left")
        ttk.Button(fr, text="?", width=1, style="Small.TButton",
                   command=lambda: self._show_param_info("test_size")).pack(side="left", padx=6)

        fr2 = ttk.Frame(lf, style="Card.TFrame")
        fr2.pack(fill="x", padx=10, pady=(4,8))
        ttk.Label(fr2, text="CD-HIT similarity threshold", style="Card.TLabel", width=28, anchor="w").pack(side="left")
        ttk.Entry(fr2, textvariable=self.cdhit_thresh, font=FONT, width=8).pack(side="left")
        ttk.Button(fr2, text="?", width=1, style="Small.TButton",
                   command=lambda: self._show_param_info("cdhit_thresh")).pack(side="left", padx=6)

    def _param_row(self, parent, label, var, lo, hi, key):
        fr = ttk.Frame(parent, style="Card.TFrame")
        fr.pack(fill="x", padx=10, pady=3)
        ttk.Label(fr, text=label, style="Card.TLabel", width=28, anchor="w").pack(side="left")
        sb = ttk.Spinbox(fr, from_=lo, to=hi, textvariable=var, width=8, font=FONT)
        sb.pack(side="left")
        ttk.Button(fr, text="?", width=1, style="Small.TButton",
                   command=lambda k=key: self._show_param_info(k)).pack(side="left", padx=6)

    # ── Dialogs ──────────────────────────────────────────────────────────────
    def _browse_file(self, key):
        path = filedialog.askopenfilename(title=f"Seleziona {key}")
        if path:
            self.file_vars[key].set(path)

    def _browse_dir(self, var):
        path = filedialog.askdirectory(title="Seleziona directory output")
        if path:
            var.set(path)

    def _show_info(self, key):
        title, desc, is_opt = FILE_EXAMPLES[key]
        header = f"--- {title.upper()} ---\n"
        opt_str = "⚡ FILE OPZIONALE\n" if is_opt else ""
        
        self._update_info_box(f"{header}{opt_str}\n{desc}")

    def _show_param_info(self, key):
        title, desc = PARAM_INFO[key]
        header = f"--- PARAMETRO: {title.upper()} ---\n"
        self._update_info_box(f"{header}\n{desc}")

    def _update_info_box(self, text):
        self.info_box.configure(state="normal")
        self.info_box.delete("1.0", "end")
        self.info_box.insert("end", text)
        
        # Color first line and optional warning
        if text.startswith("---"):
            self.info_box.tag_add("header", "1.0", "1.end")
        if "⚡" in text:
            self.info_box.tag_add("opt", "2.0", "2.end")
            
        self.info_box.configure(state="disabled")

    # ── Log helpers ──────────────────────────────────────────────────────────
    def _log(self, msg, tag=None):
        self.log.configure(state="normal")
        self.log.insert("end", msg + "\n", tag)
        self.log.see("end")
        self.log.configure(state="disabled")
        self.update_idletasks()

    # ── Redirect stdout ──────────────────────────────────────────────────────
    def _redirect_print(self):
        """Replace sys.stdout so print() from step functions appears in the log."""
        app = self
        class LogWriter:
            def write(self, msg):
                if msg.strip():
                    app._log(msg.rstrip())
            def flush(self):
                pass
        self._old_stdout = sys.stdout
        sys.stdout = LogWriter()

    def _restore_print(self):
        sys.stdout = self._old_stdout

    # ── Run pipeline ─────────────────────────────────────────────────────────
    def _on_run(self):
        if self.running:
            return

        # Validate required files
        required = ["go_obo", "uniprot_tsv", "fasta", "string_txt"]
        for k in required:
            v = self.file_vars[k].get().strip()
            if not v or not os.path.exists(v):
                messagebox.showerror("File mancante", f"Il file '{k}' non esiste:\n{v}")
                return

        self.running = True
        self.run_btn.configure(state="disabled")
        self.progress["value"] = 0
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")

        threading.Thread(target=self._run_pipeline, daemon=True).start()

    def _run_pipeline(self):
        self._redirect_print()
        out = self.output_dir.get().strip()
        os.makedirs(out, exist_ok=True)

        zones_csv   = os.path.join(out, "protein_zones.csv")
        unified_csv = os.path.join(out, "unified_protein_filter_ds.csv")
        string_hc   = os.path.join(out, "string_high_confidence.csv")
        string_fin  = os.path.join(out, "string_final_dataset.csv")
        cdhit_clstr = self.file_vars["cdhit_clstr"].get().strip()
        has_clstr   = cdhit_clstr and os.path.exists(cdhit_clstr)

        try:
            # Step 1
            self._log("\n═══ STEP 1/6: Create Protein Zones ═══", "ok")
            step1_create_protein_zones(
                self.file_vars["go_obo"].get(), self.file_vars["uniprot_tsv"].get(), out)
            self.progress["value"] = 1

            # Step 2
            self._log("\n═══ STEP 2/6: Create Unified Dataset ═══", "ok")
            step2_create_unified_dataset(
                self.file_vars["uniprot_tsv"].get(), zones_csv,
                self.file_vars["fasta"].get(), out,
                self.min_prot_len.get(), self.max_prot_len.get())
            self.progress["value"] = 2

            # Step 3
            self._log("\n═══ STEP 3/6: Cleaning STRING ═══", "ok")
            step3_cleaning_string(
                self.file_vars["string_txt"].get(), out, self.exp_thresh.get())
            self.progress["value"] = 3

            # Step 4
            self._log("\n═══ STEP 4/6: Final STRING Pos/Neg ═══", "ok")
            step4_final_string_pos_neg(
                string_hc, self.file_vars["uniprot_tsv"].get(), unified_csv,
                out, self.neg_ratio.get())
            self.progress["value"] = 4

            # Step 5
            self._log("\n═══ STEP 5/6: CD-HIT Prep ═══", "ok")
            run_cdhit = not has_clstr
            cdhit_bin = self.file_vars["cdhit_bin"].get().strip() or "cd-hit"
            if run_cdhit:
                self._log(f"  .clstr non fornito → cd-hit verra' eseguito automaticamente.", "warn")
                self._log(f"  Binario: {cdhit_bin}", "warn")
            step5_cdhit_prep(
                self.file_vars["fasta"].get(), out,
                self.min_prot_len.get(), self.max_prot_len.get(),
                run_cdhit, self.cdhit_thresh.get(), cdhit_bin)
            self.progress["value"] = 5

            # Determine clstr path for step 6
            if has_clstr:
                final_clstr = cdhit_clstr
            else:
                final_clstr = os.path.join(out, "proteins_cdhit.clstr")

            if not os.path.exists(final_clstr):
                self._log(f"\n⚠ File .clstr non trovato: {final_clstr}", "err")
                self._log("  CD-HIT potrebbe non essere installato. Step 6 saltato.", "err")
                raise FileNotFoundError(f".clstr non trovato: {final_clstr}")

            # Step 6
            self._log("\n═══ STEP 6/6: Create Dataset Kaggle ═══", "ok")
            step6_create_dataset_kaggle(
                final_clstr, string_fin, out, self.test_size.get())
            self.progress["value"] = 6

            self._log("\n✅  PIPELINE COMPLETATA CON SUCCESSO!", "ok")

        except Exception as e:
            self._log(f"\n❌ ERRORE: {e}", "err")
            import traceback
            self._log(traceback.format_exc(), "err")
        finally:
            self._restore_print()
            self.running = False
            self.run_btn.configure(state="normal")


if __name__ == "__main__":
    app = ETLApp()
    app.mainloop()
