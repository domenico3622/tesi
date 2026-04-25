import matplotlib.pyplot as plt
import numpy as np


# Dati: da ora in poi usiamo solo D-SCRIPT vs modello custom (ex config_2)
organism_data = {
    "S. cerevisiae": {
        "n_coppie": 43473,
        "dscript": {
            "AUPR": 0.4051,
            "AUROC": 0.7812,
            "Precision": 0.6925,
            "Recall": 0.2665,
            "Tempo_s": 300,
        },
        "custom": {
            "AUPR": 0.9227,
            "AUROC": 0.9884,
            "Precision": 0.8976,
            "Recall": 0.8311,
            "Tempo_s": 67,
        },
    },
    "Drosophila": {
        "n_coppie": 14899,
        "dscript": {
            "AUPR": 0.5631,
            "AUROC": 0.8260,
            "Precision": 0.7886,
            "Recall": 0.3618,
            "Tempo_s": 116,
        },
        "custom": {
            "AUPR": 0.6300,
            "AUROC": 0.8773,
            "Precision": 0.8591,
            "Recall": 0.3133,
            "Tempo_s": 24,
        },
    },
    "Candida albicans": {
        "n_coppie": None,
        "dscript": {
            "AUPR": 0.5998,
            "AUROC": 0.8611,
            "Precision": 0.8106,
            "Recall": 0.5193,
            "Tempo_s": 335,
        },
        "custom": {
            "AUPR": 0.8693,
            "AUROC": 0.9617,
            "Precision": 0.9074,
            "Recall": 0.7595,
            "Tempo_s": 68,
        },
    },
}

metrics = ["AUPR", "AUROC", "Precision", "Recall"]
organisms = list(organism_data.keys())


# ========== GRAFICO 1: Differenza metriche (custom - dscript) per organismo ==========
fig1, ax1 = plt.subplots(figsize=(13, 6))

x = np.arange(len(metrics))
width = 0.24
colors = ["#4ECDC4", "#45B7D1", "#FFA07A"]

for i, organism in enumerate(organisms):
    diffs = [
        organism_data[organism]["custom"][metric] - organism_data[organism]["dscript"][metric]
        for metric in metrics
    ]
    bars = ax1.bar(
        x + (i - 1) * width,
        diffs,
        width,
        label=organism,
        color=colors[i],
        alpha=0.9,
    )
    for bar in bars:
        value = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            value + (0.006 if value >= 0 else -0.02),
            f"{value:+.3f}",
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=9,
        )

ax1.axhline(0, color="black", linewidth=1)
ax1.set_xlabel("Metriche", fontsize=12, fontweight="bold")
ax1.set_ylabel("Differenza (custom - dscript)", fontsize=12, fontweight="bold")
ax1.set_title("Differenza metriche tra organismi: modello custom vs dscript", fontsize=14, fontweight="bold")
ax1.set_xticks(x)
ax1.set_xticklabels(metrics)
ax1.legend()
ax1.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("grafici/grafico1_diff_metriche_organismi.png", dpi=300, bbox_inches="tight")
print("✓ Grafico 1 salvato: grafico1_diff_metriche_organismi.png")


# ========== GRAFICO 2: Istogramma tempi (dscript vs custom) per organismo ==========
fig2, ax2 = plt.subplots(figsize=(12, 6))

x_org = np.arange(len(organisms))
width_time = 0.35

dscript_times = [organism_data[org]["dscript"]["Tempo_s"] for org in organisms]
custom_times = [organism_data[org]["custom"]["Tempo_s"] for org in organisms]

bars_dscript = ax2.bar(
    x_org - width_time / 2,
    dscript_times,
    width_time,
    label="dscript",
    color="#FF6B6B",
    alpha=0.85,
)
bars_custom = ax2.bar(
    x_org + width_time / 2,
    custom_times,
    width_time,
    label="modello custom",
    color="#4ECDC4",
    alpha=0.85,
)

for bar in list(bars_dscript) + list(bars_custom):
    height = bar.get_height()
    ax2.text(
        bar.get_x() + bar.get_width() / 2,
        height,
        f"{int(height)}s",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

ax2.set_ylabel("Tempo testing (secondi)", fontsize=12, fontweight="bold")
ax2.set_xlabel("Organismo", fontsize=12, fontweight="bold")
ax2.set_title("Confronto velocità: dscript vs modello custom", fontsize=14, fontweight="bold")
ax2.set_xticks(x_org)
ax2.set_xticklabels(organisms)
ax2.legend()
ax2.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("grafici/grafico2_tempo_testing_organismi.png", dpi=300, bbox_inches="tight")
print("✓ Grafico 2 salvato: grafico2_tempo_testing_organismi.png")


# ========== GRAFICO 3: Confronto metriche dscript vs custom per ogni organismo ==========
fig3, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for idx, organism in enumerate(organisms):
    ax = axes[idx]
    x_m = np.arange(len(metrics))
    width_m = 0.35

    dscript_vals = [organism_data[organism]["dscript"][metric] for metric in metrics]
    custom_vals = [organism_data[organism]["custom"][metric] for metric in metrics]

    ax.bar(
        x_m - width_m / 2,
        dscript_vals,
        width_m,
        label="dscript",
        color="#FF6B6B",
        alpha=0.85,
    )
    ax.bar(
        x_m + width_m / 2,
        custom_vals,
        width_m,
        label="modello custom",
        color="#4ECDC4",
        alpha=0.85,
    )

    ax.set_title(organism, fontsize=12, fontweight="bold")
    ax.set_xticks(x_m)
    ax.set_xticklabels(metrics, rotation=20)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim([0, 1.05])

axes[0].set_ylabel("Valore metrica", fontsize=12, fontweight="bold")
axes[1].legend(loc="upper left")
fig3.suptitle("Confronto metriche: dscript vs modello custom (per organismo)", fontsize=15, fontweight="bold")

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("grafici/grafico3_metriche_per_organismo.png", dpi=300, bbox_inches="tight")
print("✓ Grafico 3 salvato: grafico3_metriche_per_organismo.png")


# ========== Sintesi testuale ==========
print("\n" + "=" * 72)
print("📊 SINTESI MIGLIORAMENTI (modello custom vs dscript)")
print("=" * 72)

for organism in organisms:
    print(f"\n{organism}")
    for metric in metrics:
        dscript_val = organism_data[organism]["dscript"][metric]
        custom_val = organism_data[organism]["custom"][metric]
        delta_pct = ((custom_val - dscript_val) / dscript_val) * 100
        print(
            f"  {metric:9} | dscript: {dscript_val:.4f} → custom: {custom_val:.4f} "
            f"| Δ {delta_pct:+.1f}%"
        )

    d_time = organism_data[organism]["dscript"]["Tempo_s"]
    c_time = organism_data[organism]["custom"]["Tempo_s"]
    speedup = d_time / c_time
    print(f"  Tempo     | dscript: {d_time}s → custom: {c_time}s | speedup {speedup:.2f}x")

print("\nGrafici generati nella cartella: grafici/")
print("=" * 72)

plt.show()
