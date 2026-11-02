# Self-contained demo: titration curves and the pI distribution of a simulated proteome.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pi_calc import net_charge, isoelectric_point, molecular_weight

AA = list("ARNDCEQGHILKMFPSTWYV")
# Approximate amino-acid background frequencies (%) in UniProt.
BG = np.array([8.3, 5.5, 4.1, 5.5, 1.4, 6.8, 3.9, 7.1, 2.3, 5.9, 9.7, 5.8, 2.4, 3.9, 4.7, 6.6, 5.4, 1.1, 2.9, 6.9])
BG = BG / BG.sum()


def proteome(n, rng):
    out = []
    for _ in range(n):
        comp = rng.dirichlet(BG * 60)
        out.append("M" + "".join(rng.choice(AA, size=int(rng.lognormal(5.8, 0.5)), p=comp)))
    return out


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(12)
    seqs = proteome(1500, rng)
    df = pd.DataFrame({"length": [len(s) for s in seqs],
                       "mw_kda": [round(molecular_weight(s) / 1000, 2) for s in seqs],
                       "pI_EMBOSS": [round(isoelectric_point(s, "EMBOSS"), 2) for s in seqs],
                       "pI_Lehninger": [round(isoelectric_point(s, "Lehninger"), 2) for s in seqs],
                       "charge_pH7.4": [round(net_charge(s, 7.4), 2) for s in seqs]})
    df.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    ph = np.linspace(0, 14, 200)
    picks = {"most acidic": df["pI_EMBOSS"].idxmin(), "closest to pI 7": (df["pI_EMBOSS"] - 7).abs().idxmin(),
             "most basic": df["pI_EMBOSS"].idxmax()}
    for (name, i), c in zip(picks.items(), ["#c0504d", "#444444", "#1f77b4"]):
        axes[0].plot(ph, [net_charge(seqs[i], p) for p in ph], color=c, label=f"{name} (pI {df.loc[i, 'pI_EMBOSS']})")
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].axvline(7.4, color="k", ls=":", lw=1, label="cytosolic pH 7.4")
    axes[0].set_ylim(-60, 60)
    axes[0].set_xlabel("pH")
    axes[0].set_ylabel("net charge")
    axes[0].set_title("Titration curves")
    axes[0].legend(fontsize=8)

    axes[1].hist(df["pI_EMBOSS"], bins=np.arange(3, 13, 0.2), color="#2a7f62")
    axes[1].axvline(7.4, color="k", ls=":", lw=1)
    axes[1].set_xlabel("isoelectric point")
    axes[1].set_ylabel("proteins")
    axes[1].set_title("Proteome pI is bimodal, with a dip near neutral")

    axes[2].scatter(df["pI_EMBOSS"], df["pI_Lehninger"] - df["pI_EMBOSS"], s=5, alpha=0.5)
    axes[2].axhline(0, color="k", lw=0.8)
    axes[2].set_xlabel("pI (EMBOSS pKa set)")
    axes[2].set_ylabel("Lehninger minus EMBOSS")
    axes[2].set_title("The pKa table changes the answer")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(df.describe().round(2))
    print("Largest disagreement between pKa sets:", round((df["pI_Lehninger"] - df["pI_EMBOSS"]).abs().max(), 2))


if __name__ == "__main__":
    main()
