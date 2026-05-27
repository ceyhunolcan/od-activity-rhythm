import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

steps = [
    ("NHANES 2013-2014 examined participants", "n = 9,813", None),
    ("Aged ≥40\nPST eligibility threshold", "n = 3,708", "Excluded: n = 6,105"),
    ("Not pregnant at examination", "n = 3,705", "Excluded: n = 3"),
    ("Completed Pocket Smell Test\nCSXEXSTS = 1", "n = 3,094", "Excluded: n = 611"),
    ("No recent acute upper-respiratory illness\nCSQ010 = 1", "n = 2,857", "Excluded: n = 237"),
    ("Has accelerometer record\nPAXHD_H", "n = 2,857", None),
    ("Met ≥4 valid-day wear-time inclusion\n≥16 h/day wake + sleep state", "n = 2,327", "Excluded: n = 530"),
    ("Final analytic sample", "n = 2,327", None),
]

fig, ax = plt.subplots(figsize=(10.5, 12))
ax.set_xlim(0, 14)
ax.set_ylim(0, 20)
ax.axis("off")

main_x = 1.0
main_w = 7.0
main_h = 1.25

ex_x = 10.0
ex_w = 3.1
ex_h = 0.78

ys = [17.0, 14.9, 12.8, 10.7, 8.6, 6.5, 4.4, 2.3]

for i, ((label, n, excluded), y) in enumerate(zip(steps, ys)):
    main_box = FancyBboxPatch(
        (main_x, y), main_w, main_h,
        boxstyle="round,pad=0.06",
        linewidth=1.2,
        edgecolor="black",
        facecolor="#f8f9fb",
    )
    ax.add_patch(main_box)

    ax.text(
        main_x + main_w/2, y + 0.80,
        label,
        ha="center", va="center",
        fontsize=9.2,
        fontweight="bold",
        linespacing=1.15,
    )

    ax.text(
        main_x + main_w/2, y + 0.33,
        n,
        ha="center", va="center",
        fontsize=10.2,
        fontweight="bold",
    )

    if i < len(ys) - 1:
        ax.annotate(
            "",
            xy=(main_x + main_w/2, ys[i+1] + main_h + 0.04),
            xytext=(main_x + main_w/2, y - 0.04),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color="black"),
        )

    if excluded:
        ex_box = FancyBboxPatch(
            (ex_x, y + 0.24), ex_w, ex_h,
            boxstyle="round,pad=0.05",
            linewidth=0.9,
            edgecolor="#999",
            facecolor="#f6efe4",
        )
        ax.add_patch(ex_box)

        ax.text(
            ex_x + ex_w/2, y + 0.63,
            excluded,
            ha="center", va="center",
            fontsize=8.8,
            fontweight="normal",
        )

        ax.annotate(
            "",
            xy=(ex_x - 0.08, y + 0.63),
            xytext=(main_x + main_w + 0.08, y + 0.63),
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
        )

ax.text(
    7,
    19.15,
    "Figure 1. Participant flow: olfactory dysfunction × accelerometer-derived activity, NHANES 2013-2014",
    ha="center",
    fontsize=10.8,
    fontweight="normal",
)

ax.text(
    7,
    1.0,
    "Figure 1. STROBE cohort flow diagram.",
    ha="center",
    fontsize=9.5,
    fontweight="normal",
)

plt.savefig(
    "paper/figures/figure1_strobe_remade.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
)
print("Saved paper/figures/figure1_strobe_remade.png")
