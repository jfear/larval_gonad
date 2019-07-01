import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

CLUSTER_COLORS = snakemake.params.cluster_colors
LEGEND_NAMES = snakemake.params.legend_names
ONAME = snakemake.output[0]


def main():
    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(2.25, 1))

    for color, name in zip(CLUSTER_COLORS, LEGEND_NAMES):
        ax.scatter([], [], color=color, label=name)

    plt.legend(loc='center', ncol=3, frameon=False)
    ax.set_axis_off()

    fig.savefig(ONAME, bbox_inches='tight')


if __name__ == '__main__':
    main()
