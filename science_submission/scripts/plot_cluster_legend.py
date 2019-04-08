import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

cluster_colors = snakemake.params.cluster_colors
legend_names = snakemake.params.legend_names
oname = snakemake.output[0]


def main():
    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(2.25, 1))

    for color, name in zip(cluster_colors, legend_names):
        ax.scatter([], [], color=color, label=name)

    plt.legend(loc='center', ncol=3, frameon=False)
    ax.set_axis_off()

    fig.savefig(oname, bbox_inches='tight')


if __name__ == '__main__':
    main()
