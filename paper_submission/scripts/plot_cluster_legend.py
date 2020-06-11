import os
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt


def main():
    fig, ax = plt.subplots()

    for color, name in zip(snakemake.params.cluster_colors, snakemake.params.legend_names):
        ax.scatter([], [], color=color, label=name)

    plt.legend(loc='center', ncol=1, frameon=False)
    ax.set_axis_off()

    fig.savefig(snakemake.output[0])


if __name__ == '__main__':
    if os.getenv('SNAKE_DEBUG', False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config
        config = read_config("config/common.yaml")
        color_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir='paper_submission',
            params=dict(
                cluster_colors=color_config['clusters'],
                legend_names=config['legend_names']
            )
        )

    plt.style.use('../config/figure_styles.mplstyle')

    main()
