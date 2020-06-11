""" Figure 4: Chromosome locations."""
import os

from svgutils.compose import Figure, Panel, SVG, Text, Grid


def main():
    Figure("12cm", "12cm",
        Panel(
            SVG(snakemake.input.volume),
        ),
        Panel(
            SVG(snakemake.input.sphere),
        ).move(145, 0),
    ).save(snakemake.output[0])


if __name__ == '__main__':
    if os.getenv('SNAKE_DEBUG', False):
        from larval_gonad.debug import snakemake_debug
        snakemake = snakemake_debug(
            input=dict(
                volume="../../output/fish-wf/chrom_volume.svg",
                sphere="../../output/fish-wf/sphericity.svg"
            )
        )

    main()
