"""Downloads ortholog tables from GEO."""
import shutil
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


def main():
    tmpdir = TemporaryDirectory()
    gz = Path(tmpdir.name, f"{snakemake.wildcards.species}.tar.gz")
    tar = Path(tmpdir.name, f"{snakemake.wildcards.species}.tar")
    ortho = Path(tmpdir.name, f"{snakemake.wildcards.species}.ortholog.txt")

    run(["curl", "-o", gz, snakemake.params[0]])
    run(["gunzip", gz])
    run(["tar", "--directory", tmpdir.name, "-xf", tar])
    shutil.copyfile(ortho, snakemake.output[0])
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
