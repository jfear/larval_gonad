import shutil
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


def main():
    tmpdir = TemporaryDirectory()
    gz = Path(tmpdir.name, f"{snakemake.wildcards.samplename}.txt.gz")
    counts = Path(tmpdir.name, f"{snakemake.wildcards.samplename}.txt")

    run(["curl", "-o", gz, snakemake.params[0]])
    run(["gunzip", gz])
    shutil.copyfile(counts, snakemake.output[0])
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
