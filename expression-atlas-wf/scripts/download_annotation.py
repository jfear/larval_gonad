import shutil
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

URL = snakemake.params[0]
OUTPUT_FILE = snakemake.output[0]
SPECIES = snakemake.wildcards["species"]


def main():
    tmpdir = TemporaryDirectory()
    gz = Path(tmpdir.name, f"{SPECIES}.tar.gz")
    tar = Path(tmpdir.name, f"{SPECIES}.tar")
    gtf = Path(tmpdir.name, f"{SPECIES}.YO.gtf")

    run(["curl", "-o", gz, URL])
    run(["gunzip", gz])
    run(["tar", "--directory", tmpdir.name, "-xf", tar])
    shutil.copyfile(gtf, OUTPUT_FILE)
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
