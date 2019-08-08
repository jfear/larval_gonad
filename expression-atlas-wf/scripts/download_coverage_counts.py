import shutil
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

URL = snakemake.params[0]
OUTPUT_FILE = snakemake.output[0]
SAMPLENAME = snakemake.wildcards["samplename"]


def main():
    tmpdir = TemporaryDirectory()
    gz = Path(tmpdir.name, f"{SAMPLENAME}.txt.gz")
    counts = Path(tmpdir.name, f"{SAMPLENAME}.txt")

    run(["curl", "-o", gz, URL])
    run(["gunzip", gz])
    shutil.copyfile(counts, OUTPUT_FILE)
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
