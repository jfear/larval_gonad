from pathlib import Path
import shutil
from subprocess import run

RMD = Path(snakemake.input.rmd).resolve()
HTML = Path(snakemake.output.html).resolve()
ODIR = HTML.parent
TMP_RMD = Path(ODIR, RMD.name)

def main():
    shutil.copyfile(RMD, TMP_RMD)
    cmd = (
        f"cd {ODIR.as_posix()} && "
        f"""Rscript -e 'rmarkdown::render("{TMP_RMD.name}")'"""
    )
    run(cmd, shell=True)
    TMP_RMD.unlink()


if __name__ == "__main__":
    main()
