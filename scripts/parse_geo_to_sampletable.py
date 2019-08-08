"""Download metadata from GEO and make sample table.

We are using Haiwang's data from GEO (GSE99574). Here we download
sample metadtaa from GEO and create a sample table.
"""
import tarfile
import gzip
import re
from collections import namedtuple
from tempfile import NamedTemporaryFile
from xml.etree import cElementTree
from urllib.request import urlopen
from io import BytesIO

GSE = snakemake.params[0]
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# GSE = 'GSE99574'

SampleAttributes = namedtuple(
    "SampleAttributes",
    "samplename GSM BioSample SRX species FBsp dev_stage tissue sex rep plate row col",
)


def main():
    xml = download_data()
    root = cElementTree.fromstring(xml)
    samples = parse_samples(root)
    write_samples(samples)


def download_data():
    url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:5]}nnn/{GSE}/miniml/{GSE}_family.xml.tgz"
    resp = urlopen(url)
    tar = tarfile.open(fileobj=BytesIO(resp.read()), mode="r:gz")
    return tar.extractfile(tar.getmembers()[0]).read()


def parse_samples(root):
    namespace = {
        "info": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML",
        "instance": "http://www.w3.org/2001/XMLSchema-instance",
    }

    samples = []
    for sample in root.findall("info:Sample", namespace):
        samples.append(
            SampleAttributes(
                sample.find("info:Title", namespace).text,
                sample.find("info:Accession", namespace).text,
                re.search(
                    r"SAMN\d+", sample.find(".//*[@type='BioSample']").attrib["target"].strip()
                )[0],
                re.search(r"SRX\d+", sample.find(".//*[@type='SRA']").attrib["target"].strip())[0],
                sample.find(".//info:Organism", namespace).text,
                sample.find(".//*[@tag='flybase species id']").text.strip(),
                sample.find(".//*[@tag='developmental stage']").text.strip(),
                sample.find(".//*[@tag='tissue']").text.strip(),
                sample.find(".//*[@tag='Sex']").text.strip().lower(),
                sample.find(".//*[@tag='replicate']").text.strip(),
                re.match(
                    r"Plate(\d+)_\w\d+", sample.find(".//*[@tag='plate and well id']").text.strip()
                )[1],
                re.match(
                    r"Plate\d+_(\w)\d+", sample.find(".//*[@tag='plate and well id']").text.strip()
                )[1],
                re.match(
                    r"Plate\d+_\w(\d+)", sample.find(".//*[@tag='plate and well id']").text.strip()
                )[1],
            )
        )

    return samples


def write_samples(samples):
    with open(OUTPUT_FILE, "w") as fh:
        fh.write("\t".join(samples[0]._fields) + "\n")
        fh.write(
            "\n".join(
                [
                    "\t".join(s)
                    for s in samples
                    if ("ERCC" not in s.samplename) and ("leftover" not in s.samplename)
                ]
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
