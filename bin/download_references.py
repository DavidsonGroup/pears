#!/usr/bin/env python

import argparse
import urllib.request
import gzip

GENOME_LINKS = {
    "GRCh38": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
}

ANNO_LINKS = {
    "GENCODE40": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz",
    "GENCODE41": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz",
    "GENCODE42": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz",
    "GENCODE43": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz",
    "GENCODE44": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
    "GENCODE45": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz",
    "GENCODE46": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz",
    "GENCODE47": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz",
    "GENCODE48": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz",
    "GENCODE49": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
}


def main():
    parser = argparse.ArgumentParser(description="Download reference genomes for PEARS")
    parser.add_argument("--reference", required=True, help="Reference to download")
    args = parser.parse_args()

    # Reference should be specified in the format "GENOME+ANNO", e.g. "GRCh38+GENCODE49"
    if "+" not in args.reference:
        print("Error: Reference should be specified in the format 'GENOME+ANNO', e.g. 'GRCh38+GENCODE49'")
        return

    genome, anno = args.reference.split("+")
    if genome not in GENOME_LINKS:
        print(f"Error: Genome '{genome}' not found. Available genomes: {', '.join(GENOME_LINKS.keys())}")
        return

    if anno not in ANNO_LINKS:
        print(f"Error: Annotation '{anno}' not found. Available annotations: {', '.join(ANNO_LINKS.keys())}")
        return

    genome_link = GENOME_LINKS[genome]
    anno_link = ANNO_LINKS[anno]

    genome_name = genome_link.split("/")[-1]
    anno_name = anno_link.split("/")[-1]

    print(f"Downloading genome from {genome_link}...")
    urllib.request.urlretrieve(genome_link, genome_name)
    with gzip.open(genome_name, "rb") as f_in:
        with open(genome_name.replace(".gz", ""), "wb") as f_out:
            f_out.write(f_in.read())
    print("Download complete.")

    print(f"Downloading annotation from {anno_link}...")
    urllib.request.urlretrieve(anno_link, anno_name)
    with gzip.open(anno_name, "rb") as f_in:
        with open(anno_name.replace(".gz", ""), "wb") as f_out:
            f_out.write(f_in.read())
    print("Download complete.")

if __name__ == "__main__":
    main()
