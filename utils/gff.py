import os
from utils.meta import get_info
from utils.exceptions import *


def extract_gene(database, name):
    chrom_gene_annotation = {}
    gff_path = os.path.join(
        database,
        "databases",
        "gff",
        get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)["class"],
        get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)["ref"]
        + ".gff",
    )
    if not os.path.exists(gff_path):
        raise gffNotFound("Annotation file(gff) not found!")
    with open(gff_path, "r") as f2:
        lines = f2.read().strip().split("\n")
        for line in lines:
            if line[0] == "#":
                continue
            row = line.strip().split("\t")
            if row[2] == "CDS":
                if row[0] not in chrom_gene_annotation:
                    chrom_gene_annotation[row[0]] = {}
                start = row[3]
                end = row[4]
                gene = row[8].split(";")[0].replace("ID=gene-", "")
                chrom_gene_annotation[row[0]][gene] = [start, end]
    return chrom_gene_annotation
