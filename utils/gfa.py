import re
import sys
import os
from utils.common import *
from utils.odgi import *
from utils.meta import get_info, get_chromosome


def gfa_parse_link_path(database, name):
    """Parse link and path from gfa files

    Returns:
        gfa_path: paths in gfa files
        gfa_link: links in gfa files
    """
    gfa_path = {}
    gfa_link = {}
    unsigned_gfa_link = {}
    details = get_info(os.path.join(database, "metadata", "Metadata.tsv"), name)
    gfa_file = os.path.join(
        database, "databases", "gfa", details["class"], name + ".gfa"
    )
    with open(gfa_file, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "L":
                prec_node = row[1] + row[2]
                prec_node_number = row[1]
                rear_node = row[3] + row[4]
                rear_node_number = row[3]
                if prec_node in gfa_link:
                    gfa_link[prec_node].append(rear_node)
                    unsigned_gfa_link[prec_node_number].append(rear_node_number)
                else:
                    gfa_link[prec_node] = [rear_node]
                    unsigned_gfa_link[prec_node_number] = [rear_node_number]
            if row[0] == "P":
                # chrom = re.sub(r"\[.*?\]", "", row[1])
                chrom = row[1]
                gfa_path[chrom] = {}
                if "," in row[2]:
                    gfa_path[chrom]["path"] = row[2].split(",")
                    gfa_path[chrom]["overlap"] = row[3]
                else:
                    gfa_path[chrom]["path"] = [row[2]]
                    gfa_path[chrom]["overlap"] = row[3]
    return gfa_link, unsigned_gfa_link, gfa_path


def segmentStr(gfaFile):
    segment = {}
    with open(gfaFile, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "S":
                segment[row[1]] = row[2]
    return segment
