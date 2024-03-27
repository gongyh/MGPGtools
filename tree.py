import logging
import os
import pandas as pd
import toytree
import toyplot.svg
import multiprocessing
from functools import partial
from utils.odgi import *
from utils.meta import *
from utils.common import *
from utils.sequence import nucl_complement
from utils.gfa import gfa_parse_link_path, nodeStr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class Tree(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        self.outdir = options.outdir
        if options.gene is not None:
            self.gene = options.gene
        if options.genesFile is not None:
            self.genesFile = options.genesFile
        self.threads = 1 if options.threads is not None else options.threads
        self.meta = os.path.join(self.database, "metadata", "Metadata.tsv")
        self.ref = get_info(self.meta, self.name)["ref"]
        self.gfa = os.path.join(
            self.database,
            "databases",
            "gfa",
            get_info(self.meta, self.name)["class"],
            self.name + ".gfa",
        )
        self.gff = os.path.join(
            self.database,
            "databases",
            "gff",
            get_info(self.meta, self.name)["class"],
            get_info(self.meta, self.name)["ref"] + ".gff",
        )
        self.process = 16

    def drawGeneTree(self):
        refGenome = get_info(self.meta, self.name)["ref"]
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        ogBuild(self.gfa, ogFile, self.threads)
        gene_tag = self.extract_gff(list(self.gene))
        self.getGeneRecord(tmp, ogFile, gene_tag, refGenome, list(self.gene))
        with open(os.path.join(self.outdir, "tmp", self.gene + ".dnd"), "r") as file:
            treNewick = file.read()
            treNewick = treNewick.replace("\n", "")
        tre = toytree.tree(treNewick, tree_format=0)
        style = {
            "tip_labels_align": True,
            "tip_labels_style": {"font-size": "9px"},
        }
        canvas, axes, mark = tre.draw(width=400, height=300, **style)
        toyplot.svg.render(canvas, os.path.join(self.outdir, self.gene + ".svg"))
        # delete_temp_dir(os.path.join(self.outdir, "tmp"))

    def drawSpeciesTree(self):
        refGenome = get_info(self.meta, self.name)["ref"]
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        ogBuild(self.gfa, ogFile, self.threads)
        genesList = []
        with open(self.genesFile, "r") as f:
            genesList = f.read().strip().split("\n")
        gene_tag = self.extract_gff(genesList)
        pool = multiprocessing.Pool(processes=self.process)
        partial_getGeneRecord = partial(
            self.getGeneRecord,
            outdir=tmp,
            ogFile=ogFile,
            geneTag=gene_tag,
            refGenome=refGenome
        )
        pool.map(partial_getGeneRecord, genesList)
        with open(os.path.join(tmp, "total.nw"), "w") as totalnw:
            for i in genesList:
                with open(os.path.join(tmp, i + ".dnd"), "r") as singlenw:
                    fstr = singlenw.read()
                    totalnw.write(fstr)
        concatTreeCmd = [
            "astral",
            "-i",
            os.path.join(tmp, "total.nw"),
            "-o",
            os.path.join(tmp, "species.nw"),
        ]
        run(concatTreeCmd)
        with open(os.path.join(tmp, "species.nw"), "r") as fdnd:
            treNewick = fdnd.read()
            treNewick = treNewick.replace("\n", "")
        tre = toytree.tree(treNewick, tree_format=0)
        style = {
            "tip_labels_align": True,
            "tip_labels_style": {"font-size": "9px"},
        }
        canvas, axes, mark = tre.draw(width=400, height=300, **style)
        toyplot.svg.render(canvas, os.path.join(self.outdir, "speciesTree.svg"))

    def getGeneRecord(self, gene, outdir, ogFile, geneTag, refGenome):
        chrom = geneTag[gene][0]
        start = geneTag[gene][1]
        end = geneTag[gene][2]
        extractOgSortedFile = os.path.join(
            outdir, self.name + "." + gene + ".sorted.og"
        )
        geneGfa = os.path.join(outdir, self.name + "." + gene + ".gfa")
        refPath = chrom + ":" + start + "-" + end
        extractCmd = [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-r",
            refPath,
            "-d",
            "3000",
            "-o",
            "-",
        ]
        extractSortCmd = ["odgi", "sort", "-i", "-", "-o", extractOgSortedFile]
        p1 = subprocess.Popen(extractCmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(extractSortCmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p2.communicate()
        ogView(extractOgSortedFile, geneGfa, self.threads)
        gfa_node, gfa_path = gfa_parse_link_path(geneGfa)
        genome_path_list = []
        records = []
        refSeq = ""
        for chrom, p in gfa_path.items():
            genomeName = chrom.split("#")[0] + "." + chrom.split("#")[1]
            if genomeName in genome_path_list:
                continue
            seq = ""
            for i in p["path"]:
                nodeName = "node_" + i[:-1]
                if i[-1] == "+":
                    seq += gfa_node[nodeName]
                if i[-1] == "-":
                    seq += nucl_complement(gfa_node[nodeName])
            if p["tag"] == "reverse":
                seq = nucl_complement(seq)
            if genomeName == refGenome:
                refSeq = seq
            seqRecord = SeqRecord(Seq(seq), id=genomeName, description=self.name)
            records.append(seqRecord)
            genome_path_list.append(genomeName)
        SeqIO.write(records, os.path.join(outdir, gene + ".fasta"), "fasta")
        cmds = [
            "clustalw2",
            "-INFILE=" + os.path.join(outdir, gene + ".fasta"),
            "-ALIGN",
            "-TYPE=DNA",
        ]
        run(cmds)

    def extract_gff(self, genes):
        geneTag = {}
        with open(self.gff, "r") as f:
            lines = f.read().strip().split("\n")
            for line in lines:
                if line[0] == "#":
                    continue
                row = line.strip().split("\t")
                if row[2] == "CDS":
                    gene = row[8].split(";")[0].replace("ID=", "")
                    if gene in genes:
                        tag = self.ref.replace(".", "#") + "#" + row[0]
                        geneTag[gene] = [tag, row[3], row[4]]
        return geneTag
