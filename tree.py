import logging
import os
import pandas as pd
import toytree
import toyplot.svg
from utils.gff import *
from utils.odgi import *
from utils.meta import *
from utils.common import *
from utils.gfa import nodeStr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class Tree(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        self.outdir = options.outdir
        self.gene = options.gene
        self.threads = 1 if options.threads is not None else options.threads
        self.meta = os.path.join(self.database, "metadata", "Metadata.tsv")
        self.gfa = os.path.join(
            self.database,
            "databases",
            "gfa",
            get_info(self.meta, self.name)["class"],
            self.name + ".gfa",
        )

    def drawTree(self):
        refSeq, records = self.getGeneRecord()
        SeqIO.write(
            records, os.path.join(self.outdir, "tmp", self.gene + ".fasta"), "fasta"
        )
        cmds = [
            "clustalw2",
            "-INFILE=" + os.path.join(self.outdir, "tmp", self.gene + ".fasta"),
            "-ALIGN",
            "-TYPE=DNA",
        ]
        if_success, stdout, stderr = run(cmds)
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
        delete_temp_dir(os.path.join(self.outdir, "tmp"))

    def drawTreeWithGfa(self, altGfa, sampleTxt, label):
        altRecords, isAltGenome, sampleList = self.getAltGeneRecords(
            altGfa, sampleTxt, label
        )
        refSeq, records = self.getGeneRecord()
        if not isAltGenome:
            for i in sampleList:
                altRecords.append(SeqRecord(Seq(refSeq), id=i, description=self.name))
        records.extend(altRecords)
        SeqIO.write(
            records, os.path.join(self.outdir, "tmp", self.gene + ".fasta"), "fasta"
        )
        cmds = [
            "clustalw2",
            "-INFILE=" + os.path.join(self.outdir, "tmp", self.gene + ".fasta"),
            "-ALIGN",
            "-TYPE=DNA",
        ]
        if_success, stdout, stderr = run(cmds)
        with open(os.path.join(self.outdir, "tmp", self.gene + ".dnd"), "r") as file:
            treNewick = file.read()
            treNewick = treNewick.replace("\n", "")
        tre = toytree.tree(treNewick, tree_format=0)
        style = {
            "tip_labels_align": True,
            "tip_labels_style": {"font-size": "9px"},
        }
        canvas, axes, mark = tre.draw(width=400, height=300, **style)
        toyplot.svg.render(
            canvas, os.path.join(self.outdir, label + "." + self.gene + ".svg")
        )
        delete_temp_dir(os.path.join(self.outdir, "tmp"))

    def getGeneRecord(self):
        tRange, chrom = self.getGeneRange()
        check_directory(os.path.join(self.outdir, "tmp"))
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        extractOGFile = os.path.join(
            self.outdir, "tmp", self.name + "." + self.gene + ".og"
        )
        extractOgSortedFile = os.path.join(
            self.outdir, "tmp", self.name + "." + self.gene + ".sorted.og"
        )
        csvFile = os.path.join(self.outdir, "tmp", self.name + "." + self.gene + ".csv")
        geneGfa = os.path.join(self.outdir, "tmp", self.name + "." + self.gene + ".gfa")
        refGenome = get_info(self.meta, self.name)["ref"]
        genomeList = getPanTxt(self.database, self.name)
        refPath = refGenome.replace(".", "#") + "#" + chrom
        ogBuild(self.gfa, ogFile, self.threads)
        ogExtract(ogFile, extractOGFile, refPath, tRange, self.threads)
        ogSort(extractOGFile, extractOgSortedFile, self.threads)
        ogPath(extractOgSortedFile, csvFile, self.threads)
        ogView(extractOgSortedFile, geneGfa, self.threads)
        node = nodeStr(geneGfa)
        df = pd.read_csv(csvFile, delimiter="\t")
        # df_merged = df.groupby(df.iloc[:, 0].str.split('#').str[0]).apply(self.merge_rows).reset_index(drop=True)
        gene_ratio = {}
        uniqueNodes = {}
        commonNodes = {}
        records = []
        refSeq = ""
        for i, r in df.iterrows():
            uniqueCol = [column for column, value in r.iloc[3:].items() if value == 0]
            commonCol = [column for column, value in r.iloc[3:].items() if value == 1]
            tag = r["path.name"].split("#")[0] + "." + r["path.name"].split("#")[1]
            uniqueNodes[tag] = uniqueCol
            commonNodes[tag] = commonCol
        for k, v in uniqueNodes.items():
            if k == refGenome:
                continue
            length = 0
            for z in v:
                length = length + len(node[z[5:]])
            rate = length / (int(tRange[1]) - int(tRange[0]))
            gene_ratio[k] = 1 - rate
        for x, y in commonNodes.items():
            id = x
            seq = ""
            for z in y:
                seq += node[z[5:]]
            if x == refGenome:
                refSeq = seq
            seqRecord = SeqRecord(Seq(seq), id=id, description=self.name)
            records.append(seqRecord)
        otherGenome = [item for item in genomeList if item not in commonNodes.keys()]
        for a in otherGenome:
            seqRecord = SeqRecord(Seq(refSeq), id=a, description=self.name)
            records.append(seqRecord)
        return refSeq, records

    def getAltGeneRecords(self, altGfa, sampleTxt, label):
        tRange, chrom = self.getGeneRange()
        check_directory(os.path.join(self.outdir, "tmp"))
        ogFile = os.path.join(self.outdir, "tmp", label + ".sorted.og")
        extractOGFile = os.path.join(
            self.outdir, "tmp", label + "." + self.gene + ".og"
        )
        extractOgSortedFile = os.path.join(
            self.outdir, "tmp", label + "." + self.gene + ".sorted.og"
        )
        csvFile = os.path.join(self.outdir, "tmp", label + "." + self.gene + ".csv")
        geneGfa = os.path.join(self.outdir, "tmp", label + "." + self.gene + ".gfa")
        refGenome = get_info(self.meta, self.name)["ref"]
        refPath = refGenome.replace(".", "#") + "#" + chrom
        ogBuild(altGfa, ogFile, self.threads)
        ogExtract(ogFile, extractOGFile, refPath, tRange, self.threads)
        ogSort(extractOGFile, extractOgSortedFile, self.threads)
        ogPath(extractOgSortedFile, csvFile, self.threads)
        ogView(extractOgSortedFile, geneGfa, self.threads)
        sampleList = []
        with open(sampleTxt, "r") as file:
            for line in file.readlines():
                sampleList.append(line.strip())
        node = nodeStr(geneGfa)
        df = pd.read_csv(csvFile, delimiter="\t")
        # df_merged = df.groupby(df.iloc[:, 0].str.split('#').str[0]).apply(self.merge_rows).reset_index(drop=True)
        commonNodes = {}
        records = []
        refSeq = ""
        for i, r in df.iterrows():
            commonCol = [column for column, value in r.iloc[3:].items() if value == 1]
            if "#" in r["path.name"]:
                tag = r["path.name"].split("#")[0] + "." + r["path.name"].split("#")[1]
            else:
                tag = r["path.name"]
            commonNodes[tag] = commonCol
        if len(commonNodes) == 1 and next(iter(commonNodes)) == refGenome:
            sampleList.remove(refGenome)
            return records, False, sampleList
        else:
            for x, y in commonNodes.items():
                id = x
                seq = ""
                for z in y:
                    seq += node[z[5:]]
                if x == refGenome:
                    refSeq = seq
                    continue
                seqRecord = SeqRecord(Seq(seq), id=id, description=self.name)
                records.append(seqRecord)
            otherGenome = [
                item for item in sampleList if item not in commonNodes.keys()
            ]
            for a in otherGenome:
                seqRecord = SeqRecord(Seq(refSeq), id=a, description=self.name)
                records.append(seqRecord)
            return records, True, []

    def getGeneRange(self):
        G = extract_gene(self.database, self.name)
        for i, j in G.items():
            if self.gene in j:
                tRange = j[self.gene]
                chrom = i
        return tRange, chrom
