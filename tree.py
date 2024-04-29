import logging
import os
import pandas as pd
import toytree
import toyplot.svg
import multiprocessing
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from functools import partial
from collections import OrderedDict
from utils.odgi import ogBuild, ogView
from utils.meta import get_info
from utils.common import check_directory, delete_temp_dir, run
from utils.sequence import nucl_complement
from utils.gfa import gfa_parse_link_path, nodeStr
from utils.mummer import *


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
        self.fasta = options.fasta
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
        self.ref_genome = os.path.join(
            self.database,
            "databases",
            "ref_genome",
            get_info(self.meta, self.name)["class"],
            get_info(self.meta, self.name)["ref"] + "_genomic.fna",
        )
        self.process = 16
        self.ASTRAL = "/mnt/me4/zhousq/Astral/astral.5.7.8.jar"
        self.nb = 100

    def drawGeneTree(self):
        """
        Construct a phylogenetic tree for a single gene.
        """
        refGenome = get_info(self.meta, self.name)["ref"]
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        ogBuild(self.gfa, ogFile, self.threads)
        gene_tag = self.extract_gff(list(self.gene))
        self.ramxl(tmp, ogFile, gene_tag, refGenome, list(self.gene), False, "", None)
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
        """
        Construct a phylogenetic tree composed of core genes from a specific species.
        """
        refGenome = get_info(self.meta, self.name)["ref"]
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        tmp = os.path.join(self.outdir, "tmp")
        allSingleGenes_tree = os.path.join(
            self.outdir, "tmp", "allSingleGenes_tree.nwk"
        )
        allSingleGenes_bootstrap = os.path.join(
            self.outdir, "tmp", "allSingleGenes_bootstrap.txt"
        )
        astral_out = os.path.join(self.outdir, "tmp", "Astral.coalescend_out.result")
        astral_tree = os.path.join(self.outdir, "tmp", "Astral.coalescend_tree.nwk")
        check_directory(tmp)
        ogBuild(self.gfa, ogFile, self.threads)
        genesList = []
        with open(self.genesFile, "r") as f:
            genesList = f.read().strip().split("\n")
        gene_tag = self.extract_gff(genesList)
        pool = multiprocessing.Pool(processes=self.process)
        partial_ramxl = partial(
            self.ramxl,
            outdir=tmp,
            ogFile=ogFile,
            geneTag=gene_tag,
            refGenome=refGenome,
            withContig=False,
            contigName="",
            extractFastaDict=None,
        )
        pool.map(partial_ramxl, genesList)
        with open(allSingleGenes_tree, "w") as f:
            for i in genesList:
                with open(
                    os.path.join(tmp, "RAxML_bipartitions." + i + "_out.nwk"), "r"
                ) as fin:
                    f.write(fin.read())
        with open(allSingleGenes_bootstrap, "w") as f:
            for i in genesList:
                f.write(os.path.join(tmp, "RAxML_bootstrap." + i + "_out.nwk") + "\n")
        astralCmd = [
            "java",
            "-jar",
            self.ASTRAL,
            "-i",
            allSingleGenes_tree,
            "-b",
            allSingleGenes_bootstrap,
            "-r",
            str(self.nb),
            "-o",
            astral_out,
        ]
        run(astralCmd)
        os.system("tail -n 1 " + astral_out + " > " + astral_tree)
        tre = toytree.tree(astral_tree, tree_format=0)
        style = {
            "tip_labels_align": True,
            "tip_labels_style": {"font-size": "9px"},
            "node_labels": "support",
            "node_sizes": 15,
        }
        canvas, axes, mark = tre.draw(width=400, height=300, **style)
        toyplot.svg.render(canvas, os.path.join(self.outdir, "speciesTree.svg"))

    def drawTreeWithContig(self):
        refGenome = get_info(self.meta, self.name)["ref"]
        contig = self.fasta.split("/")[-1]
        contigName = contig[: contig.rfind(".")]
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        # ogBuild(self.gfa, ogFile, self.threads)
        genesList = []
        with open(self.genesFile, "r") as f:
            genesList = f.read().strip().split("\n")
        gff2fasta(self.gff, self.ref_genome, genesList, tmp)
        geneTag, geneL, genePath = self.extractGff()
        coreGeneFasta = os.path.join(tmp, "core.fasta")
        genesMapDict = filtNucmerResult(self.fasta, coreGeneFasta, geneL, tmp)
        coreGenesWithContig = []
        idConvert = {}
        bedFile = os.path.join(tmp, "core.bed")
        extractFastaFile = os.path.join(tmp, "extractCore.fasta")
        with open(bedFile, "w") as f:
            for k, v in genesMapDict.items():
                if v["assemblCore"]:
                    coreGenesWithContig.append(k)
                convertStr = (
                    ">"
                    + v["refNode"]
                    + "_"
                    + str(v["ref_start"])
                    + "-"
                    + str(v["ref_end"])
                    + ":. "
                )
                idConvert[convertStr] = k
                f.write(
                    v["refNode"]
                    + "\t"
                    + str(v["ref_start"] - 1)
                    + "\t"
                    + str(v["ref_end"])
                    + "\n"
                )
        subseqCmd = ["seqkit", "subseq", "--bed", bedFile, self.fasta]
        if_success, stdout, stderr = run(subseqCmd)
        with open(extractFastaFile, "w") as f:
            for line in stdout.strip().split("\n"):
                if ">" in line:
                    f.write(">{}\n".format(idConvert[line]))
                else:
                    f.write(line + "\n")
        extractFastaDict = SeqIO.to_dict(SeqIO.parse(extractFastaFile, "fasta"))
        for i, j in extractFastaDict.items():
            if genesMapDict[i]["orientation"] == "reverse":
                extractFastaDict[i].seq = extractFastaDict[i].seq.reverse_complement()
        gene_tag = self.extract_gff(coreGenesWithContig)
        pool = multiprocessing.Pool(processes=self.process)
        partial_ramxl = partial(
            self.ramxl,
            outdir=tmp,
            ogFile=ogFile,
            geneTag=gene_tag,
            refGenome=refGenome,
            withContig=True,
            contigName=contigName,
            extractFastaDict=extractFastaDict,
        )
        pool.map(partial_ramxl, coreGenesWithContig)
        with open(os.path.join(tmp, "total.nw"), "w") as totalnw:
            for i in coreGenesWithContig:
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

    def ramxl(
        self,
        gene,
        outdir,
        ogFile,
        geneTag,
        refGenome,
        withContig,
        contigName,
        extractFastaDict,
    ):
        """
        Handle each gene to get RAxML output.
        """
        chrom = geneTag[gene][0]
        start = str(int(geneTag[gene][1]) - 1)
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
        if withContig:
            extractFastaDict[gene].id = contigName
            records.append(extractFastaDict[gene])
        SeqIO.write(records, os.path.join(outdir, gene + ".fasta"), "fasta")
        mafftCmd = ["mafft", "--auto", os.path.join(outdir, gene + ".fasta")]
        isSuccess, stdout, stderr = run(mafftCmd)
        mafftFasta = os.path.join(outdir, gene + ".mafft.fasta")
        trimalFasta = os.path.join(outdir, gene + ".trimed.fasta")
        out_tree_file = gene + "_out.nwk"
        with open(mafftFasta, "w") as f:
            f.write(stdout)
        trimalCmd = ["trimal", "-in", mafftFasta, "-out", trimalFasta, "-automated1"]
        run(trimalCmd)
        raxmlCmd = [
            "raxmlHPC-PTHREADS",
            "-f",
            "a",
            "-N",
            str(self.nb),
            "-m",
            "PROTGAMMAJTT",
            "-x",
            "123456",
            "-p",
            "123456",
            "-s",
            trimalFasta,
            "-n",
            out_tree_file,
            "-w",
            outdir,
        ]
        isSuccess, stdout, stderr = run(raxmlCmd)

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

    def extractGff(self):
        genePath = []
        geneTag = {}
        geneLength = {}
        with open(self.gff, "r") as f:
            lines = f.read().strip().split("\n")
            # line_count = 1
            for line in lines:
                if line[0] == "#":
                    continue
                row = line.strip().split("\t")
                if row[2] == "CDS":
                    # line_count += 1
                    # if line_count > 20:
                    #     break
                    gene = row[8].split(";")[0].replace("ID=", "")
                    tag = self.ref.replace(".", "#") + "#" + row[0]
                    k = tag + ":" + row[3] + "-" + row[4]
                    genePath.append(k)
                    geneTag[k] = gene
                    l = int(row[4]) - int(row[3]) + 1
                    geneLength[gene] = l
        return geneTag, geneLength, genePath
