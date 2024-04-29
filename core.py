import logging
import os
import multiprocessing
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from functools import partial
from utils.meta import get_info, getPanTxt
from utils.common import check_directory, delete_temp_dir, run
from utils.odgi import ogBuild, ogView, ogPath
from utils.gfa import nodeLength
from utils.sequence import *
from utils.mummer import *


class Core(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        self.outdir = options.outdir
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
        self.genomesNum = int(get_info(self.meta, self.name)["genomesNum"])
        self.coreGenes = options.coreGenes
        self.process = 16

    def extract_gff(self):
        """
        Extract the start and end positions of genes, gene ID,
        and the gene path string for further analysis from the gff file.
        """
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
                    # if line_count > 10:
                    #     break
                    gene = row[8].split(";")[0].replace("ID=", "")
                    tag = self.ref.replace(".", "#") + "#" + row[0]
                    k = tag + ":" + str(int(row[3]) - 1) + "-" + row[4]
                    genePath.append(k)
                    geneTag[k] = gene
                    l = int(row[4]) - int(row[3]) + 1
                    geneLength[gene] = l
        return geneTag, geneLength, genePath

    def extractGenesOg(
        self,
        genePath,
        ogFile,
        outdir,
        geneTag,
        geneLength,
        genomeListExceptRef,
        genomeNum,
    ):
        """
        Extract subgraphs based on the positions of each gene.
        """
        geneName = geneTag[genePath]
        gene_length = geneLength[geneName]
        # ODGI file for each gene.
        extractSortedOg = os.path.join(outdir, genePath + ".sorted.og")
        extractGfa = os.path.join(outdir, genePath + ".gfa")
        extractCmd = [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-r",
            genePath,
            "-d",
            "3000",
            "-o",
            "-",
        ]
        extractSortCmd = ["odgi", "sort", "-i", "-", "-o", extractSortedOg]
        p1 = subprocess.Popen(extractCmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(extractSortCmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p2.communicate()
        # Convert the sorted og file into a gfa file.
        ogView(extractSortedOg, extractGfa, 1)
        nodeL = nodeLength(extractGfa)
        if_success, stdout, stderr = ogPath(extractSortedOg, 1)
        lines = stdout.split("\n")
        lines.pop()
        columns = lines[0].split("\t")
        data = [l.split("\t") for l in lines[1:]]
        # Matrix for each gene.
        df = pd.DataFrame(data, columns=columns)
        cols_to_extract = df.columns[3:][df.iloc[0, 3:] == "1"].tolist()
        cols_to_extract.insert(0, "path.name")
        genomeDF = df[cols_to_extract]
        # genome_df: Record the genome names present in the matrix.
        # variant_genome: Genomes with a total number of mutated bases greater than 0.2 on this gene.
        dna_80_genome = []
        for index, row in genomeDF.iterrows():
            # remove reference genome
            if genePath in str(row["path.name"]):
                continue
            # genome id
            genomeName = (
                row["path.name"].split("#")[0] + "." + row["path.name"].split("#")[1]
            )
            if genomeName in dna_80_genome:
                continue
            zeroColumns = list(genomeDF.columns[1:][row[1:] == "0"])
            length = 0
            for n in zeroColumns:
                length += nodeL[n[5:]]
            if length / gene_length <= 0.2:
                dna_80_genome.append(genomeName)
        if len(dna_80_genome) != 0:
            nodeDict = {}
            pathDict = {}
            genes = {}
            genes[geneName] = {}
            genes[geneName]["core"] = []
            genes[geneName]["prot"] = []
            with open(extractGfa, "r") as f:
                l = f.read().strip().split("\n")
                for i in l:
                    row = i.strip().split("\t")
                    if row[0] == "S":
                        nodeDict[row[1]] = row[2]
                    if row[0] == "P":
                        # genome name.
                        gName = row[1].split("#")[0] + "." + row[1].split("#")[1]
                        if gName == self.ref:
                            if "," in row[2]:
                                p = row[2].split(",")
                            else:
                                p = [row[2]]
                            seq = ""
                            for s in p:
                                if s[-1] == "+":
                                    seq += nodeDict[s[:-1]]
                                else:
                                    seq += nucl_complement(nodeDict[s[:-1]])
                            seq = get_protein_coding_regions(seq)
                            refSeq = Seq(seq)
                        else:
                            if gName in dna_80_genome:
                                if "," in row[2]:
                                    p = row[2].split(",")
                                else:
                                    p = [row[2]]
                                seq = ""
                                for s in p:
                                    if s[-1] == "+":
                                        seq += nodeDict[s[:-1]]
                                    else:
                                        seq += nucl_complement(nodeDict[s[:-1]])
                                pathDict[row[1]] = {}
                                if int(p[0][:-1]) < int(p[-1][:-1]):
                                    seq = get_protein_coding_regions(seq)
                                else:
                                    seq = get_protein_coding_regions(
                                        nucl_complement(seq)
                                    )
                                pathDict[row[1]] = Seq(seq)
            refProt = refSeq.translate(table="Bacterial")
            refProtLen = len(refProt)
            for k, v in pathDict.items():
                aligner = Align.PairwiseAligner()
                aligner.mode = "global"
                alignments = aligner.align(refProt, v.translate(table="Bacterial"))
                # Determine if the protein similarity is greater than 0.85.
                if int(alignments[0].score) / refProtLen > 0.85:
                    gName = k.split("#")[0] + "." + k.split("#")[1]
                    if gName not in genes[geneName]["core"]:
                        genes[geneName]["core"].append(gName)
                        genes[geneName]["prot"].append(
                            SeqRecord(v.translate(table="Bacterial"), id=gName)
                        )
            if len(genes[geneName]["core"]) < genomeNum:
                genes[geneName]["prot"] = []
            return genes
        else:
            return {}

    def staticCoreGene(self):
        # Check the temporary file storage directory, create if not exist.
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        # Define variables for core genes.
        coreGene = {"0-15%": [], "15-95%": [], "95-99%": [], "99-100%": [], "100%": []}
        # List containing all genomes.
        genomeList = getPanTxt(self.database, self.name)
        genomeNum = len(genomeList)
        # List containing all genomes except the reference.
        genomeList_except_ref = [x for x in genomeList if x != self.ref]
        # Utilize ODGI to construct an OG-formatted file, sort it, extract subgraphs based on a BED file,
        # display path information of the subgraphs, and convert the subgraphs to GFA format file.
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        geneTag, geneL, genePath = self.extract_gff()
        geneList = list(geneL.keys())
        totalGenesNum = len(geneList)
        ogBuild(self.gfa, ogFile, self.threads)
        ref = self.ref.split(".")[0] + self.ref.split(".")[1]
        # Missing gene dictionary.
        existGenes = {}
        result = []
        gene_99_100 = []
        # Multiple processes handle the graph for each gene, then merge the results.
        corePool = multiprocessing.Pool(processes=self.process)
        partial_extractGeneOg = partial(
            self.extractGenesOg,
            ogFile=ogFile,
            outdir=tmp,
            geneTag=geneTag,
            geneLength=geneL,
            genomeListExceptRef=genomeList_except_ref,
            genomeNum=genomeNum,
        )
        result = corePool.map(partial_extractGeneOg, genePath)
        for j in result:
            existGenes.update(j)
        corePool.close()
        corePool.join()
        # Gene matrix
        geneMatrix = pd.DataFrame(index=geneList, columns=genomeList)
        # By default, all are set to 1.
        geneMatrix[:] = 0
        geneMatrix = geneMatrix.astype(int)
        # Change the corresponding values to 0 based on the missing gene dictionary.
        for a, b in existGenes.items():
            # The number of total genomes minus the number of missing genomes,
            # gives the number of genomes that have this gene.
            l = len(b["core"])
            if 0 <= l / genomeNum < 0.15:
                coreGene["0-15%"].append(a)
            if 0.15 <= l / genomeNum < 0.95:
                coreGene["15-95%"].append(a)
            if 0.95 <= l / genomeNum < 0.99:
                coreGene["95-99%"].append(a)
            if 0.99 <= l / genomeNum < 1:
                coreGene["99-100%"].append(a)
                gene_99_100.append(a)
            if l / genomeNum == 1:
                coreGene["100%"].append(a)
                gene_99_100.append(a)
            for i in b["core"]:
                geneMatrix.at[a, i] = 1
        # Output gene matrix to gene_presence_absence.tsv.
        geneMatrix.to_csv(
            os.path.join(self.outdir, "gene_presence_absence.tsv"), sep="\t", index=True
        )
        # Write to summary_statistics.txt.
        with open(os.path.join(self.outdir, "summary_statistics.txt"), "w") as f:
            f.write("Core genes(99% <= strains <= 100%): {}\n".format(len(gene_99_100)))
            f.write(
                "Soft core genes(95% <= strains < 99%): {}\n".format(
                    len(coreGene["95-99%"])
                )
            )
            f.write(
                "Shell genes(15% <= strains < 95%): {}\n".format(
                    len(coreGene["15-95%"])
                )
            )
            f.write(
                "Cloud genes(0% <= strains < 15%): {}\n".format(len(coreGene["0-15%"]))
            )
            f.write("Total genes: {}\n".format(totalGenesNum))
        delete_temp_dir(os.path.join(self.outdir, "tmp"))

    def contigCompleteness(self):
        with open(self.coreGenes, "r") as f:
            coreGenes = f.read().strip().split("\n")
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        gff2fasta(self.gff, self.ref_genome, coreGenes, tmp)
        geneTag, geneL, genePath = self.extract_gff()
        coreGeneFasta = os.path.join(tmp, "core.fasta")
        genesMapDict = filtNucmerResult(self.fasta, coreGeneFasta, geneL, tmp)
        assemblGenes = []
        for k, v in genesMapDict.items():
            if v["assemblCore"]:
                assemblGenes.append(k)
        outputFile = os.path.join(self.outdir, "completeness.txt")
        with open(outputFile, "w") as f:
            f.write(
                "Core Gene Completeness: {}\n".format(
                    len(assemblGenes) / len(coreGenes)
                )
            )
        delete_temp_dir(os.path.join(self.outdir, "tmp"))
