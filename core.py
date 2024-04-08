import logging
import os
import multiprocessing
import pandas as pd
from Bio import Seq
from pathlib import Path
from functools import partial
from utils.meta import *
from utils.common import *
from utils.odgi import *
from utils.gfa import *


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

    # 提取gff文件中的基因起始和终止位置,基因ID, 以及用于下一步分析的基因路径字符串
    def extract_gff(self):
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
                    l = int(row[4]) - int(row[3])
                    geneLength[gene] = l
        return geneTag, geneLength, genePath

    def coreGene_fasta(self, coreGenes, outdir):
        coregff = os.path.join(outdir, "core.gff")
        corefasta = os.path.join(outdir, "core.fasta")
        with open(self.gff, "r") as fgff, open(coregff, "w") as fcoregff:
            lines = fgff.read().strip().split("\n")
            for line in lines:
                if line[0] == "#":
                    fcoregff.write(line + "\n")
                    continue
                row = line.strip().split("\t")
                if row[2] == "CDS":
                    gene = row[8].split(";")[0].replace("ID=", "")
                    if gene in coreGenes:
                        fcoregff.write(line + "\n")
        extractCoreGenesCmd = [
            "gffread",
            coregff,
            "-g",
            self.ref_genome,
            "-x",
            corefasta,
        ]
        isSuccess, stdout, stderr = run(extractCoreGenesCmd)

    def filtBlastGene(self, blastResult, geneLength):
        geneIteratorList = []
        assemblGenes = []
        with open(blastResult, "r") as f:
            lines = f.read().strip().split("\n")
            for l in lines:
                row = l.strip().split("\t")
                if row[0] in geneIteratorList:
                    continue
                else:
                    geneIteratorList.append(row[0])
                    geneName = row[0].replace("_gene", "")
                    gene_length = geneLength[geneName]
                    if (int(row[3]) * float(row[2]) / 100) / gene_length >= 0.8:
                        assemblGenes.append(geneName)
        return assemblGenes

    # 根据每个基因的位置提取子图
    def extractGenesOg(
        self, genePath, ogFile, outdir, geneTag, geneLength, genomeListExceptRef
    ):
        geneName = geneTag[genePath]
        gene_length = geneLength[geneName]
        # 每个基因的odgi文件
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
        ogView(extractSortedOg, extractGfa, 1)
        nodeL = nodeLength(extractGfa)
        if_success, stdout, stderr = ogPath(extractSortedOg, 1)
        lines = stdout.split("\n")
        lines.pop()
        columns = lines[0].split("\t")
        data = [l.split("\t") for l in lines[1:]]
        # 每个基因的矩阵
        df = pd.DataFrame(data, columns=columns)
        cols_to_extract = df.columns[3:][df.iloc[0, 3:] == "1"].tolist()
        cols_to_extract.insert(0, "path.name")
        genomeDF = df[cols_to_extract]
        # genome_df: 记录矩阵中出现的基因组名称
        # variant_genome: 这个基因上突变碱基总数大于0.2的基因组
        genome_df = []
        variant_genome = []
        for index, row in genomeDF.iterrows():
            # 去除参考基因组
            if genePath in str(row["path.name"]):
                continue
            # 基因组名称
            genomeName = (
                row["path.name"].split("#")[0] + "." + row["path.name"].split("#")[1]
            )
            if genomeName not in genome_df:
                genome_df.append(genomeName)
            else:
                if genomeName not in variant_genome:
                    continue
            zeroColumns = list(genomeDF.columns[1:][row[1:] == "0"])
            length = 0
            for n in zeroColumns:
                length += nodeL[n[5:]]
            if length / gene_length > 0.2:
                if genomeName not in variant_genome:
                    variant_genome.append(genomeName)
        absence_genome = list(set(genomeListExceptRef) - set(genome_df))
        variant_genome.extend(absence_genome)
        absenceGene = {}
        absenceGene[geneName] = variant_genome
        # delete_files(extractSortedOg)
        # delete_files(extractGfa)
        return absenceGene

    def staticCoreGene(self):
        # 检查临时文件存放目录，没有则创建
        tmp = os.path.join(self.outdir, "tmp")
        check_directory(tmp)
        # 定义核心基因的字典
        coreGene = {"0-15%": [], "15-95%": [], "95-99%": [], "99-100%": [], "100%": []}
        # 包含所有基因组的列表
        genomeList = getPanTxt(self.database, self.name)
        genomeNum = len(genomeList)
        # 包含除参考外的所有基因组的列表
        genomeList_except_ref = [x for x in genomeList if x != self.ref]
        # 利用odgi构建og格式文件，排序，根据bed文件提取子图，显示子图的paths信息，将子图转化为gfa格式文件
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        geneTag, geneL, genePath = self.extract_gff()
        geneList = list(geneL.keys())
        totalGenesNum = len(geneList)
        ogBuild(self.gfa, ogFile, self.threads)
        ref = self.ref.split(".")[0] + self.ref.split(".")[1]
        # 缺失基因字典
        absenceGene = {}
        result = []
        gene_99_100 = []
        # 多个进程处理每个基因的图, 合并结果
        corePool = multiprocessing.Pool(processes=self.process)
        partial_extractGeneOg = partial(
            self.extractGenesOg,
            ogFile=ogFile,
            outdir=tmp,
            geneTag=geneTag,
            geneLength=geneL,
            genomeListExceptRef=genomeList_except_ref,
        )
        result = corePool.map(partial_extractGeneOg, genePath)
        for j in result:
            absenceGene.update(j)
        corePool.close()
        corePool.join()
        # 基因矩阵
        geneMatrix = pd.DataFrame(index=geneList, columns=genomeList)
        # 默认全部是1
        geneMatrix[:] = 1
        geneMatrix = geneMatrix.astype(int)
        # 根据缺失基因字典将对应值改为0
        for a, b in absenceGene.items():
            # 总基因组数目减去缺失基因组数目，存在这个基因的基因组数
            l = genomeNum - len(b)
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
            for i in b:
                geneMatrix.at[a, i] = 0
        # 基因矩阵输出到gene_presence_absence.tsv
        geneMatrix.to_csv(
            os.path.join(self.outdir, "gene_presence_absence.tsv"), sep="\t", index=True
        )
        # 写入summary_statistics.txt
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
        self.coreGene_fasta(coreGenes, tmp)
        geneTag, geneL, genePath = self.extract_gff()
        blastIndex = os.path.join(tmp, "index")
        makeBlastDBCmd = [
            "makeblastdb",
            "-in",
            self.fasta,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            blastIndex,
        ]
        run(makeBlastDBCmd)
        coreGeneFasta = os.path.join(tmp, "core.fasta")
        blastResult = os.path.join(tmp, "blastResult")
        blastCmd = [
            "blastn",
            "-query",
            coreGeneFasta,
            "-db",
            blastIndex,
            "-evalue",
            "1e-6",
            "-outfmt",
            "6",
            "-num_threads",
            "6",
            "-out",
            blastResult,
        ]
        run(blastCmd)
        assemblGenes = self.filtBlastGene(blastResult, geneL)
        outputFile = os.path.join(self.outdir, "completeness.txt")
        with open(outputFile, "w") as f:
            f.write(
                "Core Gene Completeness: {}\n".format(
                    len(assemblGenes) / len(coreGenes)
                )
            )
        delete_temp_dir(os.path.join(self.outdir, "tmp"))
