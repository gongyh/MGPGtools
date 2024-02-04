import logging
import os

# import multiprocessing
import pandas as pd
from pathlib import Path

# from functools import partial
from utils.meta import *
from utils.processDF import *
from utils.common import *
from utils.odgi import *
from utils.gfa import *


class Core(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        self.outdir = options.outdir
        self.threads = 1 if options.threads is not None else options.threads
        self.meta = os.path.join(self.database, "metadata", "Metadata.tsv")
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
        self.genomesNum = int(get_info(self.meta, self.name)["genomesNum"])
        self.ref = get_info(self.meta, self.name)["ref"]

    # 提取gff文件中的基因起始和终止位置写入bed文件中，并将每个基因长度存入字典
    def extract_gff(self, bedPath):
        geneTag = {}
        geneLength = {}
        with open(self.gff, "r") as f:
            new_file = open(os.path.join(bedPath, "genes.bed"), "w")
            lines = f.read().strip().split("\n")
            # line_count = 1
            for line in lines:
                if line[0] == "#":
                    continue
                row = line.strip().split("\t")
                if row[2] == "gene":
                    # line_count += 1
                    # if line_count > 20:
                    #     break
                    gene = row[8].split(";")[0].replace("ID=gene-", "")
                    tag = self.ref.replace(".", "#") + "#" + row[0]
                    k = tag + ":" + row[3] + "-" + row[4]
                    geneTag[k] = gene
                    l = int(row[4]) - int(row[3])
                    geneLength[gene] = l
                    new_file.write(tag + "\t" + row[3] + "\t" + row[4] + "\n")
            new_file.close()
        return geneTag, geneLength

    # 更换tsv文件第一列的名称
    def changeGeneTag(self, geneTag, dir, tsv):
        filepath = os.path.join(dir, tsv)
        newFilePath = os.path.join(dir, tsv.replace(".tsv", ".New.tsv"))
        with open(filepath, "r") as file:
            lines = file.readlines()
        new_lines = []
        for line in lines:
            columns = line.strip().split("\t")
            if columns[0] in geneTag:
                columns[0] = geneTag[columns[0]]
            new_line = "\t".join(columns) + "\n"
            new_lines.append(new_line)
        with open(newFilePath, "w") as file:
            file.writelines(new_lines)

    def staticCoreGene(self):
        # 检查临时文件存放目录，没有则创建
        check_directory(os.path.join(self.outdir, "tmp"))
        # 定义核心基因的字典
        coreGene = {"0-15%": [], "15-95%": [], "95-99%": [], "99-100%": [], "100%": []}
        # 包含所有基因组的列表
        genomeList = getPanTxt(self.database, self.name)
        # 利用odgi构建og格式文件，排序，根据bed文件提取子图，显示子图的paths信息，将子图转化为gfa格式文件
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        geneTag, genel = self.extract_gff(os.path.join(self.outdir, "tmp"))
        geneList = list(genel.keys())
        ogBuild(self.gfa, ogFile, self.threads)
        extractOGFile = os.path.join(self.outdir, "tmp", "genes.og")
        extractOgSortedFile = os.path.join(self.outdir, "tmp", "genes.sorted.og")
        bedFile = os.path.join(self.outdir, "tmp", "genes.bed")
        tsvFile = os.path.join(self.outdir, "tmp", "genes.tsv")
        newTsvFile = os.path.join(self.outdir, "tmp", "genes.New.tsv")
        gfaFile = os.path.join(self.outdir, "tmp", "genes.gfa")
        ogExtractBed(ogFile, extractOGFile, bedFile, self.threads)
        ogSort(extractOGFile, extractOgSortedFile, self.threads)
        ogPath(extractOgSortedFile, tsvFile, self.threads)
        ogView(extractOgSortedFile, gfaFile, self.threads)
        # 根据子图的gfa文件(包含所有基因序列的gfa)，提取每个节点的长度
        nodel = nodeLength(gfaFile)
        self.changeGeneTag(geneTag, os.path.join(self.outdir, "tmp"), "genes.tsv")
        df = pd.read_csv(newTsvFile, delimiter="\t")
        # results字典中的键是基因的ID，值是一个字典，这个字典保存每个基因组中该基因与参考不一致的node
        # 不在results字典键中的基因是与参考100%相同的基因
        results = processDf(df, coreGene)
        absenceGene = {}
        # pool = multiprocessing.Pool(processes = 16)
        # partial_processRow = partial(processRow, df_merged)
        # results = results.update(pool.map(partial_processRow, df_merged.iterrows()))
        # results = pool.map(partial_processRow, df_merged.itertuples(index=False))
        # print(results)
        # pool.close()
        # pool.join()
        for g, v in results.items():
            currentGeneLength = genel[g]
            uniqueNum = 0
            for i, z in v.items():
                length = 0
                for q in z:
                    length = length + nodel[q[5:]]
                if length / currentGeneLength > 0.2:
                    i = i.split("#")[0]+"."+i.split("#")[1]
                    uniqueNum = uniqueNum + 1
                    if g in absenceGene:
                        absenceGene[g].append(i)
                    else:
                        absenceGene[g] = [i]
            rate = (self.genomesNum - uniqueNum) / self.genomesNum
            if rate == 1:
                coreGene["100%"].append(g)
            if 0.99 <= rate < 1:
                coreGene["99-100%"].append(g)
            if 0.95 <= rate < 0.99:
                coreGene["95-99%"].append(g)
            if 0.15 <= rate < 0.95:
                coreGene["15-95%"].append(g)
            if 0 <= rate < 0.15:
                coreGene["0-15%"].append(g)
        total_genes = len(genel)
        gene_99_100 = len(coreGene["100%"]) + len(coreGene["99-100%"])
        geneMatrix = pd.DataFrame(index=geneList, columns=genomeList)
        geneMatrix[:] = 1
        geneMatrix = geneMatrix.astype(int)
        for a, b in absenceGene.items():
            for i in b:
                geneMatrix.at[a, i] = 0
        geneMatrix.to_csv(
            os.path.join(self.outdir, "gene_presence_absence.tsv"), sep="\t", index=True
        )
        with open(os.path.join(self.outdir, "summary_statistics.txt"), "w") as f:
            f.write("Core genes(99% <= strains <= 100%): {}\n".format(gene_99_100))
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
            f.write("Total genes: {}".format(total_genes))
        delete_temp_dir(os.path.join(self.outdir, "tmp"))
