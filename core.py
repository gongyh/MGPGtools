import logging
import os
import multiprocessing
import pandas as pd
from pathlib import Path
from functools import partial
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
        self.process = 16

    # 提取gff文件中的基因起始和终止位置写入bed文件中，并将每个基因长度存入字典
    def extract_gff(self, bedPath):
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
                if row[2] == "gene":
                    # line_count += 1
                    # if line_count > 20:
                    #     break
                    gene = row[8].split(";")[0].replace("ID=", "")
                    if "GCF_001267405.1_01747" in gene:
                        tag = self.ref.replace(".", "#") + "#" + row[0]
                        k = tag + ":" + row[3] + "-" + row[4]
                        genePath.append(k)
                        geneTag[k] = gene
                        l = int(row[4]) - int(row[3])
                        geneLength[gene] = l
                        break
                    else:
                        continue
        return geneTag, geneLength, genePath

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
        geneTag, geneL, genePath = self.extract_gff(os.path.join(self.outdir, "tmp"))
        geneList = list(geneL.keys())
        total_genes = len(geneList)
        ogBuild(self.gfa, ogFile, self.threads)
        ref = self.ref.split(".")[0] + self.ref.split(".")[1]
        # 缺失基因字典
        absenceGene = {}
        result = []
        gene_99_100 = []
        # 多个进程处理每个基因的图, 合并结果
        pool = multiprocessing.Pool(processes=self.process)
        partial_extractGeneOg = partial(
            extractGenesOg,
            ogFile=ogFile,
            outdir=tmp,
            geneTag=geneTag,
            geneLength=geneL,
            genomeListExceptRef=genomeList_except_ref,
        )
        result = pool.map(partial_extractGeneOg, genePath)
        for j in result:
            absenceGene.update(j)
        pool.close()
        pool.join()
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
            f.write("Total genes: {}".format(total_genes))
        delete_temp_dir(os.path.join(self.outdir, "tmp"))