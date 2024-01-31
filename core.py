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

    def extract_gff(self, bedPath):
        geneTag = {}
        geneLength = {}
        export_bedFiles = 1
        with open(self.gff, "r") as f1:
            line_count = 0
            file_count = 1
            new_file = open(os.path.join(bedPath, str(file_count) + ".bed"), "w")
            lines = f1.read().strip().split("\n")
            for line in lines:
                if line[0] == "#":
                    continue
                row = line.strip().split("\t")
                if row[2] == "gene":
                    gene = row[8].split(";")[0].replace("ID=gene-", "")
                    tag = self.ref.replace(".", "#") + "#" + row[0]
                    k = tag + ":" + row[3] + "-" + row[4]
                    geneTag[k] = gene
                    l = int(row[4]) - int(row[3])
                    geneLength[gene] = l
                    new_file.write(tag + "\t" + row[3] + "\t" + row[4] + "\n")
                    line_count += 1
                    if line_count == 50:
                        new_file.close()
                        file_count += 1
                        line_count = 0
                        new_file = open(os.path.join(bedPath, str(file_count) + ".bed"), "w")
            new_file.close()
            export_bedFiles = file_count
        return geneTag, geneLength, export_bedFiles

    def changeGeneTag(self, geneTag, dir, csv):
        filepath = os.path.join(dir, csv)
        newFilePath = os.path.join(dir, csv.replace(".tsv", ".New.tsv"))
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
        check_directory(os.path.join(self.outdir, "tmp"))
        coreGene = {"0-15%": [], "15-95%": [], "95-99%": [], "99-100%": [], "100%": []}
        ogFile = os.path.join(self.outdir, "tmp", self.name + ".sorted.og")
        geneTag, genel, bedFiles = self.extract_gff(os.path.join(self.outdir, "tmp"))
        ogBuild(self.gfa, ogFile, self.threads)
        for i in range(1, bedFiles + 1):
            extractOGFile = os.path.join(self.outdir, "tmp", str(i) + ".og")
            extractOgSortedFile = os.path.join(self.outdir, "tmp", str(i) + ".sorted.og")
            bedFile = os.path.join(self.outdir, "tmp", str(i) + ".bed")
            tsvFile = os.path.join(self.outdir, "tmp", str(i) + ".tsv")
            newTsvFile = os.path.join(self.outdir, "tmp", str(i) + ".New.tsv")
            gfaFile = os.path.join(self.outdir, "tmp", str(i) + ".gfa")
            ogExtractBed(ogFile, extractOGFile, bedFile, self.threads)
            ogSort(extractOGFile, extractOgSortedFile, self.threads)
            ogPath(extractOgSortedFile, tsvFile, self.threads)
            ogView(extractOgSortedFile, gfaFile, self.threads)
            nodel = nodeLength(gfaFile)
            self.changeGeneTag(geneTag, os.path.join(self.outdir, "tmp"), str(i) + ".tsv")
            df = pd.read_csv(newTsvFile, delimiter="\t")
            df_merged = mergeDF(df)
            results = processDfmerged(df_merged, coreGene)
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
                for z in v.values():
                    length = 0
                    for q in z:
                        length = length + nodel[q[5:]]
                    if length / currentGeneLength < 0.2:
                        uniqueNum = uniqueNum + 1
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
            print(coreGene)
        for k, v in coreGene.items():
            with open(os.path.join(self.outdir, k + ".txt"), "w") as f:
                for i in v:
                    f.write(i + "\n")
        delete_temp_dir(os.path.join(self.outdir, "tmp"))
