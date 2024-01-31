import os
import csv
import subprocess
import re


def run(cmd):
    """Run executable program.

    Parameters
    ----------
    cmd : str
        Command to run.

    Returns
    -------
    boolean
        True if executed , else False.
    exception
        Program output if executed successfully, else Exception.
    """

    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
        )
        stdout, stderr = proc.communicate()
        return True, stdout, stderr
    except subprocess.CalledProcessError as e:
        return False, e


def parse_gfa(gfa):
    """Parse link and path from gfa files

    Returns:
        gfa_path: paths in gfa files
        gfa_link: links in gfa files
    """
    gfa_path = {}
    with open(gfa, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            row = line.strip().split("\t")
            if row[0] == "P":
                chrom = row[1]
                if "," in row[2]:
                    gfa_path[chrom] = row[2].split(",")
                else:
                    p = []
                    gfa_path[chrom] = p.append(row[2])
    return gfa_path


def get_chromosome(meta, database, genus):
    chrom_size = dict()
    info = get_info(meta, genus)
    genomic_path = os.path.join(
        database,
        "databases",
        "ref_genomic",
        info["class"],
        info["ref"] + "_genomic.fna.gz",
    )
    cmd = ["seqkit", "fx2tab", "-j", "20", "-l", "-n", "-i", "-H", genomic_path]
    if_success, stdout, stderr = run(cmd)
    a = stdout.strip().split("\n")
    del a[0]
    for i in a:
        r = i.strip().split("\t")
        chrom_size[r[0]] = r[1]
    return chrom_size


def get_info(meta, genus):
    """
    从Metadata.tsv中获取genus的信息
    """
    k = {
        "domain": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "segments": "",
        "links": "",
        "genomesNum": "",
        "ref": "",
    }
    with open(meta, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for i in reader:
            if i[5] == genus:
                results = [i[0], i[1], i[2], i[3], i[4], i[6], i[7], i[8], i[9]]
                info = dict(zip(k.keys(), results))
                break
    return info


def extractGenes(gff):
    geneLength = {}
    with open(gff, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            if line[0] == "#":
                continue
            row = line.strip().split("\t")
            gene = row[8].split(";")[0].replace("ID=", "")
            l = int(row[4]) - int(row[3])
            geneLength[gene] = l
    return geneLength


def segmentLength(gfa):
    sLen = {}
    with open(gfa, "r") as f:
        lines = f.read().strip().split("\n")
        for line in lines:
            r = []
            row = line.strip().split("\t")
            if row[0] == "S":
                sLen[row[1]] = len(row[2])
            if row[0] == "L" or row[0] == "P":
                break
    return sLen


def parse_ggff(gGff, gfa, gff, database, meta, genus, outdir):
    ref_nodes = []
    ref_chrom_list = []
    ref_chrom_path_list = []
    chrom_size = get_chromosome(meta, database, genus)
    details = get_info(meta, genus)
    ref_genome = details["ref"].replace(".", "#")
    # ref_chrom_path_list: path形式list
    # ref_chrom_list: chromosome name形式list
    for chrom in chrom_size.keys():
        ref_chrom_path_list.append(ref_genome + "#" + chrom)
        ref_chrom_list.append(chrom)
    gfa_path = parse_gfa(gfa)
    ## k: 染色体
    for k in ref_chrom_path_list:
        ref_nodes.extend(gfa_path[k])
    alt_path = {
        key: value for key, value in gfa_path.items() if key not in ref_chrom_path_list
    }
    segLength = segmentLength(gfa)
    genomesNum = int(details["genomesNum"])
    geneLength = extractGenes(gff)
    core_gene = {"0-15%": [], "15-95%": [], "95-99%": [], "99-100%": [], "100%": []}
    variant_gene_dict = {}
    if alt_path:
        with open(gGff, "r") as ggff:
            lines = ggff.read().strip().split("\n")
            for line in lines:
                row = line.strip().split("\t")
                geneID = row[5].split(";")[0].replace("ID=", "")
                if "," not in row[0]:
                    core_gene["100%"].append(geneID)
                else:
                    variant_gene_dict[geneID] = {}
                    range_p = row[0].split(",")
                    p = [re.sub(r"\[.*?\]", "", item) for item in range_p]
                    unsigned_p = [x[:-1] for x in p]
                    for k, v in alt_path.items():
                        length = 0
                        variant_gene_dict[geneID][k] = {}
                        # 负链
                        if int(v[0][:-1]) > int(v[-1][:-1]):
                            v = [
                                s.replace("+", "-").replace("-", "+")
                                for s in reversed(v)
                            ]
                        if (int(p[-1][:-1]) - int(p[0][:-1])) > (
                            int(v[-1][:-1]) - int(v[0][:-1])
                        ):
                            # if p[0] not in ref_nodes:
                            #     variant_gene_dict[k]["core"] = [0, 0]
                            #     continue
                            # else:
                            #     start_dist = range_p[0][range_p[0].find("[")+1: range_p[0].find("]")].split(":")
                            #     start_dist = int(start_dist[1]) - int(start_dist[0])
                            #     end_idx = 0
                            #     l = 0
                            #     for i in ref_nodes[ref_nodes.index(p[0])+1:-1]:
                            #         l = segLength[i[:-1]] + l
                            #         if l >= (geneLength[geneID] - start_dist):
                            #             end_idx = ref_nodes.index(i)
                            #             break
                            #     p_idx = find_index(reversed(ref_nodes[ref_nodes.index(p[0])+1:end_idx+1]),p)
                            #     p = p[0: p_idx+1]
                            variant_gene_dict[geneID][k]["core"] = [0, 0]
                            continue
                        if int(v[0][:-1]) >= int(p[-1][:-1]) or int(v[-1][:-1]) <= int(
                            p[0][:-1]
                        ):
                            variant_gene_dict[geneID][k]["core"] = [1, 1]
                            continue
                        if (
                            int(p[0][:-1])
                            < int(v[0][:-1])
                            < int(p[-1][:-1])
                            < int(v[-1][:-1])
                        ):
                            variant_end_index = find_index(reversed(p), v)
                            variant = v[0:variant_end_index]
                            variant_range = p[find_index(v, p) : -1]
                            ref_start_index = find_index(variant_range, ref_nodes)
                            ref_end_index = find_index(
                                reversed(variant_range), ref_nodes
                            )
                            ref_range = ref_nodes[ref_start_index : ref_end_index + 1]
                        if (
                            int(v[0][:-1])
                            < int(p[0][:-1])
                            < int(v[-1][:-1])
                            < int(p[-1][:-1])
                        ):
                            variant_start_index = find_index(p, v)
                            variant = v[variant_start_index:-1]
                            variant_range = p[0 : find_index(reversed(v), p) + 1]
                            ref_start_index = find_index(variant_range, ref_nodes)
                            ref_end_index = find_index(
                                reversed(variant_range), ref_nodes
                            )
                            ref_range = ref_nodes[ref_start_index : ref_end_index + 1]
                        if (
                            int(v[0][:-1])
                            <= int(p[0][:-1])
                            < int(p[-1][:-1])
                            <= int(v[-1][:-1])
                        ):
                            variant_start_index = find_index(p, v)
                            variant_end_index = find_index(reversed(p), v)
                            variant = v[variant_start_index : variant_end_index + 1]
                            ref_start_index = find_index(p, ref_nodes)
                            ref_end_index = find_index(reversed(p), ref_nodes)
                            ref_range = ref_nodes[ref_start_index : ref_end_index + 1]
                        unique_nodes = [x for x in ref_range if x not in variant]
                        unique_index = [unsigned_p.index(x[:-1]) for x in unique_nodes]
                        for i in unique_index:
                            dist = range_p[i][
                                range_p[i].find("[") + 1 : range_p[i].find("]")
                            ].split(":")
                            length = length + int(dist[1]) - int(dist[0])
                        if length / geneLength[geneID] < 0.15:
                            variant_gene_dict[geneID][k]["core"] = [
                                1,
                                1 - length / geneLength[geneID],
                            ]
                        else:
                            variant_gene_dict[geneID][k]["core"] = [
                                0,
                                1 - length / geneLength[geneID],
                            ]
        for g, geneInfo in variant_gene_dict.items():
            uniquePath = []
            uniqueNum = 0
            for chr, chrInfo in geneInfo.items():
                if chrInfo["core"][0] == 0:
                    if chr.split("#")[0] not in uniquePath:
                        uniqueNum = uniqueNum + 1
                        uniquePath.append(chr.split("#")[0])
            rate = (genomesNum - uniqueNum) / genomesNum
            if rate == 1:
                core_gene["100%"].append(g)
            if 0.99 <= rate < 1:
                core_gene["99-100%"].append(g)
            if 0.95 <= rate < 0.99:
                core_gene["95-99%"].append(g)
            if 0.15 <= rate < 0.95:
                core_gene["15-95%"].append(g)
            if 0 <= rate < 0.15:
                core_gene["0-15%"].append(g)
        for k, v in core_gene.items():
            print(k + ": " + str(len(v)))
            with open(outdir + k + ".txt", "w") as f:
                for i in v:
                    f.write(i + "\n")


def find_index(lst1, lst2):
    for item1 in lst1:
        for item2 in lst2:
            if item1[:-1] == item2[:-1]:
                idx = lst2.index(item2)
                break
        else:
            continue
        break
    return idx


if __name__ == "__main__":
    results = parse_ggff(
        "/mnt/me4/zhousq/core_gene_test/Moorella/gff/tmp/Moorella.ggff",
        "/mnt/me4/zhousq/pangenome-1.0/databases/gfa/Moorellia/Moorella.gfa",
        "/mnt/me4/zhousq/core_gene_test/Moorella/gff/tmp/GCF_001267405.1_genomic.gff",
        "/mnt/me4/zhousq/pangenome-1.0",
        "/mnt/me4/zhousq/panTools/meta/Metadata.tsv",
        "Moorella",
        "/mnt/me4/zhousq/core_gene_test/Moorella/",
    )
    for k, v in results.items():
        print(k + " : " + str(len(v)))
