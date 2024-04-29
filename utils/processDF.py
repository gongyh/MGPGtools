import os
import subprocess
import pandas as pd
from utils.odgi import *
from utils.common import delete_files
from utils.gfa import nodeLength


def processRow(row, genome_df, ref):
    """
    Process each row of gene_df for multiprocessing.
    """
    row = row[1]
    result = {}
    tmp_df = pd.DataFrame(data=row)
    # Rows with a value of 1 in column 3 and beyond correspond to nodes for that gene.
    column_indexes = (
        tmp_df[3:][tmp_df.iloc[3:] == 1].dropna(axis=0, how="all").index.tolist()
    )
    # Include the "path.name" index in the column index.
    column_indexes.insert(0, "path.name")
    # Extract columns from the genome dataframe based on the gene column.
    genome_gene_df = genome_df.loc[:, column_indexes]
    # Iterate through the path dataframe to obtain a list of nodes composed of genes that are the same as the reference gene in each genome.
    result[row["path.name"]] = {}
    for i, r in genome_gene_df.iterrows():
        if ref in str(r["path.name"]):
            continue
        columns = list(genome_gene_df.columns[2:][r[2:] == 1])
        if len(columns) == 0:
            continue
        else:
            result[row["path.name"]][r["path.name"]] = columns
    return result


def processDf(df, coreGene):
    result = {}
    # The dataframe df has genes in the first part of the rows and genome paths in the latter part of the rows.
    # Extract the rows of genome paths to form a new dataframe.
    genome_df = df[df["path.name"].str.contains("GCA|GCF") == True]
    # Extract the rows of genes to form a new dataframe.
    gene_df = df[df["path.name"].str.contains("GCA|GCF") == False]
    # Iterate through each row of the gene dataframe.
    for index, row in gene_df.iterrows():
        # tmp_df = pd.DataFrame(data=row).transpose()
        # Find all the node IDs that are 1 starting from the 3rd column (where the first two columns are path.name and path.count), 
        # which represent the nodes for that gene. Obtain the column indices formed by these nodes.
        column_indexes = row.columns[3:][row.iloc[0, 3:] == 1].tolist()
        # Add the "path.name" index to the column indices.
        column_indexes.insert(0, "path.name")
        # Extract partial columns from the genome dataframe based on the gene column. 
        # Rows with values not all equal to 0 represent mutation regions, while rows with all values equal to 0 represent non-mutation regions.
        genome_gene_df = genome_df.loc[:, column_indexes]
        genome_gene_df = genome_gene_df[(genome_gene_df.iloc[:, 1:] != 0).any(axis=1)]
        # Iterate through the genome path dataframe to obtain a list of nodes composed of genes different from the reference gene in each genome.
        for i, r in genome_gene_df.iterrows():
            columns = list(genome_gene_df.columns[2:][r[2:] == 1])
            result[row["path.name"]][r["path.name"]] = columns
    return result


def extractGenesOg(genePath, ogFile, outdir, geneTag, geneLength, genomeListExceptRef):
    """
    Extract a subgraph based on the position of each gene.
    """
    geneName = geneTag[genePath]
    gene_length = geneLength[geneName]
    # The odgi file for each gene.
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
    # The matrix for each gene.
    df = pd.DataFrame(data, columns=columns)
    cols_to_extract = df.columns[3:][df.iloc[0, 3:] == "1"].tolist()
    cols_to_extract.insert(0, "path.name")
    genomeDF = df[cols_to_extract]
    # genome_df: Record the names of the genomes that appear in the matrix.
    # variant_genome: The genomes with a total number of mutated bases on this gene greater than 0.2.
    genome_df = []
    variant_genome = []
    for index, row in genomeDF.iterrows():
        # Remove the reference genome.
        if genePath in str(row["path.name"]):
            continue
        # Genome name.
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
