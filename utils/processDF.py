import pandas as pd


# def processRow(df_merged, row):
#     result = {}
#     if "#" not in row:
#         result[row["path.name"]] = {}
#         temp_df = pd.DataFrame(data=row).transpose()
#         other_rows = df_merged[
#             df_merged["path.name"].str.contains("#")
#             == True
#         ]
#         new_df = pd.DataFrame(other_rows)
#         for column in temp_df.columns[3:]:
#             if row[column] == 0:
#                 new_df.drop(column, axis=1, inplace=True)
#         for i, r in new_df.iterrows():
#             columns = [
#                 column for column, value in r.iloc[3:].items() if value == 0
#             ]
#             result[row["path.name"]][r["path.name"]] = columns
#     return result


def processDf(df, coreGene):
    result = {}
    # df是一个dataframe，前部分行是gene，后部分行是genome paths
    # 提取genome path的行构成新的dataframe
    genome_rows = df[
                df["path.name"].str.contains("GCA|GCF") == True
            ]
    genome_df = pd.DataFrame(genome_rows)
    # 提取gene行构成新的dataframe
    gene_rows = df[
                df["path.name"].str.contains("GCA|GCF") == False
            ]
    gene_df = pd.DataFrame(gene_rows)
    # 遍历gene dataframe的每一行
    for index, row in gene_df.iterrows():
        tmp_df = pd.DataFrame(data=row).transpose()
        # 找到从第3列开始(第三列开始是node，前两列分别是path.name, path.conut)所有是1的nodeid，即是该基因的node，获得这些node构成的列索引
        column_indexes = tmp_df.columns[3:][tmp_df.iloc[0, 3:] == 1].tolist()
        # 列索引加入"path.name"索引
        column_indexes.insert(0, "path.name")
        # genome dataframe根据基因的列提取部分列，不全为0的行是突变区域，全为0的行是非突变区域
        genome_gene_df = genome_df.loc[:, column_indexes]
        genome_gene_df = genome_gene_df[(genome_gene_df.iloc[:, 1:] != 0).any(axis=1)] 
        # 如果所有的基因组path和node都是0，说明不在突变区域，是核心基因
        if (genome_gene_df.iloc[:, 1:] == 0).all().all():
            coreGene["100%"].append(row["path.name"])
            continue
        result[row["path.name"]] = {}
        # 遍历基因组path dataframe，获得每个基因组与reference gene不同的node组成的list
        for i, r in genome_gene_df.iterrows():
            columns = (list(genome_gene_df.columns[2:][r[2:] == 0]))
            result[row["path.name"]][r["path.name"]] = columns
    return result

# def mergeDF(df):
#     df_merged = (
#         df.groupby(
#             df.iloc[:, 0].str.split("#").str[:2].str.join("#")
#         )
#         .apply(merge_rows)
#         .reset_index(drop=True)
#     )
#     return df_merged


# def merge_rows(row_group):
#     merged_row = row_group.iloc[0].copy()
#     merged_row.iloc[0] = "#".join(row_group.iloc[0].iloc[0].split("#")[:2])
#     for i in range(1, len(row_group.columns)):
#         merged_row.iloc[i] = row_group.iloc[:, i].sum()
#     return merged_row
