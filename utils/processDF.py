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


def processDfmerged(df_merged, coreGene):
    result = {}
    genome_rows = df_merged[
                df_merged["path.name"].str.contains("GCA|GCF") == True
            ]
    genome_df = pd.DataFrame(genome_rows)
    gene_rows = df_merged[
                df_merged["path.name"].str.contains("GCA|GCF") == False
            ]
    gene_df = pd.DataFrame(gene_rows)
    for index, row in gene_df.iterrows():
        tmp_df = pd.DataFrame(data=row).transpose()
        column_indexes = tmp_df.columns[3:][tmp_df.iloc[0, 3:] == 1].tolist()
        column_indexes.insert(0, "path.name")
        genome_gene_df = genome_df.loc[:, column_indexes]
        if (genome_gene_df.iloc[:, 1:] == 0).all().all():
            coreGene["100%"].append(row["path.name"])
            continue
        result[row["path.name"]] = {}
        for i, r in genome_gene_df.iterrows():
            columns = (list(genome_gene_df.columns[2:][r[2:] == 0]))
            if r.iloc[1] == 0:
                columns = columns[1:]
            if r.iloc[-1] == 0:
                columns = columns[:-1]
            result[row["path.name"]][r["path.name"]] = columns
    return result


def mergeDF(df):
    df_merged = (
        df.groupby(
            df.iloc[:, 0].str.split("#").str[:2].str.join("#")
        )
        .apply(merge_rows)
        .reset_index(drop=True)
    )
    return df_merged


def merge_rows(row_group):
    merged_row = row_group.iloc[0].copy()
    merged_row.iloc[0] = "#".join(row_group.iloc[0].iloc[0].split("#")[:2])
    for i in range(1, len(row_group.columns)):
        merged_row.iloc[i] = row_group.iloc[:, i].sum()
    return merged_row
