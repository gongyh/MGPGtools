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
    for index, row in df_merged.iterrows():
        if "GCF" not in row["path.name"] and "GCA" not in row["path.name"]:
            temp_df = pd.DataFrame(data=row).transpose()
            other_rows = df_merged[
                df_merged["path.name"].str.contains("GCA|GCF") == True
            ]
            new_df = pd.DataFrame(other_rows)
            for column in temp_df.columns[3:]:
                if row[column] == 0:
                    new_df.drop(column, axis=1, inplace=True)
            all_zeros = new_df.iloc[:, 3:].apply(lambda x: (x == 0).all(), axis=1)
            if all_zeros.all():
                coreGene["100%"].append(row["path.name"])
                continue
            result[row["path.name"]] = {}
            for i, r in new_df.iterrows():
                columns = [column for column, value in r.iloc[3:].items() if value == 0]
                if r.iloc[3] == 0:
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
