# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:55:42 2022

@author: ronit
"""
import pandas as pd
import numpy as np

groupby = pd.read_table("C:\\Users\\ronit\\Downloads\\lab_rotation_2\\BFL_LVA_PMA.rbh.groupby")

orthochr = groupby.loc[groupby["count"] > 13]
orthochr = orthochr[["BFL_scaf", "LVA_scaf", "PMA_scaf"]]

orthochr_list = []
for x in orthochr.index:
    temp = orthochr.loc[x, :].values.flatten().tolist(); temp = " ".join(temp)    
    orthochr_list.append(temp)    

BFL = orthochr["BFL_scaf"].tolist(); LVA = orthochr["LVA_scaf"].tolist(); PMA = orthochr["PMA_scaf"].tolist()

groupby_filter = groupby.loc[groupby["count"] <= 13]
groupby_filter = groupby_filter.drop(["gene_group", "rbh", "alpha", "alpha_type"], axis=1)
groupby_filter = groupby_filter[["count", "BFL_scaf", "BFL_pos", "BFL_gene", "LVA_scaf", "LVA_pos", "LVA_gene", "PMA_scaf", "PMA_pos", "PMA_gene"]]

groupby_filter = groupby_filter.reset_index(drop = True)

l1 = pd.DataFrame()

iterator = 1

for (i, j, k) in zip(range(len(LVA)), range(len(PMA)), range(len(BFL))):
    groupby_filter.loc[(groupby_filter["LVA_scaf"] == LVA[i]) & (groupby_filter["PMA_scaf"] == PMA[j]), "Case"] = "1(A)"
    groupby_filter.loc[(groupby_filter["BFL_scaf"] == BFL[k]) & (groupby_filter["PMA_scaf"] == PMA[j]), "Case"] = "1(B)"
    groupby_filter.loc[(groupby_filter["BFL_scaf"] != BFL[k]) & (groupby_filter["LVA_scaf"] != LVA[i]) & (groupby_filter["PMA_scaf"] == PMA[j]), "Case"] = "2(P)"
    groupby_filter.loc[(groupby_filter["BFL_scaf"] == BFL[k]) & (groupby_filter["LVA_scaf"] == LVA[i]) & (groupby_filter["PMA_scaf"] != PMA[j]), "Case"] = "2(DP)"
    temp_list = groupby_filter["Case"].tolist()
    col_name = "Case" + str(iterator)
    l1.loc[:, col_name] = temp_list
    groupby_filter["Case"] = np.nan
    iterator += 1

l1.columns = orthochr_list

df2 = l1.stack(dropna = True)
df2 = df2.to_frame()
df2["index"] = df2.index
df2.columns = ["case", "index"]

df3 = pd.DataFrame(df2["index"].tolist()).fillna("").add_prefix("index_") 
i1 = df3["index_0"].to_list(); i2 = df3["index_1"].to_list(); i3 = df2["case"].to_list()

value_df = pd.DataFrame(list(zip(i1, i2, i3))); value_df.columns = ["row", "Ortho_Chr", "case"]

groupby_filter.reset_index(inplace = True)
value_df["count"] = np.nan
value_df['count']=value_df['row'].map(dict(zip(groupby_filter['index'],groupby_filter['count']))).fillna(value_df.count)

list_1a = []; list_1b = []

for val1 in orthochr_list:
    temp1 = value_df[value_df["Ortho_Chr"].isin([val1]) & value_df["case"].isin(["1(A)"])]
    list_1a.append(sum(temp1["count"]))

for val in orthochr_list:
    temp2 = value_df[value_df["Ortho_Chr"].isin([val]) & value_df["case"].isin(["1(B)"])]
    list_1b.append(sum(temp2["count"]))

lg_df = groupby.loc[groupby["count"] > 13]
lg_df =  lg_df.drop(["gene_group", "rbh", "alpha", "alpha_type", "BFL_pos", "BFL_gene", "LVA_pos", "LVA_gene", "PMA_pos", "PMA_gene"], axis=1)
lg_df = lg_df[["count", "BFL_scaf", "LVA_scaf", "PMA_scaf"]]
lg_df["1(A)_Count"] = list_1a
lg_df["1(B)_Count"] = list_1b
