import pandas as pd
import logomaker
import plotly.express as px
from ipywidgets import Output, VBox
from IPython.display import display, clear_output
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import logomaker
import copy
import random

def extract_sufficient_data(data_merged):
    allele_counts = pd.DataFrame({"count": data_merged.groupby(["allele_type", "allele", "peptide", "binder"]).size().groupby(["allele_type", "allele", "binder"]).size()}).reset_index()
    alleles = allele_counts["allele"].unique()
    alleles_keep = []
    for allele in alleles:
        b_cnt = list(allele_counts[(allele_counts["allele"]==allele) & (allele_counts["binder"]==1)]["count"])
        nb_cnt = list(allele_counts[(allele_counts["allele"]==allele) & (allele_counts["binder"]==0)]["count"])
        #have both binder and non-binder representatives
        if len(b_cnt)>0 and len(nb_cnt)>0:
            #have at least 1000 examples
            if b_cnt[0]+nb_cnt[0] >=1000: 
                #have at least 150 of binders/non-binders
                if b_cnt[0]>=200 and nb_cnt[0]>=200: alleles_keep.append(allele)
    return data_merged[data_merged["allele"].isin(alleles_keep)]

def drop_duplicates(data):
    sum_pairs = pd.DataFrame({"count": data.groupby(['allele','peptide']).size()}).reset_index()
    duplicates = sum_pairs[sum_pairs["count"] > 1]
    duplicates_names = (duplicates["allele"]+"-"+duplicates["peptide"]).to_list()
    return data[~(data["allele"]+"-"+data["peptide"]).isin(duplicates_names)]

def utils_logo_event(trace, points, selector, fig, labels, axs):
    allele = fig.data[0].labels[points.point_inds[0]]
    hla_peptides = []
    if allele[:3]=="HLA":
        hla = labels[labels['allele type']==allele]
        hla_b = hla[hla['binder']==1]
        hla_peptides = hla_b["peptide"].to_numpy()
    else:
        hla = labels[labels['allele']==allele]
        hla_b = hla[hla['binder']==1]
        hla_peptides = hla_b["peptide"].to_numpy()
        
    if len(hla_peptides) == 0:
        print("No known binders in dataset!")
    else:
        crp_counts_df =logomaker.alignment_to_matrix(sequences=hla_peptides, to_type='counts')
        clear_output()
        axs[0].clear()
        axs[1].clear()
        axs[2].clear()
        logomaker.Logo(crp_counts_df, ax = axs[0], color_scheme="charge")
        logomaker.Logo(crp_counts_df, ax = axs[1], color_scheme="hydrophobicity")
        logomaker.Logo(crp_counts_df, ax = axs[2])
        axs[0].set_title("Charge logo for allele: "+str(allele))
        axs[1].set_title("Hydrophobicity logo for allele: "+str(allele))
        axs[2].set_title("Default logo for allele: "+str(allele))
        display(axs[0].figure)
        
def create_split(full_dataset, pr_b, pr_nb):
    split = {"allele":[], "train_b":[], "train_nb":[], "test_b":[], "test_nb":[]}
    alleles = full_dataset["allele"].unique()
    for allele in alleles:
        allele_data = full_dataset[full_dataset["allele"]==allele]
        count_all = allele_data.shape[0]
        y = pd.DataFrame(allele_data.groupby(["allele", "binder"]).count().unstack(fill_value=0).stack().reset_index())
        nonbinder_cnt = list(y[y["binder"]==0]["peptide"])[0]
        binder_cnt = list(y[y["binder"]== 1]["peptide"])[0]

        test_size = 0.2*count_all
        binder_test_size = pr_b*binder_cnt 
        nonbinder_test_size = pr_nb*nonbinder_cnt 

        print("----------------------------------")
        print("Allele: "+allele)
        print("Number of peptides (train/test): "+str(count_all-test_size)+"/"+str(test_size))
        print("Number of binders (train/test): "+str(binder_cnt-binder_test_size)+"/"+str(binder_test_size))
        print("Number of nonbinders (train/test): "+str(nonbinder_cnt-nonbinder_test_size)+"/"+str(nonbinder_test_size))
        split["allele"].append(allele)
        split["train_b"].append(binder_cnt-binder_test_size)
        split["test_b"].append(binder_test_size)
        split["train_nb"].append(nonbinder_cnt-nonbinder_test_size)
        split["test_nb"].append(nonbinder_test_size)
        
        #reformat split
        reformat_split = {"allele":[], "cnt":[], "type":[]}

        reformat_split["allele"] = copy.deepcopy(split["allele"])
        reformat_split["allele"].extend(split["allele"])
        reformat_split["allele"].extend(split["allele"])
        reformat_split["allele"].extend(split["allele"])

        reformat_split["cnt"] = copy.deepcopy(split["train_b"])
        reformat_split["cnt"].extend(split["train_nb"])
        reformat_split["cnt"].extend(split["test_b"])
        reformat_split["cnt"].extend(split["test_nb"])

        reformat_split["type"] = ["train_b" for e in split["allele"]] 
        reformat_split["type"].extend(["train_nb" for e in split["allele"]])
        reformat_split["type"].extend(["test_b" for e in split["allele"]]) 
        reformat_split["type"].extend(["test_nb" for e in split["allele"]]) 
        
    return (split, pd.DataFrame(reformat_split))
    
    
def apply_split(full_dataset, split_df):
    alleles = full_dataset["allele"].unique()
    for i, allele in enumerate(alleles):
        allele_data = full_dataset[full_dataset["allele"]==allele]
        allele_data = allele_data.sort_values(by='ba', ascending=False)
        allele_data_b = allele_data[allele_data["binder"]==1]
        allele_data_nb = allele_data[allele_data["binder"]==0]
        allele_data_b = allele_data_b.reset_index()
        allele_data_nb = allele_data_nb.reset_index()

        allele_cnt = split_df[split_df["allele"]==allele]
        train_b_cnt = allele_cnt[allele_cnt["type"]=="train_b"]
        train_nb_cnt = allele_cnt[allele_cnt["type"]=="train_nb"]
        tbc = list(train_b_cnt["cnt"])[0]
        tnbc = list(train_nb_cnt["cnt"])[0]
        train_indexes_b = random.sample(range(len(allele_data_b)), int(tbc))
        train_indexes_nb = random.sample(range(len(allele_data_nb)), int(tnbc))

        #set(df.index) - set(blacklist)
        train_b = allele_data_b[allele_data_b.index.isin(train_indexes_b)]
        train_nb = allele_data_nb[allele_data_nb.index.isin(train_indexes_nb)]
        train_b = train_b.reset_index()
        train_nb = train_nb.reset_index()

        if i == 0:
            train_set = pd.concat([train_b, train_nb])
        else:
            train_set = pd.concat([train_set, train_b])
            train_set = pd.concat([train_set, train_nb])

        test_b = allele_data_b[~ allele_data_b.index.isin(train_indexes_b)]
        test_nb = allele_data_nb[~ allele_data_nb.index.isin(train_indexes_nb)]

        if i == 0:
            test_set = pd.concat([test_b, test_nb])
        else:
            test_set = pd.concat([test_set, test_b])
            test_set = pd.concat([test_set, test_nb])

    return train_set, test_set
    
    
def update_fold(train_set, fold_n, fold_b, fold_nb):
    pep_al_pairs = []
    for index, row in fold_b.iterrows():
        pep_al_pairs.append((row["peptide"], row["allele"]))
    for index, row in fold_nb.iterrows():
        pep_al_pairs.append((row["peptide"], row["allele"]))
    train_set.loc[train_set.apply(lambda x: (x.peptide, x.allele) in pep_al_pairs, axis=1), "fold_num"] = fold_n
    return train_set       

def reformat_stats_ds2(alleles, b_stats, nb_stats):
    # reformat
    res_df = {"allele": [], "count": [], "binder":[]}
    for allele in alleles:
        res_df["allele"].append(allele)
        res_df["allele"].append(allele)
        res_df["count"].append(list(b_stats[b_stats["allele"]==allele]["binder_cnt"])[0]) 
        if allele[0] == "C": 
            res_df["count"].append('0')
        else:
            res_df["count"].append(list(nb_stats[nb_stats["allele"]==allele]["nonbinder_cnt"])[0])
        res_df["binder"].append('1')
        res_df["binder"].append('0')
    return pd.DataFrame(res_df)  
      