import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance
import numpy as np
import Bio.PDB
import Bio.AlignIO as al
from sklearn import preprocessing


def read_pdb_starts():
    pdb_starts_dict = {}
    pdb_starts_dict["4kvm_chainsBE"] = 0
    pdb_starts_dict["4kvx_chainB"] = 0
    return pdb_starts_dict


def read_msa_fasta():
    """
    Reads multiple structure alignment from MUSTANG. It determines the structurally aligned core of the proteins.
    :return:
    """
    pdb_align_dict = {'4kvm_chainsBE': [], '4kvx_chainB': []}
    # file_path = os.path.join("../data/input/etc", "nats_alignment.afasta")
    file_path = os.path.join(".", "fasta_rmsd.afasta")
    records = al.read(open(file_path), "fasta")
    tlist = list(zip(*records))
    for i in range(0, records.get_alignment_length()):
        if '-' not in [y for y in tlist][i]:
            for rec in records:
                if not rec.id[0:4] == '4ua3':
                    ls = [i for i, e in enumerate(rec.seq) if e != '-']
                    res_cpt = ls.index(i)
                    pdb_align_dict[rec.id[:-4]].append(res_cpt + read_pdb_starts()[rec.id[:-4]])
    return pdb_align_dict


def read_msa_fasta_gaps():
    """
    Reads msa from MUSTANG in webnma3. It conserves the nan and puts them in the returned dictionary.
    :return:
    """
    pdb_align_dict = {'4kvm_chainsBE': [], '4kvx_chainB': []}
    file_path = os.path.join(".", "alignment_webnma3.afasta")
    records = al.read(open(file_path), "fasta")
    tlist = list(zip(*records))
    for i in range(0, records.get_alignment_length()):
        for rec in records:
            ls = [i for i, e in enumerate(rec.seq)]
            if not rec.seq[i] == "-":
                res_cpt = ls.index(i) + 1
                pdb_align_dict[rec.id[:-4]].append(res_cpt + read_pdb_starts()[rec.id[:-4]])
            else:
                pdb_align_dict[rec.id[:-4]].append(np.nan)
    return pdb_align_dict


# Global variables (Ugly)
ca_align_list = []
ca_align_dict = read_msa_fasta_gaps()

if __name__ == '__main__':

    resid_list = []
    resid_list_BEB = []
    fluct_list_BE = []
    fluct_list_B = []
    fluct_list_BEB = []
    train = pd.DataFrame(columns=['resid', 'fluct_score_BE', 'fluct_score_B'])

    with open('alignment_fluctuations_webnma3.txt', 'r') as f1:
        for line in f1.readlines():
            if not line.startswith("index") and not line.startswith("\n"):
                resid = int(line.split("\t")[0].strip())
                # if resid in ca_align_dict["4kvx_chainB"]:
                resid_list.append(resid)
                if math.isnan(ca_align_dict["4kvx_chainB"][resid - 1]):
                    fluct_list_B.insert(resid, math.nan)
                else:
                    fluct_list_B.append(float(line.split("\t")[1].strip().replace(",", ".")))
                fluct_list_BE.append(float(line.split("\t")[2].strip().replace(",", ".")))

    with open('4kvm_chainsBE_modesv2_fluct.dat', 'r') as f1:
        for line in f1.readlines():
            resid = int(line.split("\t")[0].strip())
            if not line.startswith("index") and not line.startswith("\n"):
                # if resid in ca_align_dict["4kvm_chainsBE"]:
                # if math.isnan(ca_align_dict["4kvm_chainsBE"][718 + resid - 1]):
                #     # fluct_list_BEB.insert(718 + resid - 1, math.nan)
                # else:
                fluct_list_BEB.append(float(line.split("\t")[1].strip().replace(",", ".")))

    train['resid'] = resid_list
    train['fluct_score_BE'] = fluct_list_BE
    train['fluct_score_B'] = fluct_list_B

    # Truncate the dataframe from when fluct_score_B starts
    train = train.truncate(before=train['fluct_score_B'].first_valid_index())

    train['fluct_score_BEB'] = fluct_list_BEB

    ## min max normalization
    # column_names_to_normalize = ['fluct_score_B', 'fluct_score_BE']
    # x = train[column_names_to_normalize].values
    # min_max_scaler = preprocessing.MinMaxScaler()
    # x_scaled = min_max_scaler.fit_transform(x)
    # df_temp = pd.DataFrame(x_scaled, columns=column_names_to_normalize, index=train.index)
    # train[column_names_to_normalize] = df_temp

    # train = (train - train.min()) / (train.max() - train.min())

    # sns.lineplot("resid", "fluct_score_BE", label="4kvm_chainsBE", data=train)
    # sns.lineplot("resid", "fluct_score_B", label="4kvx_chainsB", data=train)
    # plt.show()

    # Detect changes in alignment on chain B
    # train['C'] = train['fluct_score_B'].diff()

    ax = plt.gca()

    train.plot(kind='line', x='resid', y='fluct_score_BEB', ax=ax)
    train.plot(kind='line', x='resid', y='fluct_score_B', color='red', ax=ax)

    plt.show()
