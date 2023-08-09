import pandas as pd
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['get_adr', 'get_exo']


def get_adr(filepath):
    targets = pd.read_csv(filepath)
    targets = targets[targets.Index == "IncR"]
    targets = targets[targets.Tag == "All"]
    targets = targets[targets.Year <= 2020]
    adr = - np.mean(np.diff(np.log(targets.M)))
    return adr


def get_exo(filepath):
    targets = pd.read_csv(filepath)
    txi = targets[targets.Index == "TxI"]
    txi = txi[txi.Year <= 2020]

    txi_pub = txi[txi.Tag == 'Pub']
    txi_pri = txi[txi.Tag == 'Pri']

    txi_pub = (txi_pub.M * txi_pub.N).sum() / txi_pub.N.sum()
    txi_pri = (txi_pri.M * txi_pri.N).sum() / txi_pri.N.sum()

    csi = targets[targets.Index == "PrCSIPub"].M.mean()

    return {
        'p_txi_pub': txi_pub,
        'p_txi_eng': txi_pri,
        'p_csi_pub': csi
    }
