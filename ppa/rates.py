import numpy as np
import sims_pars as spr
from sims_pars.fitting import AbsObjectiveSimBased
import scipy.stats as sts

__author__ = 'Chu-Chang Ku'
__all__ = ['Rates', 'ObjRates', 'get_prior']


class Rates:
    def __init__(self, d, pars):
        d = dict(d)
        d.update(pars)

        self.Pars = d

        ltfu = np.array([d['TxLTFU_pub'], d['TxLTFU_eng'], d['TxLTFU_eng'] * d['rr_ltfu_pri']])
        die = np.array([d['TxDead_pub'], d['TxDead_eng'], d['TxDead_eng']])
        succ = np.array([d['TxSucc_pub'], d['TxSucc_eng'], d['TxSucc_eng']])

        self.R_Succ_Tx = np.array([d['r_succ_pub'], d['r_succ_eng'], 1 / d['dur_succ_pri']])
        self.R_LTFU_Tx = self.R_Succ_Tx * ltfu / succ
        self.R_Die_Tx = self.R_Succ_Tx * die / succ
        self.DurTx = 1 / (self.R_LTFU_Tx + self.R_Die_Tx + self.R_Succ_Tx)

        self.DetR = np.array([d['CNR_pub'], d['CNR_eng'], d['det_pri']])
        self.P_TxI = np.array([d['TxI_pub'], d['TxI_eng'], d['TxI_pri']])
        self.TxR = self.DetR * self.P_TxI

        self.PrevTx = self.TxR / (self.R_LTFU_Tx + self.R_Die_Tx + self.R_Succ_Tx)

        self.PrevUt = prev = sts.binom.rvs(n=d['Pop'], p=d['PrevUt'] * 1e-5) / d['Pop']
        self.PrevA = prev * d['Pr_Asym']
        self.PrevS = prev * d['Pr_NotAware']
        self.PrevC = prev * d['Pr_NotCS']
        self.PrevE = prev * d['Pr_NotDet']

        self.R_SelfCure = d['r_sc']
        self.R_Die_Ut = d['r_die_ut']

        self.R_Die_Asym = d['r_die_ut'] * d['rr_die_asym']
        self.R_Die_Sym = d['r_die_ut']

        r_mu_sym = self.R_SelfCure + self.R_Die_Sym

        self.R_Det = r_det = self.DetR.sum() / (prev * d['Pr_NotDet'])
        self.R_CSI = r_cs = d['Pr_NotDet'] / d['Pr_NotCS'] * (r_mu_sym + r_det)
        self.R_Aware = r_aware = d['Pr_NotCS'] / d['Pr_NotAware'] * (r_mu_sym + r_cs)
        self.R_Onset = d['Pr_NotAware'] / d['Pr_Asym'] * (r_mu_sym + r_aware)

        self.PPV = np.array([d['ppv_pub'], d['ppv_eng'], d['ppv_pri']])

        self.PrevTxAll = self.PrevTx / self.PPV
        self.ExPrevTx = sts.binom.rvs(n=d['Pop'], p=d['PrevTx'] * 1e-5) / d['Pop']

        self.PrPri = self.PrevTxAll[1:].sum() / self.PrevTxAll.sum()
        self.PrPri0 = self.PrevTxAll[1].sum() / self.PrevTxAll[:2].sum()

    def summary(self, key=0):
        su = {
            'Key': key,
            'PrevTxPub': self.PrevTx[0],
            'PrevTxEng': self.PrevTx[1],
            'PrevTxPri': self.PrevTx[2],
            'DurTxPub': self.DurTx[0],
            'DurTxEng': self.DurTx[1],
            'DurTxPri': self.DurTx[2],
            'PrPri': self.PrPri,
            'PrPri0': self.PrPri0,
            'PriVsPub': self.PrevTx[0] / self.PrevTx[1:].sum()
        }

        su.update(self.Pars)
        return su

    def to_json(self):
        return {
            'R_Succ_Tx': self.R_Succ_Tx,
            'R_LTFU_Tx': self.R_LTFU_Tx,
            'R_Die_Tx': self.R_Die_Tx,
            'R_SelfCure': self.R_SelfCure,
            'R_Die_Asym': self.R_Die_Asym,
            'R_Die_Sym': self.R_Die_Sym,
            'R_Onset': self.R_Onset,
            'R_Aware': self.R_Aware,
            'R_CSI': self.R_CSI,
            'R_Det': self.R_Det,
            'PrevUt': self.PrevUt,
            'PrevTx': self.PrevTx,
            'DetR': self.DetR,
            'P_TxI': self.P_TxI,
            'PPV': self.PPV
        }


def get_prior():
    return spr.bayes_net_from_script('''
    PCore prior {
        r_succ_pub = 2
        r_succ_eng = 2
        rr_ltfu_pri = 1.5
        ppv_pub = 0.85
        ppv_eng = 0.85
        ppv_pri ~ unif(0.2, ppv_eng)
        
        r_die_ut ~ unif(0.14, 0.18)
        r_sc ~ unif(0.1, 0.3)
        rr_die_asym ~ unif(0, 1)
        TxI_pri ~ unif(0.5, 1)
        
        dur_succ_pri ~ unif(2 / 12, 9 / 12)
        det_pri ~ unif(0.00001, 0.01)
    }
    ''')


class ObjRates(AbsObjectiveSimBased):
    def __init__(self, d):
        AbsObjectiveSimBased.__init__(self, get_prior())
        self.Data = d

    def simulate(self, pars):
        return Rates(self.Data, pars)

    def link_likelihood(self, sim):
        x = sim.PrevTxAll.sum()
        mu = sim.ExPrevTx
        return sts.norm.logpdf(x, loc=mu, scale=mu)


if __name__ == '__main__':
    import pandas as pd

    bn = get_prior()

    dat = pd.read_csv('../data/cascade/d_cascade_2019.csv').set_index('State')
    dat = {i: row.to_dict() for i, row in dat.iterrows()}
    dat = dat['Tamil Nadu']

    p0 = spr.sample(bn)

    cas = Rates(dat, p0)

    for k, v in cas.to_json().items():
        print(k, v)

    obj = ObjRates(dat)
    sim = obj.simulate(obj.sample_prior())
    print(obj.link_likelihood(sim))
