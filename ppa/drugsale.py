import scipy.stats as sts
import numpy as np
import sims_pars as spr
from sims_pars.fitting import AbsObjectiveSimBased


__all__ = ['CascadeTx', 'get_prior', 'ObjForward']


# State-space
Asym = 0
NotAware = 1
NotCS = 2
ExCS = 3
N_States = 4


class CascadeTx:
    def __init__(self, d, pars):
        pr_a, pr_s, pr_c, pr_e = d['Pr_Asym'], d['Pr_NotAware'], d['Pr_NotCS'], d['Pr_NotDet']
        self.R_Die_Asym = pars['r_die_asym']
        self.R_Die_Sym = pars['r_die_sym']
        self.R_SelfCure = pars['r_sc']

        self.R_Onset = r_onset = pars['r_onset']
        # self.IncR = pars['inc']
        self.PrevUt = d['PrevUt']
        self.PrevA = d['PrevUt'] * pr_a
        self.PrevS = self.PrevA * pr_s / pr_a
        self.PrevC = self.PrevS * pr_c / pr_s
        self.PrevE = self.PrevC * pr_e / pr_c
        self.IncR = self.PrevA * (self.R_Die_Asym + self.R_SelfCure + self.R_Onset)

        mu = self.R_Die_Sym + self.R_SelfCure
        self.R_Aware = r_aware = r_onset * pr_a / pr_s - mu
        self.R_CSI = r_csi = r_aware * pr_s / pr_c - mu
        self.R_Det = r_csi * pr_c / pr_e - mu

        # self.PrevUt = self.PrevA + self.PrevS + self.PrevC + self.PrevE

        self.P_Det = np.array([
            pars['p_det_pub'],
            (1 - pars['p_det_pub']) * pars['p_eng'],
            (1 - pars['p_det_pub']) * (1 - pars['p_eng'])
        ])

        self.DetR = self.PrevE * self.R_Det * self.P_Det

        self.P_TxI = np.array([d['TxI_pub'], d['TxI_eng'], pars['TxI_pri']])
        self.PPV = np.array([pars['ppv_pub'], pars['ppv_eng'], pars['ppv_pri']])

        self.CNR = self.DetR / self.PPV
        self.TxR = self.CNR * self.P_TxI

        LTFU = np.array([d.TxLTFU_pub, d.TxLTFU_eng, d.TxLTFU_eng])
        Succ = np.array([d.TxSucc_pub, d.TxSucc_eng, d.TxSucc_eng])
        Die = np.array([d.TxDead_pub, d.TxDead_eng, d.TxDead_eng])

        self.R_Succ_Tx = np.array([pars['r_succ_pub'], pars['r_succ_pri'], pars['r_succ_pri']])
        self.R_LTFU_Tx = self.R_Succ_Tx * LTFU / Succ
        self.R_Die_Tx = self.R_Succ_Tx * Die / Succ
        self.R_Die_Tx[2] = self.R_Die_Tx[1]
        pr_ltfu_pri = self.R_LTFU_Tx[1] / (self.R_LTFU_Tx[1] + self.R_Succ_Tx[1]) * 1
        self.R_LTFU_Tx[2] = pr_ltfu_pri / (1 - pr_ltfu_pri) * (self.R_Die_Tx[2] + self.R_Succ_Tx[2])

        self.DurTx = 1 / (self.R_LTFU_Tx + self.R_Die_Tx + self.R_Succ_Tx)

        self.PrevTx = self.DetR * self.P_TxI * self.DurTx

        self.DrugTime = self.TxR * self.DurTx

        self.PrOnPriDrug = pars['p_pridrug']
        self.OnPriDrug = self.DrugTime[1] * pars['p_pridrug'] + self.DrugTime[2]

    def __call__(self, t, y):
        dy = np.zeros_like(y)

        onset = self.R_Onset * y[Asym]
        aware = self.R_Aware * y[NotAware]
        csi = self.P_Entry * self.R_CSI * y[NotCS]
        det0 = self.P_Dx0 * csi
        fn0 = (1 - self.P_Dx0) * csi

        recsi = self.R_ReCSI * y[ExCS].reshape((-1, 1)) * self.P_Tr
        det1 = self.P_Dx1 * recsi
        fn1 = (1 - self.P_Dx1) * recsi

        det = det0 + det1.sum(0)

        dy[Asym] = det.sum() - onset
        dy[NotAware] = onset - aware
        dy[NotCS] = aware - csi.sum()
        dy[ExCS] = fn0 - recsi.sum(1) + fn1.sum(0)

        mu = np.zeros(N_States)
        mu[Asym] = self.R_Die_Asym + self.R_SelfCure
        mu[NotAware] = self.R_Die_Sym + self.R_SelfCure
        mu[NotCS] = self.R_Die_Sym + self.R_SelfCure
        mu[ExCS] = self.R_Die_Sym + self.R_SelfCure
        mu *= y

        dy -= mu
        dy[Asym] += mu.sum()

        return dy

    def get_y0(self):
        return np.array([self.PrevA, self.PrevS, self.PrevC, self.PrevE])

    def epi(self):
        return {
            'PrevUt': self.PrevUt,
            'OnPriDrug': self.OnPriDrug,
            'DurTxPri': self.DurTx[2],
            'DrugTimeEng': self.DrugTime[1],
            'DrugTimePri': self.DrugTime[2],
            'PrDetPub': self.P_Det[0],
            'PrDetEng': self.P_Det[1],
            'PrDetPri': self.P_Det[2],
            'PrOnPriDrug': self.PrOnPriDrug,
            'CnrPub': self.CNR[0],
            'CnrEng': self.CNR[1],
            'PpvPri': self.PPV[2],
            'TxiPri': self.P_TxI[2]
        }

    def json(self):
        return {
            'R_SelfCure': self.R_SelfCure,
            'R_Die_Asym': self.R_Die_Asym,
            'R_Die_Sym': self.R_Die_Sym,
            'R_Onset': self.R_Onset,
            'R_Aware': self.R_Aware,
            'R_CSI': self.R_CSI,
            'R_Det': self.R_Det,
            'P_Det': self.P_Det.tolist(),
            'PPV': self.PPV.tolist(),
            'P_TxI': self.P_TxI.tolist(),

            'R_Succ_Tx': self.R_Succ_Tx.tolist(),
            'R_LTFU_Tx': self.R_LTFU_Tx.tolist(),
            'R_Die_Tx': self.R_Die_Tx.tolist(),
            'IncR': self.IncR
        }


def get_prior():
    return spr.bayes_net_from_script('''
        PCore cascade {
            r_succ_pub = 2
            r_succ_eng ~ unif(0, 26)
            rr_ltfu_pri = 1.5
            ppv_pub = 0.85
            ppv_eng = 0.85
            ppv_pri ~ unif(0.1, ppv_eng)

            r_die_ut ~ unif(0.14, 0.18)
            r_sc ~ unif(0.1, 0.3)
            rr_die_asym ~ unif(0, 1)

            r_die_sym = r_die_ut
            r_die_asym = r_die_ut * rr_die_asym
            TxI_pri ~ unif(0.5, 1)

            dur_succ_pri ~ unif(0, 18 / 12)
            r_succ_pri = 1 / dur_succ_pri

            r_onset ~ unif(0.5, 5)
            
            p_det_pub ~ unif(0, 1)
            p_eng ~ unif(0, 1)
            p_pridrug ~ unif(0.4, 0.8)
            
            inc ~ unif(0, 0.1)
        }
        ''')


class ObjForward(AbsObjectiveSimBased):
    def __init__(self, d):
        AbsObjectiveSimBased.__init__(self, get_prior())
        self.Data = d
        self.UsingDrugSale = True

    def simulate(self, pars):
        return CascadeTx(self.Data, pars)

    def link_likelihood(self, sim):
        li = 0
        li += sts.norm.logpdf(sim.CNR[0], loc=self.Data.CNR_pub, scale=self.Data.CNR_pub)
        li += sts.norm.logpdf(sim.CNR[1], loc=self.Data.CNR_eng, scale=self.Data.CNR_eng)

        scale = (self.Data.DrugTime_U - self.Data.DrugTime_L) / 4

        li += sts.norm.logpdf(sim.DrugTime[1:].sum(), loc=self.Data.PrevTxPri, scale=self.Data.PrevTxPri / 10)

        if self.UsingDrugSale:
            li += sts.norm.logpdf(sim.OnPriDrug, loc=self.Data.DrugTime_M,
                                  scale=scale)

        return li


if __name__ == '__main__':
    import pandas as pd

    prior = get_prior()

    dat = pd.read_csv('../data/cascade/d_cascade_2019.csv')
    dat = {row['State']: row for _, row in dat.iterrows()}
    d = dat['India']

    cas = CascadeTx(d, spr.sample(prior))

    print(cas.json())

    obj = ObjForward(d)
    p0 = obj.sample_prior()
    sim0 = obj.simulate(p0)
    obj.UsingDrugSale = False
    li = obj.link_likelihood(sim0)
    print(li)
    obj.UsingDrugSale = True
    li = obj.link_likelihood(sim0)
    print(li)

