from scipy.optimize import minimize_scalar
from scipy import linalg
import numpy as np
import sims_pars as spr
from sims_pars.fitting import AbsObjectiveSimBased
from ppa.shifting import Shifting


__all__ = ['CascadeForward', 'get_prior', 'ObjForward']


# State-space
Asym = 0
NotAware = 1
NotCS = 2
ExCS_pub = 3
ExCS_eng = 4
ExCS_pri = 5
ExCS = [ExCS_pub, ExCS_eng, ExCS_pri]
N_States = 6


class CascadeForward:
    def __init__(self, d, shf: Shifting, pars):
        shf.update(pars)
        prev = d.PrevUt * 1e-5

        self.PrevA = prev * d['Pr_Asym']
        self.PrevS = prev * d['Pr_NotAware']
        self.PrevC = prev * d['Pr_NotCS']
        self.PrevE = prev * d['Pr_NotDet']

        self.R_Die_Asym = pars['r_die_asym']
        self.R_Die_Sym = pars['r_die_sym']
        self.R_SelfCure = pars['r_sc']

        self.R_Onset = r_onset = pars['r_onset']
        self.IncR = self.PrevA * (self.R_Die_Asym + self.R_SelfCure + self.R_Onset)

        mu = self.R_Die_Sym + self.R_SelfCure
        self.R_Aware = r_aware = r_onset * self.PrevA / self.PrevS - mu
        self.R_CSI = r_csi = r_aware * self.PrevS / self.PrevC - mu
        self.R_Det = r_csi * self.PrevC / self.PrevE - mu

        cs0 = r_csi * self.PrevC * shf.P_Entry
        fn0 = cs0 * (1 - shf.P_Dx0)
        det0 = cs0 * shf.P_Dx0

        self.P_Entry, self.P_Dx0 = shf.P_Entry.copy(), shf.P_Dx0.copy()
        self.P_Tr, self.P_Dx1 = trm, pdx1 = shf.P_Tr.copy(), shf.P_Dx1.copy()

        def fn(x, trm, pdx1, mu, fn0, exp_pe):
            pe = linalg.solve(x * (trm * (1 - pdx1)).T - np.eye(3) * (x + mu), -fn0)
            return (exp_pe - pe.sum()) ** 2

        opt = minimize_scalar(fn, args=(trm, pdx1, mu, fn0, self.PrevE, ), method='bounded', bounds=(0, 50))
        self.R_ReCSI = r_recsi = opt.x

        pe = linalg.solve(r_recsi * (trm * (1 - pdx1)).T - np.eye(3) * (r_recsi + mu), -fn0)

        self.PrevE_hs = pe
        self.DetR = det0 + (r_recsi * trm * pdx1 * pe).sum(0)

        self.P_TxI = np.array([d['TxI_pub'], d['TxI_eng'], pars['TxI_pri']])
        self.PPV = np.array([pars['ppv_pub'], pars['ppv_eng'], pars['ppv_pri']])

        self.CNR = self.DetR / self.PPV
        self.TxR = self.CNR * self.P_TxI

        LTFU = np.array([d.TxLTFU_pub, d.TxLTFU_eng, d.TxLTFU_eng * pars['rr_ltfu_pri']])
        Succ = np.array([d.TxSucc_pub, d.TxSucc_eng, d.TxSucc_eng])
        Die = np.array([d.TxDead_pub, d.TxDead_eng, d.TxDead_eng])

        self.R_Succ_Tx = np.array([pars['r_succ_pub'], pars['r_succ_eng'], pars['r_succ_pri']])
        self.R_LTFU_Tx = self.R_Succ_Tx * LTFU / Succ
        self.R_Die_Tx = self.R_Succ_Tx * Die / Succ
        self.DurTx = 1 / (self.R_LTFU_Tx + self.R_Die_Tx + self.R_Succ_Tx)

        self.PrevTx = self.DetR * self.P_TxI * self.DurTx

        self.DrugTime = self.TxR * self.DurTx

        self.DataCNR = np.array([d.CNR_pub, d.CNR_eng, 0])
        self.DataPrevTx = d.PrevTx * 1e-5

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
        return np.concatenate([np.array([self.PrevA, self.PrevS, self.PrevC]), self.PrevE_hs])

    def json(self):
        return {
            'R_SelfCure': self.R_SelfCure,
            'R_Die_Asym': self.R_Die_Asym,
            'R_Die_Sym': self.R_Die_Sym,
            'R_Onset': self.R_Onset,
            'R_Aware': self.R_Aware,
            'R_CSI': self.R_CSI,
            'R_ReCSI': self.R_ReCSI,

            'P_Entry': self.P_Entry.tolist(),
            'P_Tr': self.P_Tr.tolist(),
            'P_Dx0': self.P_Dx0.tolist(),
            'P_Dx1': self.P_Dx1.tolist(),

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
            r_succ_eng = 2
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

            dur_succ_pri ~ unif(2 / 12, 12 / 12)
            r_succ_pri = 1 / dur_succ_pri

            ppm ~ unif(0, 1)
            odds_pri ~ unif(0.5, 5)
            tr_pub_pub ~ unif(0, 1)
            tr_pri_pub ~ unif(0, 1)

            r_onset ~ unif(0.5, 5)
        }
        ''')


class ObjForward(AbsObjectiveSimBased):
    def __init__(self, d, shf: Shifting):
        AbsObjectiveSimBased.__init__(self, get_prior())
        self.Data = d
        self.Shifting = shf
        self.UsingPrevTx = True
        self.UsingCNR = True

    def simulate(self, pars):
        return CascadeForward(self.Data, self.Shifting, pars)

    def link_likelihood(self, sim):
        li = 0
        if self.UsingCNR:
            li -= np.power(sim.CNR[:2] / sim.DataCNR[:2] - 1, 2).sum()
        if self.UsingPrevTx:
            li -= np.power(sim.DrugTime.sum() / sim.DataPrevTx - 1, 2).sum()
        return li


if __name__ == '__main__':
    import pandas as pd

    prior = get_prior()

    dat = pd.read_csv('../data/cascade/d_cascade_2019.csv')
    dat = {row['State']: row for _, row in dat.iterrows()}
    d = dat['India']

    src_shf = pd.read_csv('../data/shifting/shifting.csv')
    shf0 = Shifting(src_shf, d.Pr_Pub_CSI)

    cas = CascadeForward(d, shf0, spr.sample(prior))

    print(cas.json())

    obj = ObjForward(d, shf0)
    p0 = obj.sample_prior()
    sim0 = obj.simulate(p0)
    li = obj.link_likelihood(sim0)
    print(li)
