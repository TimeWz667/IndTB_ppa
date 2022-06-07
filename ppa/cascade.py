from pydantic import BaseModel
import numpy as np
from scipy.optimize import minimize_scalar
from scipy import linalg
from ppa.rates import Rates
from ppa.shifting import Shifting

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade', 'bind_cascade']


class Cascade(BaseModel):
    R_SelfCure: float = 0
    R_Die_Asym: float = 0
    R_Die_Sym: float = 0
    R_Die_Tx: list = 0

    R_Onset: float
    R_Aware: float
    R_CSI: float
    R_ReCSI: float

    P_Entry: list
    P_Tr: list
    P_Dx0: list
    P_Dx1: list

    P_TxI: list

    R_Succ_Tx: list
    R_LTFU_Tx: list

    PrevUt: float
    PrevTx: list
    DetR: list

    PPV: list

    def revive(self):
        self.P_Entry = np.array(self.P_Entry)
        self.P_Tr = np.array(self.P_Tr)
        self.P_Dx0 = np.array(self.P_Dx0)
        self.P_Dx1 = np.array(self.P_Dx1)

        self.R_Succ_Tx = np.array(self.R_Succ_Tx)
        self.R_LTFU_Tx = np.array(self.R_LTFU_Tx)
        self.R_Die_Tx = np.array(self.R_Die_Tx)

        self.PrevTx = np.array(self.PrevTx)
        self.DetR = np.array(self.DetR)
        self.P_TxI = np.array(self.P_TxI)

        self.PPV = np.array(self.PPV)

    def calc_prop_det(self, tol=1e-8):
        det = self.P_Entry * self.P_Dx0
        x = self.P_Entry * (1 - self.P_Dx0)

        while x.sum() > tol:
            tr0 = self.P_Tr * x.reshape((-1, 1))
            tp, fp = tr0 * self.P_Dx1, tr0 * (1 - self.P_Dx1)

            det += tp.sum(0)

            x = fp.sum(0)

        return det.reshape(-1) / det.sum()


def find_recsi(shf: Shifting, rates: Rates):
    cs0 = rates.R_CSI * rates.PrevC * shf.P_Entry
    fn0 = cs0 * (1 - shf.P_Dx0)

    trm, pdx1 = shf.P_Tr, shf.P_Dx1
    mu = rates.R_Die_Sym + rates.R_SelfCure

    def fn(x, trm, pdx1, mu, fn0, exp_pe):
        pe = linalg.solve(x * (trm * (1 - pdx1)).T - np.eye(3) * (x + mu), -fn0)
        return (exp_pe - pe.sum()) ** 2

    opt = minimize_scalar(fn, args=(trm, pdx1, mu, fn0, rates.PrevE, ), method = 'bounded', bounds=(0, 50))
    r_recsi = opt.x
    if r_recsi < 0:
        raise ValueError

    return r_recsi[0]


def bind_cascade(shf: Shifting, rates: Rates):
    r_recsi = find_recsi(shf, rates)

    d = rates.to_json()
    for k in ['R_Succ_Tx', 'R_LTFU_Tx', 'R_Die_Tx', 'PPV', 'DetR', 'PrevTx', 'P_TxI']:
        d[k] = d[k].tolist()

    d['P_Tr'] = shf.P_Tr.tolist()
    d['P_Entry'] = shf.P_Entry.tolist()
    d['P_Dx0'] = shf.P_Dx0.tolist()
    d['P_Dx1'] = shf.P_Dx1.tolist()
    d['R_ReCSI'] = r_recsi

    return Cascade.parse_obj(d)
