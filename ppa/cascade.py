from pydantic import BaseModel, Field
import numpy as np
from scipy.optimize import minimize_scalar
from scipy import linalg
from ppa.rates import Rates
from ppa.shifting import Shifting

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade', 'bind_cascade']


class Cascade(BaseModel):
    IncR: float
    R_SelfCure: float
    R_Die_Asym: float
    R_Die_Sym: float
    R_Die_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    R_Onset: float
    R_Aware: float
    R_CSI: float
    R_ReCSI: float

    P_Entry: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    P_Tr: np.ndarray = Field(default_factory=lambda: np.zeros((3, 3)))
    P_Dx0: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    P_Dx1: np.ndarray = Field(default_factory=lambda: np.zeros((3, 3)))

    P_TxI: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    R_Succ_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    R_LTFU_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    PPV: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def parse_json(rs):
        rs = {k: (np.array(v) if isinstance(v, list) else v) for k, v in rs.items()}
        return Cascade.parse_obj(rs)

    def match_targets(self):
        p_a = self.IncR / (self.R_Onset + self.R_Die_Asym + self.R_SelfCure)
        p_s = (self.R_Onset * p_a) / (self.R_Aware + self.R_Die_Sym + self.R_SelfCure)
        p_c = (self.R_Aware * p_s) / (self.R_CSI + self.R_Die_Sym + self.R_SelfCure)

        cs0 = self.R_CSI * p_c * self.P_Entry
        fn0 = cs0 * (1 - self.P_Dx0)
        det0 = cs0 * self.P_Dx0

        mu = self.R_ReCSI + self.R_Die_Sym + self.R_SelfCure

        p_e_hs = linalg.solve(self.R_ReCSI * (self.P_Tr * (1 - self.P_Dx1)).T - np.eye(3) * mu, -fn0)
        p_e = p_e_hs.sum()

        det1 = (self.R_ReCSI * self.P_Tr * self.P_Dx1 * p_e_hs).sum(0)
        det = det0 + det1

        det_all = det / self.PPV
        tx = det_all * self.P_TxI

        dur = 1 / (self.R_LTFU_Tx + self.R_Succ_Tx + self.R_Die_Tx)

        return {
            'IncR': self.IncR,
            'PrevAsym': p_a,
            'PrevNotAware': p_s,
            'PrevNotCS': p_c,
            'PrevNotDet': p_e,
            'Prev': p_a + p_s + p_c + p_e,
            'CNR_Pub': det_all[0],
            'CNR_Eng': det_all[1],
            'Det_Pri': det[2],
            'PrTxi_pri': self.P_TxI[2],
            'Dur_Pub': dur[0],
            'Dur_Eng': dur[1]
        }


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
