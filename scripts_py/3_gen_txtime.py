import json
import pandas as pd
from pydantic import BaseModel, Field
import numpy as np
from scipy import linalg
import numpy.random as rd
import os


class CascadeRates(BaseModel):
    R_Onset: float = 0
    R_Aware: float = 0
    R_CSI: float = 0
    R_ReCSI: float = 0

    P_Entry: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    P_Tr: np.ndarray = Field(default_factory=lambda: np.zeros((3, 3)))
    P_Dx0: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    P_Dx1: np.ndarray = Field(default_factory=lambda: np.zeros((3, 3)))

    P_TxI: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    R_Succ_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))
    R_LTFU_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    R_SelfCure: float = 0
    R_Die_Asym: float = 0
    R_Die_Sym: float = 0
    R_Die_Tx: np.ndarray = Field(default_factory=lambda: np.zeros(3))

    IncR: float = 0
    PPV: np.ndarray = 0

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def parse_json(js):
        js = {k: np.array(v) if isinstance(v, list) else v for k, v in js.items()}
        return CascadeRates.parse_obj(js)

    def simulate(self):
        mu = self.R_Die_Sym + self.R_SelfCure

        prev_a = self.IncR / (self.R_Die_Asym + self.R_SelfCure + self.R_Onset)
        prev_s = (self.R_Onset * prev_a) / (self.R_Aware + mu)
        prev_c = (self.R_Aware * prev_s) / (self.R_CSI + mu)

        cs0 = self.R_CSI * prev_c * self.P_Entry
        fn0 = cs0 * (1 - self.P_Dx0)
        det0 = cs0 * self.P_Dx0
        trm, pdx1 = self.P_Tr, self.P_Dx1

        prev_e = linalg.solve(self.R_ReCSI * (trm * (1 - pdx1)).T - np.eye(3) * (self.R_ReCSI + mu), -fn0)

        det = det0 + (self.R_ReCSI * trm * pdx1 * prev_e).sum(0)

        r_to = self.R_LTFU_Tx + self.R_Succ_Tx + self.R_Die_Tx
        # r_to[2] /= 1.5
        ppv = self.PPV.copy()
        # ppv[2] /= 1.5

        cnr = det / ppv
        txi = cnr * self.P_TxI

        dur_tx = 1 / r_to

        prev_t = det * self.P_TxI * dur_tx

        prev_pos = txi * dur_tx

        return {
            'ratio_pe_p': prev_pos[2:].sum() / prev_pos[:2].sum(),
            'ratio_p_ep': prev_pos[1:].sum() / prev_pos[:1].sum(),
            'ppv_pri': ppv[2],
            'dur_tx_pri': dur_tx[2],
            'txi_pri': txi[2],
            'cnr_pub': cnr[0],
            'cnr_eng': cnr[1],
            'cnr_pri': cnr[2],
            'prev_a': prev_a,
            'prev_s': prev_s,
            'prev_c': prev_c,
            'prev_e': prev_e.sum(),
            'prev_ut': prev_a + prev_s + prev_c + prev_e.sum(),
            'prev_pos_pri': prev_pos[2],
            'prev_pos_eng': prev_pos[1],
            'prev_pos_pub': prev_pos[1],
            'prev_pos': prev_pos.sum()
        }


mss = list()

for src in os.listdir("../out/cf_pars"):
    loc = src.replace('pars_', '').replace('.json', '')
    print(loc)

    with open(f'../out/cf_pars/{src}', 'r') as f:
        pars = json.load(f)

    incs = np.array([p['IncR'] for p in pars])
    rd.shuffle(incs)
    for inc, p in zip(incs, pars):
        p['IncR'] = inc

    for i, p in enumerate(pars):
        sim = CascadeRates.parse_json(p).simulate()
        sim['Location'] = loc
        sim['Key'] = i
        sim['Target'] = 'CNR+PrevUt'
        mss.append(sim)

for src in os.listdir("../out/cft_pars"):
    loc = src.replace('pars_', '').replace('.json', '')

    with open(f'../out/cft_pars/{src}', 'r') as f:
        pars = json.load(f)

    incs = np.array([p['IncR'] for p in pars])
    rd.shuffle(incs)
    for inc, p in zip(incs, pars):
        p['IncR'] = inc

    for i, p in enumerate(pars):
        sim = CascadeRates.parse_json(p).simulate()
        sim['Location'] = loc
        sim['Key'] = i
        sim['Target'] = 'CNR+Prev'
        mss.append(sim)

mss = pd.DataFrame(mss)
mss.to_csv('../out/tx_time1.csv', index=False)
