import numpy as np
import numpy.random as rd
from scipy.integrate import solve_ivp
from sims_pars import bayes_net_from_script
from sims_pars.fit.targets import read_targets
from sims_pars.fit.base import DataModel, Particle
from ppa.inputs import *
import pandas as pd
import json

__author__ = 'Chu-Chang Ku'
__all__ = ['Obj', 'CascadeShort']


class I:
    Asym = 0
    Sym = 1
    ExCS = 2


class CascadeShort:
    def __init__(self, t0=2010, adr=0):
        self.T0Decline = t0
        self.ADR = adr

    def get_inc(self, t, pars):
        t = max(t, self.T0Decline)
        inc = pars['inc0']
        inc *= np.exp(- self.ADR * (t - self.T0Decline))
        return inc

    def __call__(self, t, y, pars):
        inc = self.get_inc(t, pars)

        onset = pars['r_onset'] * y[I.Asym]
        csi = pars['r_csi'] * y[I.Sym]
        det = pars['r_recsi'] * y[I.ExCS]

        dy = np.array([inc - onset, onset - csi, csi - det])

        drop = np.array([pars['r_die_asym'], pars['r_die_sym'], pars['r_die_sym']]) + pars['r_sc']
        dy -= y * drop
        return dy

    def measure(self, t, y, pars):
        inc = self.get_inc(t, pars)
        det = pars['r_recsi'] * y[I.ExCS]
        p_pub, ppm = pars['p_pub'], pars['ppm']

        sector = np.array([p_pub, (1 - p_pub) * ppm, (1 - p_pub) * (1 - ppm)])
        txi = np.array([pars['p_txi_pub'], pars['p_txi_eng'], pars['p_txi_eng']])
        pr = y / y.sum()
        mea = {'Year': t,
               'PrevUt': y.sum(), 'PrevAsym': y[I.Asym], 'PrevSym': y[I.Sym], 'PrevExCS': y[I.ExCS],
               'PrAsym': pr[0], 'PrSym': pr[1], 'PrExCS': pr[2],
               'IncR': inc,
               'TxiR_Pub': (det * sector * txi)[0], 'TxiR_Eng': (det * sector * txi)[1],
               'TxiR_Pri': (det * sector * txi)[2]}

        return mea


class Obj(DataModel):
    def __init__(self, data, model, bn, txo, exo):
        DataModel.__init__(self, data, bn, exo=exo)
        self.Model = model
        self.ParsTxOut = txo['pars']

    def simulate2measure(self, pars, t_eval=np.linspace(2010, 2025, 16)):
        p = rd.choice(self.ParsTxOut, 1)[0]
        p = dict(p)
        p.update(pars)

        model = self.Model

        y0 = np.zeros(3)
        ys = solve_ivp(model, [1500, max(t_eval)], y0=y0, args=(p,), dense_output=True)
        ms = pd.DataFrame([model.measure(t, ys.sol(t), p) for t in t_eval])
        ms = ms.set_index('Year')
        return ms

    def simulate(self, pars) -> Particle:
        ms = self.simulate2measure(pars, np.linspace(2015, 2020, 6))
        ms2019 = ms.loc[2019]
        ext = {
            'PrevUt_2019': ms2019.PrevUt,
            'PrA_2019': ms2019.PrAsym,
            'PrS_2019': ms2019.PrSym,
            'PrC_2019': ms2019.PrExCS,
            'txi_pub_2019': ms2019.TxiR_Pub,
            'txi_eng_2019': ms2019.TxiR_Eng,
            'txi_pri_2019': ms2019.TxiR_Pri,
        }
        return Particle(pars, ext)

    @staticmethod
    def load(root, loc='India'):
        with open(f'{root}/prior.txt', 'r') as f:
            bn = bayes_net_from_script(f.read())

        with open(f'{root}/{loc}/pars_txo.json', 'r') as f:
            txo = json.load(f)

        adr = get_adr(f'{root}/India/targets.csv')
        exo = get_exo(f'{root}/{loc}/targets.csv')

        dat = dict()

        d_ts = pd.read_csv(f'{root}/{loc}/targets.csv')
        d_ts = d_ts.rename(columns={'M': 'm', 'L': 'l', 'U': 'u'})

        dat['PrevUt_2019'] = dict(d_ts[d_ts.Index == 'PrevUt'].iloc[0, :])
        dat['PrA_2019'] = dict(d_ts[d_ts.Index == 'PrAsym'].iloc[0, :])
        dat['PrS_2019'] = dict(d_ts[d_ts.Index == 'PrSym'].iloc[0, :])
        dat['PrC_2019'] = dict(d_ts[d_ts.Index == 'PrExCS'].iloc[0, :])

        for k in ['txi_pub', 'txi_eng', 'txi_pri']:
            vs = np.array([p[k] for p in txo['pars']])
            qs = np.quantile(vs, [0.025, 0.5, 0.975])
            dat[f'{k}_2019'] = {'l': qs[0], 'm': qs[1], 'u': qs[2]}

        dat = read_targets(dat)
        for d in dat.values():
            try:
                d.Range = np.abs(d.Range)
            except AttributeError:
                pass

        model = CascadeShort(adr=adr)

        return Obj(dat, model, bn, txo, exo)


if __name__ == '__main__':
    obj0 = Obj.load('../../data', 'Delhi')

    p0 = obj0.sample_prior()
    sim0 = obj0.simulate(p0)
    for k, v in sim0.Sims.items():
        print(k, '\t', v)

    print(obj0.calc_distance(sim0))
