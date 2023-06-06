import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class I:
    Asym = 0
    Sym = 1
    ExCS = 2
    TxPub = 3
    TxPri = 4
    SelfCured = 5
    UtDeath = 6
    TxSucc = 7
    TxLTFU = 8
    TxDeath = 9

    N_States = 10


class Model:
    def get_y0(self, pars):
        y0 = np.zeros(I.N_States)
        y0[I.Asym] = 1
        return y0

    def measure(self, t, y, pars):
        return {
            'Time': t,
            'N': y.sum(),
            'S_Asym': y[I.Asym],
            'S_Sym': y[I.Sym],
            'S_ExCS': y[I.ExCS],
            'S_TxPub': y[I.TxPub],
            'S_TxPri': y[I.TxPri],
            'E_SelfCured': y[I.SelfCured],
            'E_UtDeath': y[I.UtDeath],
            'E_TxSucc': y[I.TxSucc],
            'E_TxLTFU': y[I.TxLTFU],
            'E_TxDeath': y[I.TxDeath]
        }

    def calculate(self, t, y, pars):
        calc = dict()

        calc['sc_a'] = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = pars['r_sc'] * y[I.ExCS]
        calc['die_a'] = pars['r_die_asym'] * y[I.Asym]
        calc['die_s'] = pars['r_die_sym'] * y[I.Sym]
        calc['die_c'] = pars['r_die_sym'] * y[I.ExCS]

        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']

        p_ent, p_txi = pars['p_ent'], pars['p_txi']
        pdx0 = p_ent * pars['p_dx0'] * p_txi
        pdx1 = p_ent * pars['p_dx1'] * p_txi

        calc['onset'] = pars['r_onset'] * y[I.Asym]
        calc['tp0'] = tp0 = r_csi * pdx0 * y[I.Sym]
        calc['fn0'] = r_csi * (p_ent - pdx0) * y[I.Sym]
        calc['tp1'] = tp1 = r_recsi * pdx1 * y[I.ExCS]

        calc['txi'] = np.dot(tp0 + tp1, pars['tx_alo'])

        calc['txs'] = y[[I.TxPub, I.TxPri]] * pars['r_txs']
        calc['txl'] = y[[I.TxPub, I.TxPri]] * pars['r_txl']
        calc['txd'] = y[[I.TxPub, I.TxPri]] * pars['r_txd']

        return calc

    def __call__(self, t, y, pars):
        calc = self.calculate(t, y, pars)
        dy = np.zeros_like(y)

        onset = calc['onset']
        sc_a, sc_s, sc_c = calc['sc_a'], calc['sc_s'], calc['sc_c']
        die_a, die_s, die_c = calc['die_a'], calc['die_s'], calc['die_c']

        tp0, fn0, tp1 = calc['tp0'], calc['fn0'], calc['tp1']

        txi = calc['txi']
        txs, txd, txl = calc['txs'], calc['txd'], calc['txl']

        dy[I.Asym] += - onset - sc_a - die_a
        dy[I.Sym] += onset - tp0.sum() - fn0.sum() - sc_s - die_s
        dy[I.ExCS] += fn0.sum() - tp1.sum() - sc_c - die_c
        dy[I.TxPub] += txi[0] - txs[0] - txl[0] - txd[0]
        dy[I.TxPri] += txi[1] - txs[1] - txl[1] - txd[1]

        dy[I.UtDeath] += die_a + die_s + die_c
        dy[I.SelfCured] += sc_a + sc_s + sc_c
        dy[I.TxDeath] += txd.sum()
        dy[I.TxSucc] += txs.sum()
        dy[I.TxLTFU] += txl.sum()

        return dy

    def simulate(self, pars, n_year=3):
        y0 = self.get_y0(pars)
        y_end = 20
        ys = solve_ivp(self, [0, y_end], y0=y0, args=(pars,), dense_output=True)

        ms = [self.measure(t, ys.sol(t), pars) for t in np.linspace(0, n_year, n_year * 10 + 1)]
        ms.append(self.measure(y_end, ys.sol(y_end), pars))
        ms = pd.DataFrame(ms)
        return ms


if __name__ == '__main__':
    from sims_pars import bayes_net_from_script, sample
    from ppa.pars import *

    with open('../../data/prior.txt', 'r') as f:
        prior = bayes_net_from_script(f.read())

    cr = CasRepo.load('../../data/pars_india.json')

    exo = sample(prior)
    ps = cr.prepare_pars(exo=exo)

    m = Model()
    ms = m.simulate(ps)

    print(ms)
