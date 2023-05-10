import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class I:
    NonTB = 0
    Asym = 1
    Sym = 2
    ExCS = 3
    TxPub = 4
    TxPri = 5
    FpPub = 6
    FpPri = 7

    N_States = 8


class Model:
    def __init__(self, inc_timing='A'):
        self.IncTiming = inc_timing

    def get_y0(self, pars):
        y0 = np.zeros(I.N_States)
        y0[I.Asym], y0[I.Sym], y0[I.ExCS] = pars['prev_asc']
        y0[I.NonTB] = 1 - y0.sum()
        return y0

    def measure(self, t, y, pars):
        calc = self.calculate(t, y, pars)

        return {
            'Time': t,
            'N': y.sum(),
            'PrevUt': y[I.Asym] + y[I.Sym] + y[I.ExCS],
            'PrevA': y[I.Asym],
            'PrevS': y[I.Sym],
            'PrevC': y[I.ExCS],
            'PrevTxPub': y[I.TxPub],
            'PrevTxPri': y[I.TxPri],
            'Inc': calc['inc'],
            'Tp_pub': calc['tp0'][0] + calc['tp1'][0],
            'Tp_eng': calc['tp0'][1] + calc['tp1'][1],
            'Tp_pri': calc['tp0'][2] + calc['tp1'][2],
            'Fp_pub': calc['fp'][0],
            'Fp_eng': calc['fp'][1],
            'Fp_pri': calc['fp'][2]
        }

    def calculate(self, t, y, pars):
        calc = dict()

        calc['sc_a'] = sc_a = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = sc_s = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = sc_c = pars['r_sc'] * y[I.ExCS]
        calc['die_a'] = die_a = pars['r_die_asym'] * y[I.Asym]
        calc['die_s'] = die_s = pars['r_die_sym'] * y[I.Sym]
        calc['die_c'] = die_c = pars['r_die_sym'] * y[I.ExCS]

        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']

        p_ent, p_txi = pars['p_ent'], pars['p_txi']
        pdx0 = p_ent * pars['p_dx0'] * p_txi
        pdx1 = p_ent * pars['p_dx1'] * p_txi

        calc['onset'] = onset = pars['r_onset'] * y[I.Asym]
        calc['tp0'] = tp0 = r_csi * pdx0 * y[I.Sym]
        calc['fn0'] = r_csi * (p_ent - pdx0) * y[I.Sym]
        calc['tp1'] = tp1 = r_recsi * pdx1 * y[I.ExCS]

        calc['txi'] = txi = np.dot(tp0 + tp1, pars['tx_alo'])
        calc['txe'] = y[[I.TxPub, I.TxPri]] / pars['tx_dur']

        # calc['inc'] = pars['inc']
        if self.IncTiming == 'A':
            calc['inc'] = onset + sc_a + die_a
        elif self.IncTiming == 'S':
            calc['inc'] = calc['fn0'].sum() + tp0.sum() + sc_a + die_a + sc_s + die_s
        else:
            calc['inc'] = tp0.sum() + tp1.sum() + sc_a + die_a + sc_s + die_s + sc_c + die_c

        ppv = pars['ppv']
        calc['fp'] = fp = (tp0 + tp1) * (1 - ppv) / ppv
        calc['ftxi'] = np.dot(fp, pars['tx_alo'])
        calc['ftxe'] = y[[I.FpPub, I.FpPri]] / pars['tx_dur']

        return calc

    def __call__(self, t, y, pars):
        calc = self.calculate(t, y, pars)
        dy = np.zeros_like(y)

        inc, onset = calc['inc'], calc['onset']
        sc_a, sc_s, sc_c = calc['sc_a'], calc['sc_s'], calc['sc_c']
        die_a, die_s, die_c = calc['die_a'], calc['die_s'], calc['die_c']

        tp0, fn0, tp1 = calc['tp0'], calc['fn0'], calc['tp1']

        txi, txe = calc['txi'], calc['txe']

        dy[I.Asym] += inc - onset - sc_a - die_a
        dy[I.Sym] += onset - tp0.sum() - fn0.sum() - sc_s - die_s
        dy[I.ExCS] += fn0.sum() - tp1.sum() - sc_c - die_c
        dy[I.TxPub] += txi[0] - txe[0]
        dy[I.TxPri] += txi[1] - txe[1]

        dy[I.NonTB] += txe.sum() + sc_a + die_a + sc_s + die_s + sc_c + die_c - inc

        ftxi, ftxe = calc['ftxi'], calc['ftxe']

        dy[I.FpPub] += ftxi[0] - ftxe[0]
        dy[I.FpPri] += ftxi[1] - ftxe[1]
        dy[I.NonTB] += ftxe.sum() - ftxi.sum()

        return dy

    def simulate(self, pars):
        y0 = self.get_y0(pars)
        ys = solve_ivp(self, [2000, 2025], y0=y0, args=(pars,), dense_output=True)

        ms = [self.measure(t, ys.sol(t), pars) for t in np.linspace(2010, 2025, 16)]
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

    m = Model('A')
    ms = m.simulate(ps).tail(10)

    print(ms)
