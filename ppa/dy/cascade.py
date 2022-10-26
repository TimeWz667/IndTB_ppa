from ppa.dy.pathway import AbsPathway
from ppa.dy.base import AbsModel
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['I']


class I:
    NonTB = 0
    Asym = 1
    Sym = 2

    ExCsPub = 3
    ExCsEng = 4
    ExCsPri = 5

    TpPub = 6
    TpEng = 7
    TpPri = 8

    FpPub = 9
    FpEng = 10
    FpPri = 11

    ExCs = [ExCsPub, ExCsEng, ExCsPri]
    Tp = [TpPub, TpEng, TpPri]
    Fp = [FpPub, FpEng, FpPri]

    N_States = 12


class CascadeModel(AbsModel):
    def __init__(self, pp: AbsPathway, t1=2025):
        AbsModel.__init__(self, 2015, t1, 0.5, pp)

    def get_y0(self, t, pars) -> np.ndarray:
        k = np.exp(- pars['adr'] * (t - 2020))
        prev = k * pars['prv0']

        pp = pars['pp']
        r_det_all = (pars['r_det'] * pp['p_entry'] * pp['p_txi']).sum()

        asc = np.array([
            (pars['rs'] + pars['r_aware'] - pars['adr']) / pars['r_sym'],
            1,
            pars['r_aware'] / (pars['rc'] + r_det_all - pars['adr'])
        ])
        asc /= asc.sum()

        y0 = np.zeros(I.N_States)
        y0[I.NonTB] = 1 - prev
        y0[I.Asym] = asc[0] * prev * k
        y0[I.Sym] = asc[1] * prev * k
        y0[I.ExCs] = pp['prop_cs'] * asc[2] * prev * k

        return y0

    def calc(self, t, y, pars):
        pp = pars['pp']

        calc = dict()

        calc['onset'] = onset = pars['r_sym'] * y[I.Asym]

        calc['sc_a'] = sc_a = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = pars['r_sc'] * y[I.ExCs]
        calc['die_a'] = die_a = pars['r_death_a'] * y[I.Asym]
        calc['die_s'] = pars['r_death_s'] * y[I.Sym]
        calc['die_c'] = pars['r_death_s'] * y[I.ExCs]

        p_entry, p_dx0 = pp['p_entry'], pp['p_dx0']
        p_dx1, p_tr, p_trtx = pp['p_dx1'], pp['p_tr'], pp['p_trtx']

        r_csi, r_recsi = pp['r_csi'], pp['r_recsi']

        ppv, p_txi = pp['ppv'], pp['p_txi']

        calc['tp0'] = r_csi * p_entry * p_dx0 * y[I.Sym]
        calc['txi0'] = r_csi * p_entry * p_dx0 * p_txi * y[I.Sym]
        calc['loss0'] = r_csi * p_entry * (1 - p_dx0 * p_txi) * y[I.Sym]

        calc['tp1'] = r_recsi * (p_tr * p_dx1) * y[I.ExCs].reshape((-1, 1))
        calc['txi1'] = r_recsi * (p_tr * p_dx1 * p_txi) * y[I.ExCs].reshape((-1, 1))
        calc['loss1'] = r_recsi * (p_tr * (1 - p_dx1 * p_txi)) * y[I.ExCs].reshape((-1, 1))

        det = calc['tp0'] + calc['tp1'].sum(0)

        calc['det_fp'] = det * (1 - ppv) / ppv
        calc['txi_fp'] = calc['det_fp'] * p_txi

        calc['acf_a'] = acf_a = pars['r_acf'] * pars['sens_acf'] * y[I.Asym]
        calc['acf_s'] = pars['r_acf'] * pars['sens_acf'] * y[I.Sym]
        calc['acf_c'] = pars['r_acf'] * pars['sens_acf'] * y[I.ExCs]

        calc['acf_fp'] = pars['r_acf'] * pars['spec_acf'] * y[I.NonTB]

        calc['inc'] = onset + acf_a + sc_a + die_a - pars['adr'] * y[I.Asym]

        dur = pp['dur']
        p_ts, p_tl, p_td = pp['p_ts'], pp['p_tl'], pp['p_td']

        calc['to_tp'] = to_tp = dur * y[I.Tp]
        calc['ts_tp'] = p_ts * to_tp
        calc['tl_tp'] = p_tl * to_tp
        calc['td_tp'] = p_td * to_tp

        calc['to_fp'] = to_fp = dur * y[I.Fp]
        calc['ts_fp'] = p_ts * to_fp
        calc['tl_fp'] = p_tl * to_fp
        calc['td_fp'] = p_td * to_fp
        return calc

    def __call__(self, t, y, pars):
        calc = self.calc(t, y, pars)
        dy = np.zeros_like(y)

        inc, onset = calc['inc'], calc['onset']
        sc_a, sc_s, sc_c = calc['sc_a'], calc['sc_s'], calc['sc_c']
        die_a, die_s, die_c = calc['die_a'], calc['die_s'], calc['die_c']
        acf_a, acf_s, acf_c = calc['acf_a'], calc['acf_s'], calc['acf_c']

        txi0, loss0, txi1 = calc['txi0'], calc['loss0'], calc['txi1']

        dy[I.NonTB] += sc_a + sc_s + sc_c.sum() + die_a + die_s + die_c.sum() - inc
        dy[I.Asym] += inc - onset - acf_a - sc_a - die_a
        dy[I.Sym] += onset - (txi0 + loss0).sum() - acf_s - sc_s - die_s
        dy[I.ExCs] += txi0 + loss0 - txi1.sum(1) - acf_c - sc_c - die_c

        txi = txi0 + txi1.sum(0)
        p_trtx = pars['pp']['p_trtx']

        dy[I.Tp] += (p_trtx * txi.reshape((-1, 1))).sum(0)
        dy[I.TpPub] += acf_a + acf_s + acf_c.sum()

        acf_fp = calc['acf_fp']
        txi_fp = calc['txi_fp']

        dy[I.NonTB] -= txi_fp.sum() + acf_fp

        dy[I.Fp] = (p_trtx * txi_fp.reshape((-1, 1))).sum(0)
        dy[I.FpPub] += acf_fp

        return dy

    def measure(self, t, y, pars):
        calc = self.calc(t, y, pars)
        mea = {'Time': t}
        mea['N'] = n = y.sum()
        mea['PrevA'] = y[I.Asym]
        mea['PrevS'] = y[I.Sym]
        mea['PrevC'] = y[I.ExCs].sum()
        mea['TxPub'] = y[I.TpPub] + y[I.FpPub]
        mea['TxEng'] = y[I.TpEng] + y[I.FpEng]
        mea['TxPri'] = y[I.TpPri] + y[I.FpPri]
        mea['IncR'] = calc['inc']

        return mea


if __name__ == '__main__':
    import json
    from ppa.dy.pathway import PathwayPlain
    import matplotlib.pylab as plt

    loc = 'India'

    with open(f'../../docs/pars/pars_nods_{loc}.json', 'r') as f:
        post = json.load(f)

    m = CascadeModel(PathwayPlain())

    p0 = post[0]
    print(p0.keys())
    ys, ms = m.simulate(p0)

    print('ADR pars: ', p0['adr'])
    print('ADR sims: ', - np.diff(np.log(ms.IncR))[:10] * 2)
