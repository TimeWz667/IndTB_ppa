from ppa.dy.pathway import AbsPathway
from ppa.dy.base import AbsModel
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['I', 'SteadyStateModel']


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


class SteadyStateModel(AbsModel):
    def __init__(self, pp: AbsPathway, t1=2025):
        AbsModel.__init__(self, 2015, t1, 0.5, pp)

    def get_y0(self, t, pars) -> np.ndarray:
        k = np.exp(- pars['adr'] * (t - 2020))
        prev = k * pars['prv0']

        p_det = np.array(
            [pars['p_pub'], (1 - pars['p_pub']) * pars['p_ppm'], (1 - pars['p_pub']) * (1 - pars['p_ppm'])])

        p_txi = np.array([pars['p_txi_pub'], pars['p_txi_eng'], pars['txi_pri']])

        r_det_all = pars['r_det'] * (p_det * p_txi).sum()
        r_acf = pars['r_acf'] * pars['sens_acf']
        rs = rc = pars['r_sc'] + pars['r_death_s'] + r_acf

        asc = np.array([
            (rs + pars['r_aware'] - pars['adr']) / pars['r_sym'],
            1,
            pars['r_aware'] / (rc + r_det_all - pars['adr'])
        ])
        asc /= asc.sum()

        y0 = np.zeros(I.N_States)

        y0[I.Asym] = asc[0] * prev
        y0[I.Sym] = asc[1] * prev
        y0[I.ExCs] = p_det * asc[2] * prev

        ppv = np.array([pars['ppv_pub'], pars['ppv_eng'], pars['ppv_pri']])
        dur = np.array([pars['dur_pub'], pars['dur_pub'], pars['dur_pri']])

        det = pars['r_det'] * p_det * p_txi * y0[I.ExCs].sum()
        det_fp = det * (1 - ppv) / ppv
        p_pri_on_pub = pars['p_pri_on_pub']
        det = np.array([det[0], det[1] * p_pri_on_pub, det[1] * (1 - p_pri_on_pub) + det[2]])
        det_fp = np.array([det_fp[0], det_fp[1] * p_pri_on_pub, det_fp[1] * (1 - p_pri_on_pub) + det_fp[2]])

        r = (1 / dur - pars['adr'])

        acf = r_acf * y0[[I.Asym, I.Sym] + I.ExCs].sum()
        y0[I.Tp] = det / r
        y0[I.Fp] = det_fp / r
        y0[I.TpPub] += acf / r[0]
        y0[I.FpPub] += pars['r_acf'] * (1 - pars['spec_acf']) * (1 - y0.sum()) / r[0]

        y0[I.NonTB] = 1 - y0.sum()

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

        calc['to_tp'] = to_tp = 1 / dur * y[I.Tp]
        calc['ts_tp'] = p_ts * to_tp
        calc['tl_tp'] = p_tl * to_tp
        calc['td_tp'] = p_td * to_tp

        calc['to_fp'] = to_fp = 1 / dur * y[I.Fp]
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

        txi0, loss0, txi1, loss1 = calc['txi0'], calc['loss0'], calc['txi1'], calc['loss1']

        dy[I.NonTB] += sc_a + sc_s + sc_c.sum() + die_a + die_s + die_c.sum() - inc
        dy[I.Asym] += inc - onset - acf_a - sc_a - die_a
        dy[I.Sym] += onset - (txi0 + loss0).sum() - acf_s - sc_s - die_s
        dy[I.ExCs] += txi0 + loss0 - txi1.sum(1) - loss1.sum(1) + loss1.sum(0) - acf_c - sc_c - die_c

        txi = txi0 + txi1.sum(0)
        p_trtx = pars['pp']['p_trtx']

        dy[I.Tp] += (p_trtx * txi.reshape((-1, 1))).sum(0) - calc['to_tp']
        dy[I.TpPub] += acf_a + acf_s + acf_c.sum()

        acf_fp = calc['acf_fp']
        txi_fp = calc['txi_fp']

        dy[I.Fp] = (p_trtx * txi_fp.reshape((-1, 1))).sum(0) - calc['to_fp']
        dy[I.FpPub] += acf_fp

        dy[I.NonTB] += calc['to_tp'].sum() + calc['to_fp'].sum() - txi_fp.sum() - acf_fp

        return dy

    def measure(self, t, y, pars):
        calc = self.calc(t, y, pars)
        mea = {'Time': t}
        mea['N'] = y.sum()
        mea['PrevA'] = y[I.Asym]
        mea['PrevS'] = y[I.Sym]
        mea['PrevC'] = y[I.ExCs].sum()
        mea['TpPub'] = y[I.TpPub]
        mea['TpEng'] = y[I.TpEng]
        mea['TpPri'] = y[I.TpPri]
        mea['FpPub'] = y[I.FpPub]
        mea['FpEng'] = y[I.FpEng]
        mea['FpPri'] = y[I.FpPri]
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

    m = SteadyStateModel(PathwayPlain())

    p0 = post[0]
    print(p0.keys())
    ys, ms = m.simulate(p0)

    print('ADR pars: ', p0['adr'])
    print('ADR sims: ', - np.diff(np.log(ms.TpPub)) * 2)
    print('ADR sims: ', - np.diff(np.log(ms.FpPub)) * 2)
