from ppa.dy.pathway import AbsPathway
from ppa.dy.base import AbsModel
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['I', 'CohortModel']


class I:
    Asym = 0
    Sym = 1
    ExCsPub = 2
    ExCsEng = 3
    ExCsPri = 4
    TxPub = 5
    TxEng = 6
    TxPri = 7
    Cured = 8
    SelfCured = 9
    LTFU = 10
    DeadU = 11
    DeadT = 12

    ExCs = [ExCsPub, ExCsEng, ExCsPri]
    Tx = [TxPub, TxEng, TxPri]

    N_States = 13


class CohortModel(AbsModel):
    def __init__(self, pp: AbsPathway, t0=2020):
        AbsModel.__init__(self, t0=t0, t1=t0+20, dt=0.1, pp=pp)

    def get_y0(self, t, pars) -> np.ndarray:
        y0 = np.zeros(I.N_States)
        y0[I.Asym] = 1
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

        calc['txi0'] = r_csi * p_entry * p_dx0 * p_txi * y[I.Sym]
        calc['loss0'] = r_csi * p_entry * (1 - p_dx0 * p_txi) * y[I.Sym]

        calc['txi1'] = r_recsi * (p_tr * p_dx1 * p_txi) * y[I.ExCs].reshape((-1, 1))
        calc['loss1'] = r_recsi * (p_tr * (1 - p_dx1 * p_txi)) * y[I.ExCs].reshape((-1, 1))

        calc['acf_a'] = acf_a = pars['r_acf'] * pars['sens_acf'] * y[I.Asym]
        calc['acf_s'] = pars['r_acf'] * pars['sens_acf'] * y[I.Sym]
        calc['acf_c'] = pars['r_acf'] * pars['sens_acf'] * y[I.ExCs]

        dur = pp['dur']
        p_ts, p_tl, p_td = pp['p_ts'], pp['p_tl'], pp['p_td']

        calc['to_tp'] = to_tp = 1 / dur * y[I.Tx]
        calc['ts_tp'] = p_ts * to_tp
        calc['tl_tp'] = p_tl * to_tp
        calc['td_tp'] = p_td * to_tp
        return calc

    def __call__(self, t, y, pars):
        calc = self.calc(t, y, pars)
        dy = np.zeros_like(y)

        onset = calc['onset']
        sc_a, sc_s, sc_c = calc['sc_a'], calc['sc_s'], calc['sc_c']
        die_a, die_s, die_c = calc['die_a'], calc['die_s'], calc['die_c']
        acf_a, acf_s, acf_c = calc['acf_a'], calc['acf_s'], calc['acf_c']

        txi0, loss0, txi1, loss1 = calc['txi0'], calc['loss0'], calc['txi1'], calc['loss1']

        dy[I.Asym] += - onset - acf_a - sc_a - die_a
        dy[I.Sym] += onset - (txi0 + loss0).sum() - acf_s - sc_s - die_s
        dy[I.ExCs] += txi0 + loss0 - txi1.sum(1) - loss1.sum(1) + loss1.sum(0) - acf_c - sc_c - die_c

        txi = txi0 + txi1.sum(0)
        p_trtx = pars['pp']['p_trtx']

        dy[I.Tx] += (p_trtx * txi.reshape((-1, 1))).sum(0) - calc['to_tp']
        dy[I.TxPub] += acf_a + acf_s + acf_c.sum()

        dy[I.Cured] += calc['ts_tp'].sum()
        dy[I.SelfCured] += sc_a + sc_s + sc_c.sum()
        dy[I.LTFU] += calc['tl_tp'].sum()
        dy[I.DeadU] += die_a + die_s + die_c.sum()
        dy[I.DeadT] += calc['td_tp'].sum()


        return dy

    def measure(self, t, y, pars):
        calc = self.calc(t, y, pars)
        mea = {'Time': t}
        mea['N'] = y.sum()
        mea['S_Asym'] = y[I.Asym]
        mea['S_Sym'] = y[I.Sym]
        mea['S_ExCs'] = y[I.ExCs].sum()
        mea['S_TxPub'] = y[I.TxPub]
        mea['S_TxPri'] = y[I.TxEng] + y[I.TxPri]
        mea['S_Cured'] = y[I.Cured]
        mea['S_SelfCured'] = y[I.SelfCured]
        mea['S_LTFU'] = y[I.LTFU]
        mea['S_DeadU'] = y[I.DeadU]
        mea['S_DeadT'] = y[I.DeadT]

        mea['F_Onset'] = calc['onset']
        mea['F_CSI'] = (calc['txi0'] + calc['loss0']).sum()
        mea['F_TxI'] = calc['txi0'].sum() + calc['txi1'].sum()
        mea['F_Acf'] = calc['acf_a'] + calc['acf_s'] + calc['acf_c'].sum()
        return mea


if __name__ == '__main__':
    import json
    from ppa.dy.pathway import PathwayPlain
    import matplotlib.pyplot as plt

    loc = 'India'

    with open(f'../../docs/pars/pars_nods_{loc}.json', 'r') as f:
        post = json.load(f)

    m = CohortModel(PathwayPlain())

    p0 = post[0]
    print(p0.keys())
    ys, ms = m.simulate(p0)

    ms[['S_DeadU', 'S_DeadT', 'S_Cured', 'S_SelfCured']].plot()

    # ms[['CNR_Acf', 'CNR_Pub', 'CNR_Eng']].plot()
    plt.show()


