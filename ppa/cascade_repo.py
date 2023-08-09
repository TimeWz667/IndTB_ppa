import numpy as np
import numpy.random as rd
import json

__author__ = 'Chu-Chang Ku'
__all__ = ['CasRepo']


class CasRepo:
    def __init__(self, pars, prev, txo):
        self.Pars = pars
        self.Prev = prev
        self.TxO = {
            'Succ': np.zeros(2),
            'Die': np.zeros(2),
        }

        for d in txo:
            if d['Index'] == 'TxSucc':
                if d['Tag'] == 'Pub':
                    self.TxO['Succ'][0] = d['M']
                else:
                    self.TxO['Succ'][1] = d['M']
            else:
                if d['Tag'] == 'Pub':
                    self.TxO['Die'][0] = d['M']
                else:
                    self.TxO['Die'][1] = d['M']

        self.TxO['LTFU'] = 1 - self.TxO['Succ'] - self.TxO['Die']

    def prepare_pars(self, exo, pp=None):
        if pp is None:
            pp = rd.choice(self.Pars, 1)[0]

        prev = self.Prev
        p_a, p_s, p_c = prev['PrevAsym'], prev['PrevSym'], prev['PrevExCS']

        mu_a = exo['r_die_asym'] + exo['r_sc']
        mu_s = mu_c = exo['r_die_sym'] + exo['r_sc']

        txis = np.array([pp[f'txi_{sec}'] for sec in ['pub', 'eng', 'pri']])
        txi = txis.sum()

        p_txi = np.array([exo[f'p_txi_{sec}'] for sec in ['pub', 'eng', 'pri']])

        det = txis / p_txi
        ppm = det[1] / det[1:].sum()
        p_ent_pub = exo['p_ent_pub']
        p_ent = np.array([p_ent_pub, (1 - p_ent_pub) * ppm, (1 - p_ent_pub) * (1 - ppm)])

        p_dx = det / p_ent
        p_dx *= exo['p_dx_pub'] / p_dx[0]

        p_dx_all = (p_txi * p_dx * p_ent).sum()

        r_det = txi / p_c
        r_csi = p_c * (r_det + mu_c) / p_s
        r_onset = p_s * (r_csi + mu_s) / p_a

        txi0 = p_s * r_csi * p_dx_all

        # assert txi > txi0
        txi1 = max(txi - txi0, 0)
        r_recsi = txi1 / p_c / p_dx_all

        p_dx0 = p_dx1 = p_dx

        alo = np.array([[1, 0], [pp['p_pri_on_pub'], 1 - pp['p_pri_on_pub']], [0, 1]])

        inc = (mu_a + r_onset) * p_a
        dur = np.array([0.5, pp['dur_pri']])

        ps = {
            'p_txi': p_txi,
            'p_dx': p_dx,
            'p_dx0': p_dx0,
            'p_dx1': p_dx1,
            'p_ent': p_ent,
            'r_onset': r_onset,
            'r_csi': r_csi,
            'r_recsi': r_recsi,
            'mu_a': mu_a,
            'mu_s': mu_s,
            'mu_c': mu_c,
            'tx_alo': alo,
            'tx_dur': dur,
            'r_txs': self.TxO['Succ'] / dur,
            'r_txd': self.TxO['Die'] / dur,
            'r_txl': self.TxO['LTFU'] / dur,
            'ppv': np.array([pp[f'ppv_{sec}'] for sec in ['pub', 'eng', 'pri']]),
            'inc': inc,
            'prev_ut': prev['PrevUt'],
            'prev_asc': (p_a, p_s, p_c)

        }
        ps.update(exo)

        return ps

    def reform(self, ps, rat_01=1):
        prev = self.Prev
        p_s, p_c = prev['PrevSym'], prev['PrevExCS']

        ps1 = dict(ps)

        r_csi, r_recsi, p_ent, p_txi = ps['r_csi'], ps['r_recsi'], ps['p_ent'], ps['p_txi']

        pdx = (p_ent * p_txi * ps['p_dx']).sum()

        txi = p_s * r_csi * pdx + p_c * r_recsi * pdx

        k = (p_s * r_csi * rat_01 + p_c * r_recsi * 1) / txi

        ps1['p_dx0'] = ps['p_dx'] / pdx * rat_01 / k
        ps1['p_dx1'] = ps['p_dx'] / pdx / k
        return ps1

    @staticmethod
    def load(file):
        with open(file, 'r') as f:
            js = json.load(f)

        return CasRepo(js['pars'], js['prev'], js['txo'])


if __name__ == '__main__':
    from sims_pars import bayes_net_from_script, sample

    with open('../../data/prior.txt', 'r') as f:
        prior = bayes_net_from_script(f.read())

    cr = CasRepo.load('../../data/Bihar/pars_txo.json')

    exo = sample(prior)
    ps = cr.prepare_pars(exo=exo)

    for k, v in ps.items():
        print(k, ': ', v)
