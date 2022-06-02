import numpy as np
from scipy.optimize import minimize

__author__ = 'Chu-Chang Ku'
__all__ = ['Shifting', 'find_shifting']


class Shifting:
    def __init__(self, src, pr_pub0):
        self.Pr_Pub0 = pr_pub0
        n_dx0, n_entry = 0, 0
        n_dx1, n_tr = np.zeros(2), np.zeros(2)

        for _, tr in src.iterrows():
            fr, to = tr['From'], tr['To']

            if fr == 'UC':
                if to.startswith('pub'):
                    n_entry += tr['N']
                    if to.endswith('det'):
                        n_dx0 += tr['N']
            else:
                if fr.startswith('pub'):
                    n_tr[0] += tr['N']
                    if to.endswith('det'):
                        n_dx1[0] += tr['N']
                else:
                    n_tr[1] += tr['N']
                    if to.endswith('det'):
                        n_dx1[1] += tr['N']

        self.Dx0 = n_dx0 / n_entry
        self.Dx1 = np.repeat(n_dx1 / n_tr, [2, 1])

        self.P_Entry = np.zeros(3)
        self.P_Tr = np.zeros((3, 3))
        self.P_Dx1 = np.zeros((3, 3))
        self.P_Dx0 = np.zeros(3)

    def update(self, p):
        ppm, odds_pri, tr_pub_pub, tr_pri_pub = p['ppm'], p['odds_pri'], p['tr_pub_pub'], p['tr_pri_pub']
        pr_pub = self.Pr_Pub0

        self.P_Entry = np.array([pr_pub, (1 - pr_pub) * ppm, (1 - pr_pub) * (1 - ppm)])
        dx0, dx1 = self.Dx0, self.Dx1

        dx_pri = dx0 / (1 - dx0) * odds_pri
        dx_pri = dx_pri / (1 + dx_pri)
        self.P_Dx0 = np.array([dx0, dx0, dx_pri])

        dx_pri = dx1 / (1 - dx1) * odds_pri
        dx_pri = dx_pri / (1 + dx_pri)
        self.P_Dx1 = np.array([dx1, dx1, dx_pri]).T

        self.P_Tr = p_trm = np.zeros((3, 3))
        p_trm[0] = np.array([tr_pub_pub, (1 - tr_pub_pub) * ppm, (1 - tr_pub_pub) * (1 - ppm)])
        p_trm[1:] = np.array([tr_pri_pub, (1 - tr_pri_pub) * ppm, (1 - tr_pri_pub) * (1 - ppm)])

    def info(self, tol=1e-7):
        det = self.P_Entry * self.P_Dx0
        x = self.P_Entry * (1 - self.P_Dx0)

        i = 1
        n_vis = i * det.sum()

        while x.sum() > tol:
            x = x.reshape((-1, 1))
            i += 1
            tr0 = self.P_Tr * x
            tp, fp = tr0 * self.P_Dx1, tr0 * (1 - self.P_Dx1)
            det0 = tp.sum(0)

            det += det0
            n_vis += i * det0.sum()
            x = fp.sum(0)

        return {'det': det.reshape(-1) / det.sum(), 'vis': n_vis}


def find_shifting(src, pr_pub, ex_det, ex_vis=2.75):
    shf = Shifting(src, pr_pub)

    def fn(x, shf, ex_det, ex_vis, full=False):
        p = {'ppm': x[0], 'odds_pri': x[1], 'tr_pub_pub': x[2], 'tr_pri_pub': x[3]}

        shf.update(p)
        info = shf.info()
        det, vis = info['det'], info['vis']
        if full:
            return det, vis

        return ((det - ex_det) ** 2).sum() + (vis - ex_vis) ** 2

    bds = [(0, 1), (0.1, 10), (0, 1), (0, 1)]

    x0 = np.array([0.5, 0.7, 0.5, 0.5])

    opt = minimize(fn, x0=x0, args=(shf, ex_det, ex_vis,), method='SLSQP', bounds=bds)

    shf.update({
        'ppm': opt.x[0],
        'odds_pri': opt.x[1],
        'tr_pub_pub': opt.x[2],
        'tr_pri_pub': opt.x[3]
    })
    return shf
