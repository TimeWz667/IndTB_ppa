import numpy as np
from abc import ABCMeta, abstractmethod

__author__ = 'Chu-Chang Ku'
__all__ = ['AbsPathway', 'PathwayPlain']


class AbsPathway(metaclass=ABCMeta):
    @abstractmethod
    def __call__(self, pars):
        pass


class PathwayPlain(AbsPathway):
    def __call__(self, pars):
        pp = dict()

        pp['p_dx0'] = np.zeros(3)

        pp['p_dx1'] = np.ones((3, 3))
        pp['p_tr'] = np.diag(np.ones(3))
        p_pub, p_ppm = pars['p_pub'], pars['p_ppm']

        pp['prop_cs'] = pp['p_entry'] = np.array([p_pub, (1 - p_pub) * p_ppm, (1 - p_pub) * (1 - p_ppm)])

        p_pri_on_pub = pars['p_pri_on_pub']

        pp['p_trtx'] = np.array([
            [1, 0, 0],
            [0, p_pri_on_pub, 1 - p_pri_on_pub],
            [0, 0, 1]
        ])

        pp['r_csi'] = pars['r_aware']
        pp['r_recsi'] = pars['r_det']

        pp['p_txi'] = np.array([pars['p_txi_pub'], pars['p_txi_eng'], pars['txi_pri']])
        pp['ppv'] = np.array([pars['ppv_pub'], pars['ppv_eng'], pars['ppv_pri']])
        pp['dur'] = np.array([pars['dur_pub'], pars['dur_pub'], pars['dur_pri']])
        pp['p_ts'] = np.array([pars['p_txs_pub'], pars['p_txs_eng'], pars['p_txs_eng']])
        pp['p_tl'] = np.array([pars['p_txl_pub'], pars['p_txl_eng'], pars['p_txl_eng']])
        pp['p_td'] = np.array([pars['p_txd_pub'], pars['p_txd_eng'], pars['p_txd_eng']])
        return pp
