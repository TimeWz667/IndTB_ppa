from ppa.healthcare.diagnosis import *
import numpy as np
import numpy.random as rd

__author__ = 'Chu-Chang Ku'
__all__ = ['get_system']


class Sector:
    def __init__(self, p_ent, alg):
        self.Entry = p_ent
        self.Algorithms = alg

    def seek_care(self, n_tb, n_nontb):
        res = Results()
        for i, p in enumerate(self.Entry):
            res += self.Algorithms[i].dx(n_tb * p, n_nontb * p)
        return res

    def seek_care_sto(self, n_tb, n_nontb):
        ns_tb = rd.multinomial(n_tb, self.Entry)
        ns_nontb = rd.multinomial(n_nontb, self.Entry)
        res = Results()
        for i in range(len(self.Entry)):
            res += self.Algorithms[i].dx_sto(ns_tb[i], ns_nontb[i])
        return res


class System:
    def __init__(self, p_ent, pub, eng, pri):
        self.Entry = p_ent
        self.Public = pub
        self.Engaged = eng
        self.Private = pri

    def seek_care(self, n_tb, n_nontb):
        n_tb_pub, n_tb_eng, n_tb_pri = self.Entry * n_tb
        n_nontb_pub, n_nontb_eng, n_nontb_pri = self.Entry * n_nontb
        return {
            'Public': self.Public.seek_care(n_tb_pub, n_nontb_pub),
            'Engaged': self.Engaged.seek_care(n_tb_eng, n_nontb_eng),
            'Private': self.Private.seek_care(n_tb_pri, n_nontb_pri)
        }

    def seek_care_sto(self, n_tb, n_nontb):
        n_tb_pub, n_tb_eng, n_tb_pri = rd.multinomial(n_tb, self.Entry)
        n_nontb_pub, n_nontb_eng, n_nontb_pri = rd.multinomial(n_nontb, self.Entry)
        return {
            'Public': self.Public.seek_care_sto(n_tb_pub, n_nontb_pub),
            'Engaged': self.Engaged.seek_care_sto(n_tb_eng, n_nontb_eng),
            'Private': self.Private.seek_care_sto(n_tb_pri, n_nontb_pri)
        }


def get_system(p, has_cdx=True):
    sputum = Specimen('Sputum', p['p_loss_sputum'])
    ssm = Test('SSM', p['sens_ssm'], p['spec_ssm'], sputum)
    xpert = Test('Xpert', p['sens_xpert'], p['spec_xpert'], sputum)
    xpert_ssm = Test('Xpert_ss-', p['sens_xpert_ss-'], p['spec_xpert'], sputum)
    if has_cdx:
        cdx = Test('CDx', p['sens_cdx'], p['spec_cdx'])
    else:
        cdx = None

    alg1 = Algorithm('SSM > Xpert > CDx', ssm=ssm, xpert=xpert_ssm, cdx=cdx)
    alg2 = Algorithm('SSM > CDx', ssm=ssm, cdx=cdx)
    alg3 = Algorithm('Xpert > CDx', xpert=xpert, cdx=cdx)
    alg4 = Algorithm('CDx', cdx=cdx)

    ent_pub = np.array([
        p['p_ava_ssm_pub'] * p['p_ava_xpert_pub'],
        p['p_ava_ssm_pub'] * (1 - p['p_ava_xpert_pub']),
        (1 - p['p_ava_ssm_pub']) * p['p_ava_xpert_pub'],
        (1 - p['p_ava_ssm_pub']) * (1 - p['p_ava_xpert_pub'])
    ])

    ent_eng = np.array([
        p['p_ava_xpert_eng'],
        (1 - p['p_ava_xpert_eng'])
    ])

    ent = np.array([
        p['p_csi_pub'],
        (1 - p['p_csi_pub']) * p['p_csi_ppm'],
        (1 - p['p_csi_pub']) * (1 - p['p_csi_ppm'])
    ])

    public = Sector(ent_pub, [alg1, alg2, alg3, alg4])
    engaged = Sector(ent_eng, [alg3, alg4])
    private = Sector(np.array([1]), [alg4])

    return System(ent, public, engaged, private)


if __name__ == '__main__':
    p0 = {
        'sens_ssm': 0.64,
        'spec_ssm': 0.98,
        'sens_xpert': 0.85,
        'sens_xpert_ss-': 0.64,
        'spec_xpert': 0.98,
        'p_ava_xpert_pub': 0.2,
        'p_ava_ssm_pub': 0.8,
        'p_ava_xpert_eng': 0.3,
        'p_loss_sputum': 0.15,
        'sens_cdx': 0.7,
        'spec_cdx': 0.95,
        'p_csi_pub': 0.483,
        'p_csi_ppm': 0.6,
        'dur_pri': 0.7,
        'dur_pub': 0.5,
        'p_txi_pub': 0.9,
        'p_txi_eng': 0.8,
        'p_txi_pri': 0.8,
        'p_refer_i2u': 0.3,
    }

    system = get_system(p0)

    system.Public.seek_care(1e4, 1e5).print()
    system.Public.seek_care_sto(1e4, 1e5).print()

    res = system.seek_care(1e4, 1e5)
    for k, v in res.items():
        print(k)
        v.print()

    print('---------------------------------------------')
    res = system.seek_care_sto(1e4, 1e5)
    for k, v in res.items():
        print(k)
        v.print()

    system = get_system(p0, has_cdx=False)
    res = system.seek_care(1000, 1000)

    print('---------------------------------------------')
    print('Public')
    for alg in system.Public.Algorithms:
        print('Algorithm ', alg.Key)
        res0 = alg.dx(1000, 0)
        res1 = alg.dx(0, 1000)
        print('-- TP_Bac: ', res0.TruePos / 1000)
        print('-- FN_Bac: ', res0.FalseNeg / 1000)
        print('-- T_SSM: ', res0['N_test_SSM'] / 1000)
        print('-- T_Xpert: ', (res0['N_test_Xpert_ss-'] + res0['N_test_Xpert']) / 1000)
        print('-- FP_Bac: ', res1.FalsePos / 1000)
        print('-- F_SSM: ', res1['N_test_SSM'] / 1000)
        print('-- F_Xpert: ', (res1['N_test_Xpert_ss-'] + res1['N_test_Xpert']) / 1000)
