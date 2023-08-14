import numpy as np
import json
from sims_pars import sample


__author__ = 'Chu-Chang Ku'
__all__ = ['get_adr', 'attach_pars_cas', 'restructure_pars_cas', 'calc_stats']


def get_adr(targets):
    targets = targets[targets.Index == "IncR"]
    targets = targets[targets.Tag == "All"]
    targets = targets[targets.Year <= 2020]
    adr = - np.mean(np.diff(np.log(targets.M)))
    return adr


def _attach_pars_cas(p, p_a, p_s, p_c, adr):
    mu_a = p['r_die_asym'] + p['r_sc']
    mu_s = mu_c = p['r_die_sym'] + p['r_sc']
    p['adr'] = adr

    txis = np.array([p[f'txi_{sec}'] for sec in ['pub', 'eng', 'pri']])
    txi = txis.sum()

    p_txi = np.array([p[f'p_txi_{sec}'] for sec in ['pub', 'eng', 'pri']])

    det = txis / p_txi
    ppm = det[1] / det[1:].sum()
    p_ent_pub = p['p_csi_pub']
    p_ent = np.array([p_ent_pub, (1 - p_ent_pub) * ppm, (1 - p_ent_pub) * (1 - ppm)])

    p_dx = det / p_ent
    p_dx *= p['p_dx_pub'] / p_dx[0]

    p['p_dx_pub'], p['p_dx_eng'], p['p_dx_pri'] = p_dx

    p_dx_all = (p_txi * p_dx * p_ent).sum()

    r_det = txi / p_c
    p['r_csi'] = r_csi = p_c * (r_det + mu_c - adr) / p_s
    p['r_onset'] = r_onset = p_s * (r_csi + mu_s - adr) / p_a

    txi0 = p_s * r_csi * p_dx_all
    txi1 = max(txi - txi0, 0)
    p['r_recsi'] = txi1 / p_c / p_dx_all

    p['inc'] = (mu_a + r_onset - adr) * p_a


def _validate_cas(p):
    if p['r_recsi'] <= 0:
        return False

    pdx = np.array([p[f'p_dx_{sector}'] for sector in ['pub', 'eng', 'pri']])
    if np.any(pdx > 1):
        return False
    if np.any(pdx < 0):
        return False

    return True


def attach_pars_cas(js, bn, adr, max_attempt=1):
    prev = js['prev']
    p_a, p_s, p_c = prev['PrevAsym'], prev['PrevSym'], prev['PrevExCS']
    exo = js['exo']

    for p in js['pars']:
        n_attempt = 0

        while True:
            p.update(sample(bn, exo))
            _attach_pars_cas(p, p_a, p_s, p_c, adr)
            n_attempt += 1
            is_valid = _validate_cas(p)
            if is_valid or n_attempt >= max_attempt:
                break
        p['is_valid'] = is_valid


def restructure_pars_cas(p):
    pc = {k: p[k] for k in ['r_onset', 'r_csi', 'r_recsi', 'r_sc', 'r_die_asym', 'r_die_sym']}
    alo = np.array([[1, 0], [p['p_pri_on_pub'], 1 - p['p_pri_on_pub']], [0, 1]])

    p_ent_pub, ppm = p['p_csi_pub'], p['ppm']
    p_ent = np.array([p_ent_pub, (1 - p_ent_pub) * ppm, (1 - p_ent_pub) * (1 - ppm)])

    pc.update({
        'p_ent': p_ent,
        'p_dx': np.array([p[f'p_dx_{sector}'] for sector in ['pub', 'eng', 'pri']]),
        'p_txi': np.array([p[f'p_txi_{sector}'] for sector in ['pub', 'eng', 'pri']]),
        'tx_alo': alo,
        'tx_dur': np.array([0.5, p['dur_pri']]),
        'r_txs': np.array([p[f'r_txs_{sector}'] for sector in ['pub', 'pri']]),
        'r_txd': np.array([p[f'r_txd_{sector}'] for sector in ['pub', 'pri']]),
        'r_txl': np.array([p[f'r_txl_{sector}'] for sector in ['pub', 'pri']]),
        'ppv': np.array([p[f'ppv_{sector}'] for sector in ['pub', 'eng', 'pri']]),
    })

    return pc


def calc_stats(p):
    inc = p['inc']
    adr = p['adr']
    p = restructure_pars_cas(p)

    p_dx = p['p_ent'] * p['p_dx'] * p['p_txi']
    p_dx_all = p_dx.sum()

    txi = ((p_dx / p_dx_all).reshape((3, -1)) * p['tx_alo']).sum(0)

    p_a = inc / (p['r_die_asym'] + p['r_sc'] + p['r_onset'] - adr)
    p_s = p['r_onset'] * p_a / (p['r_die_sym'] + p['r_sc'] + p['r_csi'] - adr)
    p_c = p['r_csi'] * (1 - p_dx_all) * p_s / (p['r_die_sym'] + p['r_sc'] + p['r_recsi'] * p_dx_all - adr)

    mor = p_a * p['r_die_asym'] + (p_s + p_c) * p['r_die_sym']

    dur_a = 1 / (p['r_die_asym'] + p['r_sc'] + p['r_onset'])
    dur_s = 1 / (p['r_die_sym'] + p['r_sc'] + p['r_csi'])
    dur_c = 1 / (p['r_die_sym'] + p['r_sc'] + p['r_recsi'] * p_dx_all)

    p_a_die = p['r_die_asym'] * dur_a
    p_a_sc = p['r_sc'] * dur_a
    p_a_prog = p['r_onset'] * dur_a

    p_s_die = p['r_die_sym'] * dur_s
    p_s_sc = p['r_sc'] * dur_s
    p_s_prog = p['r_csi'] * (1 - p_dx_all) * dur_s
    p_s_det = p['r_csi'] * p_dx_all * dur_s
    p_s_tu = p_s_det * txi[0] / txi.sum()
    p_s_ti = p_s_det * txi[1] / txi.sum()

    p_c_die = p['r_die_sym'] * dur_c
    p_c_sc = p['r_sc'] * dur_c
    p_c_det = p['r_recsi'] * p_dx_all * dur_c
    p_c_tu = p_c_det * txi[0] / txi.sum()
    p_c_ti = p_c_det * txi[1] / txi.sum()

    p_to_tu = np.array([p['r_txs'][0], p['r_txl'][0], p['r_txd'][0]])
    p_to_tu /= p_to_tu.sum()
    p_to_ti = np.array([p['r_txs'][1], p['r_txl'][1], p['r_txd'][1]])
    p_to_ti /= p_to_ti.sum()

    p_ts_die = (txi * p['tx_dur'] * p['r_txd'])
    p_ts_ltfu = (txi * p['tx_dur'] * p['r_txl'])
    p_ts_succ = (txi * p['tx_dur'] * p['r_txs'])
    p_t_die = p_ts_die.sum()
    p_t_ltfu = p_ts_ltfu.sum()
    p_t_succ = p_ts_succ.sum()

    drop_die_a = p_a_die
    drop_die_s = p_a_prog * p_s_die
    drop_die_c = p_a_prog * p_s_prog * p_c_die

    drop_sc_a = p_a_sc
    drop_sc_s = p_a_prog * p_s_sc
    drop_sc_c = p_a_prog * p_s_prog * p_c_sc

    cas_onset = p_a_prog
    cas_csi = cas_onset * (p_s_prog + p_s_det)
    cas_txi = cas_onset * (p_s_prog * p_c_det + p_s_det)
    cas_txs = cas_txi * p_t_succ

    drop_die_t = cas_txi * p_t_die
    drop_ltfu_t = cas_txi * p_t_ltfu

    delay_pat = dur_s
    delay_sys = dur_c * (1 - p_dx_all)

    return {
        'Inc': inc,
        'PrevA': p_a,
        'PrevS': p_s,
        'PrevC': p_c,
        'MorUt': mor,
        'DelayPat': delay_pat,
        'DelaySys': delay_sys,
        'DelayTot': (delay_pat + delay_sys),
        'DropDieA': drop_die_a,
        'DropDieS': drop_die_s,
        'DropDieC': drop_die_c,
        'DropDieT': drop_die_t,
        'DropSelfCureA': drop_sc_a,
        'DropSelfCureS': drop_sc_s,
        'DropSelfCureC': drop_sc_c,
        'DropLTFU': drop_ltfu_t,
        'p_a2die': p_a_die,
        'p_a2sc': p_a_sc,
        'p_a2s': p_a_prog,

        'p_s2die': p_s_die,
        'p_s2sc': p_s_sc,
        'p_s2e': p_s_prog,
        'p_s2tu': p_s_tu,
        'p_s2ti': p_s_ti,

        'p_c2die': p_c_die,
        'p_c2sc': p_c_sc,
        'p_c2tu': p_c_tu,
        'p_c2ti': p_c_ti,

        'p_tu2txs': p_to_tu[0],
        'p_tu2txl': p_to_tu[1],
        'p_tu2txd': p_to_tu[2],
        'p_ti2txs': p_to_ti[0],
        'p_ti2txl': p_to_ti[1],
        'p_ti2txd': p_to_ti[2],

        'CasOnset': cas_onset,
        'CasCSI': cas_csi,
        # 'CasDet': cas_csi / p_txi_all,
        'CasTxI': cas_txi,
        'CasTxS': cas_txs,
        'CDR': cas_txi
    }


if __name__ == '__main__':
    from sims_pars import bayes_net_from_script
    import pandas as pd

    with open('../data/prior.txt', 'r') as f:
        prior = bayes_net_from_script(f.read())

    adr = get_adr(pd.read_csv('../data/India/targets.csv'))

    with open('../data/Gujarat/pars_txo.json', 'r') as f:
        js0 = json.load(f)

    attach_pars_cas(js0, prior, adr)

    p0 = js0['pars'][5]
    p1 = restructure_pars_cas(p0)

    for k, v in p1.items():
        print(k, ': ', v)

    inc = p0['inc0']
    stats = calc_stats(p0)

    print(f'Inc: {inc * 1e5:.1f}')
    print(f'Mor: {stats["MorUt"] * 1e5:.1f}')

    print(f'P_A: {js0["prev"]["PrevAsym"] * 1e5:.1f}, {stats["PrevA"] * 1e5:.1f}')
    print(f'P_S: {js0["prev"]["PrevSym"] * 1e5:.1f}, {stats["PrevS"] * 1e5:.1f}')
    print(f'P_C: {js0["prev"]["PrevExCS"] * 1e5:.1f}, {stats["PrevC"] * 1e5:.1f}')

    for k, v in stats.items():
        print(k, ': ', v)
