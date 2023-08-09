import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['attach_pars_tx']

#
# def extract_txo(targets):
#     tx_die = targets[targets.Index == 'TxDie']
#     tx_succ = targets[targets.Index == 'TxSucc']
#
#     txo = dict()
#
#     sel = tx_die[tx_die.Tag == 'Pub']
#     txi = sel.N.sum()
#     txd = (sel.M * sel.N).sum()
#
#     sel = tx_succ[tx_succ.Tag == 'Pub']
#     txs = (sel.M * sel.N).sum()
#     txl = txi - txs - txd
#
#     txo['pub'] = {'txi': txi, 'txd': txd, 'txs': txs, 'txl': txl}
#
#     sel = tx_die[tx_die.Tag == 'Pri']
#     txi = sel.N.sum()
#     txd = (sel.M * sel.N).sum()
#
#     sel = tx_succ[tx_succ.Tag == 'Pri']
#     txs = (sel.M * sel.N).sum()
#     txl = txi - txs - txd
#
#     txo['eng'] = {'txi': txi, 'txd': txd, 'txs': txs, 'txl': txl}
#     return txo


def attach_pars_tx(js):
    txo = js['txo']
    ps_pub = np.array([txo['pub'][o] for o in ['txs', 'txd', 'txl']]) / txo['pub']['txi']

    for p in js['pars']:
        txi_ep = txo['eng']['txi'] * p['p_pri_on_pub']
        tx_eng_ei = np.array([txo['eng'][o] for o in ['txs', 'txd', 'txl']]) - txi_ep * ps_pub
        # tx_eng_ei[tx_eng_ei < 0] = 0
        ps_pri = tx_eng_ei / (txo['eng']['txi'] * (1 - p['p_pri_on_pub']))

        p['r_txs_pub'], p['r_txd_pub'], p['r_txl_pub'] = ps_pub / 0.5
        p['r_txs_pri'], p['r_txd_pri'], p['r_txl_pri'] = ps_pri / p['dur_pri']

        # Public is always better assumption
        if p['r_txd_pri'] < p['r_txd_pub']:
            diff = p['r_txd_pub'] - p['r_txd_pri']
            p['r_txd_pri'] = p['r_txd_pub']
            p['r_txl_pri'] -= diff

        if p['r_txs_pri'] > p['r_txs_pub']:
            diff = p['r_txs_pri'] - p['r_txs_pub']
            p['r_txs_pri'] = p['r_txs_pub']
            p['r_txl_pri'] += diff


if __name__ == '__main__':
    import json
    import pandas as pd

    loc = 'Delhi'
    js = json.load(open(f'../data/{loc}/pars_under.json', 'r'))
    targets0 = pd.read_csv(f'../data/{loc}/targets.csv')
    attach_pars_tx(js)
    print(loc, '\t ', sum([p['r_txs_pri'] > 0 for p in js['pars']]) / len(js['pars']))

    txo = js['txo']
    p0 = js['pars'][6]

    rates_pub = np.array([p0[f'r_{k}_pub'] for k in ['txs', 'txd', 'txl']])
    rates_pri = np.array([p0[f'r_{k}_pri'] for k in ['txs', 'txd', 'txl']])
    print('Rates, public', rates_pub)
    print('Rates, private', rates_pri)

    print('Validation')
    print('--Public')
    print('Expected: ', np.array([txo['pub'][o] for o in ['txs', 'txd', 'txl']]) / txo['pub']['txi'])
    print('Simulated: ', rates_pub / rates_pub.sum())

    print('--Engaged private')
    print('Expected: ', np.array([txo['eng'][o] for o in ['txs', 'txd', 'txl']]) / txo['eng']['txi'])
    s0 = rates_pri / rates_pri.sum() * (1 - p0['p_pri_on_pub'])
    s0 += rates_pub / rates_pub.sum() * p0['p_pri_on_pub']
    print('Simulated: ', s0)
