import pandas as pd
from sims_pars.fitting import ApproxBayesComSMC
from ppa import *


dat = pd.read_csv('../data/cascade/d_cascade_2019.csv')
dat = {row['State']: row for _, row in dat.iterrows()}


alg = ApproxBayesComSMC(max_round=35, n_collect=500, n_core=4, verbose=0)

src_shf = pd.read_csv('../data/shifting/shifting.csv')


if __name__ == '__main__':
    import json
    import os

    os.makedirs("../out/cft_pars", exist_ok=True)
    os.makedirs("../out/cf_pars", exist_ok=True)

    for loc, d in dat.items():
        print(loc)
        shf = Shifting(src_shf, d.Pr_Pub_CSI)

        obj = ObjForward(d, shf)

        # obj.UsingCNR = True
        # obj.UsingPrevTx = False
        # alg.fit(obj)
        # css = [obj.simulate(p).json() for p in alg.Collector.ParameterList]
        #
        # with open(f'../out/cf_pars/pars_{loc}.json', 'w') as f:
        #     json.dump(css, f)
        #
        obj.UsingCNR = True
        obj.UsingPrevTx = True
        alg.fit(obj)
        css = [obj.simulate(p).json() for p in alg.Collector.ParameterList]

        with open(f'../out/cft_pars/pars_{loc}.json', 'w') as f:
            json.dump(css, f)
