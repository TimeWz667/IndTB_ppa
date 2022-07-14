import pandas as pd
from sims_pars.fitting import ApproxBayesComSMC
from ppa.drugsale import *


dat = pd.read_csv('../data/cascade/d_cascade_2019.csv')
dat = {row['State']: row for _, row in dat.iterrows()}


alg = ApproxBayesComSMC(max_round=35, n_collect=500, n_core=4, verbose=0)


if __name__ == '__main__':
    import json
    import os

    os.makedirs("../out/cds0_pars", exist_ok=True)
    os.makedirs("../out/cds1_pars", exist_ok=True)

    locs = ['India']
    for loc in locs:
        d = dat[loc]
    #
    # for loc, d in dat.items():
        print(loc)

        obj = ObjForward(d)

        obj.UsingDrugSale = False
        alg.fit(obj)
        sss = [obj.simulate(p) for p in alg.Collector.ParameterList]
        css = [s.json() for s in sss]
        ess = [s.epi() for s in sss]

        with open(f'../out/cds0_pars/pars_{loc}.json', 'w') as f:
            json.dump(css, f)
        with open(f'../out/cds0_pars/fit_{loc}.json', 'w') as f:
            json.dump(ess, f)

        obj.UsingDrugSale = True
        alg.fit(obj)
        sss = [obj.simulate(p) for p in alg.Collector.ParameterList]
        css = [s.json() for s in sss]
        ess = [s.epi() for s in sss]

        with open(f'../out/cds1_pars/pars_{loc}.json', 'w') as f:
            json.dump(css, f)
        with open(f'../out/cds1_pars/fit_{loc}.json', 'w') as f:
            json.dump(ess, f)
