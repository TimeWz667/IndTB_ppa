import pandas as pd
from ppa import *
from sims_pars.fitting import ApproxBayesComSMC
import os

__author__ = 'Chu-Chang Ku'


dat = pd.read_csv('../data/cascade/d_cascade_2019.csv')
dat = {row['State']: row.to_dict() for _, row in dat.iterrows()}

src_shf = pd.read_csv('../data/shifting/shifting.csv')


alg = ApproxBayesComSMC(max_round=40, n_collect=300, n_core=4, verbose=2)


if __name__ == '__main__':
    import json

    cnr_year = 2019

    folder = f'../out/cascade_{cnr_year}'

    os.makedirs(folder, exist_ok=True)

    for loc, d in dat.items():

        print(loc)

        obj = ObjRates(d)
        alg.fit(obj)
        post = alg.Collector

        ms = [Rates(d, dict(p)) for p in post.ParameterList]

        cascades = []

        for rates in ms:
            ex_det = rates.DetR
            ex_det = ex_det / ex_det.sum()
            ex_vis = 2.75

            shf = find_shifting(src_shf, d['Pr_Pub_CSI'], ex_det=ex_det, ex_vis=ex_vis)

            cascade = bind_cascade(shf, rates)

            cascades.append(cascade.dict())


        with open(f'{folder}/Cascade_{loc}.json', 'w') as f:
            json.dump(cascades, f)
