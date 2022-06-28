import numpy as np
import pandas as pd
import json


src = pd.read_csv('../data/shifting/shifting.csv')

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


print(n_dx0 / n_entry)
print(n_dx1 / n_tr)

d = {
    'N_Dx0': n_dx0,
    'N_Entry0': n_entry,
    'N_Dx1': n_dx1.tolist(),
    'N_Tr': n_tr.tolist()
}
print(d)

with open('../data/cascade/pdx.json', 'w') as f:
    json.dump(d, f)
