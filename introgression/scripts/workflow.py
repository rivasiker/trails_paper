from gwf import Workflow
from os.path import exists
import numpy as np

gwf = Workflow()

mu = 2e-8
g = 25
N_AB = 10_000*2
N_ABC = 10_000*2
N_ref = N_ABC
t_1 = 60_000/g
t_2 = (600_000-t_1)/g
t_3 = (6_000_000-t_2)/g
t_A = t_1
t_B = t_1
t_C = t_1+t_2
r = 1e-8/mu

t_m = t_1-40_000/g
m = 0.05

dct = {1:8, 3:57, 5:121}

tot_lst = []
t_A = t_1
t_B = t_1
t_C = t_1+t_2

model = 'introgression_error_model'

for n_int_AB in [1]:
    for n_int_ABC in [1, 3, 5]:
        for seed in range(1, 6):
            # if exists(f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv'):
            #     pass
            tot_lst.append('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model))
            gwf.target('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model),
                inputs=['optimize_introgression.py'],
                outputs=['../results/{}_{}_{}_{}_{}.csv'.format(x, n_int_AB, n_int_ABC, seed, model) for x in ['trees', 'sim']],
                cores=n_int_ABC,
                memory='{}g'.format(n_int_ABC*4),
                walltime= '{}:00:00'.format(dct[n_int_ABC]),
                account='Primategenomes') << f"""
            python optimize_introgression.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {model} {t_m} {m}
            """

#[print(i, end = ' ') for i in tot_lst]
